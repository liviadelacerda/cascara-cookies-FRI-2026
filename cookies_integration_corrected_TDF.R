# =============================================================================
# Corrigendum - Food Research International
# "Effects of post-harvest processing on chemical composition and sensory
#  quality of coffee cascara, using clean-label cookies as a model bakery system"
#
# Multiblock integration: CA of CATA data, Multiple Factor Analysis (MFA,
# Fig. 7) and PLS-R models (Supplementary Table S6), recomputed with the
# corrected total dietary fiber (TDF) values of NAT and PARCH.
#
# Input files (same folder as this script, values only):
#   Liking_CATA_cookies.xlsx              - consumer liking and CATA data
#   Volatiles_relative_areas_top12.xlsx   - relative areas of selected volatiles
#   NonVolatiles_cascara_corrected.xlsx   - proximate composition and
#                                           non-volatiles (corrected TDF, column FA)
#   Texture_hardness.xlsx                 - cookie hardness
#
# Outputs: export_fig/ (figures; mfa-biplot-clean.png = corrected Fig. 7)
#          exports_check/ (intermediate tables)
#
# Usage: set the working directory to this folder and source the script, e.g.
#   setwd("path/to/Data_and_scripts")
#   source("cookies_integration_corrected_TDF.R")
#
# R packages: dplyr, tidyr, purrr, stringr, readxl, readr, janitor, FactoMineR,
# factoextra, ade4, mixOmics, ggplot2, ggrepel, patchwork, tibble, knitr,
# flextable, officer
# =============================================================================

## -----------------------------------------------------------------------------
#| label: setup
#| message: false
#| warning: false
if (!interactive()) grDevices::pdf(NULL)  # avoid Rplots.pdf when run with Rscript
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
  library(readxl); library(janitor)
  library(FactoMineR); library(factoextra)
  library(ade4)         # RV / Procrustes
  library(mixOmics)     # PLS
})
set.seed(123)



## -----------------------------------------------------------------------------
#| label: params
# Caminho-base (todas as planilhas estão nesta pasta)
base_dir <- getwd()  # folder containing this script and the four .xlsx data files

paths <- list(
  # 1) Sensory (CATA + liking), consumer-level data
  sensory_candidates = file.path(base_dir, "Liking_CATA_cookies.xlsx"),
  sensory_sheet_candidates = c("Sheet1"),
  # 2) Volatiles (GC-MS), relative areas of the 12 selected compounds
  volatiles_candidates = file.path(base_dir, "Volatiles_relative_areas_top12.xlsx"),
  volatiles_sheet_candidates = c("Sheet1"),
  # 3) Non-volatiles / proximate composition (corrected TDF of NAT and PARCH)
  nonvol_file  = file.path(base_dir, "NonVolatiles_cascara_corrected.xlsx"),
  nonvol_sheet = "novolatiles_bs"
)




## -----------------------------------------------------------------------------
#| label: utils
#| message: false
#| warning: false

# retorna o primeiro caminho existente em um vetor de candidatos
pick_first_existing <- function(v){
  v <- v[!is.na(v)]
  hit <- which(file.exists(v))[1]
  if (length(hit) && !is.na(hit)) v[hit] else NA_character_
}

# escolhe a primeira planilha existente dentre candidatos (por nome ou índice)
pick_first_sheet <- function(path, sheet_candidates){
  shs <- tryCatch(readxl::excel_sheets(path), error = function(e) character(0))
  for (s in sheet_candidates) {
    if ((is.character(s) && s %in% shs) || (is.numeric(s) && length(shs) >= s)) return(s)
  }
  if (length(shs)) shs[1] else 1
}

# padroniza rótulos de amostra
std_sample <- function(x) toupper(trimws(as.character(x)))

# imputação simples por mediana (numéricas)
impute_median <- function(df){
  df[] <- lapply(df, function(col){
    if (is.numeric(col)) { col[is.na(col)] <- stats::median(col, na.rm = TRUE); col } else col
  })
  df
}

# utilitário para pegar a primeira coluna existente num conjunto de candidatos
pick_col <- function(df, candidates){
  hit <- candidates[candidates %in% names(df)][1]
  if (is.na(hit)) NA_character_ else hit
}



## -----------------------------------------------------------------------------
#| label: ingest-sensory
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(knitr)
})

# -- helpers (só define se não existirem) -------------------------------------
if (!exists("pick_first_existing")) {
  pick_first_existing <- function(paths) {
    paths <- paths[!is.na(paths)]
    paths <- paths[file.exists(paths)]
    if (length(paths)) paths[1] else NA_character_
  }
}
if (!exists("pick_first_sheet")) {
  pick_first_sheet <- function(path, sheet_candidates) {
    for (s in sheet_candidates) {
      ok <- tryCatch({ readxl::read_excel(path, sheet = s, n_max = 1); TRUE }, error = function(e) FALSE)
      if (ok) return(s)
    }
    stop("Nenhuma aba válida encontrada entre os candidatos.")
  }
}
if (!exists("pick_col")) {
  pick_col <- function(d, candidates){
    nm <- names(d)
    hit <- nm[match(tolower(candidates), tolower(nm), nomatch = 0)]
    hit <- hit[hit != 0]
    if (length(hit)) hit[1] else NA_character_
  }
}
if (!exists("std_sample")) {
  std_sample <- function(x) {
    x <- trimws(as.character(x))
    toupper(x)
  }
}

# 1) arquivo e aba
sens_path  <- pick_first_existing(paths$sensory_candidates)
if (is.na(sens_path)) stop("Sensory file not found in 'paths$sensory_candidates'.")
sens_sheet <- pick_first_sheet(sens_path, paths$sensory_sheet_candidates)

raw_sens <- readxl::read_excel(sens_path, sheet = sens_sheet)

# 2) padronização mínima (consumer, sample, liking, cluster + CATA 0/1)
idc   <- pick_col(raw_sens, c("consumer","avaliador","Avaliador","ID","id"))
smpc  <- pick_col(raw_sens, c("sample","amostra","tratamento","Tratamento","Sample"))
likec <- pick_col(raw_sens, c("liking","Nota Hedônica","nota_hedonica","Hedonic","Hedonica","liking_score"))
cluc  <- pick_col(raw_sens, c("cluster","Cluster","k2","kcluster","cluster_k2","grupo","Group","segment","seg","grp"))

if (any(is.na(c(idc, smpc, likec)))) {
  stop("Colunas essenciais ausentes (consumer/sample/liking). Verifique nomes na planilha.")
}

d_sens <- raw_sens %>%
  dplyr::rename(consumer = !!idc, sample = !!smpc, liking = !!likec) %>%
  dplyr::mutate(
    consumer = as.character(consumer),
    sample   = std_sample(sample),
    liking   = suppressWarnings(as.numeric(liking)),
    cluster  = if (!is.na(cluc) && cluc %in% names(raw_sens)) as.character(.data[[cluc]]) else NA_character_
  )

# detectar atributos CATA binários (0/1)
fixed <- c("consumer","sample","liking","cluster")
cand  <- setdiff(names(d_sens), fixed)
is_bin01 <- function(x){
  v <- suppressWarnings(as.numeric(x))
  ux <- sort(unique(na.omit(v)))
  length(ux) <= 2 && all(ux %in% c(0,1))
}
attr_cols <- cand[vapply(d_sens[cand], is_bin01, logical(1))]
if (!length(attr_cols)) stop("Nenhum atributo CATA binário (0/1) foi detectado.")
d_sens[attr_cols] <- lapply(d_sens[attr_cols], function(x) as.integer(suppressWarnings(as.numeric(x))))

# 3) definir cluster se não existir (k=2 via k-means no perfil de liking)
has_two_levels <- function(x){ ux <- unique(x); ux <- ux[!is.na(ux)]; length(ux) >= 2 }
if (!has_two_levels(d_sens$cluster)) {
  wide <- d_sens %>%
    dplyr::group_by(consumer, sample) %>%
    dplyr::summarise(like = mean(liking, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = sample, values_from = like)

  num_cols <- setdiff(names(wide), "consumer")
  if (length(num_cols) >= 1 && nrow(wide) >= 2) {
    for (cc in num_cols) {
      v <- wide[[cc]]; med <- suppressWarnings(stats::median(v, na.rm = TRUE))
      wide[[cc]][is.na(v)] <- if (is.finite(med)) med else 0
    }
    set.seed(123)
    km <- stats::kmeans(scale(as.matrix(wide[, num_cols, drop = FALSE])), centers = 2, nstart = 25)
    d_sens <- d_sens %>%
      dplyr::left_join(tibble::tibble(consumer = wide$consumer, cluster = paste0("C", km$cluster)),
                       by = "consumer", suffix = c("", ".km")) %>%
      dplyr::mutate(cluster = ifelse(is.na(.data$cluster), .data$cluster.km, .data$cluster)) %>%
      dplyr::select(-cluster.km)
  }
}

d_sens <- d_sens %>%
  dplyr::mutate(
    cluster = factor(cluster),
    sample  = factor(sample),
    liking  = as.numeric(liking)
  )

# 4) saídas resumidas (unidade: amostra × cluster)
liking_summ <- d_sens %>%
  dplyr::group_by(sample, cluster) %>%
  dplyr::summarise(
    mean_liking = mean(liking, na.rm = TRUE),
    sd_liking   = sd(liking,   na.rm = TRUE),
    n           = dplyr::n(),
    .groups = "drop"
  )

cata_props <- d_sens %>%
  dplyr::select(consumer, sample, cluster, dplyr::all_of(attr_cols)) %>%
  tidyr::pivot_longer(-c(consumer, sample, cluster), names_to = "attr", values_to = "pres") %>%
  dplyr::group_by(sample, cluster, attr) %>%
  dplyr::summarise(prop = mean(pres, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = attr, values_from = prop, names_prefix = "prop_")

penalty_by_cluster <- d_sens %>%
  dplyr::group_by(cluster) %>%
  dplyr::group_modify(function(df, key){
    purrr::map_dfr(attr_cols, function(att){
      m1 <- mean(df$liking[df[[att]]==1], na.rm=TRUE)
      m0 <- mean(df$liking[df[[att]]==0], na.rm=TRUE)
      tibble::tibble(attr = att, delta = m1 - m0)
    })
  }) %>% dplyr::ungroup()

penalty_index <- cata_props %>%
  tidyr::pivot_longer(dplyr::starts_with("prop_"), names_to = "attr", values_to = "prop") %>%
  dplyr::mutate(attr = sub("^prop_", "", attr)) %>%
  dplyr::left_join(penalty_by_cluster, by = c("cluster","attr")) %>%
  dplyr::group_by(sample, cluster) %>%
  dplyr::summarise(penalty_index = sum(prop * delta, na.rm = TRUE), .groups = "drop")

# ++ novos objetos esperados por outros chunks:
dados_cata <- d_sens %>% dplyr::select(consumer, sample, dplyr::all_of(attr_cols))
dados_clusterizados <- d_sens %>% dplyr::select(consumer, sample, Cluster = cluster, liking)

# (opcional) previews
# knitr::kable(head(dados_cata), caption = "dados_cata")
# knitr::kable(head(dados_clusterizados), caption = "dados_clusterizados")
# knitr::kable(liking_summ, caption = "Médias de aceitação por amostra × cluster")



## -----------------------------------------------------------------------------
#| label: ca-cata
#| message: false
#| warning: false

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(FactoMineR); library(factoextra) })

# Seleção de colunas de CATA: filtradas se disponíveis
prop_cols <- if (exists("cata_attrs_keep")) {
  keep <- paste0("prop_", cata_attrs_keep)
  keep[keep %in% names(cata_props)]
} else {
  grep("^prop_", names(cata_props), value = TRUE)
}

cata_mat <- cata_props %>%
  dplyr::mutate(sample = std_sample(sample)) %>%
  tidyr::unite("sample_cluster", sample, cluster, sep = " :: ", remove = FALSE) %>%
  dplyr::select(sample, cluster, sample_cluster, dplyr::all_of(prop_cols))

# matriz X (linhas = sample::cluster; colunas = atributos)
Xcata <- cata_mat %>% dplyr::select(dplyr::all_of(prop_cols)) %>% as.data.frame()
# remove colunas constantes (sem variação)
if (ncol(Xcata) > 0) {
  sdv <- apply(Xcata, 2, sd, na.rm = TRUE)
  Xcata <- Xcata[, is.finite(sdv) & sdv > 0, drop = FALSE]
}
rownames(Xcata) <- cata_mat$sample_cluster

if (nrow(Xcata) >= 2 && ncol(Xcata) >= 2) {
  res_ca <- FactoMineR::CA(Xcata, graph = FALSE)

  # Plot (biplot simples)
  factoextra::fviz_ca_biplot(
    res_ca, repel = TRUE,
    col.row = "steelblue", col.col = "darkorange",
    title = "CA — CATA (samples × attributes)"
  )
} else {
  warning("CA skipped: insufficient rows/columns in CATA matrix (after filtering/constant removal).")
}



## -----------------------------------------------------------------------------
#| label: ca-cata-scores
#| message: false
#| warning: false

# Garante que o CA foi calculado
stopifnot(exists("res_ca"))

# ------------------ Scores de indivíduos (sample::cluster) ------------------
ca_scores <- as.data.frame(res_ca$row$coord) %>%
  tibble::rownames_to_column("row_id") %>%
  tidyr::separate(row_id, into = c("sample", "cluster"), sep = " :: ", remove = FALSE) %>%
  dplyr::mutate(
    sample  = factor(sample),
    cluster = factor(cluster)
  )

# ------------------ Loadings dos atributos ------------------
ca_loadings <- as.data.frame(res_ca$col$coord) %>%
  tibble::rownames_to_column("attribute")

# Exporta para uso global (necessário para o RV)
assign("ca_scores", ca_scores, envir = .GlobalEnv)
assign("ca_loadings", ca_loadings, envir = .GlobalEnv)



## -----------------------------------------------------------------------------
#| label: ingest-volatiles
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(stringr); library(janitor); library(tidyr)
})

# conversor numérico robusto (pt/eng, ND etc.)
coerce_num <- function(x){
  v <- as.character(x)
  v[tolower(v) %in% c("nd","n.d.","na","nan","--","-","")] <- NA
  v1 <- suppressWarnings(readr::parse_number(v))
  v2 <- suppressWarnings(readr::parse_number(
    v,
    locale = readr::locale(decimal_mark = ",", grouping_mark = ".")
  ))
  if (sum(!is.na(v2)) > sum(!is.na(v1))) v2 else v1
}

# 1) arquivo
vol_path <- pick_first_existing(paths$volatiles_candidates)
if (is.na(vol_path)) stop("Volatiles file not found in 'paths$volatiles_candidates'.")

# 2) escolher aba (se nada bater, usa a primeira)
sheets <- readxl::excel_sheets(vol_path)
sheet_candidates <- tryCatch(paths$volatiles_sheet_candidates, error = function(e) NULL)

sheet_ok <- NULL
if (!is.null(sheet_candidates)) {
  for (s in sheet_candidates) {
    if (is.character(s) && s %in% sheets) {
      sheet_ok <- s; break
    }
    if (is.numeric(s) && s >= 1 && s <= length(sheets)) {
      sheet_ok <- sheets[as.integer(s)]; break
    }
  }
}
if (is.null(sheet_ok)) sheet_ok <- sheets[1]

# 3) leitura — cabeçalho na PRIMEIRA linha (sem skip)
vol_raw <- readxl::read_excel(vol_path, sheet = sheet_ok, col_names = TRUE)

# precisa ter colunas: group, sample, seguida dos 12 compostos
if (!all(c("group","sample") %in% names(vol_raw))) {
  stop("A planilha de voláteis deve conter as colunas 'group' e 'sample'.")
}

# Estrutura esperada (exemplo):
# group   sample     2,3-butanediol  2-acetylpyrrole  ...  Isovaleric acid
# NAT     NAT_R1     ...

# 4) padronização: group = tratamento (NAT, PARCH...), sample_raw = NAT_R1 etc.
vol_repl <- vol_raw %>%
  dplyr::rename(group_id = group, sample_raw = sample) %>%
  dplyr::mutate(
    sample = toupper(trimws(as.character(group_id)))   # NAT, PARCH, SC4...
  ) %>%
  # nomes limpos para os compostos
  dplyr::rename_with(~ janitor::make_clean_names(.x),
                     -c(group_id, sample_raw, sample)) %>%
  # conversão numérica das colunas de compostos
  dplyr::mutate(
    dplyr::across(
      -c(group_id, sample_raw, sample),
      coerce_num
    )
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::mutate(rep = dplyr::row_number()) %>%  # cria ID de réplica (1,2,3)
  dplyr::ungroup()

# 5) remover colunas 100% NA e de variância zero
num_cols <- setdiff(names(vol_repl), c("group_id","sample_raw","sample","rep"))

# remove colunas totalmente NA
all_na <- vapply(vol_repl[num_cols], function(x) all(is.na(x)), logical(1))
num_cols <- num_cols[!all_na]

# remove colunas com variância zero (na amostra de réplicas)
num_ok <- vapply(
  vol_repl[num_cols],
  function(x) is.numeric(x) && stats::sd(x, na.rm = TRUE) > 0,
  logical(1)
)

vol_repl <- vol_repl %>%
  dplyr::select(sample, rep, dplyr::all_of(num_cols[num_ok]))

# 6) matriz agregada por amostra (média) — para JOIN posterior
vol_master_vars <- vol_repl %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(dplyr::across(-rep, ~ mean(.x, na.rm = TRUE)), .groups = "drop")

vol_varnames <- setdiff(names(vol_master_vars), "sample")


## -----------------------------------------------------------------------------
#| label: volatiles-preprocess
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(janitor); library(knitr)
})

# Garante que temos vol_repl vindo do chunk anterior
stopifnot(exists("vol_repl"))

# Mantém apenas sample, rep e compostos numéricos
num_cols <- setdiff(names(vol_repl), c("sample","rep"))
num_ok   <- vapply(
  vol_repl[num_cols],
  function(x) is.numeric(x) && stats::sd(x, na.rm = TRUE) > 0,
  logical(1)
)

vol_repl <- vol_repl %>%
  dplyr::select(sample, rep, dplyr::all_of(num_cols[num_ok]))

# Remover artefatos (se existirem nomes com esses termos)
drop_terms <- c("glycerin","phthalate","siloxane","squalene","unknown")
keep <- !Reduce(`|`, lapply(drop_terms,
                            function(rx) grepl(rx, names(vol_repl), ignore.case = TRUE)))
vol_repl <- vol_repl[, keep]

# -------- Pré-processamento: TIC -> log10 -> autoscale -------------
vol_mat <- vol_repl %>% dplyr::select(-sample, -rep)

# normalização por soma (TIC)
row_sums <- rowSums(vol_mat, na.rm = TRUE)
vol_tic  <- sweep(vol_mat, 1, ifelse(row_sums == 0, 1, row_sums), "/")

# log10(x + eps)
eps <- 1e-9
vol_log <- log10(vol_tic + eps)

# autoscale por coluna
col_means <- apply(vol_log, 2, mean, na.rm = TRUE)
col_sds   <- apply(vol_log, 2, sd,   na.rm = TRUE)
col_sds[col_sds == 0 | !is.finite(col_sds)] <- 1
vol_scaled <- sweep(sweep(vol_log, 2, col_means, "-"), 2, col_sds, "/")
vol_scaled <- cbind(vol_repl[, c("sample","rep")], as.data.frame(vol_scaled))

# filtro de prevalência (≥ 80% de valores finitos)
prev_thr <- 0.80
keep_prev <- vapply(
  vol_scaled[, setdiff(names(vol_scaled), c("sample","rep"))],
  function(v) mean(is.finite(v)) >= prev_thr,
  logical(1)
)
vol_scaled <- vol_scaled %>%
  dplyr::select(sample, rep, dplyr::all_of(names(keep_prev)[keep_prev]))

# Relatórios rápidos (opcionais)
knitr::kable(
  tibble::tibble(
    samples   = dplyr::n_distinct(vol_scaled$sample),
    replicas  = nrow(vol_scaled),
    variables = ncol(vol_scaled) - 2
  ),
  caption = "Volatiles matrix after TIC+log+autoscale and filtering"
)

knitr::kable(head(colnames(vol_scaled), 12),
             col.names = "first columns",
             caption = "Preview of column names in 'vol_scaled'")



## -----------------------------------------------------------------------------
#| label: fisher-top12-volatiles
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(stringr); library(knitr)
})

stopifnot(exists("vol_repl"))

# -------- utilitário para limpar nomes químicos --------
clean_chem_name <- function(x){
  x <- gsub("^x", "", x)              # remove prefixo x
  x <- gsub("_", " ", x)              # troca _ por espaço
  x <- gsub("([0-9]) ([a-zA-Z])", "\\1-\\2", x)  # 2 octenoic → 2-octenoic
  x <- gsub("([0-9]),([0-9])", "\\1,\\2", x)     # 2,3 pattern
  return(x)
}

# -------- detectar colunas principais --------
pick_col <- function(d, candidates){
  nm <- names(d)
  hit <- nm[match(tolower(candidates), tolower(nm), nomatch = 0)]
  hit[hit != 0][1]
}

dat <- vol_repl

sample_col <- pick_col(dat, c("Sample","Amostra","sample","amostra"))
rep_col    <- pick_col(dat, c("Rep","Replica","Replicate","rep","réplica","replicate"))

if (is.na(sample_col)) stop("Não encontrei a coluna de amostra em vol_repl.")
names(dat)[names(dat) == sample_col] <- "Sample"
if (!is.na(rep_col)) names(dat)[names(dat) == rep_col] <- "Rep"

if (!"Group" %in% names(dat)) {
  dat$Group <- dat$Sample
}

# garantir que não existam colunas list
is_list <- vapply(dat, is.list, logical(1))
if (any(is_list)) {
  for (nm in names(dat)[is_list]) {
    li <- dat[[nm]]
    ok <- all(vapply(li, function(v) length(v) %in% c(0,1), logical(1)))
    if (ok) {
      vec <- suppressWarnings(as.numeric(unlist(li)))
      if (!all(is.na(vec))) dat[[nm]] <- vec
    }
  }
}

# -------- matriz numérica --------
num_cols <- names(dat)[vapply(dat, is.numeric, logical(1))]
num_cols <- setdiff(num_cols, c("Rep"))
stopifnot(length(num_cols) >= 2)

M <- as.data.frame(dat[, num_cols, drop = FALSE])

if ("Rep" %in% names(dat)) {
  rownames(M) <- make.unique(paste0(dat$Sample, "_R", dat$Rep))
} else {
  rownames(M) <- make.unique(as.character(dat$Sample))
}

Group <- factor(dat$Group)

# -------- presença --------
presence_pct <- colMeans((as.matrix(M) > 0) & is.finite(as.matrix(M)), na.rm = TRUE) * 100

# -------- log1p + filtro variância --------
X <- log1p(as.matrix(M))
sd_all <- apply(X, 2, sd, na.rm = TRUE)
keep   <- is.finite(sd_all) & sd_all > 0
X <- X[, keep, drop = FALSE]
presence_pct <- presence_pct[colnames(X)]

# -------- Fisher ANOVA --------
fisher_F <- function(x, g){
  ok <- is.finite(x) & !is.na(g)
  x <- x[ok]; g <- droplevels(g[ok])
  k <- nlevels(g); if (k < 2) return(NA_real_)
  n_i <- as.numeric(table(g))
  mu  <- mean(x, na.rm = TRUE)
  mu_i <- tapply(x, g, mean, na.rm = TRUE)
  ssb  <- sum(n_i * (mu_i - mu)^2, na.rm = TRUE)
  ssw  <- sum(tapply(x, g, function(v) sum((v - mean(v, na.rm=TRUE))^2 )), na.rm = TRUE)
  df1 <- k - 1; df2 <- sum(n_i - 1)
  if (!is.finite(ssw) || df2 <= 0) return(NA_real_)
  (ssb/df1) / (ssw/df2)
}

Fvals <- vapply(seq_len(ncol(X)), function(j) fisher_F(X[, j], Group), numeric(1))

fisher_table <- tibble(
  var          = colnames(X),
  fisher       = as.numeric(Fvals),
  presence_pct = as.numeric(presence_pct)
) |>
  arrange(desc(fisher)) |>
  filter(!is.na(fisher))

# -------- top 12 --------
presence_min   <- 30
rare_keep_top  <- 1
top_n          <- 12
r_seq          <- seq(0.90, 0.99, by = 0.01)

pool <- fisher_table |> filter(presence_pct >= presence_min)
rare <- fisher_table |> filter(presence_pct < presence_min) |>
  slice_max(order_by = fisher, n = rare_keep_top, with_ties = FALSE)
pool <- bind_rows(pool, rare) |> distinct(var, .keep_all = TRUE)

vars_pool <- intersect(pool$var, colnames(X))
Xc <- X[, vars_pool, drop = FALSE]
C  <- suppressWarnings(cor(Xc, use = "pairwise.complete.obs"))

pick_by_r <- function(r_thr, vars_local, C_local, rank_tbl){
  if (length(vars_local) <= 1) return(vars_local)
  D <- 1 - abs(C_local[vars_local, vars_local, drop = FALSE]); D[!is.finite(D)] <- 1
  hc <- hclust(as.dist(D), method = "average")
  cl <- cutree(hc, h = 1 - r_thr)
  unlist(lapply(split(vars_local, cl), function(vs){
    vs[which.max(rank_tbl$fisher[match(vs, rank_tbl$var)])]
  }), use.names = FALSE)
}

sel <- character(0)
for (r_thr in r_seq) {
  cand <- pick_by_r(r_thr, vars_pool, C, pool)
  cand <- pool |>
    filter(var %in% cand) |>
    slice_max(order_by = fisher, n = top_n, with_ties = FALSE) |>
    pull(var)
  if (length(cand) >= top_n) { sel <- cand[1:top_n]; break }
}

if (length(sel) < top_n) {
  faltam <- top_n - length(sel)
  extras <- pool |> filter(!var %in% sel) |>
    slice_max(order_by = fisher, n = max(faltam,0)) |>
    pull(var)
  sel <- unique(c(sel, extras))[1:top_n]
}

vol_top12 <- sel
vol_top12_clean <- clean_chem_name(sel)   # <<< versão amigável

# -------- tabela final --------
knitr::kable(
  fisher_table |> slice(match(vol_top12, var)) |>
    transmute(
      compound = clean_chem_name(var),
      fisher   = sprintf("%.3f", fisher),
      presence = sprintf("%.1f", presence_pct)
    ),
  caption = "Top-12 volatiles by Fisher F (presence ≥ 30%, pruning by ≥0.90 correlation)"
)



## -----------------------------------------------------------------------------
#| label: ingest-nonvol
#| message: false
#| warning: false

# Arquivo e aba de não-voláteis / composição
nonvol_path  <- paths$nonvol_file
nonvol_sheet <- paths$nonvol_sheet
if (!file.exists(nonvol_path)) stop("Non-volatile composition file not found at 'paths$nonvol_file'.")

nonvol_raw <- readxl::read_excel(nonvol_path, sheet = nonvol_sheet)

nonvol_master <- nonvol_raw %>%
  janitor::clean_names()

# Detecta a coluna de amostra (tratamento) de forma robusta
sample_candidates <- c("tratamento","sample","amostra","cascara","cascara_id","tratamento_codigo")
sample_col <- sample_candidates[sample_candidates %in% names(nonvol_master)][1]
if (is.na(sample_col)) {
  stop(paste0("Não encontrei a coluna de amostra. Colunas disponíveis: ",
              paste(names(nonvol_master), collapse = ", ")))
}

# Padroniza nome da amostra e trata 'nd' -> NA, depois converte para numérico quando possível
nonvol_master <- nonvol_master %>%
  dplyr::rename(sample = dplyr::all_of(sample_col)) %>%
  dplyr::mutate(sample = std_sample(sample)) %>%
  dplyr::mutate(dplyr::across(where(is.character), ~ ifelse(tolower(.) %in% c("nd","n.d.","na"), NA_character_, .))) %>%
  dplyr::mutate(dplyr::across(-sample, ~ suppressWarnings(readr::parse_number(as.character(.))))) %>%
  # renomeia variáveis para nomes estáveis, se existirem
  dplyr::rename(
    protein       = dplyr::any_of("ptn"),
    lipids        = dplyr::any_of("lip"),
    tdf           = dplyr::any_of("fa"),
    caffeine      = dplyr::any_of("caf"),
    trigonelline  = dplyr::any_of("trig"),
    cga           = dplyr::any_of("5-cqa"),
    ncga          = dplyr::any_of("3-cqa"),
    tpc           = dplyr::any_of("tpc"),
    flav          = dplyr::any_of("flav"),
    anthocyanins  = dplyr::any_of("ant"),
    tannins       = dplyr::any_of("tan"),
    fructose      = dplyr::any_of("fruc"),
    glucose       = dplyr::any_of("gluc")
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), ~ mean(.x, na.rm = TRUE)), .groups = "drop")




## -----------------------------------------------------------------------------
#| label: ingest-texture
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(janitor); library(readr)
})

# caminho: mesma pasta base_dir das demais planilhas
texture_path  <- file.path(base_dir, "Texture_hardness.xlsx")
texture_sheet <- "para MFA"

if (!file.exists(texture_path)) {
  stop("Arquivo de textura não encontrado em: ", texture_path)
}

texture_raw <- readxl::read_excel(texture_path, sheet = texture_sheet)

texture_master <- texture_raw %>%
  janitor::clean_names() %>%          # -> sample, hardness
  dplyr::rename(
    sample   = sample,
    Hardness = hardness               # nome final da variável no projeto
  ) %>%
  dplyr::mutate(
    sample = std_sample(sample),
    Hardness = dplyr::case_when(
      is.numeric(Hardness) ~ Hardness,
      TRUE ~ readr::parse_number(
        as.character(Hardness),
        locale = readr::locale(decimal_mark = ",", grouping_mark = ".")
      )
    )
  ) %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    Hardness = mean(Hardness, na.rm = TRUE),
    .groups = "drop"
  )

# preview opcional
knitr::kable(texture_master, caption = "Average hardness (N) per sample")



## -----------------------------------------------------------------------------
#| label: cata-select-features
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(rstatix); library(knitr)
  library(FactoMineR)
})

# ---------------- Params ----------------
alpha_q      <- 0.05   # BH em Cochran’s Q
prev_min     <- 0.10   # 10% <= prevalência <= 90%
prev_max     <- 0.90
use_ca_cos2  <- TRUE
cos2_min     <- 0.15
# Política de seleção final:
# "min_count" (padrão), "cochran_only", "union_q_pen_cos2", "intersect_all", "q_and_cos2"
selection_policy <- "min_count"

# ------------ Checagens ------------
stopifnot(exists("dados_cata"))           # consumer, sample, ATR(0/1)...
stopifnot(exists("dados_clusterizados"))  # consumer, sample, Cluster, liking (num)

# ------------ Long + Prevalência ------------
cata_long <- dados_cata %>%
  pivot_longer(-c(consumer, sample), names_to = "attribute", values_to = "presence") %>%
  mutate(
    consumer  = as.character(consumer),
    sample    = as.character(sample),
    attribute = as.character(attribute),
    presence  = as.integer(presence)
  )

prev_tbl <- cata_long %>%
  group_by(attribute) %>%
  summarise(prevalence = mean(presence, na.rm = TRUE), .groups = "drop")

attrs_prev_ok <- prev_tbl %>%
  filter(prevalence >= prev_min, prevalence <= prev_max) %>%
  pull(attribute)

# ------------ Cochran’s Q (dentro de sujeito) ------------
cochran_res <- cata_long %>%
  filter(attribute %in% attrs_prev_ok) %>%
  group_by(attribute) %>%
  group_modify(~ {
    out <- tryCatch(cochran_qtest(data = .x, presence ~ sample | consumer),
                    error = function(e) NULL)
    if (is.null(out)) tibble(statistic = NA_real_, p = NA_real_) else
      tibble(statistic = out$statistic, p = out$p)
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p, method = "BH"))

attrs_q_ok <- cochran_res %>%
  filter(!is.na(p_adj), p_adj < alpha_q) %>%
  pull(attribute)

# ------------ Penalty (±95% CI não cruza 0 em ≥1 cluster) ------------
cata_like <- cata_long %>%
  left_join(
    dados_clusterizados %>% select(consumer, sample, Cluster, liking) %>%
      mutate(Cluster = as.factor(Cluster)),
    by = c("consumer","sample")
  )

penalty_one_attr <- function(df){
  n1 <- sum(df$presence == 1, na.rm = TRUE); n0 <- sum(df$presence == 0, na.rm = TRUE)
  m1 <- if (n1>0) mean(df$liking[df$presence==1], na.rm = TRUE) else NA_real_
  m0 <- if (n0>0) mean(df$liking[df$presence==0], na.rm = TRUE) else NA_real_
  s1 <- if (n1>1)  sd(df$liking[df$presence==1], na.rm = TRUE) else 0
  s0 <- if (n0>1)  sd(df$liking[df$presence==0], na.rm = TRUE) else 0
  if (is.na(m1) || is.na(m0)) return(tibble(diff = NA_real_, lwr = NA_real_, upr = NA_real_))
  se <- sqrt((s1^2)/pmax(n1,1) + (s0^2)/pmax(n0,1))
  diff <- m1 - m0
  tibble(diff = diff, lwr = diff - 1.96*se, upr = diff + 1.96*se)
}

penalty_sig_attrs <- cata_like %>%
  filter(attribute %in% attrs_prev_ok, !is.na(Cluster), is.finite(liking)) %>%
  group_by(attribute, Cluster) %>%
  group_modify(~ penalty_one_attr(.x)) %>%
  ungroup() %>%
  mutate(sig = !is.na(lwr) & !is.na(upr) & (lwr > 0 | upr < 0)) %>%
  group_by(attribute) %>%
  summarise(sig_any_cluster = any(sig, na.rm = TRUE), .groups = "drop") %>%
  filter(sig_any_cluster) %>%
  pull(attribute)

# ------------ CA (cos² em Dim 1–2) ------------
if (use_ca_cos2) {
  cata_table <- cata_long %>%
    filter(attribute %in% attrs_prev_ok) %>%
    group_by(sample, attribute) %>%
    summarise(freq = sum(presence, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = attribute, values_from = freq, values_fill = 0)
  rownames(cata_table) <- cata_table$sample
  cata_table$sample <- NULL
  ca_res <- FactoMineR::CA(cata_table, graph = FALSE)
  cos2_attr <- as.data.frame(ca_res$col$cos2)
  cos2_attr$attribute <- rownames(cos2_attr)
  cos2_attr$cos2_12 <- rowSums(cos2_attr[, c(1,2), drop = FALSE])
  attrs_cos_ok <- cos2_attr %>%
    filter(cos2_12 >= cos2_min) %>%
    pull(attribute)
} else {
  attrs_cos_ok <- attrs_prev_ok
}

# ------------ Estratégias de seleção ------------
strategy_sets <- list(
  cochran_only     = attrs_q_ok,
  union_q_pen_cos2 = intersect(union(attrs_q_ok, penalty_sig_attrs), attrs_cos_ok),
  intersect_all    = Reduce(intersect, list(attrs_q_ok, penalty_sig_attrs, attrs_cos_ok)),
  q_and_cos2       = intersect(attrs_q_ok, attrs_cos_ok)
)

strategy_counts <- vapply(strategy_sets, function(s) length(unique(s)), integer(1))

# Escolha final
choose_set <- function(policy, sets, counts){
  nonempty <- names(sets)[counts > 0]
  if (!length(nonempty)) return(list(name = NA_character_, attrs = character(0)))
  if (policy == "cochran_only" && counts["cochran_only"] > 0)
    return(list(name = "cochran_only", attrs = unique(sets$cochran_only)))
  if (policy == "union_q_pen_cos2" && counts["union_q_pen_cos2"] > 0)
    return(list(name = "union_q_pen_cos2", attrs = unique(sets$union_q_pen_cos2)))
  if (policy == "intersect_all" && counts["intersect_all"] > 0)
    return(list(name = "intersect_all", attrs = unique(sets$intersect_all)))
  if (policy == "q_and_cos2" && counts["q_and_cos2"] > 0)
    return(list(name = "q_and_cos2", attrs = unique(sets$q_and_cos2)))
  # default: min_count (quebra empate priorizando cochran_only > q_and_cos2 > intersect_all > union)
  if (policy == "min_count") {
    min_n <- min(counts[counts > 0])
    cands <- names(counts)[counts == min_n]
    pref  <- c("cochran_only","q_and_cos2","intersect_all","union_q_pen_cos2")
    picked <- intersect(pref, cands)[1]
    return(list(name = picked, attrs = unique(sets[[picked]])))
  }
  # fallback: menor não-nulo
  min_n <- min(counts[counts > 0])
  picked <- names(counts)[counts == min_n][1]
  list(name = picked, attrs = unique(sets[[picked]]))
}

pick <- choose_set(selection_policy, strategy_sets, strategy_counts)
cata_attrs_keep <- sort(unique(pick$attrs))

# ------------ Relatórios ------------
summary_steps <- tibble::tibble(
  Step  = c("All terms", "Prevalence 10–90%", "Cochran Q (BH<0.05)",
            "Penalty ±95% CI (≥1 cluster)", "CA cos² ≥ 0.15"),
  Count = c(n_distinct(cata_long$attribute),
            length(attrs_prev_ok),
            length(attrs_q_ok),
            length(penalty_sig_attrs),
            length(attrs_cos_ok))
)

summary_strat <- tibble::tibble(
  Strategy = names(strategy_counts),
  Count    = as.integer(strategy_counts)
) %>% arrange(Count)

knitr::kable(summary_steps, caption = "CATA feature selection — counts per step")
knitr::kable(summary_strat, caption = "CATA feature selection — strategy sizes")
knitr::kable(tibble::tibble(policy_used = selection_policy,
                            selected_strategy = pick$name,
                            n_selected = length(cata_attrs_keep)),
             caption = "Final selection policy and outcome")

if (length(cata_attrs_keep)) {
  knitr::kable(tibble::tibble(term = cata_attrs_keep),
               col.names = "Selected CATA terms", caption = NULL)
} else {
  warning("Nenhum termo selecionado — ajuste thresholds ou policy.")
}

# exporta
assign("cata_attrs_keep", cata_attrs_keep, envir = .GlobalEnv)
assign("cata_feature_report", summary_steps, envir = .GlobalEnv)
assign("cata_feature_strategies", summary_strat, envir = .GlobalEnv)




## -----------------------------------------------------------------------------
#| label: cata-props-and-master
#| message: false
#| warning: false

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(knitr) })

# --- checagens mínimas ---
stopifnot(
  exists("cata_attrs_keep"),
  exists("dados_cata"),
  exists("dados_clusterizados"),
  exists("liking_summ"),
  exists("vol_master_vars_top12") || exists("vol_master_vars"),
  exists("nonvol_master"),
  exists("texture_master")    # textura
)

# compatibilidade: se só existir vol_master_vars, cria _top12 a partir de vol_top12
if (!exists("vol_master_vars_top12") && exists("vol_master_vars") && exists("vol_top12")) {
  vol_master_vars_top12 <- vol_master_vars %>%
    dplyr::select(sample, dplyr::all_of(vol_top12)) %>%
    dplyr::rename_with(~ paste0("Vol.", .x), -sample)
}

# --- normaliza nomes de coluna 'cluster' em todos os objetos relevantes ---
norm_cluster <- function(df){
  if ("cluster" %in% names(df)) {
    df
  } else if ("Cluster" %in% names(df)) {
    dplyr::rename(df, cluster = Cluster)
  } else {
    stop("Nenhuma coluna 'cluster' ou 'Cluster' encontrada.")
  }
}

dados_clusterizados_n <- norm_cluster(dados_clusterizados)
liking_summ_n         <- norm_cluster(liking_summ)

# --- 1) Proporções CATA somente para os atributos selecionados ---
cata_long <- dados_cata %>%
  tidyr::pivot_longer(-c(consumer, sample),
                      names_to = "attribute",
                      values_to = "presence") %>%
  dplyr::filter(attribute %in% cata_attrs_keep) %>%
  dplyr::left_join(
    dados_clusterizados_n %>% dplyr::select(consumer, sample, cluster),
    by = c("consumer","sample")
  ) %>%
  dplyr::mutate(cluster = as.factor(cluster))

cata_props_filtered <- cata_long %>%
  dplyr::group_by(sample, cluster, attribute) %>%
  dplyr::summarise(prop = mean(presence, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(var = paste0("prop_", attribute)) %>%
  dplyr::select(sample, cluster, var, prop) %>%
  tidyr::pivot_wider(names_from = var, values_from = prop)

# --- 2) Voláteis top-12 por amostra (já agregados) ---
vol_master_vars_top12 <- vol_master_vars_top12 %>%
  dplyr::arrange(sample)

# --- 3) MASTER integrado (hedonic + CATA + voláteis + não-voláteis + textura) ---
master_integrado <- liking_summ_n %>%
  dplyr::left_join(cata_props_filtered,   by = c("sample","cluster")) %>%
  dplyr::left_join(vol_master_vars_top12, by = "sample") %>%
  dplyr::left_join(nonvol_master,         by = "sample") %>%
  dplyr::left_join(texture_master,        by = "sample") %>%  # Hardness
  dplyr::arrange(sample, cluster)

saveRDS(master_integrado, "master_integrado.rds")

# --- 4) Relatório final de blocos ---
nonvol_vars <- c("moisture","ash","lipids","protein","tdf","glucose","fructose",
                 "tpc","flav","anthocyanins","tannins","caffeine","trigonelline",
                 "5-cqa","3-cqa")

tbl_blocks <- tibble::tibble(
  block = c("Hedonic",
            "CATA props (filtered)",
            "Volatiles",
            "Non-volatiles",
            "Texture"),
  vars  = c(
    sum(names(master_integrado) == "mean_liking"),
    sum(grepl("^prop_", names(master_integrado))),
    sum(grepl("^Vol\\.", names(master_integrado))),
    sum(names(master_integrado) %in% nonvol_vars),
    sum(names(master_integrado) == "Hardness")
  )
)

knitr::kable(tbl_blocks, caption = "Blocks in master_integrado after CATA filtering")



## -----------------------------------------------------------------------------
#| label: mfa-biplot-clean
#| message: false
#| warning: false
#| fig.width: 10
#| fig.height: 7
#| dpi: 600
#| dev: png

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(FactoMineR); library(factoextra)
})

stopifnot(exists("master_integrado"))

# -------- helpers --------
impute_median <- function(df){
  out <- df
  for (j in seq_along(out)) if (is.numeric(out[[j]])) {
    med <- suppressWarnings(stats::median(out[[j]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    out[[j]][is.na(out[[j]])] <- med
  }
  out
}

abbr_nonvol <- c(
  moisture="MOIST", ash="ASH", lipids="LIP", protein="PTN", tdf="TDF",
  glucose="GLUC", fructose="FRUC", tpc="TPC", flav="FLAV", anthocyanins="ANT",
  tannins="TAN", caffeine="CAF", trigonelline="TRIG", cga="5-CQA", ncga="3-CQA"
)
abbr_all <- c(abbr_nonvol, Hardness = "HARD")

nonvol_vars <- names(abbr_nonvol)
texture_var <- "Hardness"

# ---- Paleta idêntica ao chunk dos voláteis ----
mfa_levels  <- c("NAT","PARCH","SC4","SC10","SC16","SC20")
mfa_palette <- scales::hue_pal()(length(mfa_levels))
names(mfa_palette) <- mfa_levels

clean_chem_name <- function(x){
  x <- gsub("^x", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("([0-9]),([0-9])", "\\1,\\2", x)
  x <- gsub("([0-9]) ([A-Za-z])", "\\1-\\2", x)
  x
}

# -------- 1) blocks --------
hed_col   <- intersect("mean_liking", names(master_integrado))
cata_cols <- grep("^prop_", names(master_integrado), value = TRUE)
vol_cols  <- grep("^Vol\\.", names(master_integrado), value = TRUE)
non_cols  <- intersect(nonvol_vars, names(master_integrado))
tex_cols  <- intersect(texture_var, names(master_integrado))

block_list <- list(
  Hedonic   = master_integrado %>% select(any_of(hed_col)),
  CATA      = master_integrado %>% select(any_of(cata_cols)),
  Volatiles = master_integrado %>% select(any_of(vol_cols)),
  NonVol    = master_integrado %>% select(any_of(non_cols)),
  Texture   = master_integrado %>% select(any_of(tex_cols))
)
block_list <- block_list[vapply(block_list, function(x) ncol(x) > 0, logical(1))]
stopifnot(length(block_list) >= 2)

X <- do.call(cbind, block_list) %>% as.data.frame()
rownames(X) <- with(master_integrado, paste(sample, cluster, sep = " :: "))
X_imp <- impute_median(X)

group_sizes <- vapply(block_list, ncol, integer(1))
mfa_res <- FactoMineR::MFA(
  X_imp,
  group      = as.numeric(group_sizes),
  type       = rep("s", length(group_sizes)),
  name.group = names(group_sizes),
  graph      = FALSE
)

# -------- 2) dados para o biplot --------
ind_df <- as.data.frame(mfa_res$ind$coord) %>%
  tibble::rownames_to_column("row_id") %>%
  tidyr::separate(row_id, into = c("sample","cluster"), sep = " :: ", remove = FALSE) %>%
  mutate(
    sample  = factor(toupper(sample), levels = mfa_levels),
    cluster = factor(cluster, levels = c("C1","C2"))
  ) %>%
  filter(!is.na(sample))

var_coord <- as.data.frame(mfa_res$quanti.var$coord) %>% tibble::rownames_to_column("var")
var_cos2  <- as.data.frame(mfa_res$quanti.var$cos2)  %>% tibble::rownames_to_column("var")

var_groups <- rep(names(block_list), times = group_sizes)
names(var_groups) <- colnames(X_imp)

var_df <- var_coord %>%
  left_join(var_cos2, by = "var", suffix = c("", "_cos2")) %>%
  mutate(
    group = var_groups[var],
    core  = sub("^[^.]+\\.", "", var),
    core  = sub("^prop_", "", core),
    core  = sub("^Vol\\.",  "", core),
    core  = case_when(
      core == "mean_liking" ~ "LIKING",
      group == "Volatiles"  ~ clean_chem_name(core),
      core %in% names(abbr_all) ~ abbr_all[core],
      TRUE ~ core
    ),
    label    = core,
    cos2_sum = coalesce(Dim.1_cos2, 0) + coalesce(Dim.2_cos2, 0)
  )

# limites por bloco
limits  <- c(Hedonic = 1L, CATA = 8L, Volatiles = 12L, NonVol = 15L, Texture = 1L)
lim_tbl <- tibble::tibble(group = names(limits), n_lim = as.integer(limits))

var_top <- var_df %>%
  left_join(lim_tbl, by = "group") %>%
  mutate(n_lim = coalesce(n_lim, 999L)) %>%
  group_by(group) %>%
  arrange(desc(cos2_sum), .by_group = TRUE) %>%
  mutate(.rank = row_number()) %>%
  filter(.rank <= n_lim) %>%
  ungroup() %>%
  select(-.rank, -n_lim)

# reescalonamento das variáveis
rng_ind1 <- range(ind_df$Dim.1, na.rm = TRUE)
rng_ind2 <- range(ind_df$Dim.2, na.rm = TRUE)
rng_var1 <- range(var_top$Dim.1, na.rm = TRUE)
rng_var2 <- range(var_top$Dim.2, na.rm = TRUE)

sx <- ifelse(diff(rng_var1) > 0, diff(rng_ind1)/diff(rng_var1), 1)
sy <- ifelse(diff(rng_var2) > 0, diff(rng_ind2)/diff(rng_var2), 1)

scale_fac <- 0.85 * min(sx, sy)

var_sc <- var_top %>%
  mutate(
    Dim.1 = Dim.1 * scale_fac,
    Dim.2 = Dim.2 * scale_fac
  )

# cores dos rótulos
var_col_map <- c(
  Hedonic   = "#223D6C",
  CATA      = "#223D6C",
  Volatiles = "black",
  NonVol    = "grey30",
  Texture   = "#8c564b"
)
var_sc$col <- var_col_map[var_sc$group]

var_sc_like  <- filter(var_sc, label == "LIKING")
var_sc_other <- filter(var_sc, label != "LIKING")

perc <- round(mfa_res$eig[1:2, 2], 1)

# -------- 3) plot + export high-res --------
p_mfa <- ggplot() +
  geom_segment(
    data = var_sc,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    colour = "grey75", linewidth = 0.3
  ) +
  geom_point(
    data = ind_df,
    aes(Dim.1, Dim.2, color = sample, shape = cluster),
    size = 3.2, alpha = 0.95
  ) +
  # ---- rótulos apenas das variáveis ----
  geom_text_repel(
    data = var_sc_other,
    aes(Dim.1, Dim.2, label = label),
    color = var_sc_other$col, size = 2.8, seed = 2,
    max.overlaps = 300, show.legend = FALSE
  ) +
  geom_text_repel(
    data = var_sc_like,
    aes(Dim.1, Dim.2, label = label),
    fontface = "bold", color = var_sc_like$col, size = 3.3, seed = 3,
    max.overlaps = 300, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey60") +

  # ---- paleta idêntica + ordem correta ----
  scale_color_manual(
    values = mfa_palette,
    breaks = mfa_levels,
    drop   = FALSE,
    name   = "Processing"
  ) +
  scale_shape_manual(
    values = c(C1 = 15, C2 = 16),
    breaks = c("C1","C2"),
    drop   = FALSE,
    name   = "Cluster"
  ) +

  labs(
    title = "MFA biplot",
    x     = paste0("Dim 1 (", perc[1], "%)"),
    y     = paste0("Dim 2 (", perc[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank()
  )

print(p_mfa)

# diretório de exportação e salvamento em alta resolução (PNG)
export_dir <- file.path(getwd(), "export_fig")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

ggplot2::ggsave(
  filename = file.path(export_dir, "mfa-biplot-clean.png"),
  plot     = p_mfa,
  width    = 10,
  height   = 7,
  units    = "in",
  dpi      = 600
)



## -----------------------------------------------------------------------------
#| label: mfa-eigen-tables
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(knitr); library(tibble)
})

# Garante que o MFA já foi ajustado
stopifnot(exists("mfa_res"))

# ---------------- Autovalores (eigenvalues) ----------------
eig_df_raw <- as.data.frame(mfa_res$eig)
colnames(eig_df_raw) <- c("eigenvalue", "percent", "cumulative")

eig_df <- eig_df_raw %>%
  rownames_to_column("Dim_label") %>%          # ignora o rótulo textual
  mutate(
    Dim        = dplyr::row_number(),          # Dim = 1, 2, 3, ...
    eigenvalue = as.numeric(eigenvalue),
    percent    = as.numeric(percent),
    cumulative = as.numeric(cumulative)
  ) %>%
  select(Dim, eigenvalue, percent, cumulative)

kable(
  eig_df,
  digits  = c(0, 3, 1, 1),
  caption = "MFA — autovalores e variância explicada por dimensão"
)

# ---------------- Autovetores (loadings) ----------------
load_df_raw <- as.data.frame(mfa_res$quanti.var$coord) %>%
  rownames_to_column("variable")

# número de dimensões a mostrar (ajuste se quiser mais)
max_dim <- min(5, ncol(load_df_raw) - 1)

load_df <- load_df_raw %>%
  select(variable, dplyr::all_of(colnames(load_df_raw)[2:(1 + max_dim)]))

kable(
  load_df,
  digits  = 3,
  caption = "MFA — autovetores (loadings) das variáveis quantitativas"
)



## -----------------------------------------------------------------------------
#| label: mfa-no-liking
#| message: false
#| warning: false
#| fig.width: 10
#| fig.height: 7

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(FactoMineR); library(factoextra)
})

stopifnot(exists("master_integrado"))

# ---------- helpers ----------
impute_median <- function(df){
  out <- df
  for (j in seq_along(out)) if (is.numeric(out[[j]])) {
    med <- suppressWarnings(stats::median(out[[j]], na.rm = TRUE))
    if (!is.finite(med)) med <- 0
    out[[j]][is.na(out[[j]])] <- med
  }
  out
}

abbr_nonvol <- c(
  moisture="MOIST", ash="ASH", lipids="LIP", protein="PTN", tdf="TDF",
  glucose="GLUC", fructose="FRUC", tpc="TPC", flav="FLAV", anthocyanins="ANT",
  tannins="TAN", caffeine="CAF", trigonelline="TRIG", cga="5-CQA", ncga="3-CQA"
)
abbr_all <- c(abbr_nonvol, Hardness = "HARD")

nonvol_vars <- names(abbr_nonvol)
texture_var <- "Hardness"

# mesma paleta do PCA (Processing)
sample_cols <- c(
  NAT   = "#ff9896",
  SC4   = "#2ca02c",
  SC10  = "#1f77b4",
  SC20  = "#e377c2",
  PARCH = "#bcbd22",
  SC8   = "#17becf",
  SC16  = "#9467bd"
)

clean_chem_name <- function(x){
  x <- gsub("^x", "", x)
  x <- gsub("_", " ", x)
  x <- gsub("([0-9]),([0-9])", "\\1,\\2", x)
  x <- gsub("([0-9]) ([A-Za-z])", "\\1-\\2", x)
  x
}

# ---------- 1) blocos (sem Hedonic) ----------
cata_cols <- grep("^prop_", names(master_integrado), value = TRUE)
vol_cols  <- grep("^Vol\\.", names(master_integrado), value = TRUE)
non_cols  <- intersect(nonvol_vars, names(master_integrado))
tex_cols  <- intersect(texture_var, names(master_integrado))

block_list <- list(
  CATA      = master_integrado %>% dplyr::select(any_of(cata_cols)),
  Volatiles = master_integrado %>% dplyr::select(any_of(vol_cols)),
  NonVol    = master_integrado %>% dplyr::select(any_of(non_cols)),
  Texture   = master_integrado %>% dplyr::select(any_of(tex_cols))
)
block_list <- block_list[vapply(block_list, function(x) ncol(x) > 0, logical(1))]
stopifnot(length(block_list) >= 2)

X <- do.call(cbind, block_list) %>% as.data.frame()
rownames(X) <- with(master_integrado, paste(sample, cluster, sep = " :: "))
X_imp <- impute_median(X)

group_sizes <- vapply(block_list, ncol, integer(1))
mfa_res_nolik <- FactoMineR::MFA(
  X_imp,
  group      = as.numeric(group_sizes),
  type       = rep("s", length(group_sizes)),
  name.group = names(group_sizes),
  graph      = FALSE
)

# ---------- 2) dados para o biplot ----------
ind_df <- as.data.frame(mfa_res_nolik$ind$coord) %>%
  tibble::rownames_to_column("row_id") %>%
  tidyr::separate(row_id, into = c("sample","cluster"), sep = " :: ", remove = FALSE) %>%
  dplyr::mutate(
    sample  = factor(toupper(sample), levels = names(sample_cols)),
    cluster = factor(cluster, levels = c("C1","C2"))
  ) %>%
  dplyr::filter(!is.na(sample))

var_coord <- as.data.frame(mfa_res_nolik$quanti.var$coord) %>%
  tibble::rownames_to_column("var")
var_cos2  <- as.data.frame(mfa_res_nolik$quanti.var$cos2) %>%
  tibble::rownames_to_column("var")

var_groups <- rep(names(block_list), times = group_sizes)
names(var_groups) <- colnames(X_imp)

var_df <- var_coord %>%
  dplyr::left_join(var_cos2, by = "var", suffix = c("", "_cos2")) %>%
  dplyr::mutate(
    group = var_groups[var],
    core  = sub("^[^.]+\\.", "", var),
    core  = sub("^prop_", "", core),
    core  = sub("^Vol\\.",  "", core),
    core  = dplyr::case_when(
      group == "Volatiles"        ~ clean_chem_name(core),
      core %in% names(abbr_all)   ~ abbr_all[core],
      TRUE                        ~ core
    ),
    label    = core,
    cos2_sum = dplyr::coalesce(Dim.1_cos2, 0) + dplyr::coalesce(Dim.2_cos2, 0)
  )

# limites por bloco (mesma ideia do MFA com liking)
limits  <- c(CATA = 8L, Volatiles = 12L, NonVol = 15L, Texture = 1L)
lim_tbl <- tibble::tibble(group = names(limits), n_lim = as.integer(limits))

var_top <- var_df %>%
  dplyr::left_join(lim_tbl, by = "group") %>%
  dplyr::mutate(n_lim = dplyr::coalesce(n_lim, 999L)) %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(dplyr::desc(cos2_sum), .by_group = TRUE) %>%
  dplyr::mutate(.rank = dplyr::row_number()) %>%
  dplyr::filter(.rank <= n_lim) %>%
  dplyr::ungroup() %>%
  dplyr::select(-.rank, -n_lim)

# reescala variáveis para caber no espaço dos indivíduos
rng_ind1 <- range(ind_df$Dim.1, na.rm = TRUE)
rng_ind2 <- range(ind_df$Dim.2, na.rm = TRUE)
rng_var1 <- range(var_top$Dim.1, na.rm = TRUE)
rng_var2 <- range(var_top$Dim.2, na.rm = TRUE)
sx <- ifelse(diff(rng_var1) > 0, diff(rng_ind1)/diff(rng_var1), 1)
sy <- ifelse(diff(rng_var2) > 0, diff(rng_ind2)/diff(rng_var2), 1)
scale_fac <- 0.85 * min(sx, sy)

var_sc <- var_top %>%
  dplyr::mutate(
    Dim.1 = Dim.1 * scale_fac,
    Dim.2 = Dim.2 * scale_fac
  )

# cores dos rótulos (CATA azul; textura marrom; demais neutros)
var_col_map <- c(
  CATA      = "#223D6C",
  Volatiles = "black",
  NonVol    = "grey30",
  Texture   = "#8c564b"
)
var_sc$col <- var_col_map[var_sc$group]

perc <- round(mfa_res_nolik$eig[1:2, 2], 1)

# ---------- 3) plot (mesmo estilo do MFA com liking, sem LIKING) ----------
ggplot() +
  # segmentos dos loadings (sem seta)
  geom_segment(
    data = var_sc,
    aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
    colour = "grey75", linewidth = 0.3
  ) +
  # indivíduos
  geom_point(
    data = ind_df,
    aes(Dim.1, Dim.2, color = sample, shape = cluster),
    size = 3.2, alpha = 0.95
  ) +
  geom_text_repel(
    data = ind_df,
    aes(Dim.1, Dim.2, label = sample, color = sample),
    size = 3.0, seed = 1, max.overlaps = 100, show.legend = FALSE
  ) +
  # variáveis
  geom_text_repel(
    data = var_sc,
    aes(Dim.1, Dim.2, label = label),
    color = var_sc$col, size = 2.8, seed = 2,
    max.overlaps = 300, show.legend = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey60") +
  scale_color_manual(values = sample_cols, drop = TRUE, name = "Processing") +
  scale_shape_manual(values = c(C1 = 15, C2 = 16), drop = TRUE, name = "Cluster") +
  labs(
    title    = "MFA biplot",
        x = paste0("Dim 1 (", perc[1], "%)"),
    y = paste0("Dim 2 (", perc[2], "%)")
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid    = element_blank()
  )



## -----------------------------------------------------------------------------
#| label: vol-pca-for-rv
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr); library(tibble)
})

# Precisamos de:
# - vol_top12 : vetor de nomes das 12 variáveis selecionadas pelo Fisher
# - vol_repl  : matriz de voláteis em nível de réplica, com coluna 'sample'
stopifnot(exists("vol_top12"), exists("vol_repl"))

# Garante que todas as variáveis de interesse existem em vol_repl
missing_vars <- setdiff(vol_top12, names(vol_repl))
if (length(missing_vars) > 0) {
  stop("As seguintes variáveis de 'vol_top12' não estão em 'vol_repl': ",
       paste(missing_vars, collapse = ", "))
}

# ---- 1) Agrega por amostra: média das 12 variáveis (top-12) ----
vol_top12_by_sample_raw <- vol_repl %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    dplyr::across(dplyr::all_of(vol_top12), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Ordena por amostra e prepara matriz numérica
vol_df <- vol_top12_by_sample_raw %>%
  dplyr::arrange(sample) %>%
  tibble::column_to_rownames("sample") %>%
  as.data.frame()

# Garante apenas colunas numéricas
vol_df <- vol_df[, vapply(vol_df, is.numeric, logical(1)), drop = FALSE]

# ---- 2) Imputação simples de NA por mediana e remoção de colunas com variância zero ----
for (j in seq_along(vol_df)) {
  v   <- vol_df[[j]]
  med <- suppressWarnings(stats::median(v, na.rm = TRUE))
  if (!is.finite(med)) med <- 0
  v[is.na(v)] <- med
  vol_df[[j]] <- v
}
sd_cols <- vapply(vol_df, stats::sd, na.rm = TRUE, FUN.VALUE = numeric(1))
vol_df  <- vol_df[, sd_cols > 0 & is.finite(sd_cols), drop = FALSE]

# ---- 3) PCA em escala (mesmo critério do biplot) ----
vol_pca <- prcomp(vol_df, center = TRUE, scale. = TRUE)

# Scores PC1 e PC2 por amostra (para o RV / acoplamento)
vol_master <- as.data.frame(vol_pca$x[, 1:2, drop = FALSE]) %>%
  tibble::rownames_to_column("sample") %>%
  dplyr::rename(PC1_vol = PC1, PC2_vol = PC2)

# Disponibiliza globalmente para outros chunks (RV, etc.)
assign("vol_pca", vol_pca, envir = .GlobalEnv)
assign("vol_master", vol_master, envir = .GlobalEnv)
assign("vol_top12_by_sample_raw", vol_top12_by_sample_raw, envir = .GlobalEnv)

# (opcional) pré-visualização:

print(vol_master)




## -----------------------------------------------------------------------------
#| label: pls-01-prepare-base
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

stopifnot(exists("master_integrado"))

# tratamentos válidos
pls_samples <- c("NAT","PARCH","SC4","SC10","SC16","SC20")

pls_base <- master_integrado %>%
  filter(sample %in% pls_samples) %>%
  arrange(sample, cluster)

# checagem
stopifnot(all(pls_samples %in% pls_base$sample))



## -----------------------------------------------------------------------------
#| label: pls-02a-Y-liking
#| message: false
#| warning: false

Y_liking <- pls_base %>%
  select(sample, cluster, mean_liking) %>%
  tidyr::pivot_wider(
    names_from  = cluster,
    values_from = mean_liking,
    names_prefix = ""
  ) %>%
  rename(
    C1_Liking = C1,
    C2_Liking = C2
  ) %>%
  mutate(
    Overall_Liking = rowMeans(select(., C1_Liking, C2_Liking), na.rm = TRUE)
  ) %>%
  select(sample, Overall_Liking, C1_Liking, C2_Liking)



## -----------------------------------------------------------------------------
#| label: diag-cata-props
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(knitr)
})

stopifnot(exists("cata_props"))

prop_cols_all <- grep("^prop_", names(cata_props), value = TRUE)

# contagens por “tipo” de termo (odr/flv e os 3 termos-chave)
odr_cols   <- grep("_odr$", prop_cols_all, value = TRUE, ignore.case = TRUE)
flv_cols   <- grep("_flv$", prop_cols_all, value = TRUE, ignore.case = TRUE)
sweet_cols <- grep("^prop_.*sweet$",  prop_cols_all, value = TRUE, ignore.case = TRUE)
bitter_cols<- grep("^prop_.*bitter$", prop_cols_all, value = TRUE, ignore.case = TRUE)
yummy_cols <- grep("^prop_.*yummy$",  prop_cols_all, value = TRUE, ignore.case = TRUE)

diag_tbl <- tibble::tibble(
  item = c("Linhas em cata_props (amostra × cluster)",
           "Total colunas prop_",
           "prop_ com sufixo _odr",
           "prop_ com sufixo _flv",
           "prop_ contendo sweet",
           "prop_ contendo bitter",
           "prop_ contendo yummy"),
  n = c(nrow(cata_props),
        length(prop_cols_all),
        length(odr_cols),
        length(flv_cols),
        length(sweet_cols),
        length(bitter_cols),
        length(yummy_cols))
)

knitr::kable(diag_tbl, caption = "Diagnóstico de cata_props (antes de filtragem)")

# Preview dos nomes (primeiros 40 para não poluir)
knitr::kable(
  tibble::tibble(prop_cols_preview = head(prop_cols_all, 40)),
  caption = "Primeiras colunas 'prop_' em cata_props"
)

# Também ajuda ver o que existe especificamente em odr/flv/sweet/bitter/yummy
knitr::kable(tibble::tibble(odr = head(odr_cols, 30)), caption = "Exemplos de prop_..._odr")
knitr::kable(tibble::tibble(flv = head(flv_cols, 30)), caption = "Exemplos de prop_..._flv")
knitr::kable(tibble::tibble(sweet = sweet_cols), caption = "Colunas com sweet (prop_...)")
knitr::kable(tibble::tibble(bitter = bitter_cols), caption = "Colunas com bitter (prop_...)")
knitr::kable(tibble::tibble(yummy = yummy_cols), caption = "Colunas com yummy (prop_...)")



## -----------------------------------------------------------------------------
#| label: export-cata-props
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(readr)
})

stopifnot(exists("cata_props"))

export_dir <- file.path(getwd(), "exports_check")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

csv_path <- file.path(export_dir, "cata_props_FULL.csv")
rds_path <- file.path(export_dir, "cata_props_FULL.rds")

readr::write_csv(cata_props, csv_path, na = "")
saveRDS(cata_props, rds_path)

cat("EXPORT OK\n")
cat("CSV:", csv_path, "\n")
cat("RDS:", rds_path, "\n")



## -----------------------------------------------------------------------------
#| label: plsC-01-build-Y-by-cluster
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tidyr)
})

stopifnot(exists("cata_props"), exists("liking_summ_n") || exists("liking_summ"))

pls_samples <- c("NAT","PARCH","SC4","SC10","SC16","SC20")

# liking por cluster (garante nomes)
lik_tbl <- (if (exists("liking_summ_n")) liking_summ_n else liking_summ) %>%
  select(sample, cluster, mean_liking) %>%
  mutate(sample = toupper(as.character(sample)),
         cluster = as.character(cluster)) %>%
  filter(sample %in% pls_samples)

# CATA por cluster (do próprio cata_props)
cata_sup_regex <- "(?i)(_odr$|_flv$|^sweet$|^bitter$|^yummy$)"

cata_tbl <- cata_props %>%
  mutate(sample = toupper(as.character(sample)),
         cluster = as.character(cluster)) %>%
  filter(sample %in% pls_samples) %>%
  select(sample, cluster, starts_with("prop_")) %>%
  rename_with(~ str_remove(.x, "^prop_"), starts_with("prop_")) %>%
  select(sample, cluster, matches(cata_sup_regex))

# junta liking + CATA (mesmo cluster)
Y_by_cluster <- lik_tbl %>%
  left_join(cata_tbl, by = c("sample","cluster")) %>%
  arrange(sample, cluster)

# separa
Y_C1 <- Y_by_cluster %>% filter(cluster == "C1") %>% select(-cluster)
Y_C2 <- Y_by_cluster %>% filter(cluster == "C2") %>% select(-cluster)

stopifnot(nrow(Y_C1) == 6, nrow(Y_C2) == 6)




## -----------------------------------------------------------------------------
#| label: export-cata-by-cluster
#| message: false
#| warning: false

suppressPackageStartupMessages({ library(readr) })

export_dir <- file.path(getwd(), "exports_check")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

readr::write_csv(Y_C1, file.path(export_dir, "Y_C1_liking_plus_CATA.csv"), na = "")
readr::write_csv(Y_C2, file.path(export_dir, "Y_C2_liking_plus_CATA.csv"), na = "")

cat("Exportados:\n",
    file.path(export_dir, "Y_C1_liking_plus_CATA.csv"), "\n",
    file.path(export_dir, "Y_C2_liking_plus_CATA.csv"), "\n")



## -----------------------------------------------------------------------------
#| label: plsC-02-build-X
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr)
})

stopifnot(exists("vol_top12"), exists("vol_master_vars"), exists("nonvol_master"))

nonvol_vars2 <- c(
  "glucose","fructose","tpc","flav","anthocyanins","tannins",
  "caffeine","trigonelline","cga","ncga"
)

# voláteis (sem prefixo Vol.)
X_vol <- vol_master_vars %>%
  mutate(sample = toupper(as.character(sample))) %>%
  filter(sample %in% pls_samples) %>%
  select(sample, all_of(vol_top12)) %>%
  arrange(sample)

# não-voláteis (fixos)
X_non <- nonvol_master %>%
  mutate(sample = toupper(as.character(sample))) %>%
  filter(sample %in% pls_samples) %>%
  select(sample, any_of(nonvol_vars2)) %>%
  arrange(sample)

# X final
X_df <- X_vol %>%
  left_join(X_non, by = "sample") %>%
  arrange(sample)

X_mat <- X_df %>% select(-sample) %>% as.matrix()
rownames(X_mat) <- X_df$sample

# imputação simples + scale
for (j in seq_len(ncol(X_mat))) {
  v <- X_mat[, j]
  med <- suppressWarnings(median(v, na.rm = TRUE))
  if (!is.finite(med)) med <- 0
  v[is.na(v)] <- med
  X_mat[, j] <- v
}
X_mat <- scale(X_mat)



## -----------------------------------------------------------------------------
#| label: chem-label-map
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# Converte nomes "R-friendly" -> rótulos legíveis (sem x inicial; sem underscores)
chem_label <- function(x){
  x <- as.character(x)
  x <- str_remove(x, "^Vol\\.")           # se ainda houver prefixo (não deveria no X final)
  x <- str_remove(x, "^x")                # remove x inicial (caso "x2_3_butanediol")
  x <- str_replace_all(x, "_", " ")       # underscores -> espaço
  x <- str_replace_all(x, " +", " ")
  x <- str_trim(x)

  # ajustes típicos para números e hifens (melhora leitura tipo "2 3" -> "2,3")
  x <- str_replace_all(x, "(\\d) (\\d)", "\\1,\\2")  # "2 3" -> "2,3" quando apropriado
  x <- str_replace_all(x, "(\\d) ([A-Za-z])", "\\1-\\2") # "2 octenoic" -> "2-octenoic"
  x
}

# Dicionário para os voláteis efetivamente usados na PLS
stopifnot(exists("vol_top12"))
vol_label_map <- tibble::tibble(
  var = vol_top12,
  label = chem_label(vol_top12)
)

# (Opcional) ver rapidamente o mapeamento
print(vol_label_map)



## -----------------------------------------------------------------------------
#| label: plsC-03-fit-C1-C2
#| message: false
#| warning: false

suppressPackageStartupMessages({ library(mixOmics) })

# --- Y matrizes ---
make_Y <- function(df){
  Y <- df %>% select(-sample) %>% as.matrix()
  rownames(Y) <- df$sample
  Y
}

Ymat_C1 <- make_Y(Y_C1)
Ymat_C2 <- make_Y(Y_C2)

# alinha linhas com X
Ymat_C1 <- Ymat_C1[rownames(X_mat), , drop = FALSE]
Ymat_C2 <- Ymat_C2[rownames(X_mat), , drop = FALSE]

stopifnot(identical(rownames(X_mat), rownames(Ymat_C1)))
stopifnot(identical(rownames(X_mat), rownames(Ymat_C2)))

set.seed(123)
pls_C1 <- pls(X = X_mat, Y = Ymat_C1, ncomp = 1, scale = FALSE)
pls_C2 <- pls(X = X_mat, Y = Ymat_C2, ncomp = 1, scale = FALSE)



## -----------------------------------------------------------------------------
#| label: plsC-04-metrics-article
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(knitr)
})

stopifnot(exists("pls_C1"), exists("pls_C2"))

perf_C1 <- mixOmics::perf(pls_C1, validation = "loo", progressBar = FALSE)
perf_C2 <- mixOmics::perf(pls_C2, validation = "loo", progressBar = FALSE)

get_val <- function(meas_list, feature_name, comp = 1){
  df <- meas_list$values
  df %>%
    filter(.data$feature == feature_name, .data$comp == comp) %>%
    slice(1) %>%
    pull(.data$value) %>%
    { if (!length(.)) NA_real_ else as.numeric(.[1]) }
}

# R2 de AJUSTE (fit) a partir de yhat do modelo (comp 1)
R2_fit <- function(pls_obj, y_name = "mean_liking", comp = 1){
  Y <- pls_obj$Y[, y_name, drop = TRUE]
  # mixOmics: fitted.values é um array n x q x ncomp
  Yhat <- pls_obj$fitted.values[, y_name, comp, drop = TRUE]
  sse <- sum((Y - Yhat)^2, na.rm = TRUE)
  sst <- sum((Y - mean(Y, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(sse) || !is.finite(sst) || sst == 0) return(NA_real_)
  1 - (sse / sst)
}

metrics_one <- function(pls_obj, perf_obj, model_name){
  tibble(
    Model     = model_name,
    Y         = "mean_liking",
    R2_fit    = R2_fit(pls_obj, "mean_liking", comp = 1),
    Q2_LOO    = get_val(perf_obj$measures$Q2, "mean_liking", comp = 1),
    RMSEP_LOO = get_val(perf_obj$measures$RMSEP, "mean_liking", comp = 1),
    Q2_total  = get_val(perf_obj$measures$Q2.total, "Y", comp = 1)
  )
}

tbl_metrics_article <- bind_rows(
  metrics_one(pls_C1, perf_C1, "PLS-C1"),
  metrics_one(pls_C2, perf_C2, "PLS-C2")
)

knitr::kable(tbl_metrics_article, digits = 3,
             caption = "PLS-R metrics for the manuscript (1 component): model-fit R² and LOO cross-validation (Q², RMSEP)")



## -----------------------------------------------------------------------------
#| label: plsC-05-permutation-manual-and-word
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(mixOmics)
  library(dplyr)
  library(tibble)
  library(officer)
  library(flextable)
})

stopifnot(exists("pls_C1"), exists("pls_C2"))

# ---------- helpers ----------
get_val <- function(meas_list, feature_name, comp = 1){
  df <- meas_list$values
  out <- df %>%
    filter(.data$feature == feature_name, .data$comp == comp) %>%
    slice(1) %>%
    pull(.data$value)
  if (!length(out)) NA_real_ else as.numeric(out[1])
}

R2_fit <- function(pls_obj, y_name = "mean_liking", comp = 1){
  Y <- pls_obj$Y[, y_name, drop = TRUE]
  Yhat <- pls_obj$fitted.values[, y_name, comp, drop = TRUE]
  sse <- sum((Y - Yhat)^2, na.rm = TRUE)
  sst <- sum((Y - mean(Y, na.rm = TRUE))^2, na.rm = TRUE)
  if (!is.finite(sse) || !is.finite(sst) || sst == 0) return(NA_real_)
  1 - (sse / sst)
}

# Q2 observado (LOO) de um modelo já ajustado
Q2_obs_from_model <- function(pls_obj){
  perf_obj <- mixOmics::perf(pls_obj, validation = "loo", progressBar = FALSE)
  get_val(perf_obj$measures$Q2, "mean_liking", comp = 1)
}

# Q2 (LOO) com Y permutado (refaz o modelo)
Q2_perm_one <- function(X, Y, ncomp = 1){
  pls_p <- mixOmics::pls(X = X, Y = Y, ncomp = ncomp, mode = "regression")
  perf_p <- mixOmics::perf(pls_p, validation = "loo", progressBar = FALSE)
  get_val(perf_p$measures$Q2, "mean_liking", comp = 1)
}

# Permutation test manual
perm_test_Q2 <- function(pls_obj, nperm = 500, seed = 123){
  set.seed(seed)

  X <- pls_obj$X
  Y <- pls_obj$Y

  # Observado
  Q2_obs <- Q2_obs_from_model(pls_obj)

  # Permutações (permuta as linhas de Y mantendo X fixo)
  Q2_perm <- numeric(nperm)
  for (b in seq_len(nperm)) {
    idx <- sample.int(nrow(Y), replace = FALSE)
    Yp  <- Y[idx, , drop = FALSE]
    Q2_perm[b] <- Q2_perm_one(X, Yp, ncomp = 1)
  }

  # p-valor unilateral: proporção de permutados >= observado
  # (com correção padrão para evitar p=0)
  p_emp <- (sum(Q2_perm >= Q2_obs, na.rm = TRUE) + 1) / (sum(is.finite(Q2_perm)) + 1)

  list(Q2_obs = Q2_obs, Q2_perm = Q2_perm, p_value = p_emp)
}

# ---------- roda para C1 e C2 ----------
res_C1 <- perm_test_Q2(pls_C1, nperm = 500, seed = 123)
res_C2 <- perm_test_Q2(pls_C2, nperm = 500, seed = 123)

# ---------- métricas “para artigo” ----------
perf_C1 <- mixOmics::perf(pls_C1, validation = "loo", progressBar = FALSE)
perf_C2 <- mixOmics::perf(pls_C2, validation = "loo", progressBar = FALSE)

tbl_article <- tibble(
  Model      = c("PLS-C1", "PLS-C2"),
  Y          = "mean_liking",
  R2_fit     = c(R2_fit(pls_C1, "mean_liking", 1), R2_fit(pls_C2, "mean_liking", 1)),
  Q2_LOO     = c(res_C1$Q2_obs, res_C2$Q2_obs),
  RMSEP_LOO  = c(get_val(perf_C1$measures$RMSEP, "mean_liking", 1),
                get_val(perf_C2$measures$RMSEP, "mean_liking", 1)),
  Q2_total   = c(get_val(perf_C1$measures$Q2.total, "Y", 1),
                get_val(perf_C2$measures$Q2.total, "Y", 1)),
  p_perm_Q2  = c(res_C1$p_value, res_C2$p_value),
  nperm      = 500
)

# ---------- exporta Word ----------
ft <- flextable(tbl_article) %>%
  colformat_double(j = c("R2_fit","Q2_LOO","RMSEP_LOO","Q2_total","p_perm_Q2"),
                   digits = 3) %>%
  set_header_labels(
    Model     = "Model",
    Y         = "Response (Y)",
    R2_fit    = "R² (fit)",
    Q2_LOO    = "Q² (LOO)",
    RMSEP_LOO = "RMSEP (LOO)",
    Q2_total  = "Q² (total)",
    p_perm_Q2 = "p (perm. Q²)",
    nperm     = "Permutations (n)"
  ) %>%
  autofit()

doc <- read_docx() %>%
  body_add_par(
    "Table X. PLS-R performance metrics and permutation test results (1 latent component).",
    style = "Normal"
  ) %>%
  body_add_flextable(ft)

out_file <- "Table_PLS_metrics_permutation_FINAL.docx"
print(doc, target = out_file)

out_file




## -----------------------------------------------------------------------------
#| label: plsC-06-coef-vip-C1-C2
#| echo: false
#| message: false
#| warning: false
#| fig-width: 14
#| fig-height: 8
#| dpi: 300
#| fig-cap: "PLS-R (1 component): regression coefficients (top) and VIP scores (bottom) for mean liking in clusters C1 and C2."

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tibble)
})

stopifnot(exists("pls_C1"), exists("pls_C2"))

if (!exists("nonvol_vars2")) {
  nonvol_vars2 <- c("glucose","fructose","tpc","flav","anthocyanins","tannins",
                    "caffeine","trigonelline","cga","ncga")
}

clean_name <- function(x){
  x %>%
    str_remove("^Vol\\.") %>%
    str_remove("^x") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}
predictor_type <- function(pred_raw){
  ifelse(pred_raw %in% nonvol_vars2, "Non-volatile", "Volatile")
}

extract_wpq <- function(pls_obj, comp = 1){
  W <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) W <- pls_obj$loading.weights$X
  if (is.null(W) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) W <- pls_obj$loadings$X
  if (is.null(W)) stop("Não encontrei pesos X (W) no objeto PLS.")

  P <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) P <- pls_obj$loadings$X
  if (is.null(P)) P <- W

  Q <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Q <- pls_obj$loadings$Y
  if (is.null(Q)) stop("Não encontrei loadings de Y (Q) no objeto PLS.")

  W <- as.matrix(W); P <- as.matrix(P); Q <- as.matrix(Q)

  if (ncol(W) < comp || ncol(P) < comp || ncol(Q) < comp)
    stop("O objeto PLS não tem o número de componentes solicitado: comp = ", comp)

  w <- W[, comp, drop = TRUE]
  p <- P[, comp, drop = TRUE]
  q <- Q[, comp, drop = TRUE]

  if (is.null(names(w)) && !is.null(rownames(W))) names(w) <- rownames(W)
  if (is.null(names(p)) && !is.null(rownames(P))) names(p) <- rownames(P)

  list(w = w, p = p, q = q)
}

coef_kcomp_singleY <- function(pls_obj, comp_max = 1){
  # b = sum_k w_k * q_k / (p_k' w_k)   (para 1 Y; aqui usamos mean_liking)
  # retorna um vetor com coeficientes por preditor
  Wmat <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) Wmat <- as.matrix(pls_obj$loading.weights$X)
  if (is.null(Wmat) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) Wmat <- as.matrix(pls_obj$loadings$X)
  if (is.null(Wmat)) stop("Não encontrei pesos X (W) no objeto PLS.")

  Pmat <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) Pmat <- as.matrix(pls_obj$loadings$X)
  if (is.null(Pmat)) Pmat <- Wmat

  Qmat <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Qmat <- as.matrix(pls_obj$loadings$Y)
  if (is.null(Qmat)) stop("Não encontrei loadings de Y (Q) no objeto PLS.")

  K <- min(comp_max, ncol(Wmat), ncol(Pmat), ncol(Qmat))
  if (K < 1) stop("Sem componentes disponíveis para calcular coeficientes.")

  b <- rep(0, nrow(Wmat))
  names(b) <- rownames(Wmat)

  for (k in seq_len(K)) {
    w <- Wmat[, k, drop = TRUE]
    p <- Pmat[, k, drop = TRUE]
    q <- Qmat[, k, drop = TRUE]
    if (length(q) > 1) q <- q[1]
    denom <- sum(p * w, na.rm = TRUE)
    if (!is.finite(denom) || denom == 0) next
    b <- b + (w * as.numeric(q) / denom)
  }
  b
}

vip_kcomp <- function(pls_obj, comp_max = 1){
  # VIP_j = sqrt(p * sum_k SSY_k * w_jk^2 / sum_j w_jk^2  / sum_k SSY_k)
  # Para 1 Y: SSY_k proporcional a ||t_k||^2 * q_k^2
  Wmat <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) Wmat <- as.matrix(pls_obj$loading.weights$X)
  if (is.null(Wmat) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) Wmat <- as.matrix(pls_obj$loadings$X)
  if (is.null(Wmat)) stop("Não encontrei pesos X (W) no objeto PLS.")

  Tmat <- NULL
  if (!is.null(pls_obj$variates) && !is.null(pls_obj$variates$X)) Tmat <- as.matrix(pls_obj$variates$X)
  if (is.null(Tmat)) stop("Não encontrei scores X (t) no objeto PLS (variates$X).")

  Qmat <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Qmat <- as.matrix(pls_obj$loadings$Y)
  if (is.null(Qmat)) stop("Não encontrei loadings de Y (Q) no objeto PLS.")

  K <- min(comp_max, ncol(Wmat), ncol(Tmat), ncol(Qmat))
  p <- nrow(Wmat)
  if (K < 1) stop("Sem componentes disponíveis para calcular VIP.")

  SSY <- numeric(K)
  for (k in seq_len(K)) {
    t <- Tmat[, k, drop = TRUE]
    q <- Qmat[, k, drop = TRUE]
    if (length(q) > 1) q <- q[1]
    SSY[k] <- sum(t^2, na.rm = TRUE) * (as.numeric(q)^2)
  }
  if (sum(SSY, na.rm = TRUE) == 0) SSY[] <- 1

  W2 <- Wmat[, seq_len(K), drop = FALSE]^2
  denom_k <- colSums(W2, na.rm = TRUE)
  denom_k[denom_k == 0 | !is.finite(denom_k)] <- 1

  vip2 <- (p * rowSums(sweep(W2, 2, denom_k, "/") * rep(SSY, each = p), na.rm = TRUE)) / sum(SSY, na.rm = TRUE)
  vip <- sqrt(vip2)
  names(vip) <- rownames(Wmat)
  vip
}

# ======= ESCOLHA DO NÚMERO DE COMPONENTES PARA COEF/VIP =======
ncomp_show <- 1  # <- mude para 2 se decidir mostrar comp 1–2

get_cluster_tables <- function(pls_obj, cluster_name, ncomp_use = 1){
  b <- coef_kcomp_singleY(pls_obj, comp_max = ncomp_use)
  v <- vip_kcomp(pls_obj, comp_max = ncomp_use)

  coef_df <- tibble(predictor_raw = names(b), coef = as.numeric(b)) %>%
    mutate(predictor = clean_name(predictor_raw),
           type = predictor_type(predictor_raw),
           Cluster = cluster_name)

  vip_df <- tibble(predictor_raw = names(v), VIP = as.numeric(v)) %>%
    mutate(predictor = clean_name(predictor_raw),
           type = predictor_type(predictor_raw),
           Cluster = cluster_name)

  list(coef = coef_df, vip = vip_df)
}

t1 <- get_cluster_tables(pls_C1, "Cluster 1", ncomp_use = ncomp_show)
t2 <- get_cluster_tables(pls_C2, "Cluster 2", ncomp_use = ncomp_show)

coef_df <- bind_rows(t1$coef, t2$coef)
vip_df  <- bind_rows(t1$vip,  t2$vip)

g_coef <- coef_df %>%
  group_by(Cluster) %>%
  arrange(coef, .by_group = TRUE) %>%
  ungroup() %>%
  ggplot(aes(x = coef, y = reorder(predictor, coef), fill = type)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Cluster, scales = "free_y") +
  scale_fill_manual(values = c("Non-volatile" = "#F8766D", "Volatile" = "#00BFC4")) +
  labs(
    x = paste0("Regression coefficient (", ncomp_show, " comp; scaled inputs)"),
    y = NULL,
    title = "PLS-R coefficients for mean liking"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "right")

g_vip <- vip_df %>%
  group_by(Cluster) %>%
  arrange(VIP, .by_group = TRUE) %>%
  ungroup() %>%
  ggplot(aes(x = VIP, y = reorder(predictor, VIP), fill = type)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_wrap(~Cluster, scales = "free_y") +
  scale_fill_manual(values = c("Non-volatile" = "#F8766D", "Volatile" = "#00BFC4")) +
  labs(
    x = paste0("VIP score (", ncomp_show, " comp)"),
    y = NULL,
    title = "VIP scores for mean liking"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "right")

p_comp <- g_coef / g_vip

print(p_comp)

# -------- export PNG (working directory) --------
export_dir <- file.path(getwd(), "export_fig")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(export_dir, paste0("Figure_PLS_coef_VIP_C1_C2_comp", ncomp_show, ".png"))

ggplot2::ggsave(
  filename = out_png,
  plot     = p_comp,
  width    = 14,
  height   = 8,
  units    = "in",
  dpi      = 300
)

message("Figura exportada: ", out_png)



## -----------------------------------------------------------------------------
#| label: plsC-06-coef-vip-C1-C2b
#| echo: false
#| message: false
#| warning: false
#| fig-width: 14
#| fig-height: 8
#| dpi: 300
#| fig-cap: "PLS-R (1 component): regression coefficients (top) and VIP scores (bottom) for mean liking in clusters C1 and C2."

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(stringr)
  library(tibble)
})

stopifnot(exists("pls_C1"), exists("pls_C2"))

if (!exists("nonvol_vars2")) {
  nonvol_vars2 <- c("glucose","fructose","tpc","flav","anthocyanins","tannins",
                    "caffeine","trigonelline","cga","ncga")
}

clean_name <- function(x){
  x %>%
    str_remove("^Vol\\.") %>%
    str_remove("^x") %>%
    str_replace_all("_", " ") %>%
    str_squish()
}

# ======= MAPA DE RÓTULOS (APENAS PARA A FIGURA) =======
# Depois do clean_name(): "cga" -> "5-CQA"; "ncga" -> "3-CQA"
label_map <- c(
  "cga"  = "5-CQA",
  "ncga" = "3-CQA"
)

predictor_type <- function(pred_raw){
  ifelse(pred_raw %in% nonvol_vars2, "Non-volatile", "Volatile")
}

extract_wpq <- function(pls_obj, comp = 1){
  W <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) W <- pls_obj$loading.weights$X
  if (is.null(W) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) W <- pls_obj$loadings$X
  if (is.null(W)) stop("Não encontrei pesos X (W) no objeto PLS.")

  P <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) P <- pls_obj$loadings$X
  if (is.null(P)) P <- W

  Q <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Q <- pls_obj$loadings$Y
  if (is.null(Q)) stop("Não encontrei loadings de Y (Q) no objeto PLS.")

  W <- as.matrix(W); P <- as.matrix(P); Q <- as.matrix(Q)

  if (ncol(W) < comp || ncol(P) < comp || ncol(Q) < comp)
    stop("O objeto PLS não tem o número de componentes solicitado: comp = ", comp)

  w <- W[, comp, drop = TRUE]
  p <- P[, comp, drop = TRUE]
  q <- Q[, comp, drop = TRUE]

  if (is.null(names(w)) && !is.null(rownames(W))) names(w) <- rownames(W)
  if (is.null(names(p)) && !is.null(rownames(P))) names(p) <- rownames(P)

  list(w = w, p = p, q = q)
}

coef_kcomp_singleY <- function(pls_obj, comp_max = 1){
  # b = sum_k w_k * q_k / (p_k' w_k)   (para 1 Y; aqui usamos mean_liking)
  # retorna um vetor com coeficientes por preditor
  Wmat <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) Wmat <- as.matrix(pls_obj$loading.weights$X)
  if (is.null(Wmat) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) Wmat <- as.matrix(pls_obj$loadings$X)
  if (is.null(Wmat)) stop("Não encontrei pesos X (W) no objeto PLS.")

  Pmat <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) Pmat <- as.matrix(pls_obj$loadings$X)
  if (is.null(Pmat)) Pmat <- Wmat

  Qmat <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Qmat <- as.matrix(pls_obj$loadings$Y)
  if (is.null(Qmat)) stop("Não encontrei loadings de Y (Q) no objeto PLS.")

  K <- min(comp_max, ncol(Wmat), ncol(Pmat), ncol(Qmat))
  if (K < 1) stop("Sem componentes disponíveis para calcular coeficientes.")

  b <- rep(0, nrow(Wmat))
  names(b) <- rownames(Wmat)

  for (k in seq_len(K)) {
    w <- Wmat[, k, drop = TRUE]
    p <- Pmat[, k, drop = TRUE]
    q <- Qmat[, k, drop = TRUE]
    if (length(q) > 1) q <- q[1]
    denom <- sum(p * w, na.rm = TRUE)
    if (!is.finite(denom) || denom == 0) next
    b <- b + (w * as.numeric(q) / denom)
  }
  b
}

vip_kcomp <- function(pls_obj, comp_max = 1){
  # VIP_j = sqrt(p * sum_k SSY_k * w_jk^2 / sum_j w_jk^2  / sum_k SSY_k)
  # Para 1 Y: SSY_k proporcional a ||t_k||^2 * q_k^2
  Wmat <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) Wmat <- as.matrix(pls_obj$loading.weights$X)
  if (is.null(Wmat) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) Wmat <- as.matrix(pls_obj$loadings$X)
  if (is.null(Wmat)) stop("Não encontrei pesos X (W) no objeto PLS.")

  Tmat <- NULL
  if (!is.null(pls_obj$variates) && !is.null(pls_obj$variates$X)) Tmat <- as.matrix(pls_obj$variates$X)
  if (is.null(Tmat)) stop("Não encontrei scores X (t) no objeto PLS (variates$X).")

  Qmat <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Qmat <- as.matrix(pls_obj$loadings$Y)
  if (is.null(Qmat)) stop("Não encontrei loadings de Y (Q) no objeto PLS.")

  K <- min(comp_max, ncol(Wmat), ncol(Tmat), ncol(Qmat))
  p <- nrow(Wmat)
  if (K < 1) stop("Sem componentes disponíveis para calcular VIP.")

  SSY <- numeric(K)
  for (k in seq_len(K)) {
    t <- Tmat[, k, drop = TRUE]
    q <- Qmat[, k, drop = TRUE]
    if (length(q) > 1) q <- q[1]
    SSY[k] <- sum(t^2, na.rm = TRUE) * (as.numeric(q)^2)
  }
  if (sum(SSY, na.rm = TRUE) == 0) SSY[] <- 1

  W2 <- Wmat[, seq_len(K), drop = FALSE]^2
  denom_k <- colSums(W2, na.rm = TRUE)
  denom_k[denom_k == 0 | !is.finite(denom_k)] <- 1

  vip2 <- (p * rowSums(sweep(W2, 2, denom_k, "/") * rep(SSY, each = p), na.rm = TRUE)) / sum(SSY, na.rm = TRUE)
  vip <- sqrt(vip2)
  names(vip) <- rownames(Wmat)
  vip
}

# ======= ESCOLHA DO NÚMERO DE COMPONENTES PARA COEF/VIP =======
ncomp_show <- 1  # <- mude para 2 se decidir mostrar comp 1–2

get_cluster_tables <- function(pls_obj, cluster_name, ncomp_use = 1){
  b <- coef_kcomp_singleY(pls_obj, comp_max = ncomp_use)
  v <- vip_kcomp(pls_obj, comp_max = ncomp_use)

  coef_df <- tibble(predictor_raw = names(b), coef = as.numeric(b)) %>%
    mutate(
      predictor = dplyr::recode(clean_name(predictor_raw), !!!label_map,
                                .default = clean_name(predictor_raw)),
      type = predictor_type(predictor_raw),
      Cluster = cluster_name
    )

  vip_df <- tibble(predictor_raw = names(v), VIP = as.numeric(v)) %>%
    mutate(
      predictor = dplyr::recode(clean_name(predictor_raw), !!!label_map,
                                .default = clean_name(predictor_raw)),
      type = predictor_type(predictor_raw),
      Cluster = cluster_name
    )

  list(coef = coef_df, vip = vip_df)
}

t1 <- get_cluster_tables(pls_C1, "Cluster 1", ncomp_use = ncomp_show)
t2 <- get_cluster_tables(pls_C2, "Cluster 2", ncomp_use = ncomp_show)

coef_df <- bind_rows(t1$coef, t2$coef)
vip_df  <- bind_rows(t1$vip,  t2$vip)

g_coef <- coef_df %>%
  group_by(Cluster) %>%
  arrange(coef, .by_group = TRUE) %>%
  ungroup() %>%
  ggplot(aes(x = coef, y = reorder(predictor, coef), fill = type)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~Cluster, scales = "free_y") +
  scale_fill_manual(values = c("Non-volatile" = "#F8766D", "Volatile" = "#00BFC4")) +
  labs(
    x = paste0("Regression coefficient (", ncomp_show, " comp; scaled inputs)"),
    y = NULL,
    title = "PLS-R coefficients for mean liking"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "right")

g_vip <- vip_df %>%
  group_by(Cluster) %>%
  arrange(VIP, .by_group = TRUE) %>%
  ungroup() %>%
  ggplot(aes(x = VIP, y = reorder(predictor, VIP), fill = type)) +
  geom_col(width = 0.75) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  facet_wrap(~Cluster, scales = "free_y") +
  scale_fill_manual(values = c("Non-volatile" = "#F8766D", "Volatile" = "#00BFC4")) +
  labs(
    x = paste0("VIP score (", ncomp_show, " comp)"),
    y = NULL,
    title = "VIP scores for mean liking"
  ) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "right")

p_comp <- g_coef / g_vip

print(p_comp)

# -------- export PNG (working directory) --------
export_dir <- file.path(getwd(), "export_fig")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

out_png <- file.path(export_dir, paste0("Figure_PLS_coef_VIP_C1_C2_compb", ncomp_show, ".png"))

ggplot2::ggsave(
  filename = out_png,
  plot     = p_comp,
  width    = 14,
  height   = 8,
  units    = "in",
  dpi      = 300
)

message("Figura exportada: ", out_png)



## -----------------------------------------------------------------------------
#| label: plsC-07a-refit-2comp-for-biplot
#| echo: false
#| message: false
#| warning: false

suppressPackageStartupMessages(library(mixOmics))

stopifnot(exists("pls_C1"), exists("pls_C2"))

ensure_2comp <- function(pls_obj){
  X <- pls_obj$X
  Y <- pls_obj$Y
  if (is.null(X) || is.null(Y)) {
    stop("O objeto PLS não contém X/Y internamente (pls_obj$X ou pls_obj$Y é NULL). ",
         "Nesse caso, precisamos refazer com os objetos X e Y usados no ajuste.")
  }
  if (!is.null(pls_obj$variates$X) && ncol(as.matrix(pls_obj$variates$X)) >= 2) return(pls_obj)

  mixOmics::pls(X = X, Y = Y, ncomp = 2, mode = "regression", scale = TRUE)
}

pls_C1_bi <- ensure_2comp(pls_C1)
pls_C2_bi <- ensure_2comp(pls_C2)

assign("pls_C1_bi", pls_C1_bi, envir = .GlobalEnv)
assign("pls_C2_bi", pls_C2_bi, envir = .GlobalEnv)

message("OK: pls_C1_bi e pls_C2_bi prontos para biplot (2 componentes).")



## -----------------------------------------------------------------------------
#| label: plsC-07b-biplots-C1-C2
#| echo: false
#| message: false
#| warning: false
#| fig-width: 12
#| fig-height: 7
#| dpi: 300

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tibble)
  library(ggplot2); library(ggrepel); library(patchwork)
})

stopifnot(exists("pls_C1_bi"), exists("pls_C2_bi"))

treat_levels  <- c("NAT","PARCH","SC4","SC10","SC16","SC20")
max_x_labels  <- 12
max_cata_lab  <- 12  # nº máx de CATA a rotular (liking sempre entra)
max_y_labels  <- max_cata_lab + 1

if (!exists("nonvol_vars2")) {
  nonvol_vars2 <- c("glucose","fructose","tpc","flav","anthocyanins","tannins",
                    "caffeine","trigonelline","cga","ncga")
}

# ========= NOMEAÇÃO (vetorizada) =========
# - SWT_odr -> SWEET_odr (mantém CATA odor separado do SWEET taste)
# - SWEET (taste) fica "SWEET" (sem virar "Sweet")
# - CATA mantém underscore (_odr/_flv)
# - químicos/voláteis viram "nome com espaço"
clean_name <- function(x){

  x0  <- as.character(x)

  # normaliza separadores para detectar casos
  x1  <- x0 %>%
    str_replace_all("[-\\.\\s]+", "_") %>%
    str_replace_all("_+", "_")

  low <- tolower(x1)

  out <- x0

  # --- casos críticos (vetorizado) ---
  out[low %in% c("swt_odr","swt_odor","sweet_odr")] <- "SWEET_odr"
  out[low %in% c("sweet","sweet_taste")]            <- "SWEET"

  # --- limpeza padrão ---
  out <- out %>%
    str_remove("^Vol\\.") %>%
    str_remove("^prop_") %>%
    str_remove("^x") %>%
    str_squish()

  # se for CATA (termina em _odr ou _flv, ou for BITTER/YUMMY/SWEET), mantém underscore
  is_cata <- str_detect(tolower(out), "(_odr$|_flv$)") | toupper(out) %in% c("BITTER","YUMMY","SWEET","SWEET_ODR")

  out[ is_cata] <- out[is_cata] %>%
    str_replace_all("\\s+", "_") %>%  # garante underscore
    str_replace_all("_+", "_")

  out[!is_cata] <- out[!is_cata] %>%
    str_replace_all("_", " ") %>%
    str_squish()

  out
}

x_type <- function(pred_raw){
  ifelse(pred_raw %in% nonvol_vars2, "Non-volatile", "Volatile")
}

# ========= getters robustos =========
get_T <- function(pls_obj){
  T <- NULL
  if (!is.null(pls_obj$variates) && !is.null(pls_obj$variates$X)) T <- pls_obj$variates$X
  if (is.null(T)) stop("Não encontrei scores X em pls_obj$variates$X.")
  T <- as.matrix(T)
  if (ncol(T) < 2) stop("Objeto PLS com <2 componentes. Use pls_*_bi (2 comps) no biplot.")
  T
}

get_W <- function(pls_obj){
  W <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) W <- pls_obj$loading.weights$X
  if (is.null(W) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) W <- pls_obj$loadings$X
  if (is.null(W)) stop("Não encontrei pesos/loadings de X em loading.weights$X nem em loadings$X.")
  W <- as.matrix(W)

  if (is.null(rownames(W)) || anyNA(rownames(W))) {
    if (!is.null(pls_obj$X) && !is.null(colnames(pls_obj$X)) && nrow(W) == ncol(pls_obj$X)) {
      rownames(W) <- colnames(pls_obj$X)
    } else {
      rownames(W) <- paste0("X", seq_len(nrow(W)))
    }
  }
  W
}

get_Q <- function(pls_obj){
  Q <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Q <- pls_obj$loadings$Y
  if (is.null(Q) && !is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$Y)) Q <- pls_obj$loading.weights$Y
  if (is.null(Q)) stop("Não encontrei loadings/pesos de Y em loadings$Y nem em loading.weights$Y.")
  Q <- as.matrix(Q)
  if (is.null(rownames(Q)) || anyNA(rownames(Q))) rownames(Q) <- paste0("Y", seq_len(nrow(Q)))
  Q
}

# ========= VIP (seguro) =========
vip_kcomp_safe <- function(pls_obj, comp_max = 2){
  tryCatch({
    W <- get_W(pls_obj)
    T <- get_T(pls_obj)
    Q <- get_Q(pls_obj)

    K <- min(comp_max, ncol(W), ncol(T), ncol(Q))
    p <- nrow(W)

    SSY <- numeric(K)
    for (k in seq_len(K)) {
      t <- T[, k, drop = TRUE]
      q <- Q[, k, drop = TRUE]
      SSY[k] <- sum(t^2, na.rm = TRUE) * sum(q^2, na.rm = TRUE)
    }
    if (!is.finite(sum(SSY)) || sum(SSY) == 0) SSY[] <- 1

    W2 <- W[, seq_len(K), drop = FALSE]^2
    denom_k <- colSums(W2, na.rm = TRUE)
    denom_k[denom_k == 0 | !is.finite(denom_k)] <- 1

    # VIP clássico (por variável X)
    wnorm <- sweep(W2, 2, denom_k, "/")
    vip2  <- (p * rowSums(wnorm %*% matrix(SSY, ncol = 1), na.rm = TRUE)) / sum(SSY, na.rm = TRUE)
    vip   <- sqrt(as.numeric(vip2))
    names(vip) <- rownames(W)
    vip
  }, error = function(e) NULL)
}

# ========= extrator =========
extract_biplot_df <- function(pls_obj, cluster_label){

  # scores
  T <- get_T(pls_obj)
  ind_df <- as.data.frame(T[, 1:2, drop = FALSE]) %>%
    setNames(c("Comp1","Comp2")) %>%
    tibble::rownames_to_column("row_id") %>%
    mutate(
      sample = str_extract(.data$row_id, paste(treat_levels, collapse = "|")),
      sample = factor(toupper(as.character(.data$sample)), levels = treat_levels),
      Cluster = cluster_label
    ) %>%
    filter(!is.na(.data$sample))

  # X
  W <- get_W(pls_obj)
  x_df <- as.data.frame(W[, 1:2, drop = FALSE]) %>%
    setNames(c("Comp1","Comp2")) %>%
    tibble::rownames_to_column("var_raw") %>%
    mutate(
      type  = x_type(.data$var_raw),
      label = clean_name(.data$var_raw),
      mag2d = sqrt(.data$Comp1^2 + .data$Comp2^2),
      Cluster = cluster_label
    )

  # Y (CATA + liking)
  Q <- get_Q(pls_obj)
  y_df <- as.data.frame(Q[, 1:2, drop = FALSE]) %>%
    setNames(c("Comp1","Comp2")) %>%
    tibble::rownames_to_column("var_raw") %>%
    mutate(
      type  = ifelse(toupper(.data$var_raw) %in% c("MEAN_LIKING","LIKING","OVERALL_LIKING"), "Liking", "CATA"),
      label = clean_name(.data$var_raw),
      mag2d = sqrt(.data$Comp1^2 + .data$Comp2^2),
      Cluster = cluster_label
    )

  # ---- seleção X por VIP (fallback: mag2d) ----
  vip_vec <- vip_kcomp_safe(pls_obj, comp_max = 2)

  if (!is.null(vip_vec)) {
    x_keep <- x_df %>%
      mutate(VIP = unname(vip_vec[.data$var_raw])) %>%
      filter(is.finite(.data$VIP), .data$VIP >= 1) %>%
      arrange(desc(.data$VIP))

    if (nrow(x_keep) == 0) x_keep <- x_df %>% slice_max(.data$mag2d, n = max_x_labels)
    if (nrow(x_keep) > max_x_labels) x_keep <- x_keep %>% slice_head(n = max_x_labels)
  } else {
    x_keep <- x_df %>% slice_max(.data$mag2d, n = max_x_labels)
  }

  # ---- seleção Y: liking sempre; CATA filtra (top mag2d) ----
  y_like <- y_df %>% filter(.data$type == "Liking")
  y_cata <- y_df %>% filter(.data$type == "CATA") %>% slice_max(.data$mag2d, n = max_cata_lab)
  y_keep <- bind_rows(y_like, y_cata)

  list(ind = ind_df, x = x_keep, y = y_keep)
}

# ========= reescala variáveis no espaço dos scores =========
rescale_vars <- function(ind_df, var_df, scale_fac = 0.85){
  rng_i1 <- range(ind_df$Comp1, na.rm = TRUE); rng_i2 <- range(ind_df$Comp2, na.rm = TRUE)
  rng_v1 <- range(var_df$Comp1, na.rm = TRUE); rng_v2 <- range(var_df$Comp2, na.rm = TRUE)
  sx <- ifelse(diff(rng_v1) > 0, diff(rng_i1)/diff(rng_v1), 1)
  sy <- ifelse(diff(rng_v2) > 0, diff(rng_i2)/diff(rng_v2), 1)
  sc <- scale_fac * min(sx, sy)
  var_df %>% mutate(Comp1 = .data$Comp1 * sc, Comp2 = .data$Comp2 * sc)
}

# ========= biplot (SEM legenda; eixos principais mais espessos) =========
make_biplot <- function(ind_df, x_df, y_df, title_txt, shift_tbl = NULL, shift_samples = NULL){

  x_sc <- rescale_vars(ind_df, x_df, scale_fac = 0.95)
  y_sc <- rescale_vars(ind_df, y_df, scale_fac = 0.95)

  col_map <- c(
    "Volatile"     = "#00BFC4",
    "Non-volatile" = "#F8766D",
    "CATA"         = "#7CAE00",
    "Liking"       = "#C77CFF"
  )

  # --- scores: rótulos com nudges manuais quando necessário ---
  ind2 <- ind_df %>%
    mutate(label = as.character(.data$sample),
           label_key = str_to_lower(str_squish(.data$label)),
           x_lab = .data$Comp1,
           y_lab = .data$Comp2)

  if (!is.null(shift_samples) && nrow(shift_samples) > 0) {
    ss2 <- shift_samples %>%
      transmute(label_key = str_to_lower(str_squish(.data$label)),
                nudge_x = .data$nudge_x, nudge_y = .data$nudge_y)

    ind2 <- ind2 %>%
      left_join(ss2, by = "label_key") %>%
      mutate(
        nudge_x = dplyr::coalesce(.data$nudge_x, 0),
        nudge_y = dplyr::coalesce(.data$nudge_y, 0),
        x_lab = .data$Comp1 + .data$nudge_x,
        y_lab = .data$Comp2 + .data$nudge_y
      )
  }

  # --- offsets para variáveis (X + Y) ---
  var_all <- bind_rows(
    x_sc %>% mutate(group = "X"),
    y_sc %>% mutate(group = "Y")
  ) %>%
    mutate(label_key = str_to_lower(str_squish(.data$label)),
           x_lab = .data$Comp1, y_lab = .data$Comp2)

  if (!is.null(shift_tbl) && nrow(shift_tbl) > 0) {
    st2 <- shift_tbl %>%
      transmute(label_key = str_to_lower(str_squish(.data$label)),
                nudge_x = .data$nudge_x, nudge_y = .data$nudge_y)

    var_all <- var_all %>%
      left_join(st2, by = "label_key") %>%
      mutate(
        nudge_x = dplyr::coalesce(.data$nudge_x, 0),
        nudge_y = dplyr::coalesce(.data$nudge_y, 0),
        x_lab = .data$Comp1 + .data$nudge_x,
        y_lab = .data$Comp2 + .data$nudge_y
      )
  }

  ggplot() +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.55),
      panel.grid.minor = element_blank(),
      axis.text  = element_text(color = "grey55"),
      axis.title = element_text(color = "grey55"),
      plot.title = element_text(color = "grey25", face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 1.0) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 1.0) +

    # scores
    geom_point(data = ind_df, aes(.data$Comp1, .data$Comp2), color = "grey45", size = 3.2) +
    geom_text_repel(
      data = ind2,
      aes(.data$x_lab, .data$y_lab, label = .data$label),
      color = "grey45", size = 4.2, seed = 1, max.overlaps = 200
    ) +

    # variáveis (X e Y já filtradas) — pontos + rótulos
    geom_point(data = var_all, aes(.data$Comp1, .data$Comp2, color = .data$type),
               size = 3.2, alpha = 0.95) +
    geom_text_repel(
      data = var_all,
      aes(.data$x_lab, .data$y_lab, label = .data$label, color = .data$type),
      size = 3.9,
      fontface = ifelse(var_all$type == "Liking", "bold", "plain"),
      seed = 2, max.overlaps = 300
    ) +

    scale_color_manual(values = col_map) +
    labs(title = title_txt, x = "Component 1", y = "Component 2")
}

# ===== Extrai dados (com VIP p/ X e filtro p/ CATA) =====
dC1 <- extract_biplot_df(pls_C1_bi, "Cluster 1")
dC2 <- extract_biplot_df(pls_C2_bi, "Cluster 2")

# ===== NUDGES (baseados no seu diagnóstico de colisões) =====
# Cluster 1 (6 colisões: 5 var-var + caffeine vs SC16)
shift_C1 <- tibble::tribble(
  ~label,           ~nudge_x, ~nudge_y,
  "tpc",              -0.25,     0.10,
  "cga",              -0.25,    -0.10,
  "SWEET",             0.20,    -0.20,
  "flav",             -0.25,     0.12,
  "benzaldehyde",      0.20,    -0.25,
  "SWEET_odr",         0.30,    -0.00
  )

# Cluster 2 (4 colisões var-var)
shift_C2 <- tibble::tribble(
  ~label,           ~nudge_x, ~nudge_y,
  "flav",              0.25,    -0.10,
  "tpc",              -0.22,     0.10,
  "trigonelline",     -0.25,     0.15
)

# desloca SC16 para não “sumir” sob caffeine
shift_samples_C1 <- tibble::tribble(
  ~label, ~nudge_x, ~nudge_y,
  "SC16",   0.35,     0.25
)

shift_samples_C2 <- tibble::tribble(
  ~label, ~nudge_x, ~nudge_y,
  "SC16",   0.35,    -0.25
)

# ===== Plota =====
pC1 <- make_biplot(
  dC1$ind, dC1$x, dC1$y,
  "PLS-R biplot — Cluster 1",
  shift_tbl = shift_C1,
  shift_samples = shift_samples_C1
)

pC2 <- make_biplot(
  dC2$ind, dC2$x, dC2$y,
  "PLS-R biplot — Cluster 2",
  shift_tbl = shift_C2,
  shift_samples = shift_samples_C2
)

p_all <- pC1 | pC2
print(p_all)

# ===== Export (robusto) =====
# Em Quarto/knitr, esta é a pasta do .qmd
doc_dir <- getwd()  # (reanálise) current_input(dir=TRUE) devolvia o caminho do arquivo
if (is.null(doc_dir) || is.na(doc_dir)) doc_dir <- getwd()

export_dir <- file.path(doc_dir, "export_fig")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

f1 <- file.path(export_dir, "Figure_PLSR_biplot_C1_comp1_comp2.png")
f2 <- file.path(export_dir, "Figure_PLSR_biplot_C2_comp1_comp2.png")
f3 <- file.path(export_dir, "Figure_PLSR_biplot_C1_C2_comp1_comp2.png")

ggsave(f1, plot = pC1, width = 10.5, height = 7.5, units = "in", dpi = 300)
ggsave(f2, plot = pC2, width = 10.5, height = 7.5, units = "in", dpi = 300)
ggsave(f3, plot = p_all, width = 18,  height = 7.5, units = "in", dpi = 300)

message("Biplots exportados em (ABS): ", normalizePath(export_dir, winslash = "/", mustWork = FALSE))
message("Arquivos criados? ",
        "\n - C1: ", file.exists(f1),
        "\n - C2: ", file.exists(f2),
        "\n - C1C2: ", file.exists(f3))



## -----------------------------------------------------------------------------
#| label: plsC-07b-biplots-C1-C2b
#| echo: false
#| message: false
#| warning: false
#| fig-width: 12
#| fig-height: 7
#| dpi: 300

suppressPackageStartupMessages({
  library(dplyr); library(stringr); library(tibble)
  library(ggplot2); library(ggrepel); library(patchwork)
})

stopifnot(exists("pls_C1_bi"), exists("pls_C2_bi"))

treat_levels  <- c("NAT","PARCH","SC4","SC10","SC16","SC20")
max_x_labels  <- 12
max_cata_lab  <- 12  # nº máx de CATA a rotular (liking sempre entra)
max_y_labels  <- max_cata_lab + 1

if (!exists("nonvol_vars2")) {
  nonvol_vars2 <- c("glucose","fructose","tpc","flav","anthocyanins","tannins",
                    "caffeine","trigonelline","cga","ncga")
}

# ======= MAPA DE RÓTULOS (APENAS PARA A FIGURA) =======
# Depois do clean_name(): "cga" -> "5-CQA"; "ncga" -> "3-CQA"
label_map <- c(
  "cga"  = "5-CQA",
  "ncga" = "3-CQA"
)

# ========= NOMEAÇÃO (vetorizada) =========
# - SWT_odr -> SWEET_odr (mantém CATA odor separado do SWEET taste)
# - SWEET (taste) fica "SWEET" (sem virar "Sweet")
# - CATA mantém underscore (_odr/_flv)
# - químicos/voláteis viram "nome com espaço"
clean_name <- function(x){

  x0  <- as.character(x)

  # normaliza separadores para detectar casos
  x1  <- x0 %>%
    str_replace_all("[-\\.\\s]+", "_") %>%
    str_replace_all("_+", "_")

  low <- tolower(x1)

  out <- x0

  # --- casos críticos (vetorizado) ---
  out[low %in% c("swt_odr","swt_odor","sweet_odr")] <- "SWEET_odr"
  out[low %in% c("sweet","sweet_taste")]            <- "SWEET"

  # --- limpeza padrão ---
  out <- out %>%
    str_remove("^Vol\\.") %>%
    str_remove("^prop_") %>%
    str_remove("^x") %>%
    str_squish()

  # se for CATA (termina em _odr ou _flv, ou for BITTER/YUMMY/SWEET), mantém underscore
  is_cata <- str_detect(tolower(out), "(_odr$|_flv$)") | toupper(out) %in% c("BITTER","YUMMY","SWEET","SWEET_ODR")

  out[ is_cata] <- out[is_cata] %>%
    str_replace_all("\\s+", "_") %>%  # garante underscore
    str_replace_all("_+", "_")

  out[!is_cata] <- out[!is_cata] %>%
    str_replace_all("_", " ") %>%
    str_squish()

  # ===== aplica recode somente no rótulo final (apenas para figura) =====
  out <- dplyr::recode(out, !!!label_map, .default = out)

  out
}

x_type <- function(pred_raw){
  ifelse(pred_raw %in% nonvol_vars2, "Non-volatile", "Volatile")
}

# ========= getters robustos =========
get_T <- function(pls_obj){
  T <- NULL
  if (!is.null(pls_obj$variates) && !is.null(pls_obj$variates$X)) T <- pls_obj$variates$X
  if (is.null(T)) stop("Não encontrei scores X em pls_obj$variates$X.")
  T <- as.matrix(T)
  if (ncol(T) < 2) stop("Objeto PLS com <2 componentes. Use pls_*_bi (2 comps) no biplot.")
  T
}

get_W <- function(pls_obj){
  W <- NULL
  if (!is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$X)) W <- pls_obj$loading.weights$X
  if (is.null(W) && !is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$X)) W <- pls_obj$loadings$X
  if (is.null(W)) stop("Não encontrei pesos/loadings de X em loading.weights$X nem em loadings$X.")
  W <- as.matrix(W)

  if (is.null(rownames(W)) || anyNA(rownames(W))) {
    if (!is.null(pls_obj$X) && !is.null(colnames(pls_obj$X)) && nrow(W) == ncol(pls_obj$X)) {
      rownames(W) <- colnames(pls_obj$X)
    } else {
      rownames(W) <- paste0("X", seq_len(nrow(W)))
    }
  }
  W
}

get_Q <- function(pls_obj){
  Q <- NULL
  if (!is.null(pls_obj$loadings) && !is.null(pls_obj$loadings$Y)) Q <- pls_obj$loadings$Y
  if (is.null(Q) && !is.null(pls_obj$loading.weights) && !is.null(pls_obj$loading.weights$Y)) Q <- pls_obj$loading.weights$Y
  if (is.null(Q)) stop("Não encontrei loadings/pesos de Y em loadings$Y nem em loading.weights$Y.")
  Q <- as.matrix(Q)
  if (is.null(rownames(Q)) || anyNA(rownames(Q))) rownames(Q) <- paste0("Y", seq_len(nrow(Q)))
  Q
}

# ========= VIP (seguro) =========
vip_kcomp_safe <- function(pls_obj, comp_max = 2){
  tryCatch({
    W <- get_W(pls_obj)
    T <- get_T(pls_obj)
    Q <- get_Q(pls_obj)

    K <- min(comp_max, ncol(W), ncol(T), ncol(Q))
    p <- nrow(W)

    SSY <- numeric(K)
    for (k in seq_len(K)) {
      t <- T[, k, drop = TRUE]
      q <- Q[, k, drop = TRUE]
      SSY[k] <- sum(t^2, na.rm = TRUE) * sum(q^2, na.rm = TRUE)
    }
    if (!is.finite(sum(SSY)) || sum(SSY) == 0) SSY[] <- 1

    W2 <- W[, seq_len(K), drop = FALSE]^2
    denom_k <- colSums(W2, na.rm = TRUE)
    denom_k[denom_k == 0 | !is.finite(denom_k)] <- 1

    # VIP clássico (por variável X)
    wnorm <- sweep(W2, 2, denom_k, "/")
    vip2  <- (p * rowSums(wnorm %*% matrix(SSY, ncol = 1), na.rm = TRUE)) / sum(SSY, na.rm = TRUE)
    vip   <- sqrt(as.numeric(vip2))
    names(vip) <- rownames(W)
    vip
  }, error = function(e) NULL)
}

# ========= extrator =========
extract_biplot_df <- function(pls_obj, cluster_label){

  # scores
  T <- get_T(pls_obj)
  ind_df <- as.data.frame(T[, 1:2, drop = FALSE]) %>%
    setNames(c("Comp1","Comp2")) %>%
    tibble::rownames_to_column("row_id") %>%
    mutate(
      sample = str_extract(.data$row_id, paste(treat_levels, collapse = "|")),
      sample = factor(toupper(as.character(.data$sample)), levels = treat_levels),
      Cluster = cluster_label
    ) %>%
    filter(!is.na(.data$sample))

  # X
  W <- get_W(pls_obj)
  x_df <- as.data.frame(W[, 1:2, drop = FALSE]) %>%
    setNames(c("Comp1","Comp2")) %>%
    tibble::rownames_to_column("var_raw") %>%
    mutate(
      type  = x_type(.data$var_raw),
      label = clean_name(.data$var_raw),
      mag2d = sqrt(.data$Comp1^2 + .data$Comp2^2),
      Cluster = cluster_label
    )

  # Y (CATA + liking)
  Q <- get_Q(pls_obj)
  y_df <- as.data.frame(Q[, 1:2, drop = FALSE]) %>%
    setNames(c("Comp1","Comp2")) %>%
    tibble::rownames_to_column("var_raw") %>%
    mutate(
      type  = ifelse(toupper(.data$var_raw) %in% c("MEAN_LIKING","LIKING","OVERALL_LIKING"), "Liking", "CATA"),
      label = clean_name(.data$var_raw),
      mag2d = sqrt(.data$Comp1^2 + .data$Comp2^2),
      Cluster = cluster_label
    )

  # ---- seleção X por VIP (fallback: mag2d) ----
  vip_vec <- vip_kcomp_safe(pls_obj, comp_max = 2)

  if (!is.null(vip_vec)) {
    x_keep <- x_df %>%
      mutate(VIP = unname(vip_vec[.data$var_raw])) %>%
      filter(is.finite(.data$VIP), .data$VIP >= 1) %>%
      arrange(desc(.data$VIP))

    if (nrow(x_keep) == 0) x_keep <- x_df %>% slice_max(.data$mag2d, n = max_x_labels)
    if (nrow(x_keep) > max_x_labels) x_keep <- x_keep %>% slice_head(n = max_x_labels)
  } else {
    x_keep <- x_df %>% slice_max(.data$mag2d, n = max_x_labels)
  }

  # ---- seleção Y: liking sempre; CATA filtra (top mag2d) ----
  y_like <- y_df %>% filter(.data$type == "Liking")
  y_cata <- y_df %>% filter(.data$type == "CATA") %>% slice_max(.data$mag2d, n = max_cata_lab)
  y_keep <- bind_rows(y_like, y_cata)

  list(ind = ind_df, x = x_keep, y = y_keep)
}

# ========= reescala variáveis no espaço dos scores =========
rescale_vars <- function(ind_df, var_df, scale_fac = 0.85){
  rng_i1 <- range(ind_df$Comp1, na.rm = TRUE); rng_i2 <- range(ind_df$Comp2, na.rm = TRUE)
  rng_v1 <- range(var_df$Comp1, na.rm = TRUE); rng_v2 <- range(var_df$Comp2, na.rm = TRUE)
  sx <- ifelse(diff(rng_v1) > 0, diff(rng_i1)/diff(rng_v1), 1)
  sy <- ifelse(diff(rng_v2) > 0, diff(rng_i2)/diff(rng_v2), 1)
  sc <- scale_fac * min(sx, sy)
  var_df %>% mutate(Comp1 = .data$Comp1 * sc, Comp2 = .data$Comp2 * sc)
}

# ========= biplot (SEM legenda; eixos principais mais espessos) =========
make_biplot <- function(ind_df, x_df, y_df, title_txt, shift_tbl = NULL, shift_samples = NULL){

  x_sc <- rescale_vars(ind_df, x_df, scale_fac = 0.95)
  y_sc <- rescale_vars(ind_df, y_df, scale_fac = 0.95)

  col_map <- c(
    "Volatile"     = "#00BFC4",
    "Non-volatile" = "#F8766D",
    "CATA"         = "#7CAE00",
    "Liking"       = "#C77CFF"
  )

  # --- scores: rótulos com nudges manuais quando necessário ---
  ind2 <- ind_df %>%
    mutate(label = as.character(.data$sample),
           label_key = str_to_lower(str_squish(.data$label)),
           x_lab = .data$Comp1,
           y_lab = .data$Comp2)

  if (!is.null(shift_samples) && nrow(shift_samples) > 0) {
    ss2 <- shift_samples %>%
      transmute(label_key = str_to_lower(str_squish(.data$label)),
                nudge_x = .data$nudge_x, nudge_y = .data$nudge_y)

    ind2 <- ind2 %>%
      left_join(ss2, by = "label_key") %>%
      mutate(
        nudge_x = dplyr::coalesce(.data$nudge_x, 0),
        nudge_y = dplyr::coalesce(.data$nudge_y, 0),
        x_lab = .data$Comp1 + .data$nudge_x,
        y_lab = .data$Comp2 + .data$nudge_y
      )
  }

  # --- offsets para variáveis (X + Y) ---
  var_all <- bind_rows(
    x_sc %>% mutate(group = "X"),
    y_sc %>% mutate(group = "Y")
  ) %>%
    mutate(label_key = str_to_lower(str_squish(.data$label)),
           x_lab = .data$Comp1, y_lab = .data$Comp2)

  if (!is.null(shift_tbl) && nrow(shift_tbl) > 0) {
    st2 <- shift_tbl %>%
      transmute(label_key = str_to_lower(str_squish(.data$label)),
                nudge_x = .data$nudge_x, nudge_y = .data$nudge_y)

    var_all <- var_all %>%
      left_join(st2, by = "label_key") %>%
      mutate(
        nudge_x = dplyr::coalesce(.data$nudge_x, 0),
        nudge_y = dplyr::coalesce(.data$nudge_y, 0),
        x_lab = .data$Comp1 + .data$nudge_x,
        y_lab = .data$Comp2 + .data$nudge_y
      )
  }

  ggplot() +
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.55),
      panel.grid.minor = element_blank(),
      axis.text  = element_text(color = "grey55"),
      axis.title = element_text(color = "grey55"),
      plot.title = element_text(color = "grey25", face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    geom_hline(yintercept = 0, color = "grey65", linewidth = 1.0) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 1.0) +

    # scores
    geom_point(data = ind_df, aes(.data$Comp1, .data$Comp2), color = "grey45", size = 3.2) +
    geom_text_repel(
      data = ind2,
      aes(.data$x_lab, .data$y_lab, label = .data$label),
      color = "grey45", size = 4.2, seed = 1, max.overlaps = 200
    ) +

    # variáveis (X e Y já filtradas) — pontos + rótulos
    geom_point(data = var_all, aes(.data$Comp1, .data$Comp2, color = .data$type),
               size = 3.2, alpha = 0.95) +
    geom_text_repel(
      data = var_all,
      aes(.data$x_lab, .data$y_lab, label = .data$label, color = .data$type),
      size = 3.9,
      fontface = ifelse(var_all$type == "Liking", "bold", "plain"),
      seed = 2, max.overlaps = 300
    ) +

    scale_color_manual(values = col_map) +
    labs(title = title_txt, x = "Component 1", y = "Component 2")
}

# ===== Extrai dados (com VIP p/ X e filtro p/ CATA) =====
dC1 <- extract_biplot_df(pls_C1_bi, "Cluster 1")
dC2 <- extract_biplot_df(pls_C2_bi, "Cluster 2")

# ===== NUDGES (baseados no seu diagnóstico de colisões) =====
# Cluster 1 (6 colisões: 5 var-var + caffeine vs SC16)
shift_C1 <- tibble::tribble(
  ~label,           ~nudge_x, ~nudge_y,
  "tpc",              -0.25,     0.10,
  "cga",              -0.25,    -0.10,
  "SWEET",             0.20,    -0.20,
  "flav",             -0.25,     0.12,
  "benzaldehyde",      0.20,    -0.25,
  "SWEET_odr",         0.30,    -0.00
)

# Cluster 2 (4 colisões var-var)
shift_C2 <- tibble::tribble(
  ~label,           ~nudge_x, ~nudge_y,
  "flav",              0.25,    -0.10,
  "tpc",              -0.22,     0.10,
  "trigonelline",     -0.25,     0.15
)

# desloca SC16 para não “sumir” sob caffeine
shift_samples_C1 <- tibble::tribble(
  ~label, ~nudge_x, ~nudge_y,
  "SC16",   0.35,     0.25
)

shift_samples_C2 <- tibble::tribble(
  ~label, ~nudge_x, ~nudge_y,
  "SC16",   0.35,    -0.25
)

# ===== Plota =====
pC1 <- make_biplot(
  dC1$ind, dC1$x, dC1$y,
  "PLS-R biplot — Cluster 1",
  shift_tbl = shift_C1,
  shift_samples = shift_samples_C1
)

pC2 <- make_biplot(
  dC2$ind, dC2$x, dC2$y,
  "PLS-R biplot — Cluster 2",
  shift_tbl = shift_C2,
  shift_samples = shift_samples_C2
)

p_all <- pC1 | pC2
print(p_all)

# ===== Export (robusto) =====
# Em Quarto/knitr, esta é a pasta do .qmd
doc_dir <- getwd()  # (reanálise) current_input(dir=TRUE) devolvia o caminho do arquivo
if (is.null(doc_dir) || is.na(doc_dir)) doc_dir <- getwd()

export_dir <- file.path(doc_dir, "export_fig")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

f1 <- file.path(export_dir, "Figure_PLSR_biplot_C1_comp1_comp2b.png")
f2 <- file.path(export_dir, "Figure_PLSR_biplot_C2_comp1_comp2b.png")
f3 <- file.path(export_dir, "Figure_PLSR_biplot_C1_C2_comp1_comp2b.png")

ggsave(f1, plot = pC1, width = 10.5, height = 7.5, units = "in", dpi = 300)
ggsave(f2, plot = pC2, width = 10.5, height = 7.5, units = "in", dpi = 300)
ggsave(f3, plot = p_all, width = 18,  height = 7.5, units = "in", dpi = 300)

message("Biplots exportados em (ABS): ", normalizePath(export_dir, winslash = "/", mustWork = FALSE))
message("Arquivos criados? ",
        "\n - C1: ", file.exists(f1),
        "\n - C2: ", file.exists(f2),
        "\n - C1C2: ", file.exists(f3))



## -----------------------------------------------------------------------------
#| label: plsC-07c-diagnose-label-collisions-v2
#| echo: false
#| message: false
#| warning: false

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

stopifnot(exists("dC1"), exists("dC2"))

find_collisions_v2 <- function(ind_df, x_df, y_df, thresh_vv = 0.22, thresh_vs = 0.28){

  rescale_vars <- function(ind_df, var_df, scale_fac = 0.85){
    rng_i1 <- range(ind_df$Comp1, na.rm = TRUE); rng_i2 <- range(ind_df$Comp2, na.rm = TRUE)
    rng_v1 <- range(var_df$Comp1, na.rm = TRUE); rng_v2 <- range(var_df$Comp2, na.rm = TRUE)
    sx <- ifelse(diff(rng_v1) > 0, diff(rng_i1)/diff(rng_v1), 1)
    sy <- ifelse(diff(rng_v2) > 0, diff(rng_i2)/diff(rng_v2), 1)
    sc <- scale_fac * min(sx, sy)
    var_df %>% mutate(Comp1 = .data$Comp1 * sc, Comp2 = .data$Comp2 * sc)
  }

  # variáveis (X+Y) já escaladas no mesmo espaço do score
  x_sc <- rescale_vars(ind_df, x_df) %>%
    transmute(label = as.character(label), type = as.character(type),
              Comp1 = as.numeric(Comp1), Comp2 = as.numeric(Comp2),
              group = "X")

  y_sc <- rescale_vars(ind_df, y_df) %>%
    transmute(label = as.character(label), type = as.character(type),
              Comp1 = as.numeric(Comp1), Comp2 = as.numeric(Comp2),
              group = "Y")

  v <- bind_rows(x_sc, y_sc) %>%
    filter(is.finite(Comp1), is.finite(Comp2)) %>%
    mutate(label_key = str_to_lower(str_squish(label)))

  s <- ind_df %>%
    transmute(label = as.character(sample), type = "Score",
              Comp1 = as.numeric(Comp1), Comp2 = as.numeric(Comp2),
              group = "S") %>%
    filter(is.finite(Comp1), is.finite(Comp2)) %>%
    mutate(label_key = str_to_lower(str_squish(label)))

  # (A) variável x variável
  collisions_vv <- tibble()
  if (nrow(v) >= 2){
    idx <- combn(seq_len(nrow(v)), 2)
    a <- v[idx[1,], , drop = FALSE]
    b <- v[idx[2,], , drop = FALSE]
    dist <- sqrt((a$Comp1 - b$Comp1)^2 + (a$Comp2 - b$Comp2)^2)

    collisions_vv <- tibble(
      kind  = "var-var",
      label_1 = a$label, type_1 = a$type, group_1 = a$group,
      label_2 = b$label, type_2 = b$type, group_2 = b$group,
      dist = dist
    ) %>%
      filter(.data$dist <= thresh_vv) %>%
      arrange(.data$dist)
  }

  # (B) variável x score
  collisions_vs <- tibble()
  if (nrow(v) >= 1 && nrow(s) >= 1){
    # produto cartesiano v x s
    a <- v[rep(seq_len(nrow(v)), each = nrow(s)), , drop = FALSE]
    b <- s[rep(seq_len(nrow(s)), times = nrow(v)), , drop = FALSE]
    dist <- sqrt((a$Comp1 - b$Comp1)^2 + (a$Comp2 - b$Comp2)^2)

    collisions_vs <- tibble(
      kind  = "var-score",
      label_1 = a$label, type_1 = a$type, group_1 = a$group,
      label_2 = b$label, type_2 = b$type, group_2 = b$group,
      dist = dist
    ) %>%
      filter(.data$dist <= thresh_vs) %>%
      arrange(.data$dist)
  }

  bind_rows(collisions_vv, collisions_vs)
}

collisions_C1 <- find_collisions_v2(dC1$ind, dC1$x, dC1$y, thresh_vv = 0.22, thresh_vs = 0.33)
collisions_C2 <- find_collisions_v2(dC2$ind, dC2$x, dC2$y, thresh_vv = 0.22, thresh_vs = 0.33)

message("=== COLISÕES (Cluster 1) — var-var + var-score ===")
print(collisions_C1, n = min(40, nrow(collisions_C1)))

message("=== COLISÕES (Cluster 2) — var-var + var-score ===")
print(collisions_C2, n = min(40, nrow(collisions_C2)))




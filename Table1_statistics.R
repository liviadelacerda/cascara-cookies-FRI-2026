# =============================================================================
# Corrigendum - Food Research International 231 (2026) 118673
# Table 1: proximate composition, sugars and non-volatile compounds of the
# cascara flours - univariate statistics.
#
# Procedures:
#   - Shapiro-Wilk test on the residuals of the one-way ANOVA model
#   - Levene test (median-centred) for homoscedasticity of the residuals
#   - one-way ANOVA followed by Tukey's HSD test (alpha = 0.05), or
#   - Kruskal-Wallis test followed by Dunn's multiple pairwise comparisons
#     (two-tailed, alpha = 0.05)
#   Letters: column letter_order gives the direction used in Table 1
#   ("ascending": a = lowest mean; "descending": a = highest mean/mean rank).
#   Glucose and fructose were analysed only for the treatments in which they
#   were detected (nd values excluded).
#
# Input : NonVolatiles_cascara_corrected.xlsx (sheet novolatiles_bs, values only)
# Output: Table1_statistics_results.xlsx (Table 1 layout and detailed tests)
#
# Usage:
#   setwd("path/to/Data_and_scripts")
#   source("Table1_statistics.R")
#
# R packages: readxl, writexl (base R for all statistics)
# =============================================================================

if (!exists("input_file")) input_file <- "NonVolatiles_cascara_corrected.xlsx"
if (!exists("output_file")) output_file <- "Table1_statistics_results.xlsx"

alpha <- 0.05
trt_order <- c("PARCH", "NAT", "SC4", "SC8", "SC10", "SC16", "SC20")

# Variables, column names in the data file and test used in Table 1
vars <- data.frame(
  label  = c("Moisture", "Ashes", "Lipids", "Proteins", "TDF", "GLUC", "FRUC",
             "TPC", "FLAV", "ANT", "TAN", "5-CQA", "3-CQA", "CAF", "TRIG"),
  column = c("Moisture", "Ash", "LIP", "Ptn", "FA", "GLUC", "FRUC",
             "TPC", "FLAV", "ANT", "TAN", "CGA", "NCGA", "CAF", "TRIG"),
  test   = c("ANOVA", "KW", "KW", "ANOVA", "KW", "ANOVA", "ANOVA",
             "ANOVA", "ANOVA", "ANOVA", "ANOVA", "KW", "KW", "ANOVA", "ANOVA"),
  detected_only = c(rep(FALSE, 5), TRUE, TRUE, rep(FALSE, 8)),
  letter_order = c("ascending", "descending", "descending", "ascending",
                   "descending", "descending", "descending", "ascending",
                   "ascending", "ascending", "ascending", "descending",
                   "descending", "ascending", "ascending"),
  stringsAsFactors = FALSE
)

# ---- data -------------------------------------------------------------------
if (grepl("\\.csv$", input_file)) {
  dat <- read.csv(input_file, check.names = FALSE)
} else {
  dat <- as.data.frame(readxl::read_excel(input_file, sheet = "novolatiles_bs"))
}
names(dat)[1] <- "TRATAMENTO"
dat$TRATAMENTO <- trimws(dat$TRATAMENTO)

# ---- helpers ----------------------------------------------------------------
# Compact letter display: groups ordered by the statistic, letters assigned to
# maximal runs of consecutive groups that do not differ (underlining method).
cld_letters <- function(groups, nonsig) {
  k <- length(groups); runs <- list()
  for (i in seq_len(k)) {
    j <- i
    while (j < k && all(nonsig[groups[i:(j + 1)], groups[i:(j + 1)]])) j <- j + 1
    if (!any(vapply(runs, function(r) i >= r[1] && j <= r[2], logical(1))))
      runs[[length(runs) + 1]] <- c(i, j)
  }
  out <- setNames(rep("", k), groups)
  for (l in seq_along(runs)) {
    idx <- runs[[l]][1]:runs[[l]][2]
    out[groups[idx]] <- paste0(out[groups[idx]], letters[l])
  }
  out
}

# Levene test (median-centred, Brown-Forsythe form) on the ANOVA residuals
levene_median <- function(y, g) {
  dev <- abs(y - ave(y, g, FUN = median))
  a <- anova(lm(dev ~ g))
  c(F = a$`F value`[1], p = a$`Pr(>F)`[1])
}

dunn_pvalues <- function(y, g) {
  N <- length(y); r <- rank(y)
  Rbar <- tapply(r, g, mean); n <- tapply(r, g, length)
  ties <- table(y); C <- sum(ties^3 - ties) / (12 * (N - 1))
  lv <- names(Rbar)
  P <- matrix(1, length(lv), length(lv), dimnames = list(lv, lv))
  for (a in lv) for (b in lv) if (a != b) {
    z <- (Rbar[a] - Rbar[b]) / sqrt((N * (N + 1) / 12 - C) * (1 / n[a] + 1 / n[b]))
    P[a, b] <- 2 * pnorm(-abs(z))
  }
  list(p = P, mean_rank = Rbar)
}

# Tukey HSD pairwise p-values (pooled variance, studentized range distribution)
tukey_pvalues <- function(y, g) {
  tk <- TukeyHSD(aov(y ~ g))$g; lv <- levels(g)
  P <- matrix(1, length(lv), length(lv), dimnames = list(lv, lv))
  for (pr in rownames(tk)) {
    ab <- strsplit(pr, "-")[[1]]
    P[ab[1], ab[2]] <- P[ab[2], ab[1]] <- tk[pr, "p adj"]
  }
  P
}

fmt_p <- function(p) ifelse(p < 0.0001, "<0.0001",
                            ifelse(p < 0.001, sprintf("%.4f", p), sprintf("%.3f", p)))
# rounding half away from zero (as in Excel)
r2 <- function(x) sprintf("%.2f", sign(x) * floor(abs(x) * 100 + 0.5 + 1e-9) / 100)

# ---- analysis ---------------------------------------------------------------
details <- list(); letters_tab <- list(); tests_tab <- list()
for (v in seq_len(nrow(vars))) {
  lab <- vars$label[v]; col <- vars$column[v]
  d <- dat[, c("TRATAMENTO", col)]; names(d) <- c("g", "y")
  d$y <- suppressWarnings(as.numeric(d$y))
  if (vars$detected_only[v]) {
    det <- tapply(d$y, d$g, function(x) any(x > 0, na.rm = TRUE))
    d <- d[d$g %in% names(det)[det], ]
  }
  d <- d[!is.na(d$y), ]
  d$g <- factor(d$g, levels = intersect(trt_order, unique(d$g)))

  fit <- aov(y ~ g, data = d)
  sw  <- shapiro.test(residuals(fit))
  lev <- levene_median(residuals(fit), d$g)
  an  <- summary(fit)[[1]]
  p_anova <- an$`Pr(>F)`[1]
  kw <- kruskal.test(y ~ g, data = d)

  if (vars$test[v] == "ANOVA") {
    NS <- tukey_pvalues(d$y, d$g) >= alpha
    m <- tapply(d$y, d$g, mean)
    ord <- names(sort(m, decreasing = vars$letter_order[v] == "descending"))
    let <- cld_letters(ord, NS)
  } else {
    dn <- dunn_pvalues(d$y, d$g)
    NS <- dn$p >= alpha
    ord <- names(sort(dn$mean_rank, decreasing = vars$letter_order[v] == "descending"))
    let <- cld_letters(ord, NS)
  }

  mean_sd <- tapply(d$y, d$g, function(x) c(mean(x), sd(x)))
  for (t in trt_order) {
    val <- if (t %in% names(mean_sd)) {
      ms <- mean_sd[[t]]
      paste0(r2(ms[1]), let[t], " \u00b1 ", r2(ms[2]))
    } else "nd"
    letters_tab[[length(letters_tab) + 1]] <-
      data.frame(Cascara = t, Variable = lab, Value = val, stringsAsFactors = FALSE)
  }
  tests_tab[[lab]] <- data.frame(
    Variable = lab, Test = ifelse(vars$test[v] == "ANOVA", "ANOVA + Tukey HSD",
                                  "Kruskal-Wallis + Dunn"),
    n = nrow(d), Shapiro_W = sw$statistic, p_Shapiro = sw$p.value,
    Levene_F = lev["F"], p_Levene = lev["p"],
    ANOVA_F = an$`F value`[1], p_ANOVA = p_anova,
    KW_H = unname(kw$statistic), p_KW = kw$p.value,
    Letters = paste(paste0(names(let), "=", let), collapse = "; "),
    stringsAsFactors = FALSE)
}

long <- do.call(rbind, letters_tab)
table1 <- reshape(long, idvar = "Cascara", timevar = "Variable", direction = "wide")
names(table1) <- sub("^Value\\.", "", names(table1))
tests <- do.call(rbind, tests_tab); rownames(tests) <- NULL

foot <- data.frame(Cascara = c("p Shapiro-Wilk", "p Levene", "p ANOVA", "p Kruskal-Wallis"))
for (lab in vars$label) {
  tt <- tests[tests$Variable == lab, ]
  foot[[lab]] <- c(fmt_p(tt$p_Shapiro), fmt_p(tt$p_Levene),
                   ifelse(tt$Test == "ANOVA + Tukey HSD", fmt_p(tt$p_ANOVA), "-"),
                   ifelse(tt$Test == "ANOVA + Tukey HSD", "-", fmt_p(tt$p_KW)))
}
table1_out <- rbind(table1[, c("Cascara", vars$label)], foot)

print(table1_out, row.names = FALSE)
writexl::write_xlsx(list(Table1 = table1_out, Tests = tests), output_file)
cat("\nResults written to", output_file, "\n")

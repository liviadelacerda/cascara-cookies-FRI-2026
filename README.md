# Cascara cookies – data and R scripts (Food Research International, 2026)

Data and R scripts supporting the corrigendum to:

de Oliveira, L. L., de Sousa, J. F., da Costa, A. B., Alves, L. P., Resende, G. N., Oliveira, N. L., Lima, H. B. V., & Alves, E. A. (2026). Effects of post-harvest processing on chemical composition and sensory quality of coffee cascara, using clean-label cookies as a model bakery system. *Food Research International*, 231, 118673. https://doi.org/10.1016/j.foodres.2026.118673

## Correction

The total dietary fiber (TDF) values of the naturally dried (NAT) and parchment (PARCH) cascara samples were interchanged in the data set used for Table 1 and the Multiple Factor Analysis (MFA, Fig. 7) of the published article. The corrected values are 55.50 ± 0.19 g/100 g for NAT and 87.10 ± 0.66 g/100 g for PARCH (dry basis). The files in this repository contain the corrected data; the corrected TDF cells are highlighted in red in `NonVolatiles_cascara_corrected.xlsx`, which also lists the values used in the published article.

## Files

All spreadsheets contain values only (no formulas).

| File | Content |
|---|---|
| `Liking_CATA_cookies.xlsx` | Consumer liking and CATA responses (sheet `Sheet1`) |
| `Volatiles_relative_areas_top12.xlsx` | Relative peak areas of the 12 selected volatile compounds (sheet `Sheet1`) |
| `NonVolatiles_cascara_corrected.xlsx` | Proximate composition and non-volatile compounds of the cascara flours, dry basis (sheet `novolatiles_bs`; column `FA` = TDF) |
| `Texture_hardness.xlsx` | Cookie hardness (sheet `para MFA`) |
| `Table1_statistics.R` | Univariate statistics of Table 1: Shapiro–Wilk, Levene, one-way ANOVA + Tukey HSD or Kruskal–Wallis + Dunn |
| `cookies_integration_corrected_TDF.R` | Correspondence analysis of CATA data, MFA (Fig. 7) and PLS-R models (Supplementary Table S6) |

## Usage

Place the scripts and the spreadsheets in the same folder, set it as the working directory and run each script:

```r
setwd("path/to/cascara-cookies-FRI-2026")
source("Table1_statistics.R")                  # -> Table1_statistics_results.xlsx
source("cookies_integration_corrected_TDF.R")  # -> export_fig/ and exports_check/
```

`export_fig/mfa-biplot-clean.png` is the corrected Fig. 7.

R packages: readxl, writexl, dplyr, tidyr, purrr, stringr, readr, tibble, janitor, FactoMineR, factoextra, ade4, mixOmics (Bioconductor), ggplot2, ggrepel, patchwork, knitr, flextable, officer.

## License

Code: MIT License (see `LICENSE`). Data: Creative Commons Attribution 4.0 International (CC BY 4.0).

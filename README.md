# Hippo–YAP state in CD34+ CML cells treated with imatinib

This repository contains the R code used to generate the
in silico Hippo–YAP analyses reported in the manuscript:

> [Akbari S. *et al.*] Time-resolved Hippo–YAP and microRNA responses to imatinib in CML K562 cells with in silico validation in CD34⁺ CML progenitors. *PLOS ONE* (submitted).

The analyses quantify Hippo kinase and YAP/TEAD output scores in CD34+ CML
cells before vs 7 days of imatinib (GSE12211) and in baseline imatinib
responders vs non-responders (GSE14671).

## Repository structure

- `code/R/`
  - `00_setup.R` – package loading, helper functions, gene sets.
  - `01_prepare_GSE12211_GSE14671.R` – download GEO data, collapse probes to genes, compute Hippo/YAP scores.
  - `02_plot_GSE12211.R` – figures for paired pre vs day 7 analyses (GSE12211).
  - `03_plot_GSE14671.R` – figures for baseline responder vs non-responder analyses (GSE14671).
  - `04_stats_hippo_scores.R` – summary statistics and effect sizes (tables S3–S4).


## Reproducing the analyses

1. Set the working directory to the root of this repository in R.

2. Run:

```r
source("code/R/01_prepare_GSE12211_GSE14671.R")
source("code/R/02_plot_GSE12211.R")
source("code/R/03_plot_GSE14671.R")
source("code/R/04_stats_hippo_scores.R")

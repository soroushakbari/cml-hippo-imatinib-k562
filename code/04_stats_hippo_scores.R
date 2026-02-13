library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(rlang)

dir.create(file.path("results", "tables"), recursive = TRUE, showWarnings = FALSE)


df12211 <- readRDS(file.path("results", "rds", "GSE12211_scores.rds")) %>%
  dplyr::filter(!is.na(timepoint), !is.na(patient_id))

df12211$timepoint <- factor(df12211$timepoint, levels = c("pre", "day7"))

cohen_d_paired <- function(pre, post) {
  diff <- post - pre
  mean(diff) / sd(diff)
}

make_stats_12211 <- function(score_col, label) {
  wide <- df12211 %>%
    dplyr::select(patient_id, timepoint, !!score_col) %>%
    tidyr::pivot_wider(
      names_from  = timepoint,
      values_from = !!score_col
    )
  
  pre  <- wide$pre
  post <- wide$day7
  
  t_res <- t.test(post, pre, paired = TRUE)
  
  tibble::tibble(
    score_type = label,
    n_pairs    = length(pre),
    mean_pre   = mean(pre),
    sd_pre     = sd(pre),
    mean_day7  = mean(post),
    sd_day7    = sd(post),
    mean_delta = mean(post - pre),
    sd_delta   = sd(post - pre),
    t_stat     = unname(t_res$statistic),
    df         = unname(t_res$parameter),
    p_value    = unname(t_res$p.value),
    cohen_d    = cohen_d_paired(pre, post)
  )
}

stats_12211_hippo <- make_stats_12211(sym("hippo_score"), "Hippo kinase score")
stats_12211_yap   <- make_stats_12211(sym("yap_score"),   "YAP/TEAD output score")

stats_12211 <- dplyr::bind_rows(stats_12211_hippo, stats_12211_yap)
message("[06] GSE12211 summary stats:")
print(stats_12211)

readr::write_tsv(
  stats_12211,
  file.path("results", "tables", "Table_GSE12211_HippoYAP_scores.tsv")
)


df14671 <- readRDS(file.path("results", "rds", "GSE14671_scores.rds"))
df14671$response_group <- factor(df14671$response_group, levels = c("NR", "R"))

cohen_d_indep <- function(x, g) {
  g  <- droplevels(g)
  x1 <- x[g == levels(g)[1]]
  x2 <- x[g == levels(g)[2]]
  n1 <- length(x1)
  n2 <- length(x2)
  m1 <- mean(x1)
  m2 <- mean(x2)
  s1 <- stats::var(x1)
  s2 <- stats::var(x2)
  spooled <- sqrt(((n1 - 1) * s1 + (n2 - 1) * s2) / (n1 + n2 - 2))
  (m2 - m1) / spooled   # d > 0 means R > NR
}

make_stats_14671 <- function(score_col, label) {
  d <- df14671 %>%
    dplyr::select(response_group, !!score_col)
  
  grp_stats <- d %>%
    dplyr::group_by(response_group) %>%
    dplyr::summarise(
      n    = dplyr::n(),
      mean = mean(!!score_col),
      sd   = sd(!!score_col),
      .groups = "drop"
    )
  
  p_wil <- wilcox.test(
    d[[as_name(score_col)]] ~ d$response_group
  )$p.value
  
  d_cohen <- cohen_d_indep(
    d[[as_name(score_col)]],
    d$response_group
  )
  
  tibble::tibble(
    score_type = label,
    n_NR       = grp_stats$n[grp_stats$response_group == "NR"],
    mean_NR    = grp_stats$mean[grp_stats$response_group == "NR"],
    sd_NR      = grp_stats$sd[grp_stats$response_group == "NR"],
    n_R        = grp_stats$n[grp_stats$response_group == "R"],
    mean_R     = grp_stats$mean[grp_stats$response_group == "R"],
    sd_R       = grp_stats$sd[grp_stats$response_group == "R"],
    mean_diff  = mean_R - mean_NR,
    p_wilcox   = p_wil,
    cohen_d    = d_cohen
  )
}

stats_14671_hippo <- make_stats_14671(sym("hippo_score"), "Hippo kinase score")
stats_14671_yap   <- make_stats_14671(sym("yap_score"),   "YAP/TEAD output score")

stats_14671 <- dplyr::bind_rows(stats_14671_hippo, stats_14671_yap)
message("[06] GSE14671 summary stats:")
print(stats_14671)

readr::write_tsv(
  stats_14671,
  file.path("results", "tables", "Table_GSE14671_HippoYAP_scores.tsv")
)

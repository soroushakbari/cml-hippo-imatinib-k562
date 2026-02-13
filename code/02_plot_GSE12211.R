library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)

df12211_scores <- readRDS(file.path("results", "rds", "GSE12211_scores.rds")) %>%
  dplyr::filter(!is.na(timepoint), !is.na(patient_id))

df12211_scores$timepoint <- factor(df12211_scores$timepoint,
                                   levels = c("pre", "day7"))



df_long <- df12211_scores %>%
  dplyr::select(patient_id, timepoint, hippo_score, yap_score) %>%
  tidyr::pivot_longer(
    cols      = c(hippo_score, yap_score),
    names_to  = "score_type",
    values_to = "score"
  )

df_long$score_type <- factor(
  df_long$score_type,
  levels = c("hippo_score", "yap_score"),
  labels = c("Hippo kinase score", "YAP/TEAD output score")
)


calc_paired_p <- function(d) {
  wide <- tidyr::pivot_wider(
    d,
    names_from  = timepoint,
    values_from = score
  )
  t.test(wide$day7, wide$pre, paired = TRUE)$p.value
}

paired_list <- split(df_long, df_long$score_type)
paired_pvals <- tibble::tibble(
  score_type = names(paired_list),
  t_pvalue   = vapply(paired_list, calc_paired_p, numeric(1))
)
message("[02] Paired t-test p-values (GSE12211):")
print(paired_pvals)

p1 <- ggplot(df_long,
             aes(x = timepoint, y = score, group = patient_id)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2) +
  facet_wrap(~ score_type, scales = "free_y") +
  labs(
    x = "Timepoint (imatinib 400 mg/day)",
    y = "Mean normalized expression",
    title = "Hippo–YAP state in CD34+ CML cells\nbefore vs 7 days of imatinib (GSE12211)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text      = element_text(face = "bold")
  )

ggsave(
  filename = file.path("results", "figures", "Fig_GSE12211_paired_scores.png"),
  plot     = p1,
  width    = 7,
  height   = 4.5,
  dpi      = 600
)

ggsave(
  filename = file.path("results", "figures", "Fig_GSE12211_paired_scores.pdf"),
  plot     = p1,
  width    = 7,
  height   = 4.5
)


seg <- df12211_scores %>%
  dplyr::select(patient_id, timepoint, hippo_score, yap_score) %>%
  tidyr::pivot_wider(
    names_from  = timepoint,
    values_from = c(hippo_score, yap_score)
  ) %>%
  dplyr::filter(
    !is.na(hippo_score_pre),
    !is.na(hippo_score_day7),
    !is.na(yap_score_pre),
    !is.na(yap_score_day7)
  )

p2 <- ggplot(seg,
             aes(x = hippo_score_pre, y = yap_score_pre)) +
  geom_point(size = 2) +
  geom_segment(
    aes(xend = hippo_score_day7,
        yend = yap_score_day7),
    arrow = arrow(length = unit(0.15, "cm")),
    alpha = 0.8
  ) +
  geom_point(
    aes(x = hippo_score_day7,
        y = yap_score_day7),
    size = 2
  ) +
  labs(
    x = "Hippo kinase score (mean normalized expr)",
    y = "YAP/TEAD output score (mean normalized expr)",
    title = "Imatinib shifts Hippo–YAP state of CD34+ CML cells\n(GSE12211, pre \u2192 day 7)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  filename = file.path("results", "figures", "Fig_GSE12211_state_space.png"),
  plot     = p2,
  width    = 5.5,
  height   = 5,
  dpi      = 600
)

ggsave(
  filename = file.path("results", "figures", "Fig_GSE12211_state_space.pdf"),
  plot     = p2,
  width    = 5.5,
  height   = 5
)

message("[02] Figures for GSE12211 saved in results/figures/")

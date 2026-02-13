library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)


dir.create(file.path("results", "figures"), recursive = TRUE, showWarnings = FALSE)

df14671_scores <- readRDS(file.path("results", "rds", "GSE14671_scores.rds"))

df14671_scores$response_group <- factor(
  df14671_scores$response_group,
  levels = c("NR", "R")
)

df14671_long <- df14671_scores %>%
  dplyr::select(sample, response_group, hippo_score, yap_score) %>%
  tidyr::pivot_longer(
    cols      = c(hippo_score, yap_score),
    names_to  = "score_type",
    values_to = "score"
  )

df14671_long$score_type <- factor(
  df14671_long$score_type,
  levels = c("hippo_score", "yap_score"),
  labels = c("Hippo kinase score", "YAP/TEAD output score")
)


compute_cohen_d <- function(x, g) {
  g <- droplevels(g)
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

stats_df <- df14671_long %>%
  dplyr::group_by(score_type) %>%
  dplyr::summarise(
    p_value = wilcox.test(score ~ response_group)$p.value,
    cohen_d = compute_cohen_d(score, response_group),
    y_pos   = max(score, na.rm = TRUE) +
      0.03 * diff(range(score, na.rm = TRUE)),
    .groups = "drop"
  )

message("[03] GSE14671 Wilcoxon p-values and Cohen's d:")
print(stats_df)


p3 <- ggplot(df14671_long,
             aes(x = response_group, y = score)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.7, size = 1.8) +
  facet_wrap(~ score_type, scales = "free_y") +
  labs(
    x = "Baseline imatinib response group",
    y = "Mean normalized expression",
    title = "Baseline Hippo–YAP signatures in responders vs non-responders\n(GSE14671)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title      = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    strip.text      = element_text(face = "bold")
  ) +
  geom_text(
    data        = stats_df,
    aes(
      x     = 1.5,
      y     = y_pos,
      label = sprintf("Wilcoxon p = %.2f\nCohen's d = %.2f", p_value, cohen_d)
    ),
    size        = 3,
    inherit.aes = FALSE
  )

ggsave(
  filename = file.path("results", "figures", "Fig_GSE14671_box_scores.png"),
  plot     = p3,
  width    = 7,
  height   = 4.5,
  dpi      = 600
)

ggsave(
  filename = file.path("results", "figures", "Fig_GSE14671_box_scores.pdf"),
  plot     = p3,
  width    = 7,
  height   = 4.5
)

message("[03] Figure for GSE14671 saved in results/figures/")

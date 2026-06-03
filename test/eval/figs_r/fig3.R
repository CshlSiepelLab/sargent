#!/usr/bin/env Rscript
## fig3.R — Cross-ancestry portability (manuscript Figure 3).
##
## Two panels (more can be added later toward the 4-5-panel target):
##   A: target-injection axis at F_ST=0.16, typed=0.70.  Grouped bars over
##      mix_B ∈ {0, 0.05, 0.20} for 4 arms.
##   B: F_ST ladder at typed=0.70, mix_B=0.20 (the predictive regime).
##      Lines for 4 arms with light error bars.
##
## Input: test/eval/block1_1_v2_results.csv  (committed in ad348d5).
## Output: ~/newwriting/arg-mapping/figs2/fig3.{pdf,png}

args_all  <- commandArgs(trailingOnly = FALSE)
file_arg  <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

csv_path <- file.path(script_dir, "..", "block1_1_v2_results.csv")
df <- read.csv(csv_path, check.names = FALSE)

## Long-format with one row per (arm, cell).
arms_df <- bind_rows(
  data.frame(df, arm = "bayesR",  r2c = df$bayesR_r2c),
  data.frame(df, arm = "gblup",   r2c = df$gblup_r2c),
  data.frame(df, arm = "argunif", r2c = df$argunif_r2c),
  data.frame(df, arm = "sargent", r2c = df$sargent_r2c)
)
arms_df$arm_lab <- arm_factor(arms_df$arm)

mean_se <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(mean = NA_real_, se = NA_real_))
  c(mean = mean(x),
    se   = if (length(x) > 1) sd(x) / sqrt(length(x)) else 0)
}

## ---- Panel A: target-injection axis at F_ST=0.16, typed=0.70 ----
panel_a <- arms_df %>%
  filter(abs(fst - 0.16) < 1e-6, abs(typed - 0.70) < 1e-6,
         mix %in% c(0.0, 0.05, 0.20)) %>%
  group_by(arm_lab, mix) %>%
  summarise(mean = mean(r2c, na.rm = TRUE),
            se   = if (n() > 1) sd(r2c, na.rm = TRUE) / sqrt(n()) else 0,
            .groups = "drop") %>%
  mutate(mix_lab = factor(mix,
                          levels = c(0.0, 0.05, 0.20),
                          labels = c("mix[B]==0",
                                     "mix[B]==0.05",
                                     "mix[B]==0.20")))

pA <- ggplot(panel_a, aes(x = mix_lab, y = mean, fill = arm_lab, color = arm_lab)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72,
           color = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.78), width = 0.18,
                linewidth = 0.45, alpha = 0.55, show.legend = FALSE) +
  scale_x_discrete(labels = scales::parse_format()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     limits = c(0, 0.33),
                     breaks = seq(0, 0.30, by = 0.10)) +
  scale_color_arm() +
  scale_fill_arm() +
  labs(x = NULL,
       y = expression(r^2),
       tag = "A") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

## ---- Panel B: F_ST ladder at typed=0.70, mix_B=0.20 ----
panel_b <- arms_df %>%
  filter(abs(typed - 0.70) < 1e-6, abs(mix - 0.20) < 1e-6) %>%
  group_by(arm_lab, fst) %>%
  summarise(mean = mean(r2c, na.rm = TRUE),
            se   = if (n() > 1) sd(r2c, na.rm = TRUE) / sqrt(n()) else 0,
            .groups = "drop")

## Small jitter ONLY on the error bars so they can be distinguished;
## points and lines remain on the un-jittered F_ST values.
## Order chosen so the two pairs that often coincide in y
## (SARGENT-U/GBLUP green-blue, and BayesR/SARGENT red-purple) take
## interleaved ranks 1/3 and 2/4, giving equal within-pair spacing.
arm_dodge_order <- c("argunif", "bayesR", "gblup", "sargent")
arm_dodge <- setNames(seq(-0.004, 0.004, length.out = length(arm_dodge_order)),
                      arm_labels[arm_dodge_order])
panel_b$fst_jit <- panel_b$fst + arm_dodge[as.character(panel_b$arm_lab)]

pB <- ggplot(panel_b, aes(x = fst, y = mean, color = arm_lab, shape = arm_lab)) +
  geom_errorbar(aes(x = fst_jit, ymin = mean - se, ymax = mean + se),
                width = 0.006, linewidth = 0.45, alpha = 0.55) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4, stroke = 0.4, fill = "white") +
  scale_color_arm() +
  scale_shape_arm() +
  scale_x_continuous(breaks = sort(unique(panel_b$fst))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     limits = c(0, 0.33),
                     breaks = seq(0, 0.30, by = 0.10)) +
  labs(x = expression(F[ST]),
       y = expression(r^2),
       tag = "B") +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

fig3 <- pA + pB + plot_layout(widths = c(1, 1.05))

save_fig(fig3, "fig3", width = 6.5, height = 3.3)

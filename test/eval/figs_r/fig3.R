#!/usr/bin/env Rscript
## fig3.R — Cross-ancestry portability (manuscript Figure 3).
##
## Two panels:
##   A: Target-injection axis at F_ST=0.16, typed=0.70.  Grouped bars over
##      mix_B ∈ {0, 0.05, 0.20} for 4 arms.  Shows that target injection
##      is what unlocks cross-ancestry prediction.
##   B: F_ST ladder at typed=0.70, mix_B=0.20 (the predictive regime).
##      Lines for 5 arms — BayesR, GBLUP, SARGENT-U, SARGENT-W, SARGENT —
##      with light error bars.  Story: at the predictive regime,
##      SARGENT-U ≈ SARGENT-W ≫ SARGENT-full; β is the cross-ancestry
##      liability, per-tree weighting is incidental once β is dropped.
##
## Inputs:
##   test/eval/block1_1_v2_results.csv           (Block 1.1 v2 baseline grid)
##   test/eval/block1_1_v2_variants_results.csv  (SARGENT-W rescue test)
## Output: ~/newwriting/arg-mapping/figs2/fig3.{pdf,png}

args_all  <- commandArgs(trailingOnly = FALSE)
file_arg  <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

csv_path <- file.path(script_dir, "..", "block1_1_v2_results.csv")
df <- read.csv(csv_path, check.names = FALSE)

variants_csv <- file.path(script_dir, "..", "block1_1_v2_variants_results.csv")
have_W <- file.exists(variants_csv)
if (have_W) {
  vdf <- read.csv(variants_csv, check.names = FALSE)
  df <- merge(df, vdf[, c("fst","mix","typed","seed","sargentW_r2c")],
              by = c("fst","mix","typed","seed"))
}

## Long-format with one row per (arm, cell)
arms_list <- list(
  data.frame(df, arm = "bayesR",  r2c = df$bayesR_r2c),
  data.frame(df, arm = "gblup",   r2c = df$gblup_r2c),
  data.frame(df, arm = "argunif", r2c = df$argunif_r2c),
  data.frame(df, arm = "sargent", r2c = df$sargent_r2c)
)
if (have_W) {
  arms_list[[5]] <- data.frame(df, arm = "sargentW", r2c = df$sargentW_r2c)
}
arms_df <- bind_rows(arms_list)

## Extended labels/colors that include SARGENT-W
arm_labels_ext  <- c(arm_labels,  sargentW = "SARGENT-W")
arm_palette_ext <- c(arm_palette, sargentW = "#76B7B2")   # teal (vine palette)
arm_shapes_ext  <- c(arm_shapes,  sargentW = 18)          # diamond
arm_order_ext   <- c(arm_order,  "sargentW")

arms_df$arm_lab <- factor(arms_df$arm, levels = arm_order_ext,
                          labels = arm_labels_ext[arm_order_ext])

scale_color_arm_ext <- function(...) {
  scale_color_manual(name = NULL,
                     values = setNames(arm_palette_ext, arm_labels_ext[names(arm_palette_ext)]), ...)
}
scale_fill_arm_ext <- function(...) {
  scale_fill_manual(name = NULL,
                    values = setNames(arm_palette_ext, arm_labels_ext[names(arm_palette_ext)]), ...)
}
scale_shape_arm_ext <- function(...) {
  scale_shape_manual(name = NULL,
                     values = setNames(arm_shapes_ext, arm_labels_ext[names(arm_shapes_ext)]), ...)
}

## ---- Panel A: target-injection axis at F_ST=0.16, typed=0.70 (5 arms) ----
## Bar order: BayesR, GBLUP, SARGENT, SARGENT-U, SARGENT-W.
panel_a_order <- c("bayesR", "gblup", "sargent", "argunif", "sargentW")
panel_a <- arms_df %>%
  filter(abs(fst - 0.16) < 1e-6, abs(typed - 0.70) < 1e-6,
         mix %in% c(0.0, 0.05, 0.20)) %>%
  mutate(arm_lab = factor(arm,
                          levels = panel_a_order,
                          labels = arm_labels_ext[panel_a_order])) %>%
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
  scale_color_arm_ext() +
  scale_fill_arm_ext() +
  annotate("text", x = Inf, y = Inf,
           label = "F[ST] == 0.16",
           parse = TRUE, hjust = 1.15, vjust = 1.6,
           size = 3, color = "gray25") +
  labs(x = NULL, y = expression(r^2), tag = "A") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

## ---- Panel B: F_ST ladder at typed=0.70, mix_B=0.20 (5 arms) ----
panel_b <- arms_df %>%
  filter(abs(typed - 0.70) < 1e-6, abs(mix - 0.20) < 1e-6) %>%
  group_by(arm_lab, fst) %>%
  summarise(mean = mean(r2c, na.rm = TRUE),
            se   = if (n() > 1) sd(r2c, na.rm = TRUE) / sqrt(n()) else 0,
            .groups = "drop")

## Error-bar dodge so the often-overlapping pairs (BayesR/SARGENT-full;
## SARGENT-U/SARGENT-W; GBLUP) are distinguishable.
b_arm_order <- c("argunif", "sargentW", "bayesR", "gblup", "sargent")
b_dodge <- setNames(seq(-0.006, 0.006, length.out = length(b_arm_order)),
                    arm_labels_ext[b_arm_order])
panel_b$fst_jit <- panel_b$fst + b_dodge[as.character(panel_b$arm_lab)]

pB <- ggplot(panel_b, aes(x = fst, y = mean, color = arm_lab, shape = arm_lab)) +
  geom_errorbar(aes(x = fst_jit, ymin = mean - se, ymax = mean + se),
                width = 0.006, linewidth = 0.45, alpha = 0.55) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4, stroke = 0.4, fill = "white") +
  scale_color_arm_ext() +
  scale_shape_arm_ext() +
  scale_x_continuous(breaks = sort(unique(panel_b$fst))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     limits = c(0, 0.33),
                     breaks = seq(0, 0.30, by = 0.10)) +
  labs(x = expression(F[ST]), y = expression(r^2), tag = "B") +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

fig3 <- pA + pB + plot_layout(widths = c(1, 1.05))
save_fig(fig3, "fig3", width = 6.5, height = 3.3)

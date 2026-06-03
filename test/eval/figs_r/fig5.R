#!/usr/bin/env Rscript
## fig5.R — Genetic architecture characterization (manuscript Figure 5).
##
## Panels rendered here (WHEN + HOW axes):
##   B: Variance-capture curves stratified by allele-age bin (WHEN axis).
##      Data: Block 2.1 v2, n=1000, 3 seeds; age bins 0-5k / 5k-20k / >20k gen.
##   C: Variance-capture curves across architectures (HOW axis).
##      Data: Block 2.5.2 polygenicity sweep, bf ∈ {0.005, 0.05, 0.5}, 30 cells.
##      Estimator curves fit without fused-lasso (pre-Block 2.6); truth curves
##      demonstrate separation; estimator separation requires fused-lasso default.
##
## Panel A (server-side) dropped in once block2_2_v1_results.csv lands.
##
## Output: ~/newwriting/arg-mapping/figs2/fig5_pBC.{pdf,png}   (Panels B+C)
##         ~/newwriting/arg-mapping/figs2/fig5.{pdf,png}        (all 3 panels, once A lands)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

args_all   <- commandArgs(trailingOnly = FALSE)
file_arg   <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

cache <- file.path(script_dir, "_cache")
dB    <- read.table(file.path(cache, "fig5_panelB.tsv"), header = TRUE, sep = "\t")
dC    <- read.table(file.path(cache, "fig5_panelC.tsv"), header = TRUE, sep = "\t")


## ============================================================
## Panel B — WHEN axis: age-stratified variance-capture curves
## ============================================================

## Age bin ordering and labels
bin_order  <- c("recent", "middle", "ancient")
bin_labels <- c(
  recent  = "Recent\n(≤5k gen)",
  middle  = "Middle\n(5k–20k gen)",
  ancient = "Ancient\n(>20k gen)"
)
dB$age_bin <- factor(dB$age_bin, levels = bin_order, labels = bin_labels[bin_order])

## Annotation: S_bin range per age bin (across seeds)
bin_n <- dB %>%
  filter(kind == "est") %>%
  group_by(age_bin) %>%
  summarise(S_min = min(S_bin), S_max = max(S_bin), .groups = "drop")

## Interpolate each seed curve onto a common k_frac grid for mean computation
common_k <- seq(0, 1, length.out = 300)

interp_seed <- function(df) {
  approx(x = df$k_frac, y = df$captured_frac,
         xout = common_k, yleft = 0, yright = 1)$y
}

dB_mean <- dB %>%
  group_by(age_bin, kind, seed) %>%
  arrange(k_frac) %>%
  do(data.frame(k_frac = common_k,
                captured_frac = interp_seed(.))) %>%
  group_by(age_bin, kind, k_frac) %>%
  summarise(captured_frac = mean(captured_frac, na.rm = TRUE),
            .groups = "drop")

## Per-seed thin lines for uncertainty band
dB$kind_f <- factor(dB$kind, levels = c("truth", "est"))

kind_colors <- c(truth = "#F28E2B", est = "#4E79A7")
kind_lty    <- c(truth = "solid",   est = "dashed")
kind_labels <- c(truth = "Truth", est = "Estimator")

pB <- ggplot() +
  ## Per-seed thin lines
  geom_line(data = dB,
            aes(x = k_frac, y = captured_frac,
                color = kind, linetype = kind,
                group = interaction(kind, seed)),
            linewidth = 0.25, alpha = 0.35) +
  ## Mean thick lines
  geom_line(data = dB_mean,
            aes(x = k_frac, y = captured_frac,
                color = kind, linetype = kind),
            linewidth = 0.9) +
  ## S_bin annotation in strip
  scale_color_manual(name = NULL, values = kind_colors, labels = kind_labels) +
  scale_linetype_manual(name = NULL, values = kind_lty,    labels = kind_labels) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
  facet_wrap(~ age_bin, nrow = 1) +
  labs(x = "Top-ranked fraction of trees (within bin)",
       y = "Fraction of variance captured",
       tag = "B") +
  theme(legend.position  = "bottom",
        legend.key.width = unit(0.6, "cm"),
        strip.text       = element_text(size = 7.5)) +
  guides(color = guide_legend(override.aes = list(linewidth = 0.9)))


## ============================================================
## Panel C — HOW axis: variance-capture by architecture
## ============================================================

bf_order  <- c("0.005", "0.05", "0.5")
bf_labels <- c(
  "0.005" = "Oligogenic\n(bf = 0.5%)",
  "0.05"  = "Intermediate\n(bf = 5%)",
  "0.5"   = "Polygenic\n(bf = 50%)"
)
bf_colors <- c(
  "0.005" = "#E15759",
  "0.05"  = "#F28E2B",
  "0.5"   = "#4E79A7"
)

dC$bf <- as.character(dC$bf)
dC$bf_f <- factor(dC$bf, levels = bf_order, labels = bf_labels[bf_order])

## Interpolate each (bf, kind, seed, demog) curve onto common log-spaced grid
common_k_log <- 10^seq(-3, 0, length.out = 300)

interp_cell <- function(df) {
  approx(x = df$k_frac, y = df$captured_frac,
         xout = common_k_log, yleft = 0, yright = 1, rule = 2)$y
}

dC_mean <- dC %>%
  group_by(bf, bf_f, kind, seed, demog) %>%
  arrange(k_frac) %>%
  do(data.frame(k_frac = common_k_log,
                captured_frac = interp_cell(.))) %>%
  group_by(bf, bf_f, kind, k_frac) %>%
  summarise(captured_frac = mean(captured_frac, na.rm = TRUE),
            .groups = "drop")

## Per-cell thin lines (seed × demog combos)
dC_thin <- dC %>%
  group_by(bf, bf_f, kind, seed, demog) %>%
  arrange(k_frac) %>%
  do(data.frame(k_frac = common_k_log,
                captured_frac = interp_cell(.))) %>%
  ungroup()

pC <- ggplot() +
  geom_line(data = dC_thin,
            aes(x = k_frac, y = captured_frac,
                color = bf_f, linetype = kind,
                group = interaction(bf_f, kind, seed, demog)),
            linewidth = 0.20, alpha = 0.25) +
  geom_line(data = dC_mean,
            aes(x = k_frac, y = captured_frac,
                color = bf_f, linetype = kind),
            linewidth = 0.85) +
  scale_color_manual(name = "Architecture",
                     values = setNames(bf_colors[bf_order], bf_labels[bf_order])) +
  scale_linetype_manual(name = NULL,
                        values = c(truth = "solid", est = "dashed"),
                        labels = c(truth = "Truth", est = "Estimator†")) +
  scale_x_log10(limits = c(0.001, 1),
                breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
  labs(x = "Top-ranked fraction of trees",
       y = "Fraction of variance captured",
       tag = "C",
       caption = "†Estimator fit without fused-lasso regularization;\n  architectural separation requires --ws2-l2 ≥ 0.1.") +
  theme(legend.position  = "bottom",
        legend.box       = "vertical",
        legend.key.width = unit(0.6, "cm"),
        plot.caption     = element_text(size = 6, hjust = 0, color = "gray40"),
        legend.title     = element_text(size = 7.5),
        legend.text      = element_text(size = 6.5)) +
  guides(color    = guide_legend(nrow = 1,
                                 override.aes = list(linewidth = 0.9, linetype = "solid")))


## ============================================================
## Combine B + C (interim, no Panel A yet)
## ============================================================

pBC <- pB + pC +
  plot_layout(ncol = 2, widths = c(1, 1))

save_fig(pBC, "fig5_pBC", width = 6.5, height = 3.5)


## ============================================================
## Full Fig 5 (placeholder until Panel A from server lands)
## ============================================================
## When test/eval/block2_2_v1_results.csv exists, drop Panel A in here:
##
## panelA_csv <- file.path(script_dir, "..", "block2_2_v1_results.csv")
## if (file.exists(panelA_csv)) {
##   dA <- read.table(panelA_csv, header=TRUE, sep="\t")
##   # ... build pA ...
##   pfull <- pA / (pB + pC) + plot_layout(heights = c(1, 1.5))
##   save_fig(pfull, "fig5", width = 6.5, height = 6.0)
## }

message("Done. Rendered fig5_pBC.{pdf,png}")

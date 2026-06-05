#!/usr/bin/env Rscript
## fig5.R — Genetic architecture characterization (manuscript Figure 5).
##
## Panel B (WHEN axis): 3×3 facet grid — rows = true architecture condition
##   (age_recent / age_middle / age_ancient), cols = observed age bin
##   (recent ≤5k / middle 5k-20k / ancient >20k gen). Diagonal facets show
##   steep capture curves; off-diagonal facets show flat (no signal). Data
##   from fig5_panelB sweep, n=500, 5 seeds.
##
## Panel C (HOW axis): variance-capture by polygenicity architecture.
##   bf ∈ {0.005, 0.05, 0.5}, 30 cells. Shows truth + FL estimator + no-FL
##   estimator to demonstrate that fused-lasso regularization (--ws2-l2 0.5)
##   is needed for architectural separation.
##
## Output: ~/newwriting/arg-mapping/figs2/fig5_pBC.{pdf,png}

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
## Panel B — WHEN axis: 3×3 architecture × age-bin grid
## ============================================================

bin_order  <- c("recent", "middle", "ancient")
arch_order <- c("age_recent", "age_middle", "age_ancient")

bin_labels <- c(
  recent  = "Observed: Recent\n(≤5k gen)",
  middle  = "Observed: Middle\n(5k–20k gen)",
  ancient = "Observed: Ancient\n(>20k gen)"
)
arch_labels <- c(
  age_recent  = "True: Recent",
  age_middle  = "True: Middle",
  age_ancient = "True: Ancient"
)

dB$age_bin        <- factor(dB$age_bin, levels = bin_order,  labels = bin_labels[bin_order])
dB$arch_condition <- factor(dB$arch,   levels = arch_order,  labels = arch_labels[arch_order])

## Interpolate each (arch, age_bin, kind, seed) curve onto a common linear grid
common_k <- seq(0, 1, length.out = 300)

interp_seed <- function(df) {
  approx(x = df$k_frac, y = df$captured_frac,
         xout = common_k, yleft = 0, yright = 1, rule = 2)$y
}

dB_mean <- dB %>%
  group_by(arch_condition, age_bin, kind, seed) %>%
  arrange(k_frac) %>%
  do(data.frame(k_frac = common_k,
                captured_frac = interp_seed(.))) %>%
  group_by(arch_condition, age_bin, kind, k_frac) %>%
  summarise(captured_frac = mean(captured_frac, na.rm = TRUE), .groups = "drop")

kind_colors_B <- c(truth = "#F28E2B", est = "#4E79A7")
kind_lty_B    <- c(truth = "solid",   est = "dashed")
kind_labels_B <- c(truth = "Truth",   est = "Estimator")

pB <- ggplot() +
  geom_line(data = dB,
            aes(x = k_frac, y = captured_frac,
                color = kind, linetype = kind,
                group = interaction(kind, seed)),
            linewidth = 0.20, alpha = 0.30) +
  geom_line(data = dB_mean,
            aes(x = k_frac, y = captured_frac,
                color = kind, linetype = kind),
            linewidth = 0.85) +
  scale_color_manual(name = NULL, values = kind_colors_B, labels = kind_labels_B) +
  scale_linetype_manual(name = NULL, values = kind_lty_B, labels = kind_labels_B) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", "0.5", "1")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
  facet_grid(arch_condition ~ age_bin) +
  labs(x = "Top-ranked fraction of trees (within bin)",
       y = "Fraction of variance captured",
       tag = "B") +
  theme(legend.position  = "bottom",
        legend.key.width = unit(0.6, "cm"),
        strip.text       = element_text(size = 6.5)) +
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

dC$bf   <- as.character(dC$bf)
dC$bf_f <- factor(dC$bf, levels = bf_order, labels = bf_labels[bf_order])

## Map kind to display label
kind_lty_C    <- c(truth = "solid", est_fl = "dashed", est = "dotted")
kind_labels_C <- c(truth = "Truth", est_fl = "Estimator (FL)", est = "Estimator (no FL)")

## Only include kinds present in the data
present_kinds <- intersect(c("truth", "est_fl", "est"), unique(dC$kind))
kind_lty_C    <- kind_lty_C[present_kinds]
kind_labels_C <- kind_labels_C[present_kinds]

## Interpolate each (bf, kind, seed, demog) curve onto common log-spaced grid
common_k_log <- 10^seq(-3, 0, length.out = 300)

interp_cell <- function(df) {
  approx(x = df$k_frac, y = df$captured_frac,
         xout = common_k_log, yleft = 0, yright = 1, rule = 2)$y
}

dC_mean <- dC %>%
  filter(kind %in% present_kinds) %>%
  group_by(bf, bf_f, kind, seed, demog) %>%
  arrange(k_frac) %>%
  do(data.frame(k_frac = common_k_log,
                captured_frac = interp_cell(.))) %>%
  group_by(bf, bf_f, kind, k_frac) %>%
  summarise(captured_frac = mean(captured_frac, na.rm = TRUE), .groups = "drop")

dC_thin <- dC %>%
  filter(kind %in% present_kinds) %>%
  group_by(bf, bf_f, kind, seed, demog) %>%
  arrange(k_frac) %>%
  do(data.frame(k_frac = common_k_log,
                captured_frac = interp_cell(.))) %>%
  ungroup()

## Weight: FL and truth are primary (thicker); no-FL is reference (thinner)
dC_mean$lw <- ifelse(dC_mean$kind == "est", 0.45, 0.85)
dC_thin$lw <- ifelse(dC_thin$kind == "est", 0.15, 0.20)
dC_thin$al <- ifelse(dC_thin$kind == "est", 0.15, 0.25)

pC <- ggplot() +
  geom_line(data = dC_thin,
            aes(x = k_frac, y = captured_frac,
                color = bf_f, linetype = kind,
                group = interaction(bf_f, kind, seed, demog),
                linewidth = I(lw), alpha = I(al))) +
  geom_line(data = dC_mean,
            aes(x = k_frac, y = captured_frac,
                color = bf_f, linetype = kind,
                linewidth = I(lw))) +
  scale_color_manual(name = "Architecture",
                     values = setNames(bf_colors[bf_order], bf_labels[bf_order])) +
  scale_linetype_manual(name = NULL,
                        values = kind_lty_C,
                        labels = kind_labels_C) +
  scale_x_log10(limits = c(0.001, 1),
                breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = c(0, 0.5, 1)) +
  labs(x = "Top-ranked fraction of trees",
       y = "Fraction of variance captured",
       tag = "C") +
  theme(legend.position  = "bottom",
        legend.box       = "vertical",
        legend.key.width = unit(0.6, "cm"),
        legend.title     = element_text(size = 7.5),
        legend.text      = element_text(size = 6.5)) +
  guides(color = guide_legend(nrow = 1,
                              override.aes = list(linewidth = 0.9, linetype = "solid")),
         linetype = guide_legend(override.aes = list(linewidth = 0.7)))


## ============================================================
## Combine B + C
## ============================================================

pBC <- pB + pC +
  plot_layout(ncol = 2, widths = c(1.5, 1))

save_fig(pBC, "fig5_pBC", width = 8.5, height = 5.0)

message("Done. Rendered fig5_pBC.{pdf,png}")

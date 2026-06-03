#!/usr/bin/env Rscript
## fig4.R — Localization headline (manuscript Figure 4).
##
## Four panels:
##   A (top, full width): Manhattan-style ŵ²_s along the genome on a
##     representative Block 1.2 cell; gold shading marks true active blocks.
##   B: Variance-capture curves (sparse vs diffuse arch., truth vs est).
##   C: Localization AUC across (contrast × demography), Block 1.2 v1.
##   D: Fused-lasso opt-in AUC contrast (polygenic vs oligogenic).
##
## Pre-computed inputs from preprocess_fig4.py in _cache/.
## Output: ~/newwriting/arg-mapping/figs2/fig4.{pdf,png}

args_all  <- commandArgs(trailingOnly = FALSE)
file_arg  <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

cache <- file.path(script_dir, "_cache")
dA    <- read.table(file.path(cache, "fig4_panelA.tsv"), header = TRUE, sep = "\t")
dB    <- read.table(file.path(cache, "fig4_panelB.tsv"), header = TRUE, sep = "\t")
dC    <- read.table(file.path(cache, "fig4_panelC.tsv"), header = TRUE, sep = "\t")
dD    <- read.table(file.path(cache, "fig4_panelD.tsv"), header = TRUE, sep = "\t")
has_fl <- any(dD$mode == "fused_lasso")
if (!has_fl) {
  message("WARNING: no fused_lasso rows in Panel D data — rendering default bars only.")
}


## ---- Panel A: Manhattan of estimated ŵ²_s ----

## Extract contiguous true-active blocks for shading
blocks <- with(dA, {
  is_act <- active == 1
  if (!any(is_act)) return(data.frame(xmin = numeric(0), xmax = numeric(0)))
  rl <- rle(is_act)
  ends   <- cumsum(rl$lengths)
  starts <- c(1, head(ends, -1) + 1)
  keep   <- rl$values
  data.frame(xmin = pos_mid_kb[starts[keep]],
             xmax = pos_mid_kb[ends[keep]])
})

pA <- ggplot(dA, aes(x = pos_mid_kb))
if (nrow(blocks) > 0) {
  pA <- pA +
    geom_rect(data = blocks,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = "#fee08b", alpha = 0.40,
              inherit.aes = FALSE)
}
## Mark the top-K trees by estimator, where K = number of true active trees.
## This emphasizes the signal where the estimator agrees with truth.
n_active <- sum(dA$active == 1)
top_thresh <- if (n_active > 0)
  sort(dA$ws2_est, decreasing = TRUE)[n_active] else Inf
dA$est_tier <- ifelse(dA$ws2_est >= top_thresh, "top", "floor")

## Rescale to display the same magnitudes in scientific notation on the
## y-axis ticks.  Linear y-axis collapses the inference noise floor below
## visibility while preserving the relative magnitudes of active trees.

pA <- pA +
  ## Truth: step line so active vs floor levels are obvious
  geom_step(aes(y = ws2_truth, color = "truth"),
            linewidth = 0.5, alpha = 0.9) +
  ## Estimator: line + small dots on linear scale
  geom_line(aes(y = ws2_est, color = "est"),
            linewidth = 0.3, alpha = 0.75) +
  geom_point(aes(y = ws2_est, color = "est"),
             size = 0.6, alpha = 0.85) +
  scale_color_manual(name = NULL,
                     breaks = c("truth", "est"),
                     values = c(truth = "#F28E2B",
                                est   = "#4E79A7"),
                     labels = c(truth = expression(True ~ italic(w)[s]^2),
                                est   = expression(Estimated ~ hat(italic(w))[s]^2))) +
  scale_y_continuous(labels = scales::label_scientific(),
                     expand = expansion(mult = c(0, 0.05))) +
  scale_x_continuous(expand = expansion(mult = c(0.005, 0.005))) +
  labs(x = "Genomic position (kb)",
       y = expression(italic(w)[s]^2),
       tag = "A") +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.key.height = unit(0.32, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))


## ---- Panel B: variance-capture curves ----

## Mean curve per (arch, kind) across cells, plus a thin band of individual cells.
## To keep the panel clean, interpolate every cell's (k, captured) onto a
## common grid and average.
common_k <- 10^seq(log10(0.001), 0, length.out = 200)

interp_one <- function(df) {
  approx(x = df$k_frac, y = df$captured_frac,
         xout = common_k, rule = 2)$y
}

dB_mean <- dB %>%
  group_by(arch, kind, cell) %>%
  do(data.frame(k_frac = common_k,
                captured_frac = interp_one(.))) %>%
  group_by(arch, kind, k_frac) %>%
  summarise(captured_frac = mean(captured_frac, na.rm = TRUE),
            .groups = "drop")

dB_mean$arch_kind <- factor(paste0(dB_mean$arch, " (", dB_mean$kind, ")"),
                            levels = c("sparse (truth)", "sparse (est)",
                                       "diffuse (truth)", "diffuse (est)"))

cap_palette <- c(
  "sparse (truth)"  = "#F28E2B",
  "sparse (est)"    = "#F28E2B",
  "diffuse (truth)" = "#4E79A7",
  "diffuse (est)"   = "#4E79A7"
)
cap_lty <- c(
  "sparse (truth)"  = "solid",
  "sparse (est)"    = "dashed",
  "diffuse (truth)" = "solid",
  "diffuse (est)"   = "dashed"
)

pB <- ggplot(dB_mean,
             aes(x = k_frac, y = captured_frac,
                 color = arch_kind, linetype = arch_kind)) +
  geom_line(linewidth = 0.7) +
  scale_color_manual(name = NULL, values = cap_palette) +
  scale_linetype_manual(name = NULL, values = cap_lty) +
  scale_x_log10(limits = c(0.001, 1),
                breaks = c(0.001, 0.01, 0.1, 1),
                labels = c("0.001", "0.01", "0.1", "1")) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25)) +
  labs(x = "Top-ranked fraction of trees",
       y = "Fraction of variance captured",
       tag = "B") +
  theme(legend.position = c(0.98, 0.02),
        legend.justification = c(1, 0),
        legend.key.height = unit(0.32, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))


## ---- Panel C: Localization AUC across contrast × demography ----

dC$contrast_f <- factor(as.character(as.integer(dC$contrast)),
                        levels = c("1000", "100", "10"))
dC$demog      <- factor(dC$demog, levels = c("admixed", "panmictic"))

dC_sum <- dC %>%
  group_by(contrast_f, demog) %>%
  summarise(auc_mean    = mean(auc, na.rm = TRUE),
            auc_se      = if (n() > 1) sd(auc, na.rm = TRUE) / sqrt(n()) else 0,
            ceil_mean   = mean(ceiling, na.rm = TRUE),
            .groups = "drop")

demog_palette <- c(admixed = "#B07AA1", panmictic = "#59A14F")

pC <- ggplot(dC_sum, aes(x = contrast_f, y = auc_mean,
                         fill = demog, color = demog)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72,
           color = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = auc_mean - auc_se, ymax = auc_mean + auc_se),
                position = position_dodge(width = 0.78), width = 0.18,
                linewidth = 0.45, alpha = 0.55, show.legend = FALSE) +
  geom_point(aes(y = ceil_mean),
             position = position_dodge(width = 0.78),
             shape = 23, size = 2.6, stroke = 0.7,
             fill = "white", show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             color = "gray40", linewidth = 0.4) +
  scale_fill_manual(name = NULL, values = demog_palette,
                    labels = c(admixed = "Admixed", panmictic = "Panmictic")) +
  scale_color_manual(name = NULL, values = demog_palette,
                     labels = c(admixed = "Admixed", panmictic = "Panmictic"),
                     guide = "none") +
  scale_y_continuous(limits = c(0.4, 1.02),
                     breaks = seq(0.4, 1.0, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = "Active-region variance enrichment",
       y = "Localization AUC",
       tag = "C") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))


## ---- Panel D: Fused-lasso opt-in (polygenic vs oligogenic) ----

dD$bf_label <- factor(dD$bf,
                      levels = c("0.5", "0.005"),
                      labels = c("polygenic\n(bf = 0.5)",
                                 "oligogenic\n(bf = 0.005)"))
dD$mode_lab <- factor(dD$mode,
                      levels = c("default", "fused_lasso"),
                      labels = c("Default", "Fused-lasso"))

dD_sum <- dD %>%
  group_by(bf_label, mode_lab) %>%
  summarise(auc_mean = mean(auc, na.rm = TRUE),
            auc_se   = if (n() > 1) sd(auc, na.rm = TRUE) / sqrt(n()) else 0,
            .groups = "drop")

## Compute gains for annotation when fused_lasso data is present
if (has_fl) {
  dD_raw <- dD %>%
    group_by(bf_label = factor(bf, levels = c("0.5", "0.005"),
                               labels = c("polygenic\n(bf = 0.5)",
                                          "oligogenic\n(bf = 0.005)")),
             mode) %>%
    summarise(auc_mean = mean(auc, na.rm = TRUE),
              auc_se   = if (n() > 1) sd(auc, na.rm = TRUE) / sqrt(n()) else 0,
              .groups = "drop")
  dD_gain <- dD_raw %>%
    tidyr::pivot_wider(names_from = mode, values_from = c(auc_mean, auc_se)) %>%
    mutate(gain  = auc_mean_fused_lasso - auc_mean_default,
           y_top = pmax(auc_mean_default, auc_mean_fused_lasso) +
                    pmax(auc_se_default, auc_se_fused_lasso, na.rm = TRUE) + 0.04)
} else {
  dD_gain <- data.frame(bf_label = character(0), y_top = numeric(0), gain = numeric(0))
}

mode_palette <- c(Default = "#999999", `Fused-lasso` = "#B07AA1")

pD <- ggplot(dD_sum, aes(x = bf_label, y = auc_mean,
                         fill = mode_lab, color = mode_lab)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7,
           color = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = auc_mean - auc_se, ymax = auc_mean + auc_se),
                position = position_dodge(width = 0.78), width = 0.18,
                linewidth = 0.45, alpha = 0.55, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             color = "gray40", linewidth = 0.4) +
  geom_text(data = dD_gain,
            aes(x = bf_label, y = y_top, label = sprintf("%+0.02f", gain)),
            inherit.aes = FALSE,
            color = "#B07AA1", size = 3, fontface = "bold") +
  scale_fill_manual(name = NULL, values = mode_palette) +
  scale_color_manual(name = NULL, values = mode_palette, guide = "none") +
  scale_y_continuous(limits = c(0.4, 1.02),
                     breaks = seq(0.4, 1.0, 0.2),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Localization AUC", tag = "D") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))


## ---- Compose ----

fig4 <- pA / (pB | pC | pD) + plot_layout(heights = c(1, 1.15))

save_fig(fig4, "fig4", width = 7.0, height = 5.4)

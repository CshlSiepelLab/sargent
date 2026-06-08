#!/usr/bin/env Rscript
# Fig 5 Panel B (per-edge architectures) review plot — shows BOTH:
#   left:  truth-side per-bin Δh² (clean diagonal by construction)
#   right: estimator-side per-bin Δh² posterior mean (collapses to ancient)
# Lets Adam compare truth vs estimator before deciding on direction.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readr)
})

script_dir <- "/Users/acs/git-src/sargent/test/eval/figs_r"
cache      <- file.path(script_dir, "_cache")
out_dir    <- file.path(script_dir, "_out")
dir.create(out_dir, showWarnings = FALSE)

est   <- read_tsv(file.path(cache, "fig5_panelB_edges.tsv"),       show_col_types = FALSE)
truth <- read_tsv(file.path(cache, "fig5_panelB_edges_truth.tsv"), show_col_types = FALSE)

arch_labels <- c(
  age_recent_edges  = "recent (young tertile)",
  age_middle_edges  = "middle tertile",
  age_ancient_edges = "ancient (deep tertile)"
)
bin_levels  <- c("recent", "middle", "ancient")

prep <- function(df, value_col, side_label) {
  df %>%
    mutate(true_arch    = factor(true_arch,
                                 levels = names(arch_labels),
                                 labels = arch_labels),
           observed_bin = factor(observed_bin, levels = bin_levels)) %>%
    group_by(true_arch, observed_bin) %>%
    summarise(mean_v = mean(.data[[value_col]]),
              se_v   = sd(.data[[value_col]]) / sqrt(n()),
              .groups = "drop") %>%
    mutate(side = side_label)
}

d_truth <- prep(truth, "truth_norm",      "Truth: Σ ws²_e · ℓ_e per bin")
d_est   <- prep(est,   "dh2_total_norm",  "Estimator: per-edge posterior Δh² per bin")

base_plot <- function(d, title) {
  ggplot(d, aes(x = observed_bin, y = mean_v, fill = observed_bin)) +
    geom_col(width = 0.7, colour = "grey25") +
    geom_errorbar(aes(ymin = pmax(mean_v - se_v, 0),
                      ymax = pmin(mean_v + se_v, 1)),
                  width = 0.2, linewidth = 0.3) +
    facet_wrap(~ true_arch, ncol = 3) +
    scale_y_continuous("Normalised Δh² per bin",
                       limits = c(0, 1.05),
                       expand = c(0, 0)) +
    scale_fill_manual(values = c(recent = "#fde725",
                                 middle = "#21918c",
                                 ancient = "#440154")) +
    labs(x = "Observed age bin (per-cell tertile)",
         title = title) +
    theme_classic(base_size = 11) +
    theme(legend.position = "none",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))
}

p_truth <- base_plot(d_truth, "(A) Truth side — per-cell tertile binning recovers diagonal")
p_est   <- base_plot(d_est,   "(B) Estimator side — posterior mean collapses to ancient bin")

p <- p_truth / p_est +
  plot_annotation(
    title    = "Fig 5 Panel B review (per-edge architectures, n=300, L=200kb, 3 seeds/arch)",
    subtitle = paste("Truth bin labels: 'recent' / 'middle' / 'ancient' are per-cell tertiles of",
                     "edge-midpoint age.\nEstimator bias: Δh²_e ∝ (ŵ²·ℓ_e·S_z)² · Σ(cᵢᵘ)² favours",
                     "long, many-descendant ancient edges by construction."),
    theme    = theme(plot.title    = element_text(face = "bold", size = 13),
                     plot.subtitle = element_text(size = 9, lineheight = 1.1))
  )

out_png <- file.path(out_dir, "fig5_panelB_edges_review.png")
ggsave(out_png, p, width = 9, height = 6.5, dpi = 140)
cat("wrote:", out_png, "\n")

#!/usr/bin/env Rscript
## fig_xancestry_age.R — "Cross-ancestry PRS failure has an allele-age signature".
##
## Four panels, all at F_ST≈0.30 (the deeper-split condition; F_ST≈0.15 is
## reported in source data for the supplementary):
##   A: Mechanism — cross-ancestry kinship by age stratum.  Bar at
##      ancient = 1.0; recent / middle = 0.  Cross-ancestry relatedness
##      lives only in deep shared branches.
##   B: THE MONEY PANEL — heatmap of the fraction of A-variants also
##      carried by B, by (age × MAF) tertile.  Horizontal stripes by age,
##      flat across frequency.  Annotates: partial corr(age|freq) z=15,
##      partial corr(freq|age) ~0 — age, not frequency, sets transfer.
##   C: Prediction payoff — ARG-BLUP cross-ancestry portability (cross-R² /
##      within-R²) by age (old vs young) and by frequency (common vs rare).
##      Old ports ~4× young; common > rare; ARG reproduces frequency and
##      adds age.
##   D: The cap — within-R² vs cross-R² connecting standard (h²=0.5, n=1000)
##      to high-power (h²=0.8, n=2800) for each F_ST.  Within climbs to
##      ~0.5; cross stays pinned at ~0.02.  Portability is bounded by the
##      shared fraction.
##
## Inputs (already committed):
##   test/eval/figs_r/xanc_data/xanc_panelA_crosskinship.csv
##   test/eval/figs_r/xanc_data/xanc_panelB_heatmap.csv
##   test/eval/figs_r/xanc_data/xanc_panelC_portability.csv
##   test/eval/figs_r/xanc_data/xanc_panelD_cap.csv
##   test/eval/figs_r/xanc_data/xanc_regression_summary.csv
##
## Output: ~/newwriting/arg-mapping/figs2/fig_xancestry_age.{pdf,png}

args_all  <- commandArgs(trailingOnly = FALSE)
file_arg  <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

suppressMessages({
  library(scales)
})

data_dir <- file.path(script_dir, "xanc_data")
dA <- read.csv(file.path(data_dir, "xanc_panelA_crosskinship.csv"))
dB <- read.csv(file.path(data_dir, "xanc_panelB_heatmap.csv"))
dC <- read.csv(file.path(data_dir, "xanc_panelC_portability.csv"))
dD <- read.csv(file.path(data_dir, "xanc_panelD_cap.csv"))
dR <- read.csv(file.path(data_dir, "xanc_regression_summary.csv"))

## Focus on the F_ST≈0.30 cell for the main panels.
fst_main <- 0.30

## Stratum ordering: ancient is the headline, put it on the right / top.
strat_levels <- c("recent", "middle", "ancient")
strat_labels <- c(recent = "Recent",
                  middle = "Middle",
                  ancient = "Ancient\n(deep / shared)")
strat_palette <- c(recent  = "#bdbdbd",
                   middle  = "#9ecae1",
                   ancient = "#3182bd")

maf_levels <- c("rare", "mid", "common")
maf_labels <- c(rare = "Rare\n(MAF < 0.05)",
                mid  = "Intermediate",
                common = "Common\n(MAF ≥ 0.25)")


## ---- Panel A — cross-kinship by age stratum (mechanism) ----
panel_a <- dA %>%
  filter(abs(fst - fst_main) < 1e-6) %>%
  mutate(stratum = factor(stratum, levels = strat_levels, labels = strat_labels[strat_levels]))

pA <- ggplot(panel_a, aes(x = stratum, y = cross_frac, fill = stratum)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.0f%%", 100 * cross_frac)),
            vjust = -0.6, size = 2.8, color = "gray20") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                     limits = c(0, 1.1),
                     breaks = seq(0, 1, 0.25),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = setNames(strat_palette, strat_labels[strat_levels]),
                    guide = "none") +
  labs(x = NULL,
       y = "Fraction of cross-ancestry\nkinship from stratum",
       tag = "A",
       subtitle = "Cross-ancestry relatedness lives only in deep / shared branches") +
  theme(plot.subtitle = element_text(size = 7.5, color = "gray25",
                                     margin = margin(b = 4)),
        axis.text.x = element_text(size = 7, lineheight = 1.0))


## ---- Panel B — age × MAF heatmap (the money panel) ----
panel_b <- dB %>%
  filter(abs(fst - fst_main) < 1e-6) %>%
  mutate(age_tertile = factor(age_tertile, levels = strat_levels,
                              labels = c("Recent", "Middle", "Ancient")),
         maf_tertile = factor(maf_tertile, levels = maf_levels,
                              labels = c("Rare", "Intermediate", "Common"))) %>%
  ## Reverse age so ancient sits at the TOP of the heatmap.
  mutate(age_tertile = factor(age_tertile,
                              levels = rev(c("Recent", "Middle", "Ancient"))))

regr <- dR %>% filter(abs(fst - fst_main) < 1e-6)
panel_b_subtitle <- sprintf(
  "Sharing stripes by age — flat across frequency\npartial corr(age | freq) z = %.1f;   partial corr(freq | age) z = %.1f",
  regr$z_age, regr$z_freq
)

pB <- ggplot(panel_b, aes(x = maf_tertile, y = age_tertile, fill = shared_pct)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(shared_pct >= 1,
                               sprintf("%.0f%%", shared_pct),
                               sprintf("%.2f%%", shared_pct)),
                color = shared_pct > 20),
            size = 2.7, fontface = "bold") +
  scale_fill_viridis_c(name = "% of A-variants\nalso in B",
                       trans = "sqrt",
                       limits = c(0, 70),
                       breaks = c(0, 1, 5, 20, 50),
                       option = "D",
                       direction = 1) +
  scale_color_manual(values = c(`FALSE` = "gray15", `TRUE` = "white"),
                     guide = "none") +
  coord_equal() +
  labs(x = "Allele frequency tertile",
       y = "Allele age tertile",
       tag = "B",
       subtitle = panel_b_subtitle) +
  theme(plot.subtitle = element_text(size = 7.0, color = "gray25",
                                     lineheight = 1.1,
                                     margin = margin(b = 4)),
        axis.text = element_text(size = 7),
        legend.position = "right",
        legend.key.height = unit(0.45, "cm"),
        legend.title = element_text(size = 7, lineheight = 1.0),
        panel.grid = element_blank())


## ---- Panel C — portability by age and by frequency ----
panel_c_levels <- c("baseline", "old", "young", "common", "rare")
panel_c_labels <- c(baseline = "Baseline\n(all)",
                    old      = "Old\n(ancient)",
                    young    = "Young\n(recent)",
                    common   = "Common",
                    rare     = "Rare")
panel_c_axis <- c(baseline = "Overall",
                  old      = "Age",
                  young    = "Age",
                  common   = "Frequency",
                  rare     = "Frequency")
panel_c_palette <- c(Overall   = "gray45",
                     Age       = "#3182bd",
                     Frequency = "#f0a040")

panel_c <- dC %>%
  filter(abs(fst - fst_main) < 1e-6,
         stratum %in% panel_c_levels) %>%
  mutate(stratum = factor(stratum, levels = panel_c_levels,
                          labels = panel_c_labels[panel_c_levels]),
         axis    = factor(panel_c_axis[as.character(
           setNames(panel_c_levels, panel_c_labels[panel_c_levels])[
             as.character(stratum)])],
           levels = c("Overall", "Age", "Frequency")))

pC <- ggplot(panel_c, aes(x = stratum, y = port_mean, fill = axis)) +
  geom_col(width = 0.65, color = "white", linewidth = 0.3) +
  geom_errorbar(aes(ymin = pmax(port_mean - port_se, 0),
                    ymax = port_mean + port_se),
                width = 0.18, linewidth = 0.4, alpha = 0.7,
                show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10)),
                     limits = c(0, max(panel_c$port_mean + panel_c$port_se) * 1.15),
                     labels = scales::label_number(accuracy = 0.01)) +
  scale_fill_manual(name = NULL, values = panel_c_palette) +
  labs(x = NULL,
       y = expression(Portability == frac(R[A %->% B]^2, R[A %->% A]^2)),
       tag = "C",
       subtitle = "ARG reproduces the frequency baseline and adds age on top") +
  theme(plot.subtitle = element_text(size = 7.5, color = "gray25",
                                     margin = margin(b = 4)),
        axis.text.x = element_text(size = 7, lineheight = 1.0),
        legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.32, "cm"))


## ---- Panel D — the cap: within-R² vs cross-R² ----
make_fst_lab <- function(x) paste0("F[ST]==", sprintf("%.2f", x))
dD$fst_lab <- factor(make_fst_lab(dD$fst),
                     levels = make_fst_lab(sort(unique(dD$fst))))

## Build segments: standard (h²=0.5) → high-power (h²=0.8) per F_ST.
dD_seg <- dD %>%
  select(fst, fst_lab, regime, within_R2, cross_R2) %>%
  tidyr::pivot_wider(names_from = regime,
                     values_from = c(within_R2, cross_R2))

panel_d_palette <- setNames(c("#4E79A7", "#E15759"),
                            make_fst_lab(c(0.15, 0.30)))

pD <- ggplot(dD, aes(x = within_R2, y = cross_R2,
                     color = fst_lab, shape = regime)) +
  geom_segment(data = dD_seg,
               aes(x = within_R2_standard_h2_0.5,
                   xend = within_R2_highpower_h2_0.8,
                   y = cross_R2_standard_h2_0.5,
                   yend = cross_R2_highpower_h2_0.8,
                   color = fst_lab),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               linewidth = 0.55, alpha = 0.8,
               inherit.aes = FALSE) +
  geom_point(size = 2.6, stroke = 0.5, fill = "white") +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted",
              color = "gray70", linewidth = 0.3) +
  scale_x_continuous(limits = c(0, 0.55),
                     breaks = seq(0, 0.5, 0.1),
                     expand = expansion(mult = c(0, 0.04))) +
  scale_y_continuous(limits = c(0, 0.10),
                     breaks = seq(0, 0.10, 0.02),
                     expand = expansion(mult = c(0, 0.04))) +
  scale_color_manual(name = NULL,
                     values = panel_d_palette,
                     labels = function(x) parse(text = x)) +
  scale_shape_manual(name = NULL,
                     values = c(standard_h2_0.5 = 16,
                                highpower_h2_0.8 = 17),
                     labels = c(standard_h2_0.5 = "Standard (h² = 0.5, n = 1000)",
                                highpower_h2_0.8 = "High power (h² = 0.8, n = 2800)")) +
  labs(x = expression(Within * "-ancestry " * R^2),
       y = expression(Cross * "-ancestry " * R^2),
       tag = "D",
       subtitle = "Within climbs to ~0.5; cross stays pinned at ~0.02") +
  theme(plot.subtitle = element_text(size = 7.5, color = "gray25",
                                     margin = margin(b = 4)),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        legend.key.size = unit(0.32, "cm"),
        legend.spacing.y = unit(0, "cm"),
        legend.text = element_text(size = 6.6))


## ---- Assemble ----
fig <- (pA + pB) / (pC + pD) +
  plot_annotation(
    title = "Cross-ancestry PRS failure has an allele-age signature",
    theme = theme(plot.title = element_text(size = 11, face = "bold",
                                            margin = margin(b = 5)))
  )

save_fig(fig, "fig_xancestry_age", width = 7.3, height = 6.5)

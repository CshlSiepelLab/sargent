#!/usr/bin/env Rscript
## fig3_supp.R — Supplementary to Figure 3: cross-ancestry without
## target injection.
##
## Single panel: F_ST ladder at typed=0.70, mix_B=0 (no target-population
## individuals in training).  Same 5 arms as Fig 3 Panel B.  Story: in
## the *no-injection* regime, cross-ancestry prediction is dismal for
## every arm; SARGENT-U remains best but the absolute R² collapse is
## severe; SARGENT-W now consistently below SARGENT-U (the "partial
## rescue" verdict from the rescue test).
##
## Companion to main Fig 3 (which uses mix_B = 0.20, the predictive
## regime where SARGENT-W ≈ SARGENT-U).

args_all   <- commandArgs(trailingOnly = FALSE)
file_arg   <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

csv_path <- file.path(script_dir, "..", "block1_1_v2_results.csv")
df <- read.csv(csv_path, check.names = FALSE)
variants_csv <- file.path(script_dir, "..", "block1_1_v2_variants_results.csv")
if (file.exists(variants_csv)) {
  vdf <- read.csv(variants_csv, check.names = FALSE)
  df <- merge(df, vdf[, c("fst","mix","typed","seed","sargentW_r2c")],
              by = c("fst","mix","typed","seed"))
}

arms_list <- list(
  data.frame(df, arm = "bayesR",   r2c = df$bayesR_r2c),
  data.frame(df, arm = "gblup",    r2c = df$gblup_r2c),
  data.frame(df, arm = "argunif",  r2c = df$argunif_r2c),
  data.frame(df, arm = "sargent",  r2c = df$sargent_r2c)
)
if ("sargentW_r2c" %in% names(df)) {
  arms_list[[5]] <- data.frame(df, arm = "sargentW", r2c = df$sargentW_r2c)
}
arms_df <- bind_rows(arms_list)

arm_labels_ext  <- c(arm_labels,  sargentW = "SARGENT-W")
arm_palette_ext <- c(arm_palette, sargentW = "#76B7B2")
arm_shapes_ext  <- c(arm_shapes,  sargentW = 18)
arm_order_ext   <- c(arm_order,  "sargentW")
arms_df$arm_lab <- factor(arms_df$arm, levels = arm_order_ext,
                          labels = arm_labels_ext[arm_order_ext])

panel <- arms_df %>%
  filter(abs(typed - 0.70) < 1e-6, abs(mix - 0.0) < 1e-6) %>%
  group_by(arm_lab, fst) %>%
  summarise(mean = mean(r2c, na.rm = TRUE),
            se   = if (n() > 1) sd(r2c, na.rm = TRUE) / sqrt(n()) else 0,
            .groups = "drop")

b_arm_order <- c("argunif", "sargentW", "bayesR", "gblup", "sargent")
b_dodge <- setNames(seq(-0.006, 0.006, length.out = length(b_arm_order)),
                    arm_labels_ext[b_arm_order])
panel$fst_jit <- panel$fst + b_dodge[as.character(panel$arm_lab)]

p <- ggplot(panel, aes(x = fst, y = mean, color = arm_lab, shape = arm_lab)) +
  geom_errorbar(aes(x = fst_jit, ymin = mean - se, ymax = mean + se),
                width = 0.006, linewidth = 0.45, alpha = 0.55) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.4, stroke = 0.4, fill = "white") +
  scale_color_manual(name = NULL,
                     values = setNames(arm_palette_ext, arm_labels_ext[names(arm_palette_ext)])) +
  scale_shape_manual(name = NULL,
                     values = setNames(arm_shapes_ext, arm_labels_ext[names(arm_shapes_ext)])) +
  scale_x_continuous(breaks = sort(unique(panel$fst))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     limits = c(0, 0.20),
                     breaks = seq(0, 0.20, by = 0.05)) +
  labs(x = expression(F[ST]),
       y = expression(r^2),
       title = expression(Without~target~injection*":"~mix[B]==0)) +
  theme(legend.position = c(0.98, 0.98),
        legend.justification = c(1, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
        plot.title = element_text(size = 9, face = "bold"))

save_fig(p, "fig3_supp", width = 4.0, height = 3.3)

#!/usr/bin/env Rscript
## fig3.R — Cross-ancestry portability under REALISTIC demographies.
##
## Block 1.1 v3 data: stdpopsim Gravel OOA + Browning AmericanAdmixture
## + Ragsdale archaic.  Supersedes the earlier v2 toy demography figure
## (history preserved in git).
##
## Two panels, 4 arms (BayesR / GBLUP / SARGENT-U / SARGENT):
##   A: Target-injection effect at OOA CEU→YRI, typed=0.70.  Grouped bars
##      over mix_B ∈ {0.00, 0.05}.  Shows the target-injection lift, with
##      SARGENT-U opening the largest gap at mix=0.05.
##   B: Cross-demography comparison at mix=0.05, typed=0.70.  Grouped bars
##      across 4 (model × pop-pair) scenarios sorted by F_ST (ADM4B18
##      EUR→ADMIX F_ST≈0.022, OOA3G09 CEU→CHB ≈0.068, CEU→YRI ≈0.076,
##      YRI→CHB ≈0.096).  Continuous F_ST axis was tried first and
##      rejected — the OOA pop-pairs cluster within ~0.03 of each other,
##      collapsing the line plot.  Story: SARGENT-U > typed-SNP methods
##      across every demography; SARGENT-full (with β) sits below
##      SARGENT-U — direct evidence that β is the cross-ancestry
##      liability (within-method gap captures the β contribution; the
##      v3 sweep did not include the explicit SARGENT-W rescue arm).
##
## ARC5R19 CEU→Neanderthal (F_ST≈0.44) is excluded — all methods collapse
## there (control, not a discriminating result).  Future supplementary
## figure can carry the archaic ladder + typed-coverage axis.
##
## Inputs:
##   test/eval/block1_1_v3_results.csv
##     (produced by figs_r/preprocess_fig3_v3.py from $HOME/test-sargent/
##      block1_1_v3/*_s*/est.predict.json + case.params.json)
## Output: ~/newwriting/arg-mapping/figs2/fig3.{pdf,png}

args_all  <- commandArgs(trailingOnly = FALSE)
file_arg  <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

csv_path <- file.path(script_dir, "..", "block1_1_v3_results.csv")
df <- read.csv(csv_path, check.names = FALSE)

## Long-format with one row per (arm, cell).
arms_df <- bind_rows(
  data.frame(df, arm = "bayesR",  r2c = df$bayesR_r2c),
  data.frame(df, arm = "gblup",   r2c = df$gblup_r2c),
  data.frame(df, arm = "argunif", r2c = df$argunif_r2c),
  data.frame(df, arm = "sargent", r2c = df$sargent_r2c)
)
arms_df$arm_lab <- arm_factor(arms_df$arm)

pop_pair_lab <- function(model, train_pop, test_pop) {
  sprintf("%s\n%s→%s", model, train_pop, test_pop)
}
arms_df$pair <- pop_pair_lab(arms_df$model, arms_df$train_pop, arms_df$test_pop)

## Drop archaic cells from the main figure (control, not main message).
arms_main <- arms_df %>% filter(model != "ARC5R19")

## ---- Panel A: target injection at CEU→YRI, typed=0.70 ----
panel_a <- arms_main %>%
  filter(model == "OOA3G09",
         train_pop == "CEU", test_pop == "YRI",
         abs(typed - 0.70) < 1e-6,
         mix %in% c(0.00, 0.05)) %>%
  group_by(arm_lab, mix) %>%
  summarise(mean = mean(r2c, na.rm = TRUE),
            se   = if (n() > 1) sd(r2c, na.rm = TRUE) / sqrt(n()) else 0,
            n    = sum(!is.na(r2c)),
            .groups = "drop") %>%
  mutate(mix_lab = factor(mix,
                          levels = c(0.00, 0.05),
                          labels = c("mix[B]==0",
                                     "mix[B]==0.05")))

pA <- ggplot(panel_a, aes(x = mix_lab, y = mean, fill = arm_lab, color = arm_lab)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72,
           color = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se),
                position = position_dodge(width = 0.78), width = 0.18,
                linewidth = 0.45, alpha = 0.55, show.legend = FALSE) +
  scale_x_discrete(labels = scales::parse_format()) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     limits = c(0, 0.30),
                     breaks = seq(0, 0.30, by = 0.10)) +
  scale_color_arm() +
  scale_fill_arm() +
  annotate("text", x = Inf, y = Inf,
           label = "OOA3G09 ~ CEU %->% YRI",
           parse = TRUE, hjust = 1.04, vjust = 1.5,
           size = 2.6, color = "gray25") +
  labs(x = NULL, y = expression(r^2), tag = "A") +
  theme(legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

## ---- Panel B: cross-demography comparison at mix=0.05, typed=0.70 ----
## Grouped bars over 4 (model × pop-pair) scenarios, sorted by F_ST so the
## gradient is preserved on the x-axis without forcing a continuous scale
## (CEU→CHB F_ST=0.068 and CEU→YRI F_ST=0.076 are visually indistinguishable
## under a continuous F_ST axis).  4 arms per group.
panel_b <- arms_main %>%
  filter(abs(typed - 0.70) < 1e-6,
         abs(mix - 0.05) < 1e-6) %>%
  group_by(model, train_pop, test_pop, arm_lab) %>%
  summarise(fst  = mean(fst, na.rm = TRUE),
            mean = mean(r2c, na.rm = TRUE),
            se   = if (n() > 1) sd(r2c, na.rm = TRUE) / sqrt(n()) else 0,
            n    = sum(!is.na(r2c)),
            .groups = "drop") %>%
  mutate(short = case_when(
    model == "ADM4B18" ~ "EUR→ADMIX",
    TRUE              ~ sprintf("%s→%s", train_pop, test_pop)
  ))

## Sort pop-pairs by F_ST so the x-axis preserves the gradient.
pair_order <- panel_b %>%
  group_by(short) %>%
  summarise(fst = mean(fst), .groups = "drop") %>%
  arrange(fst)

x_labels <- setNames(
  sprintf("%s\n(F_ST = %.3f)", pair_order$short, pair_order$fst),
  pair_order$short
)

panel_b$short <- factor(panel_b$short, levels = pair_order$short)

pB <- ggplot(panel_b, aes(x = short, y = mean, fill = arm_lab, color = arm_lab)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72,
           color = "white", linewidth = 0.25) +
  geom_errorbar(aes(ymin = pmax(mean - se, 0), ymax = mean + se),
                position = position_dodge(width = 0.78), width = 0.18,
                linewidth = 0.45, alpha = 0.55, show.legend = FALSE) +
  scale_x_discrete(labels = x_labels) +
  scale_color_arm() +
  scale_fill_arm() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08)),
                     limits = c(0, 0.30),
                     breaks = seq(0, 0.30, by = 0.10)) +
  labs(x = NULL, y = expression(r^2), tag = "B") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 7, lineheight = 1.1))

fig3 <- pA + pB + plot_layout(widths = c(1, 1.05))
save_fig(fig3, "fig3", width = 6.5, height = 3.3)

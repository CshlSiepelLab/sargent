#!/usr/bin/env Rscript
## fig2.R — Within-population: h² recovery + held-out prediction (Fig 2).
##
## Two panels, same 50 cells (n=600, sparse, panmictic, h2_total=0.6,
## L=300kb, 5 cov × 10 seeds):
##   A: ĥ² vs typed-SNP coverage, 4 arms
##   B: r²_calib vs typed-SNP coverage, 4 arms
##
## Story: SARGENT-U is coverage-invariant (flat across the axis);
## BayesR/GBLUP/SARGENT (full) collapse at low coverage.
##
## Output: ~/newwriting/arg-mapping/figs2/fig2.{pdf,png}

args_all  <- commandArgs(trailingOnly = FALSE)
file_arg  <- sub("^--file=", "", args_all[grep("^--file=", args_all)])
script_dir <- if (length(file_arg) > 0) dirname(normalizePath(file_arg)) else getwd()
source(file.path(script_dir, "theme.R"))

suppressMessages(library(jsonlite))

base <- Sys.getenv("FIG2_BASE", path.expand("~/test-sargent/fig2_n600"))

load_cells <- function(base) {
  cells <- list.dirs(base, recursive = FALSE)
  cells <- cells[grepl("asparse_cv.*_dpanmictic_s.*$", cells)]
  out <- list()
  for (cd in cells) {
    m <- regmatches(basename(cd), regexec("asparse_cv([0-9.]+)_dpanmictic_s([0-9]+)$",
                                          basename(cd)))[[1]]
    if (length(m) < 3) next
    cov  <- as.numeric(m[2])
    seed <- as.integer(m[3])
    pjf  <- file.path(cd, "est.predict.json")
    sjf  <- file.path(cd, "est.summary.json")
    if (!file.exists(pjf) || !file.exists(sjf)) next
    p <- jsonlite::fromJSON(pjf)
    s <- jsonlite::fromJSON(sjf)
    out[[length(out) + 1]] <- data.frame(
      cov         = cov,
      seed        = seed,
      bayesR_r2   = p$bayesR$r2_calib,
      gblup_r2    = p$gblup$r2_calib,
      argunif_r2  = p$arg_blup_uniform$r2_calib,
      sargent_r2  = p$sargent$r2_calib,
      bayesR_h2   = s$h2_snp,
      gblup_h2    = p$gblup_h2,
      argunif_h2  = p$arg_blup_uniform_h2,
      sargent_h2  = s$h2
    )
  }
  do.call(rbind, out)
}

df <- load_cells(base)
cat(sprintf("Loaded %d cells from %s\n", nrow(df), base))

## Long form by arm
arms_long <- function(df, metric_prefix) {
  bind_rows(
    data.frame(cov = df$cov, seed = df$seed, arm = "bayesR",
               value = df[[paste0("bayesR_", metric_prefix)]]),
    data.frame(cov = df$cov, seed = df$seed, arm = "gblup",
               value = df[[paste0("gblup_", metric_prefix)]]),
    data.frame(cov = df$cov, seed = df$seed, arm = "argunif",
               value = df[[paste0("argunif_", metric_prefix)]]),
    data.frame(cov = df$cov, seed = df$seed, arm = "sargent",
               value = df[[paste0("sargent_", metric_prefix)]])
  )
}

panel <- function(df, metric_prefix, ylab, title_text, tag,
                  ylim_max = NULL, truth_line = NULL,
                  legend_pos = c(0.02, 0.02), legend_just = c(0, 0)) {
  d <- arms_long(df, metric_prefix)
  d$arm_lab <- arm_factor(d$arm)

  s <- d %>%
    group_by(arm_lab, cov) %>%
    summarise(mean = mean(value, na.rm = TRUE),
              se   = if (n() > 1) sd(value, na.rm = TRUE) / sqrt(n()) else 0,
              .groups = "drop")

  ## Same paired-jitter trick from Fig 3 for error bars only
  arm_dodge_order <- c("argunif", "bayesR", "gblup", "sargent")
  arm_dodge <- setNames(seq(-0.012, 0.012, length.out = length(arm_dodge_order)),
                        arm_labels[arm_dodge_order])
  s$cov_jit <- s$cov + arm_dodge[as.character(s$arm_lab)]

  p <- ggplot(s, aes(x = cov, y = mean, color = arm_lab, shape = arm_lab)) +
    {if (!is.null(truth_line))
      geom_hline(yintercept = truth_line, linetype = "dashed",
                 color = "gray35", linewidth = 0.45)} +
    geom_errorbar(aes(x = cov_jit, ymin = mean - se, ymax = mean + se),
                  width = 0.018, linewidth = 0.45, alpha = 0.55) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 2.4, stroke = 0.4, fill = "white") +
    scale_color_arm() +
    scale_shape_arm() +
    scale_x_continuous(breaks = sort(unique(s$cov)),
                       trans = "reverse",
                       expand = expansion(mult = c(0.04, 0.04))) +
    scale_y_continuous(limits = c(0, ylim_max %||% NA),
                       expand = expansion(mult = c(0, 0.04))) +
    labs(x = expression(Typed-SNP~coverage~ italic(f)[typed]),
         y = ylab, tag = tag) +
    theme(legend.position = legend_pos,
          legend.justification = legend_just,
          legend.background = element_rect(fill = alpha("white", 0.85), color = NA))
  p
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

pA <- panel(df, "h2",
            ylab = expression(hat(italic(h))^2),
            title_text = "h² recovery",
            tag = "A",
            ylim_max = 0.85, truth_line = 0.60,
            legend_pos = c(0.02, 0.02), legend_just = c(0, 0))

pB <- panel(df, "r2",
            ylab = expression(r^2),
            title_text = "Prediction",
            tag = "B",
            ylim_max = 0.65,
            legend_pos = c(0.02, 0.02), legend_just = c(0, 0))

fig2 <- pA + pB + plot_layout(widths = c(1, 1))

save_fig(fig2, "fig2", width = 6.5, height = 3.3)

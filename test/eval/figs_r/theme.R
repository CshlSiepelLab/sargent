## theme.R — shared ggplot theme + palette for SARGENT manuscript figures.
## Modeled on the vine paper's makeGraphs3.R conventions (Tableau-10
## palette, base_size 8, thin gray grids, cairo_pdf output).

suppressMessages({
  library(ggplot2)
  library(scales)
  library(patchwork)
  library(dplyr)
})

## ---- Theme ----
sargent_theme <- function(base_size = 8) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title         = element_text(size = 9, face = "bold"),
      plot.tag           = element_text(size = 11, face = "bold"),
      plot.tag.location  = "plot",
      axis.title         = element_text(size = 9),
      axis.text          = element_text(size = 7),
      legend.title       = element_text(size = 8),
      legend.text        = element_text(size = 7),
      legend.background  = element_rect(fill = "white", color = NA),
      legend.key.size    = unit(0.4, "cm"),
      strip.text         = element_text(size = 8, face = "plain"),
      panel.grid.major   = element_line(color = "gray60", linewidth = 0.3),
      panel.grid.minor   = element_line(color = "gray80", linewidth = 0.2),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA)
    )
}
theme_set(sargent_theme())

## ---- Arm palette ----
## SARGENT (full): purple — the framework's headline arm within-pop.
## SARGENT-U: green     — sister framework arm; the cross-pop representative.
## GBLUP: blue          — typed-SNP infinitesimal baseline.
## BayesR: red          — typed-SNP sparse baseline.
arm_palette <- c(
  sargent  = "#B07AA1",
  argunif  = "#59A14F",
  gblup    = "#4E79A7",
  bayesR   = "#E15759"
)
arm_labels <- c(
  sargent  = "SARGENT",
  argunif  = "SARGENT-U",
  gblup    = "GBLUP",
  bayesR   = "BayesR"
)
arm_shapes <- c(
  sargent  = 16,  # filled circle
  argunif  = 15,  # filled square
  gblup    = 17,  # filled triangle
  bayesR   = 18   # filled diamond
)
arm_order <- c("sargent", "argunif", "gblup", "bayesR")

## Factor + scale helpers
arm_factor <- function(x) factor(x, levels = arm_order, labels = arm_labels[arm_order])

scale_color_arm <- function(...) {
  scale_color_manual(name = NULL, values = setNames(arm_palette, arm_labels[names(arm_palette)]), ...)
}
scale_fill_arm <- function(...) {
  scale_fill_manual(name = NULL, values = setNames(arm_palette, arm_labels[names(arm_palette)]), ...)
}
scale_shape_arm <- function(...) {
  scale_shape_manual(name = NULL, values = setNames(arm_shapes, arm_labels[names(arm_shapes)]), ...)
}

## ---- Save helper ----
## Default fig width = 6.5 in (letter, 1-in margins). Heights chosen per
## figure. Cairo PDF for embedded fonts; matching PNG for review.
save_fig <- function(plot, name, width = 6.5, height = 3.5,
                     out_dir = "~/newwriting/arg-mapping/figs2") {
  out_dir <- path.expand(out_dir)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(out_dir, paste0(name, ".pdf"))
  png_path <- file.path(out_dir, paste0(name, ".png"))
  if (capabilities("cairo")) {
    ggsave(pdf_path, plot = plot, width = width, height = height,
           units = "in", device = cairo_pdf)
  } else {
    ggsave(pdf_path, plot = plot, width = width, height = height,
           units = "in", family = "Helvetica")
  }
  ggsave(png_path, plot = plot, width = width, height = height,
         units = "in", dpi = 300, bg = "white")
  cat(sprintf("Wrote %s\nWrote %s\n", pdf_path, png_path))
}

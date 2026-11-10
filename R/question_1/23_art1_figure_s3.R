# Article 1 - Figure S3: effects of the other patient characteristics on
# dental neglect, one model each.
#
# Run from anywhere inside the project:  Rscript R/question_1/23_art1_figure_s3.R
#
# Forest plot in plain journal style (Arial, white background), as in
# 22_art1_forest_figures.R, on the probability-more-neglected scale. One row
# per characteristic, grouped by section; intervals that include 50% are grey,
# as the paths in Figure 2.
#
# Input: results/question_1/rq1_candidates_effects.csv
#   (14_rq1_candidate_effects.R; nothing is refitted)
# Outputs: results/art1_MIMIC/art1_figure_s3.pdf (fonts embedded) and .png.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

res_dir <- here("results", "question_1")
out_dir <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Figure size: full page width.
width_mm  <- 180
height_mm <- 150

# Colours, as in 22_art1_forest_figures.R; grey as in 18_art1_figure2_path.R.
surface      <- "#ffffff"
ink          <- "#0b0b0b"
ink_soft     <- "#52514e"
ink_mute     <- "#898781"
grid_col     <- "#e6e5df"
not_credible <- "#8f8d85"

theme_fig <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.background = element_rect(fill = surface, colour = NA),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft, size = 7),
    axis.line        = element_line(colour = ink_soft, linewidth = 0.3),
    axis.ticks       = element_line(colour = ink_soft, linewidth = 0.3)
  )

# ---- Figure S3 -------------------------------------------------------------

sections <- c("Chronic disease", "Smoking", "Cancer", "Cancer treatment", "Other")
effects <- read_csv(file.path(res_dir, "rq1_candidates_effects.csv"), show_col_types = FALSE) |>
  filter(scale == "probability more neglected (%)") |>
  mutate(
    section = factor(section, sections),
    # Binary: patients with it / patients where it is defined; continuous: observed values.
    term = if_else(type == "binary", sprintf("%s (%d/%d)", term, n_yes, n_defined), sprintf("%s (n = %d)", term, n_defined)),
    term = factor(term, rev(unique(term))),
    credible = if_else(lower_95 > 50 | upper_95 < 50, "yes", "no"),
    label = sprintf("%.0f [%.0f-%.0f]", median, lower_95, upper_95)
  )
stopifnot(!anyNA(effects$section), nrow(effects) == 24)

# The values sit right of the panel, at x_text; clipping is off so they show.
x_limits <- c(floor(min(effects$lower_95) / 5) * 5, ceiling(max(effects$upper_95) / 5) * 5)
x_text <- x_limits[2] + 4

figure_s3 <- ggplot(effects, aes(y = term, colour = credible)) +
  geom_vline(xintercept = 50, colour = ink_mute, linewidth = 0.3) +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95), linewidth = 0.4) +
  geom_linerange(aes(xmin = lower_50, xmax = upper_50), linewidth = 1.3) +
  geom_point(aes(x = median), fill = surface, shape = 21, size = 1.5, stroke = 0.7) +
  geom_text(aes(x = x_text, label = label), hjust = 0, size = 6.5 / .pt, family = "Arial") +
  scale_colour_manual(values = c(yes = ink, no = not_credible), guide = "none") +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 90), expand = c(0, 0)) +
  # Half a row above and below each section, so the one-row "Other" panel has height.
  scale_y_discrete(expand = expansion(add = 0.6)) +
  facet_wrap(~section, ncol = 1, scales = "free_y", space = "free_y") +
  coord_cartesian(xlim = x_limits, clip = "off") +
  labs(x = "Probability more neglected (%)", y = NULL) +
  theme_fig +
  theme(
    panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.5, face = "bold", hjust = 0),
    strip.clip = "off",
    panel.spacing.y = unit(3, "mm"),
    plot.margin = margin(4, 22, 4, 4, "mm")
  )

ggsave(file.path(out_dir, "art1_figure_s3.pdf"), figure_s3, width = width_mm, height = height_mm, units = "mm",
  device = cairo_pdf, bg = surface)
ggsave(file.path(out_dir, "art1_figure_s3.png"), figure_s3, width = width_mm, height = height_mm, units = "mm",
  dpi = 600, device = ragg::agg_png, bg = surface)
cat("Written:", file.path(out_dir, "art1_figure_s3.pdf"), "and .png\n")

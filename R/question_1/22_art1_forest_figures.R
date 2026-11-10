# Article 1 - Figure 3: effects on dental neglect, adjusted and adjusted for
# visit year only.
#
# Run from anywhere inside the project:  Rscript R/question_1/22_art1_forest_figures.R
#
# Forest plot in plain journal style (Arial, white background), as in
# 17_art1_figure_s2.R, on the probability-more-neglected scale. Top: the four
# risk type x timing contrasts from the reported model and from the model
# adjusted for visit year only. Bottom: the covariates of the reported model
# (visit year from both).
#
# Input: results/question_1/rq1_unadjusted_comparison.csv
#   (20_rq1_neglect_unadjusted.R; nothing is refitted)
# Outputs: results/art1_MIMIC/art1_figure3.pdf (fonts embedded) and .png.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(ggplot2)
})

res_dir <- here("results", "question_1")
out_dir <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Figure size: 1.5-column width.
width_mm  <- 140
height_mm <- 85

# Colours, as in 17_art1_figure_s2.R; the reported model is ink, the model
# adjusted for visit year only is slot 3 of the dataviz reference palette.
surface  <- "#ffffff"
ink      <- "#0b0b0b"
ink_soft <- "#52514e"
ink_mute <- "#898781"
grid_col <- "#e6e5df"
model_cols <- c(reported = ink, unadjusted = "#eb6834")
model_names <- c(reported = "Adjusted (reported model)", unadjusted = "Adjusted for visit year only")

theme_fig <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.background = element_rect(fill = surface, colour = NA),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft, size = 7),
    axis.line        = element_line(colour = ink_soft, linewidth = 0.3),
    axis.ticks       = element_line(colour = ink_soft, linewidth = 0.3),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(3.5, "mm")
  )

# ---- Figure 3 --------------------------------------------------------------

term_order <- c(
  "RTH vs BF, seen before treatment", "RTH vs BF, seen after treatment",
  "After vs before, BF", "After vs before, RTH",
  "Age, per 10 years", "Male vs female", "Former vs never smoker", "Current vs never smoker",
  "Seen in 2025 vs 2024", "Seen in 2026 vs 2024"
)
effects <- read_csv(file.path(res_dir, "rq1_unadjusted_comparison.csv"), show_col_types = FALSE) |>
  filter(scale == "probability more neglected (%)", section %in% c("Risk type × timing", "Covariates")) |>
  mutate(
    section = factor(section, c("Risk type × timing", "Covariates")),
    term = factor(term, rev(term_order)),
    model = factor(model, names(model_cols)),
    label = sprintf("%.0f [%.0f-%.0f]", median, lower_95, upper_95)
  )
stopifnot(!anyNA(effects$term), nrow(effects) == 2 * 4 + 6 + 2)

# The values sit right of the panel, at x_text; clipping is off so they show.
x_limits <- c(floor(min(effects$lower_95) / 5) * 5, ceiling(max(effects$upper_95) / 5) * 5)
x_text <- x_limits[2] + 4
dodge <- position_dodge(width = 0.7, orientation = "y", reverse = TRUE)

figure3 <- ggplot(effects, aes(y = term, colour = model, group = model)) +
  geom_vline(xintercept = 50, colour = ink_mute, linewidth = 0.3) +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95), linewidth = 0.4, position = dodge) +
  geom_linerange(aes(xmin = lower_50, xmax = upper_50), linewidth = 1.3, position = dodge) +
  geom_point(aes(x = median), fill = surface, shape = 21, size = 1.5, stroke = 0.7, position = dodge) +
  geom_text(aes(x = x_text, label = label), hjust = 0, size = 6.5 / .pt, family = "Arial",
    position = dodge, show.legend = FALSE) +
  scale_colour_manual(values = model_cols, labels = model_names, name = NULL) +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 90)) +
  facet_wrap(~section, ncol = 1, scales = "free_y", space = "free_y") +
  coord_cartesian(xlim = x_limits, expand = FALSE, clip = "off") +
  labs(x = "Probability more neglected (%)", y = NULL) +
  theme_fig +
  theme(
    panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.5, face = "bold", hjust = 0),
    strip.clip = "off",
    panel.spacing.y = unit(3, "mm"),
    plot.margin = margin(4, 22, 4, 4, "mm"),
    legend.position = "bottom",
    legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0)
  )

ggsave(file.path(out_dir, "art1_figure3.pdf"), figure3, width = width_mm, height = height_mm, units = "mm",
  device = cairo_pdf, bg = surface)
ggsave(file.path(out_dir, "art1_figure3.png"), figure3, width = width_mm, height = height_mm, units = "mm",
  dpi = 600, device = ragg::agg_png, bg = surface)
cat("Written:", file.path(out_dir, "art1_figure3.pdf"), "and .png\n")

# Article 1 - Figure 1: selection of the indicators of dental neglect.
#
# Run from anywhere inside the project:  Rscript R/question_1/19_art1_figure1.R
#
# Three panels, in plain journal style (Arial, white background), as in
# 17_art1_figure_s2.R:
#   A  correlations of the 7 candidate indicators, with patients per pair;
#   B  scree plots against parallel analysis, 7 candidates and 4 retained;
#   C  one-factor loadings: EFA of 7, EFA of 4, CFA of 4 with the
#      visit-X-ray residual correlation.
#
# Inputs, all written by 04_rq1_efa.R (nothing is refitted):
#   results/question_1/rq1_efa_correlations.csv
#   results/question_1/rq1_efa_scree.csv
#   results/question_1/rq1_efa_loadings.csv
# Outputs: results/art1_MIMIC/art1_figure1.pdf (fonts embedded) and .png.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

res_dir <- here("results", "question_1")
out_dir <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Figure size: double-column width.
width_mm  <- 180
height_mm <- 120

# Colours, as in 17_art1_figure_s2.R. Correlations run from blue (negative)
# through white to red (positive); the three solutions in C use slots 1-3 of
# the dataviz reference palette.
surface  <- "#ffffff"
ink      <- "#0b0b0b"
ink_soft <- "#52514e"
ink_mute <- "#898781"
grid_col <- "#e6e5df"
cor_low  <- "#2a78d6"
cor_high <- "#e34948"
solution_cols <- c("#2a78d6", "#eb6834", "#199e70")

theme_fig <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.background = element_rect(fill = surface, colour = NA),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft, size = 7),
    axis.line        = element_line(colour = ink_soft, linewidth = 0.3),
    axis.ticks       = element_line(colour = ink_soft, linewidth = 0.3),
    plot.title       = element_text(size = 8, face = "bold", margin = margin(b = 2)),
    plot.tag         = element_text(size = 10, face = "bold"),
    plot.title.position = "plot",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(3.5, "mm")
  )

# Indicator labels as in Figure 2, the 4 retained first. Every indicator is
# coded so that higher = more neglect.
labels <- c(
  hygiene_poor = "Poor hygiene",
  last_visit   = "Years since dentist",
  xray_old     = "Old / no X-ray",
  extract_log  = "Teeth to extract",
  pain         = "Visit with pain",
  treat_log    = "Teeth to treat",
  cbct         = "CBCT needed"
)
n_core <- 4

# ---- A: correlations of the 7 candidates -----------------------------------

# Lower triangle: columns are the first of each pair, rows the second.
cors <- read_csv(file.path(res_dir, "rq1_efa_correlations.csv"), show_col_types = FALSE) |>
  mutate(
    var1 = factor(labels[var1], labels[-length(labels)]),
    var2 = factor(labels[var2], rev(labels[-1]))
  )
stopifnot(nrow(cors) == 21, !anyNA(cors$var1), !anyNA(cors$var2))

# The dashed line separates the pairs of retained indicators (above) from the
# pairs that involve a dropped candidate (below).
sep_y <- length(labels) - n_core + 0.5

panel_a <- ggplot(cors, aes(var1, var2)) +
  geom_tile(aes(fill = r), colour = surface, linewidth = 0.6) +
  geom_text(aes(label = sub("-", "−", sprintf("%.2f", r))), size = 7 / .pt, family = "Arial",
    colour = ink, nudge_y = 0.13) +
  geom_text(aes(label = sprintf("n = %d", n)), size = 5.5 / .pt, family = "Arial",
    colour = ink_soft, nudge_y = -0.2) +
  annotate("segment", x = 0.5, xend = length(labels) - 0.5, y = sep_y, yend = sep_y,
    colour = ink_soft, linewidth = 0.35, linetype = "22") +
  scale_fill_gradient2(low = cor_low, mid = surface, high = cor_high, limits = c(-1, 1),
    breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("−1", "−0.5", "0", "0.5", "1"), name = "Correlation") +
  scale_x_discrete(position = "bottom", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(tag = "A", title = "Correlations of the 7 candidate indicators", x = NULL, y = NULL) +
  theme_fig +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 40, hjust = 1),
    legend.position = "bottom",
    legend.title = element_text(size = 7, vjust = 0.85),
    legend.key.width = unit(7, "mm"),
    legend.key.height = unit(2.2, "mm")
  )

# ---- B: scree plots with parallel analysis ---------------------------------

set_names <- c(
  "EFA, 7 candidate indicators" = "7 candidates",
  "EFA, 4 core indicators" = "4 retained indicators"
)
scree <- read_csv(file.path(res_dir, "rq1_efa_scree.csv"), show_col_types = FALSE) |>
  mutate(
    set = factor(set_names[set], set_names),
    series = recode(series, "Observed data" = "Observed", "Simulated random data" = "Random data (parallel analysis)")
  )
stopifnot(!anyNA(scree$set))

panel_b <- ggplot(scree, aes(factor, eigenvalue, linetype = series, shape = series)) +
  geom_hline(yintercept = 0, colour = ink_mute, linewidth = 0.3) +
  geom_line(colour = ink, linewidth = 0.4) +
  geom_point(colour = ink, fill = surface, size = 1.3, stroke = 0.5) +
  scale_linetype_manual(values = c("solid", "22"), name = NULL) +
  scale_shape_manual(values = c(16, 21), name = NULL) +
  scale_x_continuous(breaks = 1:7) +
  facet_wrap(~set, scales = "free_x") +
  labs(tag = "B", title = "Scree plots with parallel analysis", x = "Factor", y = "Eigenvalue") +
  theme_fig +
  theme(
    panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.25),
    strip.background = element_blank(),
    strip.text = element_text(size = 7, hjust = 0),
    legend.position = "bottom",
    legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.width = unit(6, "mm")
  )

# ---- C: one-factor loadings ------------------------------------------------

solution_names <- c(
  loading_all7 = "EFA, 7 candidates",
  loading_core = "EFA, 4 retained",
  loading_cfa  = "CFA, 4 retained with residual correlation"
)
loadings <- read_csv(file.path(res_dir, "rq1_efa_loadings.csv"), show_col_types = FALSE) |>
  select(indicator, all_of(names(solution_names))) |>
  pivot_longer(-indicator, names_to = "solution", values_to = "loading") |>
  filter(!is.na(loading)) |>
  mutate(
    indicator = factor(labels[indicator], rev(labels)),
    solution = factor(solution_names[solution], solution_names)
  )
stopifnot(nrow(loadings) == 7 + 4 + 4, !anyNA(loadings$indicator))

dodge <- position_dodge(width = 0.6, reverse = TRUE)
panel_c <- ggplot(loadings, aes(loading, indicator, colour = solution, shape = solution)) +
  geom_vline(xintercept = 0, colour = ink_mute, linewidth = 0.3) +
  geom_vline(xintercept = 0.4, colour = ink_mute, linewidth = 0.3, linetype = "22") +
  geom_point(size = 1.5, stroke = 0.5, fill = surface, position = dodge) +
  scale_colour_manual(values = solution_cols, name = NULL) +
  scale_shape_manual(values = c(21, 22, 24), name = NULL) +
  scale_x_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.4, 1),
    labels = c("−1", "−0.5", "0", "0.4", "1")) +
  guides(colour = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
  labs(tag = "C", title = "One-factor loadings", x = "Standardised loading", y = NULL) +
  theme_fig +
  theme(
    panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.25),
    legend.position = "bottom",
    legend.justification = "left",
    legend.margin = margin(0, 0, 0, 0),
    legend.key.spacing.y = unit(0.3, "mm")
  )

# ---- Assemble --------------------------------------------------------------

figure <- wrap_elements(full = panel_a) + (panel_b / panel_c + plot_layout(heights = c(1, 1.25))) +
  plot_layout(widths = c(1.05, 1))

ggsave(file.path(out_dir, "art1_figure1.pdf"), figure, width = width_mm, height = height_mm, units = "mm",
  device = cairo_pdf, bg = surface)
ggsave(file.path(out_dir, "art1_figure1.png"), figure, width = width_mm, height = height_mm, units = "mm",
  dpi = 600, device = ragg::agg_png, bg = surface)
cat("Written:", file.path(out_dir, "art1_figure1.pdf"), "and .png\n")

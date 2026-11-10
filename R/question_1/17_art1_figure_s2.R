# Article 1 - Figure S2: model checks of the reported MIMIC model.
#
# Run from anywhere inside the project:  Rscript R/question_1/17_art1_figure_s2.R
#
# Four panels, in plain journal style (Arial, white background):
#   A  posterior predictive check of each indicator: observed category counts
#      vs the median and 95% interval of the replicated counts;
#   B  posterior predictive check of the six pairwise Spearman correlations;
#   C  effects on neglect, reported model vs the model without hygiene;
#   D  effects on neglect under narrower, reported and wider priors.
# C and D share the effect axis, on the probability-more-neglected scale.
#
# Inputs, all written by earlier scripts (nothing is refitted):
#   results/question_1/rq1_neglect_fit.rds              05_rq1_neglect.R
#   results/question_1/rq1_neglect_ppc_correlations.csv 05_rq1_neglect.R
#   results/question_1/rq1_no_hygiene_comparison.csv    12_rq1_neglect_no_hygiene.R
#   results/question_1/rq1_prior_sens_comparison.csv    15_rq1_prior_sensitivity.R
# Outputs: results/art1_MIMIC/art1_figure_s2.pdf (fonts embedded) and .png.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(rstan)
})

in_file  <- here("data", "processed", "clean_data.rds")
res_dir  <- here("results", "question_1")
out_dir  <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Figure size: one full journal page, double-column width.
width_mm  <- 180
height_mm <- 225

# Colours. Series hues are slots 1-3 of the dataviz reference palette (aqua
# at its darker step for >= 3:1 contrast on white), validated as a set on
# white; the reported model is neutral grey in both C and D.
surface   <- "#ffffff"
ink       <- "#0b0b0b"
ink_soft  <- "#52514e"
ink_mute  <- "#898781"
grid_col  <- "#e6e5df"
obs_fill  <- "#d3d2cc"
model_cols <- c(
  reported   = ink_soft,
  no_hygiene = "#199e70",
  narrower   = "#2a78d6",
  wider      = "#eb6834"
)
model_names <- c(
  reported   = "Reported model",
  no_hygiene = "Without hygiene",
  narrower   = "Narrower priors (SDs halved)",
  wider      = "Wider priors (SDs doubled)"
)

theme_fig <- theme_classic(base_size = 8, base_family = "Arial") +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.background = element_rect(fill = surface, colour = NA),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft, size = 7),
    axis.line        = element_line(colour = ink_soft, linewidth = 0.3),
    axis.ticks       = element_line(colour = ink_soft, linewidth = 0.3),
    plot.title       = element_text(size = 8, face = "bold", margin = margin(b = 2)),
    plot.subtitle    = element_text(size = 7, colour = ink_soft, margin = margin(b = 4)),
    plot.tag         = element_text(size = 10, face = "bold"),
    plot.title.position = "plot",
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(3.5, "mm")
  )

# ---- Data (as in 05_rq1_neglect.R) -----------------------------------------

md <- readRDS(in_file) |>
  mutate(n_groups = group_bf_before + group_bf_after + group_rth_before + group_rth_after + group_other) |>
  filter(n_groups == 1, group_other == 0) |>
  transmute(
    extract = as.integer(teeth_to_extract),
    hygiene = 4L - as.integer(oral_hygiene),   # 1 good .. 3 poor
    visit = as.integer(last_dentist_visit),    # 1 <1 year .. 4 never
    xray = as.integer(dental_xray),            # 1 current .. 3 never
    age_at_visit, sex, smoking_status, visit_year
  ) |>
  drop_na(extract, hygiene, age_at_visit, sex, smoking_status, visit_year)
stopifnot(nrow(md) == 175)

fit <- readRDS(file.path(res_dir, "rq1_neglect_fit.rds"))

# ---- A: indicators, observed vs replicated counts --------------------------

indicators <- list(
  extract = list(title = "Teeth to extract", levels = 0:25, labels = NULL, x = "Number of teeth"),
  hygiene = list(title = "Oral hygiene", levels = 1:3, labels = c("Good", "Average", "Poor"), x = NULL),
  visit   = list(title = "Last dentist visit", levels = 1:4, labels = c("<1 y", "1-5 y", "5-10 y", "Never"), x = NULL),
  xray    = list(title = "Dental X-ray", levels = 1:3, labels = c("Current", ">5 y", "Never"), x = NULL)
)

# Replicated indicators come from new draws of both latent terms (see the
# generated quantities of neglect_mimic.stan), over the patients with the
# indicator observed.
ppc_counts <- imap(indicators, \(ind, name) {
  ok <- !is.na(md[[name]])
  rep <- as.matrix(fit, pars = paste0(name, "_rep"))[, ok]
  map(ind$levels, \(k) {
    r <- rowSums(rep == k)
    tibble(
      level = k, observed = sum(md[[name]][ok] == k), median = median(r),
      lower_95 = quantile(r, 0.025, names = FALSE), upper_95 = quantile(r, 0.975, names = FALSE)
    )
  }) |>
    list_rbind() |>
    mutate(indicator = name, n = sum(ok))
})

inside <- list_rbind(ppc_counts) |> mutate(inside = observed >= lower_95 & observed <= upper_95)
cat(sprintf("Category counts inside their 95%% replicated interval: %d of %d\n", sum(inside$inside), nrow(inside)))

ppc_panel <- function(name) {
  ind <- indicators[[name]]
  d <- ppc_counts[[name]]
  p <- ggplot(d, aes(x = level)) +
    geom_col(aes(y = observed), fill = obs_fill, width = 0.8) +
    geom_linerange(aes(ymin = lower_95, ymax = upper_95), colour = ink, linewidth = 0.4) +
    geom_point(aes(y = median), colour = ink, size = 0.9) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = sprintf("%s (n = %d)", ind$title, d$n[1]), x = ind$x, y = "Patients") +
    theme_fig +
    theme(panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.25))
  if (is.null(ind$labels)) {
    p + scale_x_continuous(breaks = seq(0, 25, 5), expand = expansion(add = 0.6))
  } else {
    p + scale_x_continuous(breaks = ind$levels, labels = ind$labels)
  }
}
panel_a <- (ppc_panel("extract") + labs(tag = "A")) + ppc_panel("hygiene") + ppc_panel("visit") + ppc_panel("xray") +
  plot_layout(ncol = 2)

# ---- B: correlations between indicators ------------------------------------

short_names <- c(
  "Teeth to extract" = "Teeth to extract", "Poor hygiene" = "hygiene",
  "Years since dentist" = "last visit", "Old / no X-ray" = "X-ray"
)
cor_ppc <- read_csv(file.path(res_dir, "rq1_neglect_ppc_correlations.csv"), show_col_types = FALSE) |>
  filter(model == "correlated") |>
  separate_wider_delim(pair, " - ", names = c("first", "second")) |>
  mutate(
    first = short_names[first], second = short_names[second],
    first = if_else(first == "Teeth to extract", "Extract", paste0(toupper(substr(first, 1, 1)), substring(first, 2))),
    pair = sprintf("%s - %s", first, second),
    pair = factor(pair, rev(pair))
  )
stopifnot(nrow(cor_ppc) == 6, !anyNA(cor_ppc$pair))

panel_b <- ggplot(cor_ppc, aes(y = pair)) +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95), colour = ink, linewidth = 0.4) +
  geom_point(aes(x = median), colour = ink, size = 0.9) +
  geom_point(aes(x = observed), shape = 23, fill = obs_fill, colour = ink, size = 2, stroke = 0.3) +
  scale_x_continuous(limits = c(0, 0.8), breaks = seq(0, 0.8, 0.2)) +
  labs(tag = "B", title = "Correlations between indicators", x = "Spearman correlation", y = NULL) +
  theme_fig +
  theme(panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25))

# ---- C, D: effects on neglect under the model variants ---------------------

# The four risk type x timing contrasts, the one group contrast the text of
# Methods S1 cites (RTH seen after vs BF seen before), and the covariates.
read_effects <- function(file, model_col) {
  read_csv(file.path(res_dir, file), show_col_types = FALSE) |>
    rename(model = all_of(model_col)) |>
    filter(
      scale == "probability more neglected (%)",
      section != "Group vs BF, seen before" | term == "RTH, seen after"
    ) |>
    mutate(
      term = recode(term, "RTH, seen after" = "RTH after vs BF before"),
      section = if_else(section == "Group vs BF, seen before", "Risk type × timing", section)
    ) |>
    select(model, section, term, median, lower_50, upper_50, lower_95, upper_95)
}
no_hygiene <- read_effects("rq1_no_hygiene_comparison.csv", "model")
priors     <- read_effects("rq1_prior_sens_comparison.csv", "priors")

# Both files carry the reported model; they must agree.
stopifnot(isTRUE(all.equal(
  arrange(filter(no_hygiene, model == "reported"), term),
  arrange(filter(priors, model == "reported"), term)
)))
term_order <- unique(filter(priors, model == "reported")$term)
stopifnot(length(term_order) == 11)

x_limits <- range(c(no_hygiene$lower_95, priors$lower_95, no_hygiene$upper_95, priors$upper_95))
x_limits <- c(floor(x_limits[1] / 5) * 5, ceiling(x_limits[2] / 5) * 5)

effects_panel <- function(d, models, tag, title) {
  d <- d |>
    mutate(
      model = factor(model, models),
      term = factor(term, rev(term_order)),
      section = factor(section, c("Risk type × timing", "Covariates"))
    )
  dodge <- position_dodge(width = 0.65, orientation = "y", reverse = TRUE)
  ggplot(d, aes(y = term, colour = model, group = model)) +
    geom_vline(xintercept = 50, colour = ink_mute, linewidth = 0.3) +
    geom_linerange(aes(xmin = lower_95, xmax = upper_95), linewidth = 0.35, position = dodge) +
    geom_linerange(aes(xmin = lower_50, xmax = upper_50), linewidth = 1.1, position = dodge) +
    geom_point(aes(x = median), fill = surface, shape = 21, size = 1.3, stroke = 0.6, position = dodge) +
    scale_colour_manual(values = model_cols, labels = model_names, limits = names(model_cols), drop = FALSE, name = NULL) +
    scale_x_continuous(limits = x_limits, breaks = c(10, 25, 50, 75, 90)) +
    facet_grid(section ~ ., scales = "free_y", space = "free_y") +
    labs(tag = tag, title = title, x = "Probability more neglected (%)", y = NULL) +
    theme_fig +
    theme(
      panel.grid.major.x = element_line(colour = grid_col, linewidth = 0.25),
      strip.background = element_blank(),
      strip.text.y = element_text(size = 7, face = "bold"),
      panel.spacing.y = unit(3, "mm"),
      legend.position = "bottom"
    )
}
panel_c <- effects_panel(no_hygiene, c("reported", "no_hygiene"), "C", "Without hygiene") +
  theme(strip.text.y = element_blank(), legend.position = "none")
panel_d <- effects_panel(priors, c("narrower", "reported", "wider"), "D", "Prior sensitivity") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")

# One legend for C and D, taken from a panel that holds all four series:
# patchwork cannot merge the two panels' legends, which differ in their keys.
# cowplot measures the legend on its own null device; ragg's knows Arial.
cowplot::set_null_device("agg")
shared_legend <- cowplot::get_legend(
  effects_panel(bind_rows(no_hygiene, filter(priors, model != "reported")), names(model_cols), NULL, NULL)
)

# ---- Assemble --------------------------------------------------------------

# The top and bottom halves are wrapped so that panel A does not line up its
# axes with the long effect labels of panel C.
top <- wrap_elements(full = panel_a) + wrap_elements(full = panel_b) + plot_layout(widths = c(2, 1))
bottom <- wrap_elements(full = (panel_c | panel_d) / wrap_elements(shared_legend) + plot_layout(heights = c(1, 0.05)))
figure <- top / bottom + plot_layout(heights = c(1, 1.1))

ggsave(file.path(out_dir, "art1_figure_s2.pdf"), figure, width = width_mm, height = height_mm, units = "mm",
  device = cairo_pdf, bg = surface)
ggsave(file.path(out_dir, "art1_figure_s2.png"), figure, width = width_mm, height = height_mm, units = "mm",
  dpi = 600, device = ragg::agg_png, bg = surface)
cat("Written:", file.path(out_dir, "art1_figure_s2.pdf"), "and .png\n")

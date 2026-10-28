# RQ1 - sensitivity analysis: the latent neglect model without poor hygiene.
#
# Run from anywhere inside the project:  Rscript R/question_1/12_rq1_neglect_no_hygiene.R
# after R/question_1/05_rq1_neglect.R, whose coefficients it compares against.
#
# Hygiene has by far the largest loading in the reported model, so its
# effects could be the effects on hygiene alone. This refits the reported
# model (same patients, covariates, priors, visit-X-ray residual correlation
# and sampler) with hygiene left out of the likelihood, so neglect is what
# teeth to extract, last dentist visit and X-ray share. With three
# indicators the residual correlation is separated from neglect only through
# the covariate effects, so its posterior is reported next to its prior.
#
# Outputs go to results/question_1/ with the prefix rq1_no_hygiene_.
# reports/question_1/rq1_neglect.qmd documents the comparison.

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
  library(posterior)
  library(bayesplot)
})

in_file   <- here("data", "processed", "clean_data.rds")
stan_file <- here("R", "question_1", "stan", "neglect_mimic.stan")
out_dir   <- here("results", "question_1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Sampler settings, as in 05_rq1_neglect.R
chains        <- 4
iter          <- 3000
warmup        <- 1000
max_treedepth <- 10
seed          <- 2026

# Plot colours (dataviz reference palette); the two models are validated as a
# categorical pair against the surface.
surface    <- "#fcfcfb"
ink        <- "#0b0b0b"
ink_soft   <- "#52514e"
ink_mute   <- "#898781"
grid_col   <- "#e1e0d9"
blue       <- "#2a78d6"
orange     <- "#eb6834"
chain_cols <- c("#2a78d6", "#eda100", "#1baf7a", "#e34948")
blue_ramp  <- c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#104281")

theme_rq1 <- theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft),
    strip.text       = element_text(colour = ink, hjust = 0, face = "bold")
  )
color_scheme_set(blue_ramp)

# ---- Data (as in 05_rq1_neglect.R) -----------------------------------------

md <- readRDS(in_file) |>
  mutate(n_groups = group_bf_before + group_bf_after + group_rth_before + group_rth_after + group_other) |>
  filter(n_groups == 1, group_other == 0) |>
  transmute(
    extract = as.integer(teeth_to_extract),
    hygiene = 4L - as.integer(oral_hygiene),   # 1 good .. 3 poor
    visit = as.integer(last_dentist_visit),    # 1 <1 year .. 4 never
    xray = as.integer(dental_xray),            # 1 current .. 3 never
    rth = group_rth_before + group_rth_after,
    after = group_bf_after + group_rth_after,
    age_10y = (age_at_visit - 60) / 10,
    male = as.integer(sex == "M"),
    smoking_status,
    visit_year
  ) |>
  drop_na(extract, hygiene, age_10y, male, smoking_status, visit_year)

beta_labels <- c(
  rth = "RTH vs BF, seen before treatment",
  after = "After vs before, BF",
  rth_after = "Interaction: RTH × after",
  age_10y = "Age, per 10 years",
  male = "Male vs female",
  former = "Former vs never smoker",
  current = "Current vs never smoker",
  year_2025 = "Seen in 2025 vs 2024",
  year_2026 = "Seen in 2026 vs 2024"
)
X <- with(md, cbind(
  rth = rth, after = after, rth_after = rth * after, age_10y = age_10y, male = male,
  former = as.integer(smoking_status == "former"), current = as.integer(smoking_status == "current"),
  year_2025 = as.integer(visit_year == 2025), year_2026 = as.integer(visit_year == 2026)
))[, names(beta_labels)]

indicator_labels <- c(extract = "Teeth to extract", visit = "Years since dentist", xray = "Old / no X-ray")

# No LOO here (the outcomes differ from the reported model's), so the
# per-patient log-likelihood is skipped with a 1-node rule.
stan_data <- list(
  N = nrow(md), K = ncol(X), X = X,
  extract = md$extract, hygiene = md$hygiene,
  visit = coalesce(md$visit, 0L), xray = coalesce(md$xray, 0L),
  method_factor = 1L, use_hygiene = 0L,
  Q = 1L, gh_x = as.array(0), gh_logw = as.array(0),
  K_hs = 0L, hs_scale_global = 1, hs_scale_slab = 1, hs_df_slab = 4
)

# ---- Fit -------------------------------------------------------------------

# Inits from (-0.5, 0.5) for the reason given in 05_rq1_neglect.R.
model <- stan_model(stan_file)
fit <- sampling(model,
  data = stan_data,
  chains = chains, iter = iter, warmup = warmup, cores = chains, seed = seed, init_r = 0.5,
  control = list(max_treedepth = max_treedepth),
  pars = c("eta", "log_lik", "hygiene_rep"), include = FALSE, refresh = 0
)

# ---- Diagnostics -----------------------------------------------------------

param_labels <- c(
  set_names(paste("beta:", beta_labels), paste0("beta[", seq_along(beta_labels), "]")),
  alpha = "extract: intercept",
  phi = "extract: NB shape",
  lambda_extract = "loading: teeth to extract",
  lambda_visit = "loading: years since dentist",
  lambda_xray = "loading: old / no X-ray",
  `delta[1]` = "visit-X-ray shared loading",
  `cut_visit[1]` = "cutpoint: visit <1 y|1-5 y",
  `cut_visit[2]` = "cutpoint: visit 1-5 y|5-10 y",
  `cut_visit[3]` = "cutpoint: visit 5-10 y|never",
  `cut_xray[1]` = "cutpoint: X-ray current|>5 y",
  `cut_xray[2]` = "cutpoint: X-ray >5 y|never"
)
draws_main <- as_draws_array(as.array(fit, pars = c(
  "beta", "alpha", "phi", "lambda_extract", "lambda_visit", "lambda_xray", "delta", "cut_visit", "cut_xray"
)))
variables(draws_main) <- unname(param_labels[variables(draws_main)])

latent <- summarise_draws(as_draws_array(as.array(fit, pars = c("z", "u"))), "rhat", "ess_bulk", "ess_tail")
convergence <- bind_rows(
  summarise_draws(draws_main, "mean", "sd", "median", ~ quantile2(.x, probs = c(0.025, 0.975)),
    default_convergence_measures(), default_mcse_measures()),
  tibble(
    variable = sprintf("latent scores (%d), worst", nrow(latent)),
    rhat = max(latent$rhat), ess_bulk = min(latent$ess_bulk), ess_tail = min(latent$ess_tail)
  )
) |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write_csv(convergence, file.path(out_dir, "rq1_no_hygiene_convergence.csv"))

sampler <- imap(get_sampler_params(fit, inc_warmup = FALSE), \(s, chain) {
  e <- s[, "energy__"]
  tibble(
    chain = chain,
    divergent = sum(s[, "divergent__"]),
    max_treedepth_hits = sum(s[, "treedepth__"] >= max_treedepth),
    mean_accept_stat = mean(s[, "accept_stat__"]),
    stepsize = mean(s[, "stepsize__"]),
    ebfmi = sum(diff(e)^2) / sum((e - mean(e))^2)
  )
}) |>
  list_rbind() |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(sampler, file.path(out_dir, "rq1_no_hygiene_sampler.csv"))

trace_plot <- suppressMessages(
  mcmc_trace(draws_main, np = nuts_params(fit), facet_args = list(ncol = 4)) +
    scale_colour_manual(values = chain_cols) +
    theme_rq1 + theme(strip.text = element_text(face = "plain", size = 8))
)
ggsave(file.path(out_dir, "rq1_no_hygiene_trace.png"), trace_plot, width = 15, height = 11, dpi = 300, bg = surface)

rank_plot <- suppressMessages(
  mcmc_rank_overlay(draws_main, facet_args = list(ncol = 4)) +
    scale_colour_manual(values = chain_cols) +
    theme_rq1 + theme(strip.text = element_text(face = "plain", size = 8))
)
ggsave(file.path(out_dir, "rq1_no_hygiene_rank.png"), rank_plot, width = 15, height = 11, dpi = 300, bg = surface)

# ---- Posterior predictive check --------------------------------------------

obs <- list(extract = md$extract, visit = md$visit, xray = md$xray)
level_labels <- list(visit = c("<1 y", "1-5 y", "5-10 y", "never"), xray = c("current", ">5 y", "never"))
ppc_panel <- function(ind) {
  ok <- !is.na(obs[[ind]])
  p <- ppc_bars(obs[[ind]][ok], as.matrix(fit, pars = paste0(ind, "_rep"))[, ok], prob = 0.95) +
    labs(title = sprintf("%s (n = %d)", indicator_labels[[ind]], sum(ok))) +
    theme_rq1
  if (ind == "extract") {
    p + coord_cartesian(xlim = c(-0.5, max(obs$extract) + 0.5))
  } else {
    suppressMessages(p + scale_x_continuous(breaks = seq_along(level_labels[[ind]]), labels = level_labels[[ind]]))
  }
}
ppc_indicators <- wrap_plots(map(names(obs), ppc_panel), ncol = 3, guides = "collect") +
  plot_annotation(title = "Without hygiene: observed vs replicated indicators (median and 95% interval)", theme = theme_rq1)
ggsave(file.path(out_dir, "rq1_no_hygiene_ppc_indicators.png"), ppc_indicators, width = 14, height = 4.5, dpi = 300, bg = surface)

# ---- Effects on neglect, compared with the reported model ------------------

summarise_post <- function(x, null = 0) {
  tibble(
    median = median(x),
    lower_95 = quantile(x, 0.025, names = FALSE), upper_95 = quantile(x, 0.975, names = FALSE),
    lower_50 = quantile(x, 0.25, names = FALSE), upper_50 = quantile(x, 0.75, names = FALSE),
    p_above_null = mean(x > null)
  )
}
prob_more <- function(x) 100 * pnorm(x / sqrt(2))

b <- as.matrix(fit, pars = "beta")
colnames(b) <- names(beta_labels)
structural <- list(
  "Risk type × timing" = list(
    "RTH vs BF, seen before treatment" = b[, "rth"],
    "RTH vs BF, seen after treatment" = b[, "rth"] + b[, "rth_after"],
    "After vs before, BF" = b[, "after"],
    "After vs before, RTH" = b[, "after"] + b[, "rth_after"],
    "Interaction: RTH × after" = b[, "rth_after"]
  ),
  "Group vs BF, seen before" = list(
    "BF, seen after" = b[, "after"],
    "RTH, seen before" = b[, "rth"],
    "RTH, seen after" = b[, "rth"] + b[, "after"] + b[, "rth_after"]
  ),
  "Covariates" = as.list(as.data.frame(b[, c("age_10y", "male", "former", "current", "year_2025", "year_2026")])) |>
    set_names(beta_labels[c("age_10y", "male", "former", "current", "year_2025", "year_2026")])
)
effects <- imap(structural, \(terms, section) {
  imap(terms, \(x, term) {
    bind_rows(
      summarise_post(x) |> mutate(scale = "difference in latent neglect (SD)"),
      if (!startsWith(term, "Interaction")) {
        summarise_post(prob_more(x), null = 50) |> mutate(scale = "probability more neglected (%)")
      }
    ) |>
      mutate(section = section, term = term)
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(model = "no_hygiene")

reported <- read_csv(file.path(out_dir, "rq1_neglect_coefficients.csv"), show_col_types = FALSE) |>
  filter(model == "correlated", section != "Measurement") |>
  mutate(model = "reported")
comparison <- bind_rows(reported, effects) |>
  relocate(model, section, term, scale) |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(comparison, file.path(out_dir, "rq1_no_hygiene_comparison.csv"))

# ---- Measurement -----------------------------------------------------------

# Loadings per SD of neglect and standardised loadings, computed as in
# 05_rq1_neglect.R for the three remaining indicators.
l <- as.matrix(fit, pars = c("lambda_extract", "lambda_visit", "lambda_xray"))
alpha <- as.vector(as.matrix(fit, pars = "alpha"))
phi <- as.vector(as.matrix(fit, pars = "phi"))
delta2 <- as.vector(as.matrix(fit, pars = "delta"))^2
mu <- b %*% t(X)
var_eta <- apply(mu, 1, var) + 1
mean_count <- rowMeans(exp(alpha + l[, "lambda_extract"] * mu)) * exp(l[, "lambda_extract"]^2 / 2)
resid <- cbind(log(1 + 1 / mean_count + 1 / phi), pi^2 / 3 + delta2, pi^2 / 3 + delta2)
std <- l * sqrt(var_eta) / sqrt(l^2 * var_eta + resid)

# Prior of the residual correlation implied by delta ~ Normal+(0, 2.5).
set.seed(seed)
delta_prior <- abs(rnorm(1e5, 0, 2.5))

measurement <- bind_rows(
  summarise_post(exp(l[, "lambda_extract"])) |> mutate(term = "Teeth to extract", scale = "rate ratio per SD"),
  summarise_post(exp(l[, "lambda_visit"])) |> mutate(term = "Years since dentist", scale = "odds ratio per SD"),
  summarise_post(exp(l[, "lambda_xray"])) |> mutate(term = "Old / no X-ray", scale = "odds ratio per SD"),
  map2(unname(indicator_labels), seq_len(ncol(std)), \(term, k) {
    summarise_post(std[, k]) |> mutate(term = term, scale = "standardised loading")
  }) |>
    list_rbind(),
  summarise_post(as.vector(as.matrix(fit, pars = "resid_cor_visit_xray"))) |>
    mutate(term = "Residual correlation, years since dentist - X-ray", scale = "correlation"),
  summarise_post(delta_prior^2 / (delta_prior^2 + pi^2 / 3)) |>
    mutate(term = "Residual correlation, years since dentist - X-ray", scale = "correlation, prior")
) |>
  select(-p_above_null) |>
  relocate(term, scale) |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(measurement, file.path(out_dir, "rq1_no_hygiene_measurement.csv"))

# ---- Comparison graph ------------------------------------------------------

model_names <- c(reported = "Reported model (4 indicators)", no_hygiene = "Without hygiene (3 indicators)")
fmt_ci <- function(m, lo, hi) sprintf("%.0f [%.0f, %.0f]", m, lo, hi)
dodge <- position_dodge(width = 0.6, orientation = "y", reverse = TRUE)
plot_data <- comparison |>
  filter(scale == "probability more neglected (%)", section != "Group vs BF, seen before") |>
  mutate(
    section = factor(section, unique(section)),
    term = factor(term, rev(unique(term))),
    model = factor(model_names[model], model_names),
    label = fmt_ci(median, lower_95, upper_95)
  )
comparison_plot <- ggplot(plot_data, aes(y = term, colour = model, group = model)) +
  geom_vline(xintercept = 50, colour = ink_mute, linewidth = 0.4) +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95), linewidth = 0.6, position = dodge) +
  geom_linerange(aes(xmin = lower_50, xmax = upper_50), linewidth = 2, position = dodge) +
  geom_point(aes(x = median), fill = surface, shape = 21, size = 2.4, stroke = 1.1, position = dodge) +
  geom_text(aes(x = Inf, label = label), hjust = 1, size = 2.8, colour = ink_soft, family = "Arial",
    position = dodge, show.legend = FALSE) +
  scale_colour_manual(values = set_names(c(blue, orange), model_names), name = NULL) +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 90), expand = expansion(mult = c(0.05, 0.3))) +
  facet_grid(section ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Probability of being more neglected than a comparable patient (%)", y = NULL,
    title = "Effects on latent dental neglect, with and without the hygiene indicator",
    subtitle = sprintf("n = %d. Points: posterior medians; thick bars: 50%% and thin bars: 95%% credible intervals.", nrow(md))
  ) +
  theme_rq1 +
  theme(panel.grid.major.y = element_blank(), plot.title.position = "plot", legend.position = "top",
    legend.justification = "left")
ggsave(file.path(out_dir, "rq1_no_hygiene_comparison.png"), comparison_plot, width = 10, height = 7, dpi = 300, bg = surface)

# ---- Console summary -------------------------------------------------------

cat(sprintf("\nn = %d\n\nSampler diagnostics per chain:\n", nrow(md)))
print(as.data.frame(sampler), row.names = FALSE)
cat(sprintf("\nMax R-hat %.3f, min bulk ESS %.0f, min tail ESS %.0f\n",
  max(convergence$rhat), min(convergence$ess_bulk), min(convergence$ess_tail)))
cat("\nProbability more neglected (%), reported vs without hygiene:\n")
comparison |>
  filter(scale == "probability more neglected (%)") |>
  transmute(model, section, term, value = sprintf("%.0f [%.0f, %.0f]", median, lower_95, upper_95)) |>
  pivot_wider(names_from = model, values_from = value) |>
  as.data.frame() |>
  print(row.names = FALSE)
cat("\nMeasurement, without hygiene:\n")
print(as.data.frame(measurement), row.names = FALSE)

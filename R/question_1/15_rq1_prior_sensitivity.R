# RQ1 - prior sensitivity analysis of the latent neglect model (BARG).
#
# Run from anywhere inside the project:  Rscript R/question_1/15_rq1_prior_sensitivity.R
# after R/question_1/05_rq1_neglect.R, whose coefficients it compares against.
#
# Refits the reported model (same patients, covariates, visit-X-ray residual
# correlation and sampler) twice: with every normal prior SD halved
# ("narrower") and doubled ("wider"). The NB shape keeps its Gamma(0.01, 0.01)
# prior in all three. The model is R/question_1/stan/neglect_mimic_priors.stan,
# which takes the prior SDs as data and is otherwise neglect_mimic.stan.
#
# Outputs go to results/question_1/ with the prefix rq1_prior_sens_.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(rstan)
  library(posterior)
})

in_file   <- here("data", "processed", "clean_data.rds")
stan_file <- here("R", "question_1", "stan", "neglect_mimic_priors.stan")
out_dir   <- here("results", "question_1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Sampler settings, as in 05_rq1_neglect.R
chains        <- 4
iter          <- 3000
warmup        <- 1000
max_treedepth <- 10
seed          <- 2026

# Prior SDs of the reported model, and the two alternatives.
reported_sd <- c(sd_beta = 1, sd_alpha = 2.5, sd_lambda_extract = 1, sd_lambda_ordinal = 2.5, sd_cut = 3)
prior_sets <- list(narrower = reported_sd / 2, wider = reported_sd * 2)

# Plot colours (dataviz reference palette).
surface  <- "#fcfcfb"
ink      <- "#0b0b0b"
ink_soft <- "#52514e"
ink_mute <- "#898781"
grid_col <- "#e1e0d9"
blue     <- "#2a78d6"
orange   <- "#eb6834"

theme_rq1 <- theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft),
    strip.text       = element_text(colour = ink, hjust = 0, face = "bold")
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

indicator_labels <- c(
  extract = "Teeth to extract", hygiene = "Poor hygiene",
  visit = "Years since dentist", xray = "Old / no X-ray"
)

# No LOO here, so the per-patient log-likelihood is skipped with a 1-node rule.
stan_data <- list(
  N = nrow(md), K = ncol(X), X = X,
  extract = md$extract, hygiene = md$hygiene,
  visit = coalesce(md$visit, 0L), xray = coalesce(md$xray, 0L),
  method_factor = 1L, use_hygiene = 1L,
  Q = 1L, gh_x = as.array(0), gh_logw = as.array(0),
  K_hs = 0L, hs_scale_global = 1, hs_scale_slab = 1, hs_df_slab = 4
)

# ---- Fit -------------------------------------------------------------------

# Inits from (-0.5, 0.5) for the reason given in 05_rq1_neglect.R.
model <- stan_model(stan_file)
fits <- map(prior_sets, \(sds) {
  sampling(model,
    data = c(stan_data, as.list(sds)),
    chains = chains, iter = iter, warmup = warmup, cores = chains, seed = seed, init_r = 0.5,
    control = list(max_treedepth = max_treedepth),
    pars = c("eta", "log_lik", "extract_rep", "hygiene_rep", "visit_rep", "xray_rep"), include = FALSE,
    refresh = 0
  )
})

# ---- Diagnostics -----------------------------------------------------------

main_pars <- c("beta", "alpha", "phi", "lambda_extract", "lambda_hygiene", "lambda_visit", "lambda_xray",
  "delta", "cut_hygiene", "cut_visit", "cut_xray")
convergence <- imap(fits, \(f, priors) {
  latent <- summarise_draws(as_draws_array(as.array(f, pars = c("z", "u"))), "rhat", "ess_bulk", "ess_tail")
  bind_rows(
    summarise_draws(as_draws_array(as.array(f, pars = main_pars)), "rhat", "ess_bulk", "ess_tail"),
    tibble(
      variable = sprintf("latent scores (%d), worst", nrow(latent)),
      rhat = max(latent$rhat), ess_bulk = min(latent$ess_bulk), ess_tail = min(latent$ess_tail)
    )
  ) |>
    mutate(priors = priors, .before = 1)
}) |>
  list_rbind() |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(convergence, file.path(out_dir, "rq1_prior_sens_convergence.csv"))

sampler <- imap(fits, \(f, priors) {
  imap(get_sampler_params(f, inc_warmup = FALSE), \(s, chain) {
    e <- s[, "energy__"]
    tibble(
      priors = priors,
      chain = chain,
      divergent = sum(s[, "divergent__"]),
      max_treedepth_hits = sum(s[, "treedepth__"] >= max_treedepth),
      mean_accept_stat = mean(s[, "accept_stat__"]),
      stepsize = mean(s[, "stepsize__"]),
      ebfmi = sum(diff(e)^2) / sum((e - mean(e))^2)
    )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(sampler, file.path(out_dir, "rq1_prior_sens_sampler.csv"))

# ---- Effects and loadings, compared with the reported model ----------------

summarise_post <- function(x, null = 0) {
  tibble(
    median = median(x),
    lower_95 = quantile(x, 0.025, names = FALSE), upper_95 = quantile(x, 0.975, names = FALSE),
    lower_50 = quantile(x, 0.25, names = FALSE), upper_50 = quantile(x, 0.75, names = FALSE),
    p_above_null = mean(x > null)
  )
}
prob_more <- function(x) 100 * pnorm(x / sqrt(2))

# Effects on neglect, as in 05_rq1_neglect.R.
effects_of <- function(f) {
  b <- as.matrix(f, pars = "beta")
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
  imap(structural, \(terms, section) {
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
    list_rbind()
}

# Loadings per SD of neglect, the residual correlation and the standardised
# loadings, as in 05_rq1_neglect.R.
measurement_of <- function(f) {
  b <- as.matrix(f, pars = "beta")
  l <- as.matrix(f, pars = c("lambda_extract", "lambda_hygiene", "lambda_visit", "lambda_xray"))
  alpha <- as.vector(as.matrix(f, pars = "alpha"))
  phi <- as.vector(as.matrix(f, pars = "phi"))
  delta2 <- as.vector(as.matrix(f, pars = "delta"))^2
  mu <- b %*% t(X)
  var_eta <- apply(mu, 1, var) + 1
  mean_count <- rowMeans(exp(alpha + l[, "lambda_extract"] * mu)) * exp(l[, "lambda_extract"]^2 / 2)
  resid <- cbind(log(1 + 1 / mean_count + 1 / phi), pi^2 / 3, pi^2 / 3 + delta2, pi^2 / 3 + delta2)
  std <- l * sqrt(var_eta) / sqrt(l^2 * var_eta + resid)
  bind_rows(
    summarise_post(exp(l[, "lambda_extract"])) |> mutate(term = "Teeth to extract", scale = "rate ratio per SD"),
    summarise_post(exp(l[, "lambda_hygiene"])) |> mutate(term = "Poor hygiene", scale = "odds ratio per SD"),
    summarise_post(exp(l[, "lambda_visit"])) |> mutate(term = "Years since dentist", scale = "odds ratio per SD"),
    summarise_post(exp(l[, "lambda_xray"])) |> mutate(term = "Old / no X-ray", scale = "odds ratio per SD"),
    summarise_post(as.vector(as.matrix(f, pars = "resid_cor_visit_xray"))) |>
      mutate(term = "Residual correlation, years since dentist - X-ray", scale = "correlation"),
    map2(unname(indicator_labels), seq_len(ncol(std)), \(term, k) {
      summarise_post(std[, k]) |> mutate(term = term, scale = "standardised loading")
    }) |>
      list_rbind()
  ) |>
    mutate(section = "Measurement") |>
    select(-p_above_null)
}

refits <- imap(fits, \(f, priors) bind_rows(effects_of(f), measurement_of(f)) |> mutate(priors = priors)) |>
  list_rbind()
reported <- read_csv(file.path(out_dir, "rq1_neglect_coefficients.csv"), show_col_types = FALSE) |>
  filter(model == "correlated") |>
  select(-model) |>
  mutate(priors = "reported")
comparison <- bind_rows(reported, refits) |>
  mutate(priors = factor(priors, c("narrower", "reported", "wider"))) |>
  arrange(priors) |>
  relocate(priors, section, term, scale) |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(comparison, file.path(out_dir, "rq1_prior_sens_comparison.csv"))

# ---- Comparison graph ------------------------------------------------------

prior_names <- c(
  narrower = "Narrower priors (SDs halved)",
  reported = "Reported priors",
  wider = "Wider priors (SDs doubled)"
)
fmt_ci <- function(m, lo, hi) sprintf("%.0f [%.0f, %.0f]", m, lo, hi)
dodge <- position_dodge(width = 0.7, orientation = "y", reverse = TRUE)
plot_data <- comparison |>
  filter(scale == "probability more neglected (%)", section != "Group vs BF, seen before") |>
  mutate(
    section = factor(section, unique(section)),
    term = factor(term, rev(unique(term))),
    priors = factor(prior_names[as.character(priors)], prior_names),
    label = fmt_ci(median, lower_95, upper_95)
  )
comparison_plot <- ggplot(plot_data, aes(y = term, colour = priors, group = priors)) +
  geom_vline(xintercept = 50, colour = ink_mute, linewidth = 0.4) +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95), linewidth = 0.6, position = dodge) +
  geom_linerange(aes(xmin = lower_50, xmax = upper_50), linewidth = 2, position = dodge) +
  geom_point(aes(x = median), fill = surface, shape = 21, size = 2.2, stroke = 1.1, position = dodge) +
  geom_text(aes(x = Inf, label = label), hjust = 1, size = 2.6, colour = ink_soft, family = "Arial",
    position = dodge, show.legend = FALSE) +
  scale_colour_manual(values = set_names(c(blue, ink_soft, orange), prior_names), name = NULL) +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 90), expand = expansion(mult = c(0.05, 0.3))) +
  facet_grid(section ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Probability of being more neglected than a comparable patient (%)", y = NULL,
    title = "Effects on latent dental neglect under narrower, reported and wider priors",
    subtitle = sprintf("n = %d. Points: posterior medians; thick bars: 50%% and thin bars: 95%% credible intervals.", nrow(md))
  ) +
  theme_rq1 +
  theme(panel.grid.major.y = element_blank(), plot.title.position = "plot", legend.position = "top",
    legend.justification = "left")
ggsave(file.path(out_dir, "rq1_prior_sens_comparison.png"), comparison_plot, width = 10, height = 8, dpi = 300, bg = surface)

# ---- Console summary -------------------------------------------------------

cat(sprintf("\nn = %d\n\nSampler diagnostics per chain:\n", nrow(md)))
print(as.data.frame(sampler), row.names = FALSE)
convergence |>
  group_by(priors) |>
  summarise(max_rhat = max(rhat), min_ess_bulk = min(ess_bulk), min_ess_tail = min(ess_tail)) |>
  as.data.frame() |>
  print(row.names = FALSE)
cat("\nProbability more neglected (%) and P(beta > 0) under each prior set:\n")
comparison |>
  filter(scale %in% c("probability more neglected (%)", "difference in latent neglect (SD)"),
    section != "Group vs BF, seen before") |>
  group_by(priors, section, term) |>
  summarise(
    value = paste(
      if (any(scale == "probability more neglected (%)")) {
        with(filter(pick(everything()), scale == "probability more neglected (%)"), fmt_ci(median, lower_95, upper_95))
      } else "",
      sprintf("P>0 %.2f", p_above_null[scale == "difference in latent neglect (SD)"])
    ),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = priors, values_from = value) |>
  as.data.frame() |>
  print(row.names = FALSE)
cat("\nStandardised loadings:\n")
comparison |>
  filter(scale == "standardised loading") |>
  transmute(priors, term, value = sprintf("%.2f [%.2f, %.2f]", median, lower_95, upper_95)) |>
  pivot_wider(names_from = priors, values_from = value) |>
  as.data.frame() |>
  print(row.names = FALSE)

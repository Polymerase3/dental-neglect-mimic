# RQ1 - latent dental-neglect model (MIMIC).
#
# Run from anywhere inside the project:  Rscript R/question_1/05_rq1_neglect.R
#
# Fits R/question_1/stan/neglect_mimic.stan to the patients in exactly one of
# groups 1-4 with complete covariates, twice:
# with the residual correlation between last dentist visit and X-ray (the
# reported model, "correlated") and without it ("independent"), for the
# LOO comparison. The indicators come from R/question_1/04_rq1_efa.R.
#
# Missing last-visit / X-ray values are left out of the likelihood. That is
# valid if they are missing at random given the model's variables; visit
# year, which drives the missingness, is a covariate for that reason.
#
# Outputs go to results/ with the prefix rq1_neglect_; the reported fit is
# saved as results/question_1/rq1_neglect_fit.rds.

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
  library(loo)
})

in_file   <- here("data", "processed", "clean_data.rds")
stan_file <- here("R", "question_1", "stan", "neglect_mimic.stan")
out_dir   <- here("results", "question_1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Sampler settings
chains        <- 4
iter          <- 3000
warmup        <- 1000
max_treedepth <- 10
seed          <- 2026
gh_nodes      <- 41   # Gauss-Hermite nodes per latent dimension (LOO); within 0.01 per patient of 51 nodes
n_ppc_draws   <- 1000 # draws used for the correlation check

# Plot colours (dataviz reference palette). Chain colours are validated as a
# set for overlapping lines; the legend is the secondary encoding.
surface    <- "#fcfcfb"
ink        <- "#0b0b0b"
ink_soft   <- "#52514e"
ink_mute   <- "#898781"
grid_col   <- "#e1e0d9"
blue       <- "#2a78d6"
blue_lt    <- "#9ec5f4"
orange     <- "#eb6834"
chain_cols <- c("#2a78d6", "#eda100", "#1baf7a", "#e34948")
blue_ramp  <- c("#cde2fb", "#9ec5f4", "#5598e7", "#2a78d6", "#1c5cab", "#104281")

theme_rq1 <- theme_minimal(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft),
    strip.text       = element_text(colour = ink, hjust = 0, face = "bold")
  )
color_scheme_set(blue_ramp)

# ---- Data ------------------------------------------------------------------

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

gh <- statmod::gauss.quad.prob(gh_nodes, dist = "normal")
stan_data <- list(
  N = nrow(md), K = ncol(X), X = X,
  extract = md$extract, hygiene = md$hygiene,
  visit = coalesce(md$visit, 0L), xray = coalesce(md$xray, 0L),
  Q = gh_nodes, gh_x = gh$nodes, gh_logw = log(gh$weights), use_hygiene = 1L,
  K_hs = 0L, hs_scale_global = 1, hs_scale_slab = 1, hs_df_slab = 4  # no horseshoe block here
)

# ---- Fit -------------------------------------------------------------------

# rstan also evaluates generated quantities at the initial values. With its
# default random inits in (-2, 2) on the unconstrained scale, a replicated
# extraction count could overflow the NB random number generator there, so
# inits are drawn from (-0.5, 0.5).
model <- stan_model(stan_file)
fits <- map(c(correlated = 1L, independent = 0L), \(method_factor) {
  sampling(model,
    data = c(stan_data, method_factor = method_factor),
    chains = chains, iter = iter, warmup = warmup, cores = chains, seed = seed, init_r = 0.5,
    control = list(max_treedepth = max_treedepth), pars = "eta", include = FALSE, refresh = 0
  )
})
saveRDS(fits$correlated, file.path(out_dir, "rq1_neglect_fit.rds"))

# ---- Diagnostics -----------------------------------------------------------

param_labels <- c(
  set_names(paste("beta:", beta_labels), paste0("beta[", seq_along(beta_labels), "]")),
  alpha = "extract: intercept",
  phi = "extract: NB shape",
  lambda_extract = "loading: teeth to extract",
  lambda_hygiene = "loading: poor hygiene",
  lambda_visit = "loading: years since dentist",
  lambda_xray = "loading: old / no X-ray",
  `delta[1]` = "visit-X-ray shared loading",
  `cut_hygiene[1]` = "cutpoint: hygiene good|average",
  `cut_hygiene[2]` = "cutpoint: hygiene average|poor",
  `cut_visit[1]` = "cutpoint: visit <1 y|1-5 y",
  `cut_visit[2]` = "cutpoint: visit 1-5 y|5-10 y",
  `cut_visit[3]` = "cutpoint: visit 5-10 y|never",
  `cut_xray[1]` = "cutpoint: X-ray current|>5 y",
  `cut_xray[2]` = "cutpoint: X-ray >5 y|never"
)
main_pars <- c(
  "beta", "alpha", "phi", "lambda_extract", "lambda_hygiene", "lambda_visit", "lambda_xray",
  "cut_hygiene", "cut_visit", "cut_xray"
)
pars_of <- function(m) c(main_pars, if (m == "correlated") "delta")
latent_of <- function(m) c("z", if (m == "correlated") "u")

main_draws <- function(m) {
  a <- as_draws_array(as.array(fits[[m]], pars = pars_of(m)))
  variables(a) <- unname(param_labels[variables(a)])
  a
}

# Convergence per parameter, plus the worst R-hat / ESS over the latent
# scores of every patient.
convergence <- imap(fits, \(f, m) {
  latent <- summarise_draws(as_draws_array(as.array(f, pars = latent_of(m))), "rhat", "ess_bulk", "ess_tail")
  bind_rows(
    summarise_draws(main_draws(m), "mean", "sd", "median", ~ quantile2(.x, probs = c(0.025, 0.975)),
      default_convergence_measures(), default_mcse_measures()),
    tibble(
      variable = sprintf("latent scores (%d), worst", nrow(latent)),
      rhat = max(latent$rhat), ess_bulk = min(latent$ess_bulk), ess_tail = min(latent$ess_tail)
    )
  ) |>
    mutate(model = m, .before = 1)
}) |>
  list_rbind() |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write_csv(convergence, file.path(out_dir, "rq1_neglect_convergence.csv"))

# Sampler: divergences, max-treedepth hits and E-BFMI per chain.
sampler <- imap(fits, \(f, m) {
  imap(get_sampler_params(f, inc_warmup = FALSE), \(s, chain) {
    e <- s[, "energy__"]
    tibble(
      model = m, chain = chain,
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
write_csv(sampler, file.path(out_dir, "rq1_neglect_sampler.csv"))

# Trace and rank plots of the reported model; divergences under the traces.
draws_main <- main_draws("correlated")
trace_plot <- suppressMessages(
  mcmc_trace(draws_main, np = nuts_params(fits$correlated), facet_args = list(ncol = 4)) +
    scale_colour_manual(values = chain_cols) +
    theme_rq1 + theme(strip.text = element_text(face = "plain", size = 8))
)
ggsave(file.path(out_dir, "rq1_neglect_trace.png"), trace_plot, width = 15, height = 13, dpi = 300, bg = surface)

rank_plot <- suppressMessages(
  mcmc_rank_overlay(draws_main, facet_args = list(ncol = 4)) +
    scale_colour_manual(values = chain_cols) +
    theme_rq1 + theme(strip.text = element_text(face = "plain", size = 8))
)
ggsave(file.path(out_dir, "rq1_neglect_rank.png"), rank_plot, width = 15, height = 13, dpi = 300, bg = surface)

# ---- Posterior predictive checks -------------------------------------------

# Replicated indicators come from new draws of the latent variables, so the
# checks test the model's structure, not just each patient's own fit.
obs <- list(extract = md$extract, hygiene = md$hygiene, visit = md$visit, xray = md$xray)
rep_of <- function(m, ind) as.matrix(fits[[m]], pars = paste0(ind, "_rep"))

level_labels <- list(
  hygiene = c("good", "average", "poor"),
  visit = c("<1 y", "1-5 y", "5-10 y", "never"),
  xray = c("current", ">5 y", "never")
)
ppc_panel <- function(ind) {
  ok <- !is.na(obs[[ind]])
  p <- ppc_bars(obs[[ind]][ok], rep_of("correlated", ind)[, ok], prob = 0.95) +
    labs(title = sprintf("%s (n = %d)", indicator_labels[[ind]], sum(ok))) +
    theme_rq1
  if (ind == "extract") {
    p + coord_cartesian(xlim = c(-0.5, max(obs$extract) + 0.5))
  } else {
    suppressMessages(p + scale_x_continuous(breaks = seq_along(level_labels[[ind]]), labels = level_labels[[ind]]))
  }
}
ppc_indicators <- wrap_plots(map(names(obs), ppc_panel), ncol = 2, guides = "collect") +
  plot_annotation(title = "Observed vs replicated indicators (median and 95% interval)")
ggsave(file.path(out_dir, "rq1_neglect_ppc_indicators.png"), ppc_indicators, width = 11, height = 8, dpi = 300, bg = surface)

# Pairwise Spearman correlations among the indicators: observed vs replicated,
# over the patients with both values.
set.seed(seed)
draw_idx <- sample(chains * (iter - warmup), n_ppc_draws)
pairs <- combn(names(obs), 2, simplify = FALSE)
cor_ppc <- map(names(fits), \(m) {
  reps <- map(set_names(names(obs)), \(ind) rep_of(m, ind)[draw_idx, ])
  map(pairs, \(p) {
    ok <- !is.na(obs[[p[1]]]) & !is.na(obs[[p[2]]])
    observed <- cor(obs[[p[1]]][ok], obs[[p[2]]][ok], method = "spearman")
    r_rep <- map_dbl(seq_along(draw_idx), \(s) cor(reps[[p[1]]][s, ok], reps[[p[2]]][s, ok], method = "spearman"))
    tibble(
      model = m, pair = paste(indicator_labels[[p[1]]], "-", indicator_labels[[p[2]]]), n = sum(ok),
      observed = observed, median = median(r_rep),
      lower_95 = quantile(r_rep, 0.025, names = FALSE), upper_95 = quantile(r_rep, 0.975, names = FALSE),
      p_rep_above_observed = mean(r_rep > observed)
    )
  }) |>
    list_rbind()
}) |>
  list_rbind()
write_csv(mutate(cor_ppc, across(where(is.double), \(x) round(x, 3))), file.path(out_dir, "rq1_neglect_ppc_correlations.csv"))

model_names <- c(correlated = "With visit-X-ray residual correlation", independent = "Without")
cor_plot_data <- cor_ppc |>
  mutate(model = factor(model_names[model], model_names), pair = factor(pair, rev(unique(pair))))
cor_plot <- ggplot(cor_plot_data, aes(y = pair)) +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95, colour = model), linewidth = 1.2,
    position = position_dodge(width = 0.5)) +
  geom_point(aes(x = median, colour = model), size = 2.2, position = position_dodge(width = 0.5)) +
  geom_point(aes(x = observed, shape = "Observed"), data = distinct(cor_plot_data, pair, observed),
    size = 3.2, colour = ink) +
  scale_colour_manual(values = set_names(c(blue, orange), model_names), name = "Replicated (median, 95%)") +
  scale_shape_manual(values = c(Observed = 18), name = NULL) +
  labs(x = "Spearman correlation", y = NULL, title = "Correlations between indicators: observed vs replicated") +
  theme_rq1 +
  theme(legend.position = "top", legend.box = "horizontal")
ggsave(file.path(out_dir, "rq1_neglect_ppc_correlations.png"), cor_plot, width = 9, height = 5, dpi = 300, bg = surface)

# ---- LOO -------------------------------------------------------------------

# log_lik integrates each patient's latent variables out, so PSIS-LOO
# leaves out whole patients.
loos <- map(fits, \(f) {
  ll <- extract_log_lik(f, merge_chains = FALSE)
  loo(ll, r_eff = relative_eff(exp(ll)))
})
loo_cmp <- loo_compare(loos)
loo_table <- imap(loos, \(l, m) {
  tibble(
    model = m,
    elpd_loo = l$estimates["elpd_loo", "Estimate"], se_elpd_loo = l$estimates["elpd_loo", "SE"],
    p_loo = l$estimates["p_loo", "Estimate"], looic = l$estimates["looic", "Estimate"],
    elpd_diff = loo_cmp[m, "elpd_diff"], se_diff = loo_cmp[m, "se_diff"],
    pareto_k_above_0.7 = sum(pareto_k_values(l) > 0.7), max_pareto_k = max(pareto_k_values(l))
  )
}) |>
  list_rbind() |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(loo_table, file.path(out_dir, "rq1_neglect_loo.csv"))
writeLines(c(
  "== With visit-X-ray residual correlation", capture.output(print(loos$correlated)), "",
  "== Without", capture.output(print(loos$independent)), "",
  "== Comparison", capture.output(print(loo_cmp))
), file.path(out_dir, "rq1_neglect_loo.txt"))

# ---- Coefficients ----------------------------------------------------------

summarise_post <- function(x, null = 0) {
  tibble(
    median = median(x),
    lower_95 = quantile(x, 0.025, names = FALSE), upper_95 = quantile(x, 0.975, names = FALSE),
    lower_50 = quantile(x, 0.25, names = FALSE), upper_50 = quantile(x, 0.75, names = FALSE),
    p_above_null = mean(x > null)
  )
}

# Probability more neglected: the chance that a patient with the
# characteristic is more neglected than an otherwise identical patient
# without it. Their neglect differs by beta + (z1 - z0) ~ N(beta, 2), so the
# probability is Phi(beta / sqrt(2)); 50% = no difference.
prob_more <- function(x) 100 * pnorm(x / sqrt(2))

# Effects on neglect, as residual-SD differences and as the probability more
# neglected (not for the interaction, a difference of differences). The
# group contrasts are combinations of rth, after and rth_after; the second
# set compares each group with one reference, BF seen before.
coefficients_of <- function(m) {
  b <- as.matrix(fits[[m]], pars = "beta")
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
    list_rbind() |>
    mutate(model = m, .before = 1)
}

# Standardised loadings (0-1): correlation between neglect and each
# indicator on its latent-response scale. Neglect's total variance is
# var(X beta) + 1. The residual variance is pi^2 / 3 for the ordered logits
# (plus delta^2 for visit and X-ray, the shared term) and, for the NB
# indicator, the log-normal approximation log(1 + 1 / mean count + 1 / phi)
# (Nakagawa, Johnson & Schielzeth 2017).
std_loadings_of <- function(m) {
  f <- fits[[m]]
  b <- as.matrix(f, pars = "beta")
  l <- as.matrix(f, pars = c("lambda_extract", "lambda_hygiene", "lambda_visit", "lambda_xray"))
  alpha <- as.vector(as.matrix(f, pars = "alpha"))
  phi <- as.vector(as.matrix(f, pars = "phi"))
  delta2 <- if (m == "correlated") as.vector(as.matrix(f, pars = "delta"))^2 else 0
  mu <- b %*% t(X)
  var_eta <- apply(mu, 1, var) + 1
  mean_count <- rowMeans(exp(alpha + l[, "lambda_extract"] * mu)) * exp(l[, "lambda_extract"]^2 / 2)
  resid <- cbind(log(1 + 1 / mean_count + 1 / phi), pi^2 / 3, pi^2 / 3 + delta2, pi^2 / 3 + delta2)
  std <- l * sqrt(var_eta) / sqrt(l^2 * var_eta + resid)
  map2(unname(indicator_labels), seq_len(ncol(std)), \(term, k) {
    summarise_post(std[, k]) |> mutate(term = term, scale = "standardised loading")
  }) |>
    list_rbind() |>
    mutate(model = m, section = "Measurement", .before = 1) |>
    select(-p_above_null)
}

# Loadings as the change per 1 SD of neglect: rate ratio for extractions,
# cumulative odds ratio for the ordinal indicators.
measurement_of <- function(m) {
  l <- as.matrix(fits[[m]], pars = c("lambda_extract", "lambda_hygiene", "lambda_visit", "lambda_xray"))
  bind_rows(
    summarise_post(exp(l[, "lambda_extract"])) |> mutate(term = "Teeth to extract", scale = "rate ratio per SD"),
    summarise_post(exp(l[, "lambda_hygiene"])) |> mutate(term = "Poor hygiene", scale = "odds ratio per SD"),
    summarise_post(exp(l[, "lambda_visit"])) |> mutate(term = "Years since dentist", scale = "odds ratio per SD"),
    summarise_post(exp(l[, "lambda_xray"])) |> mutate(term = "Old / no X-ray", scale = "odds ratio per SD"),
    if (m == "correlated") {
      summarise_post(as.vector(as.matrix(fits[[m]], pars = "resid_cor_visit_xray"))) |>
        mutate(term = "Residual correlation, years since dentist - X-ray", scale = "correlation")
    }
  ) |>
    mutate(model = m, section = "Measurement", .before = 1) |>
    select(-p_above_null)
}

coefficients <- map(names(fits), \(m) bind_rows(coefficients_of(m), measurement_of(m), std_loadings_of(m))) |>
  list_rbind() |>
  relocate(model, section, term, scale) |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(coefficients, file.path(out_dir, "rq1_neglect_coefficients.csv"))

# Coefficient graph of the reported model.
fmt_ci <- function(m, lo, hi, digits = 2) sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), m, lo, hi)
# Forest plot; the right-hand margin holds the "median [95% CrI]" labels.
forest <- function(data, xlab, ref, log_x = FALSE, label_room = 0.45, digits = 2, breaks = waiver()) {
  data <- mutate(data, term = factor(term, rev(unique(term))), label = fmt_ci(median, lower_95, upper_95, digits))
  x_scale <- if (log_x) scale_x_log10 else scale_x_continuous
  ggplot(data, aes(y = term)) +
    geom_vline(xintercept = ref, colour = ink_mute, linewidth = 0.4) +
    geom_linerange(aes(xmin = lower_95, xmax = upper_95), colour = blue, linewidth = 0.6) +
    geom_linerange(aes(xmin = lower_50, xmax = upper_50), colour = blue, linewidth = 2) +
    geom_point(aes(x = median), colour = blue, fill = surface, shape = 21, size = 2.6, stroke = 1.1) +
    geom_text(aes(x = Inf, label = label), hjust = 1, size = 3, colour = ink_soft) +
    x_scale(breaks = breaks, expand = expansion(mult = c(0.05, label_room))) +
    labs(x = xlab, y = NULL) +
    theme_rq1 +
    theme(panel.grid.major.y = element_blank(), plot.title.position = "plot")
}
coef_main <- filter(coefficients, model == "correlated")
structural_plot <- forest(
  filter(coef_main, scale == "probability more neglected (%)", section != "Group vs BF, seen before") |>
    mutate(section = factor(section, unique(section))),
  "Probability of being more neglected than a comparable patient (%)", 50, label_room = 0.35, digits = 0,
  breaks = c(10, 25, 50, 75, 90)
) +
  facet_grid(section ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Effects on latent dental neglect")
resid_cor <- filter(coef_main, scale == "correlation")
measurement_plot <- forest(
  filter(coef_main, scale %in% c("rate ratio per SD", "odds ratio per SD")),
  "Change per 1 SD more neglect (rate ratio for extractions, cumulative odds ratio otherwise; log scale)", 1,
  log_x = TRUE, label_room = 0.6
) +
  labs(
    title = "Measurement: how strongly each indicator reflects neglect",
    subtitle = sprintf("Residual correlation, years since dentist - X-ray: %s",
      fmt_ci(resid_cor$median, resid_cor$lower_95, resid_cor$upper_95))
  )
coef_plot <- structural_plot / measurement_plot +
  plot_layout(heights = c(11, 4)) +
  plot_annotation(
    title = "Latent dental-neglect model",
    subtitle = sprintf("n = %d. Points: posterior medians; thick bars: 50%% and thin bars: 95%% credible intervals.", nrow(md)),
    theme = theme(plot.background = element_rect(fill = surface, colour = NA))
  )
ggsave(file.path(out_dir, "rq1_neglect_coefficients.png"), coef_plot, width = 10, height = 9, dpi = 300, bg = surface)

# ---- Console summary -------------------------------------------------------

cat(sprintf("\nn = %d; years since dentist observed for %d, X-ray for %d\n\n",
  nrow(md), sum(!is.na(md$visit)), sum(!is.na(md$xray))))
cat("Sampler diagnostics per chain:\n")
print(as.data.frame(sampler), row.names = FALSE)
cat("\nWorst convergence per model:\n")
print(as.data.frame(summarise(group_by(convergence, model), max_rhat = max(rhat), min_ess_bulk = min(ess_bulk), min_ess_tail = min(ess_tail))), row.names = FALSE)
cat("\nLOO:\n")
print(as.data.frame(loo_table), row.names = FALSE)
cat("\nCoefficients (reported model):\n")
print(as.data.frame(select(coef_main, section, term, scale, median, lower_95, upper_95, p_above_null)), row.names = FALSE)

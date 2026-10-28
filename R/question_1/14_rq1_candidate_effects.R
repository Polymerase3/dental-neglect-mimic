# RQ1 - exploratory: does neglect differ by other patient characteristics?
#
# Run from anywhere inside the project:
#   Rscript R/question_1/14_rq1_candidate_effects.R
#
# A candidate is any characteristic in the clean data that was recorded
# before the visit, is not in the reported MIMIC (05_rq1_neglect.R) and is
# not itself a dental finding. Each candidate is added on its own to the
# reported model: same patients, covariates, priors (Normal(0, 1) for the
# new coefficient too), visit-X-ray residual correlation and sampler. Its
# effect then compares patients who differ in it but share patient group,
# sex, age, smoking status and visit year. One model per candidate rather
# than one joint model: the comorbidities overlap, and a joint model would
# estimate each one beyond the others.
#
# Candidates defined in one group only (e.g. months on BF, within BF after)
# are centred within that group and 0 elsewhere. The group terms absorb the
# other patients, so the effect is a within-group contrast. Years since
# quitting works the same way within former smokers. Patients missing a
# candidate are left out of that candidate's model only.
#
# A candidate is fitted only if the data can inform it beyond patient group
# and sex. A binary one needs at least 10 informative patients: the smaller
# of its yes and no counts, summed over the group x sex cells. A continuous
# one needs at least 15 observed values where it is defined.
# rq1_candidates_screening.csv lists every variable of the clean data with
# its fate.
#
# Exploratory: two dozen estimates, no multiplicity adjustment. Read the
# pattern, not single intervals.
#
# Outputs go to results/question_1/ with the prefix rq1_candidates_; the
# trace, rank and predictive-check plots of each model go to
# results/question_1/rq1_candidates_diagnostics/. The models run in
# parallel, one core each.

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
diag_dir  <- file.path(out_dir, "rq1_candidates_diagnostics")
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

# Sampler settings, as in 05_rq1_neglect.R
chains        <- 4
iter          <- 3000
warmup        <- 1000
max_treedepth <- 10
seed          <- 2026
gh_nodes      <- 41   # Gauss-Hermite nodes per latent dimension (LOO)

# Screening thresholds
min_informative <- 10  # binary candidates
min_observed    <- 15  # continuous candidates

# Plot colours (dataviz reference palette).
surface    <- "#fcfcfb"
ink        <- "#0b0b0b"
ink_soft   <- "#52514e"
ink_mute   <- "#898781"
grid_col   <- "#e1e0d9"
blue       <- "#2a78d6"
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

# ---- Data (the patients of 05_rq1_neglect.R) -------------------------------

md <- readRDS(in_file) |>
  mutate(n_groups = group_bf_before + group_bf_after + group_rth_before + group_rth_after + group_other) |>
  filter(n_groups == 1, group_other == 0) |>
  mutate(
    extract = as.integer(teeth_to_extract),
    hygiene = 4L - as.integer(oral_hygiene),   # 1 good .. 3 poor
    visit = as.integer(last_dentist_visit),    # 1 <1 year .. 4 never
    xray = as.integer(dental_xray),            # 1 current .. 3 never
    rth = group_rth_before + group_rth_after,
    after = group_bf_after + group_rth_after,
    age_10y = (age_at_visit - 60) / 10,
    male = as.integer(sex == "M")
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
X_base <- with(md, cbind(
  rth = rth, after = after, rth_after = rth * after, age_10y = age_10y, male = male,
  former = as.integer(smoking_status == "former"), current = as.integer(smoking_status == "current"),
  year_2025 = as.integer(visit_year == 2025), year_2026 = as.integer(visit_year == 2026)
))[, names(beta_labels)]

# ---- Candidates ------------------------------------------------------------

# Patients in whom each candidate is defined.
where <- with(md, list(
  all = rep(TRUE, nrow(md)),
  smokers = smoking_status != "never",
  former = smoking_status == "former",
  rth = rth == 1,
  bf_before = group_bf_before == 1,
  bf_after = group_bf_after == 1,
  rth_after = group_rth_after == 1
))

# Centred within the patients where it is defined, in `unit`s; 0 elsewhere.
centre_within <- function(x, defined, unit) {
  if_else(defined, (x - mean(x[defined], na.rm = TRUE)) / unit, 0)
}

n_comorb <- with(md, comorb_diabetes + comorb_hypertension + comorb_thyroid + comorb_cardiovascular)
values <- md |>
  transmute(
    diabetes = comorb_diabetes,
    hypertension = comorb_hypertension,
    thyroid = comorb_thyroid,
    cardiovascular = comorb_cardiovascular,
    any_comorb = as.integer(n_comorb > 0),
    n_comorb = n_comorb,
    pack_years_10 = smoking_pack_years / 10,
    years_quit_10 = centre_within(smoking_years_since_quit, where$former, 10),
    other_cancer = cancer_other,
    breast = cancer_breast,
    multiple_cancers = as.integer(cancer_head_neck + cancer_breast + cancer_other >= 2),
    oral_cavity = hn_region_oral_cavity,
    pharynx = hn_region_pharynx,
    larynx = hn_region_larynx,
    neck = hn_region_neck,
    maxilla = hn_region_maxilla,
    chemotherapy = prior_chemotherapy,
    immunotherapy = prior_immunotherapy,
    targeted = prior_targeted_therapy,
    hormone = prior_hormone_therapy,
    radiotherapy = prior_radiotherapy,
    n_systemic = prior_chemotherapy + prior_immunotherapy + prior_hormone_therapy + prior_targeted_therapy,
    months_onco_log2 = log2(months_since_onco_tx + 1),
    bf_therapeutic = if_else(where$bf_before, bf_before_therapeutic, 0L),
    bf_osteoporosis = bf_osteoporosis,
    bf_doses_5 = centre_within(bf_doses, where$bf_after, 5),
    bf_months_12 = centre_within(bf_months, where$bf_after, 12),
    rth_months_12 = centre_within(hn_rth_after_months, where$rth_after, 12),
    uck = as.integer(referring_centre == "UCK")
  )

candidates <- tribble(
  ~id, ~section, ~label, ~type, ~defined, ~source,
  "diabetes", "Chronic disease", "Diabetes", "binary", "all", "comorb_diabetes",
  "hypertension", "Chronic disease", "Hypertension", "binary", "all", "comorb_hypertension",
  "thyroid", "Chronic disease", "Thyroid disease", "binary", "all", "comorb_thyroid",
  "cardiovascular", "Chronic disease", "Cardiovascular disease", "binary", "all", "comorb_cardiovascular",
  "any_comorb", "Chronic disease", "Any of these four", "binary", "all", "comorb_*",
  "n_comorb", "Chronic disease", "Number of these four, per disease", "continuous", "all", "comorb_*",
  "pack_years_10", "Smoking", "Pack-years among smokers, per 10", "continuous", "smokers", "smoking_pack_years",
  "years_quit_10", "Smoking", "Years since quitting, former smokers, per 10", "continuous", "former", "smoking_years_since_quit",
  "other_cancer", "Cancer", "Other cancer (not breast, not head and neck)", "binary", "all", "cancer_other",
  "breast", "Cancer", "Breast cancer", "binary", "all", "cancer_breast",
  "multiple_cancers", "Cancer", "Cancers of two or more types", "binary", "all", "cancer_*",
  "oral_cavity", "Cancer", "Oral cavity, RTH groups", "binary", "rth", "hn_region_oral_cavity",
  "pharynx", "Cancer", "Pharynx, RTH groups", "binary", "rth", "hn_region_pharynx",
  "larynx", "Cancer", "Larynx, RTH groups", "binary", "rth", "hn_region_larynx",
  "neck", "Cancer", "Neck, RTH groups", "binary", "rth", "hn_region_neck",
  "maxilla", "Cancer", "Maxilla, RTH groups", "binary", "rth", "hn_region_maxilla",
  "chemotherapy", "Cancer treatment", "Prior chemotherapy", "binary", "all", "prior_chemotherapy",
  "immunotherapy", "Cancer treatment", "Prior immunotherapy", "binary", "all", "prior_immunotherapy",
  "targeted", "Cancer treatment", "Prior targeted therapy", "binary", "all", "prior_targeted_therapy",
  "hormone", "Cancer treatment", "Prior hormone therapy", "binary", "all", "prior_hormone_therapy",
  "radiotherapy", "Cancer treatment", "Prior radiotherapy, any site", "binary", "all", "prior_radiotherapy",
  "n_systemic", "Cancer treatment", "Systemic treatment types received, per type", "continuous", "all",
  "prior_chemotherapy, prior_immunotherapy, prior_hormone_therapy, prior_targeted_therapy",
  "months_onco_log2", "Cancer treatment", "Months since oncological treatment, per doubling", "continuous", "all", "months_since_onco_tx",
  "bf_therapeutic", "Cancer treatment", "BF for metastases vs prophylactic, BF before", "binary", "bf_before", "bf_before_therapeutic",
  "bf_osteoporosis", "Cancer treatment", "BF for osteoporosis", "binary", "all", "bf_osteoporosis",
  "bf_doses_5", "Cancer treatment", "BF doses received, BF after, per 5", "continuous", "bf_after", "bf_doses",
  "bf_months_12", "Cancer treatment", "Months on BF, BF after, per 12", "continuous", "bf_after", "bf_months",
  "rth_months_12", "Cancer treatment", "Months since head-and-neck RTH, RTH after, per 12", "continuous", "rth_after", "hn_rth_after_months",
  "uck", "Other", "Referred by UCK vs other centres", "binary", "all", "referring_centre"
)

# Support beyond patient group and sex: for a binary candidate, the smaller
# of its yes and no counts in each group x sex cell, summed; for a
# continuous one, the observed values where it is defined.
cells <- interaction(md$rth, md$after, md$male)
screening <- candidates |>
  mutate(
    n = map_int(id, \(v) sum(!is.na(values[[v]]))),
    n_defined = map2_int(id, defined, \(v, w) sum(!is.na(values[[v]]) & where[[w]])),
    n_yes = map2_int(id, type, \(v, t) if (t == "binary") sum(values[[v]] == 1, na.rm = TRUE) else NA_integer_),
    support = pmap_int(list(id, type, n_defined), \(v, t, nd) {
      if (t == "continuous") return(nd)
      x <- values[[v]]
      ok <- !is.na(x)
      as.integer(sum(map_dbl(split(x[ok], cells[ok]), \(y) min(sum(y == 1), sum(y == 0)))))
    }),
    kept = support >= if_else(type == "binary", min_informative, min_observed)
  )

not_candidates <- tribble(
  ~source, ~reason,
  "patient_id", "Identifier",
  "age_at_visit, sex, smoking_status, visit_year, group_*", "Already in the reported model",
  "teeth_to_extract, oral_hygiene, last_dentist_visit, dental_xray", "Indicators of neglect",
  "teeth_to_treat, visit_with_pain, oral_lesions, upper_denture, upper_denture_years, denture_night_use, cbct_needed, referred_maxillofacial, extraction_during_onco_tx",
  "Dental findings at or after the visit: consequences or parts of neglect",
  "dental_assessment_before_tx", "Earlier dental care, which the last-visit and X-ray indicators measure; the exposure of RQ3",
  "history_routine", "Describes how the history was recorded, not the patient; covered by visit year",
  "first_consult_before_tx, first_consult_during_tx, first_consult_after_tx, hn_rth_before, hn_rth_after",
  "Coincide with patient group",
  "cancer_head_neck", "Coincides with the RTH groups",
  "prior_bf_denosumab", "Coincides with BF after",
  "bf_before_prophylactic", "Complement of bf_before_therapeutic within BF before",
  "hn_rth_before_days", "Set after the visit: the dental findings can postpone radiotherapy",
  "hn_site_*", "Finer coding of the head-and-neck region, which is used instead",
  "comorb_none", "Not consistent with the disease codes (DECISIONS.md); any of the four coded diseases is used instead"
)
screening_table <- bind_rows(
  screening |>
    transmute(source, label, status = if_else(kept, "fitted", "too little support"),
      n, n_defined, n_yes, support, reason = NA_character_),
  not_candidates |> mutate(status = "not a candidate")
)
write_csv(screening_table, file.path(out_dir, "rq1_candidates_screening.csv"), na = "")

fit_ids <- screening$id[screening$kept]

# ---- Fits ------------------------------------------------------------------

param_labels <- c(
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
  "beta", "alpha", "phi", "lambda_extract", "lambda_hygiene", "lambda_visit", "lambda_xray", "delta",
  "cut_hygiene", "cut_visit", "cut_xray"
)
indicator_labels <- c(
  extract = "Teeth to extract", hygiene = "Poor hygiene",
  visit = "Years since dentist", xray = "Old / no X-ray"
)
level_labels <- list(
  hygiene = c("good", "average", "poor"),
  visit = c("<1 y", "1-5 y", "5-10 y", "never"),
  xray = c("current", ">5 y", "never")
)
gh <- statmod::gauss.quad.prob(gh_nodes, dist = "normal")

# rstan also evaluates generated quantities at the initial values; inits
# from (-0.5, 0.5) for the reason given in 05_rq1_neglect.R.
model <- stan_model(stan_file)

# Fits one candidate's model, writes its diagnostic plots and returns the
# candidate's coefficient draws with the diagnostics tables.
fit_candidate <- function(cid) {
  info <- filter(candidates, id == cid)
  keep <- !is.na(values[[cid]])
  d <- md[keep, ]
  X <- cbind(X_base[keep, ], values[[cid]][keep])
  fit <- sampling(model,
    data = list(
      N = nrow(d), K = ncol(X), X = X,
      extract = d$extract, hygiene = d$hygiene,
      visit = coalesce(d$visit, 0L), xray = coalesce(d$xray, 0L),
      method_factor = 1L, use_hygiene = 1L,
      Q = gh_nodes, gh_x = gh$nodes, gh_logw = log(gh$weights),
      K_hs = 0L, hs_scale_global = 1, hs_scale_slab = 1, hs_df_slab = 4
    ),
    chains = chains, iter = iter, warmup = warmup, cores = 1, seed = seed, init_r = 0.5,
    control = list(max_treedepth = max_treedepth), pars = "eta", include = FALSE, refresh = 0
  )

  labels <- c(set_names(paste("beta:", c(beta_labels, info$label)), paste0("beta[", seq_len(ncol(X)), "]")), param_labels)
  draws_main <- as_draws_array(as.array(fit, pars = main_pars))
  variables(draws_main) <- unname(labels[variables(draws_main)])

  latent <- summarise_draws(as_draws_array(as.array(fit, pars = c("z", "u"))), "rhat", "ess_bulk", "ess_tail")
  convergence <- bind_rows(
    summarise_draws(draws_main, "mean", "sd", "median", ~ quantile2(.x, probs = c(0.025, 0.975)),
      default_convergence_measures(), default_mcse_measures()),
    tibble(
      variable = sprintf("latent scores (%d), worst", nrow(latent)),
      rhat = max(latent$rhat), ess_bulk = min(latent$ess_bulk), ess_tail = min(latent$ess_tail)
    )
  ) |>
    mutate(candidate = cid, .before = 1)

  sampler <- imap(get_sampler_params(fit, inc_warmup = FALSE), \(s, chain) {
    e <- s[, "energy__"]
    tibble(
      candidate = cid, chain = chain,
      divergent = sum(s[, "divergent__"]),
      max_treedepth_hits = sum(s[, "treedepth__"] >= max_treedepth),
      mean_accept_stat = mean(s[, "accept_stat__"]),
      stepsize = mean(s[, "stepsize__"]),
      ebfmi = sum(diff(e)^2) / sum((e - mean(e))^2)
    )
  }) |>
    list_rbind()

  ll <- extract_log_lik(fit, merge_chains = FALSE)
  l <- loo(ll, r_eff = relative_eff(exp(ll)))
  loo_row <- tibble(
    elpd_loo = l$estimates["elpd_loo", "Estimate"], p_loo = l$estimates["p_loo", "Estimate"],
    pareto_k_above_0.7 = sum(pareto_k_values(l) > 0.7), max_pareto_k = max(pareto_k_values(l))
  )

  diag_file <- \(what) file.path(diag_dir, sprintf("%s_%s.png", cid, what))
  trace_plot <- suppressMessages(
    mcmc_trace(draws_main, np = nuts_params(fit), facet_args = list(ncol = 4)) +
      scale_colour_manual(values = chain_cols) +
      labs(title = info$label) +
      theme_rq1 + theme(strip.text = element_text(face = "plain", size = 8))
  )
  ggsave(diag_file("trace"), trace_plot, width = 15, height = 13, dpi = 300, bg = surface)
  rank_plot <- suppressMessages(
    mcmc_rank_overlay(draws_main, facet_args = list(ncol = 4)) +
      scale_colour_manual(values = chain_cols) +
      labs(title = info$label) +
      theme_rq1 + theme(strip.text = element_text(face = "plain", size = 8))
  )
  ggsave(diag_file("rank"), rank_plot, width = 15, height = 13, dpi = 300, bg = surface)

  # Replicated indicators come from new draws of the latent variables, as in 05.
  obs <- list(extract = d$extract, hygiene = d$hygiene, visit = d$visit, xray = d$xray)
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
  ppc_plot <- wrap_plots(map(names(obs), ppc_panel), ncol = 2, guides = "collect") +
    plot_annotation(title = sprintf("%s: observed vs replicated indicators (median and 95%% interval)", info$label),
      theme = theme_rq1)
  ggsave(diag_file("ppc"), ppc_plot, width = 11, height = 8, dpi = 300, bg = surface)

  list(
    candidate = cid, n = nrow(d), beta = as.matrix(fit, pars = "beta")[, ncol(X)],
    convergence = convergence, sampler = sampler, loo = loo_row
  )
}

cat(sprintf("n = %d; fitting %d candidate models in parallel\n", nrow(md), length(fit_ids)))
results <- parallel::mclapply(fit_ids, fit_candidate, mc.cores = min(length(fit_ids), parallel::detectCores()))
failed <- map_lgl(results, \(r) inherits(r, "try-error"))
if (any(failed)) stop("fits failed for ", paste(fit_ids[failed], collapse = ", "), ":\n", results[[which(failed)[1]]])
names(results) <- fit_ids
saveRDS(map(results, "beta"), file.path(out_dir, "rq1_candidates_draws.rds"))

# ---- Diagnostics -----------------------------------------------------------

convergence <- map(results, "convergence") |>
  list_rbind() |>
  mutate(across(where(is.numeric), \(x) round(x, 3)))
write_csv(convergence, file.path(out_dir, "rq1_candidates_convergence.csv"))

sampler <- map(results, "sampler") |>
  list_rbind() |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(sampler, file.path(out_dir, "rq1_candidates_sampler.csv"))

diagnostics <- map(results, \(r) {
  bind_cols(
    tibble(
      candidate = r$candidate, n = r$n,
      divergent = sum(r$sampler$divergent), max_treedepth_hits = sum(r$sampler$max_treedepth_hits),
      min_ebfmi = min(r$sampler$ebfmi),
      max_rhat = max(r$convergence$rhat, na.rm = TRUE),
      min_ess_bulk = min(r$convergence$ess_bulk, na.rm = TRUE),
      min_ess_tail = min(r$convergence$ess_tail, na.rm = TRUE)
    ),
    r$loo
  )
}) |>
  list_rbind() |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(diagnostics, file.path(out_dir, "rq1_candidates_diagnostics.csv"))

# ---- Effects ---------------------------------------------------------------

summarise_post <- function(x, null = 0) {
  tibble(
    median = median(x),
    lower_95 = quantile(x, 0.025, names = FALSE), upper_95 = quantile(x, 0.975, names = FALSE),
    lower_50 = quantile(x, 0.25, names = FALSE), upper_50 = quantile(x, 0.75, names = FALSE),
    p_above_null = mean(x > null)
  )
}
# Probability more neglected, as in 05_rq1_neglect.R: Phi(beta / sqrt(2)).
prob_more <- function(x) 100 * pnorm(x / sqrt(2))

effects <- map(results, \(r) {
  bind_rows(
    summarise_post(r$beta) |> mutate(scale = "difference in latent neglect (SD)"),
    summarise_post(prob_more(r$beta), null = 50) |> mutate(scale = "probability more neglected (%)")
  ) |>
    mutate(candidate = r$candidate, .before = 1)
}) |>
  list_rbind() |>
  left_join(select(screening, candidate = id, section, term = label, type, n, n_defined, n_yes, support), by = "candidate") |>
  relocate(candidate, section, term, type, n, n_defined, n_yes, support, scale)
write_csv(mutate(effects, across(where(is.double), \(x) round(x, 3))), file.path(out_dir, "rq1_candidates_effects.csv"), na = "")

# ---- Forest plot -----------------------------------------------------------

fmt_ci <- function(m, lo, hi) sprintf("%.0f [%.0f, %.0f]", m, lo, hi)
plot_data <- effects |>
  filter(scale == "probability more neglected (%)") |>
  mutate(
    section = factor(section, unique(candidates$section)),
    term = if_else(type == "binary", sprintf("%s (%d/%d)", term, n_yes, n_defined), sprintf("%s (n = %d)", term, n_defined)),
    term = factor(term, rev(unique(term))),
    label = fmt_ci(median, lower_95, upper_95)
  )
forest_plot <- ggplot(plot_data, aes(y = term)) +
  geom_vline(xintercept = 50, colour = ink_mute, linewidth = 0.4) +
  geom_linerange(aes(xmin = lower_95, xmax = upper_95), colour = blue, linewidth = 0.6) +
  geom_linerange(aes(xmin = lower_50, xmax = upper_50), colour = blue, linewidth = 2) +
  geom_point(aes(x = median), colour = blue, fill = surface, shape = 21, size = 2.6, stroke = 1.1) +
  geom_text(aes(x = Inf, label = label), hjust = 1, size = 3, colour = ink_soft, family = "Arial") +
  scale_x_continuous(breaks = c(10, 25, 50, 75, 90), expand = expansion(mult = c(0.05, 0.3))) +
  facet_grid(section ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Probability of being more neglected than a comparable patient (%)", y = NULL,
    title = "Patient characteristics and latent dental neglect",
    subtitle = paste0(
      "Each characteristic added on its own to the reported model, so patients are compared at the same patient group, sex,\n",
      "age, smoking status and visit year. Brackets: patients with the characteristic / patients where it is defined.\n",
      "Points: posterior medians; thick bars: 50% and thin bars: 95% credible intervals."
    )
  ) +
  theme_rq1 +
  theme(panel.grid.major.y = element_blank(), plot.title.position = "plot")
ggsave(file.path(out_dir, "rq1_candidates_effects.png"), forest_plot, width = 10, height = 9.5, dpi = 300, bg = surface)

# ---- Console summary -------------------------------------------------------

cat("\nScreening:\n")
print(as.data.frame(select(screening, id, type, n, n_defined, n_yes, support, kept)), row.names = FALSE)
cat("\nDiagnostics per model:\n")
print(as.data.frame(diagnostics), row.names = FALSE)
cat("\nProbability more neglected (%):\n")
effects |>
  filter(scale == "probability more neglected (%)") |>
  transmute(section, term, n, value = fmt_ci(median, lower_95, upper_95), p_above_50 = round(p_above_null, 3)) |>
  as.data.frame() |>
  print(row.names = FALSE)
n_excl <- sum(with(filter(effects, scale == "probability more neglected (%)"), lower_95 > 50 | upper_95 < 50))
cat(sprintf("\n%d of %d 95%% intervals exclude 50%%\n", n_excl, length(fit_ids)))

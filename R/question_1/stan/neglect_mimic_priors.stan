// Latent dental-neglect model (MIMIC) with the prior SDs passed in as data,
// for the prior sensitivity analysis in R/question_1/15_rq1_prior_sensitivity.R.
// Everything else is identical to neglect_mimic.stan, whose header describes
// the model. With sd_beta = 1, sd_alpha = 2.5, sd_lambda_extract = 1,
// sd_lambda_ordinal = 2.5 and sd_cut = 3 it is the reported model.

data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  array[N] int<lower=0> extract;
  array[N] int<lower=1, upper=3> hygiene;
  array[N] int<lower=0, upper=4> visit;
  array[N] int<lower=0, upper=3> xray;
  int<lower=0, upper=1> method_factor;
  int<lower=0, upper=1> use_hygiene;
  int<lower=0, upper=K - 1> K_hs;
  real<lower=0> hs_scale_global;
  real<lower=0> hs_scale_slab;
  real<lower=0> hs_df_slab;
  // Gauss-Hermite rule for E[f(x)], x ~ N(0, 1), used to integrate the
  // latent variables out of the per-patient log-likelihood (for LOO).
  int<lower=1> Q;
  vector[Q] gh_x;
  vector[Q] gh_logw;
  // Prior SDs. sd_lambda_ordinal covers the hygiene, visit and X-ray loadings
  // and delta; sd_cut covers all cutpoints.
  real<lower=0> sd_beta;
  real<lower=0> sd_alpha;
  real<lower=0> sd_lambda_extract;
  real<lower=0> sd_lambda_ordinal;
  real<lower=0> sd_cut;
}
transformed data {
  int n_visit = 0;
  int n_xray = 0;
  for (n in 1:N) {
    n_visit += visit[n] > 0;
    n_xray += xray[n] > 0;
  }
  array[n_visit] int idx_visit;
  array[n_xray] int idx_xray;
  {
    int i = 1;
    int j = 1;
    for (n in 1:N) {
      if (visit[n] > 0) { idx_visit[i] = n; i += 1; }
      if (xray[n] > 0) { idx_xray[j] = n; j += 1; }
    }
  }
}
parameters {
  vector[K - K_hs] beta_fixed;
  vector[N] z;
  vector[method_factor ? N : 0] u;
  real alpha;
  real<lower=0> phi;
  real<lower=0> lambda_extract;
  real<lower=0> lambda_hygiene;
  real<lower=0> lambda_visit;
  real<lower=0> lambda_xray;
  vector<lower=0>[method_factor] delta;
  ordered[2] cut_hygiene;
  ordered[3] cut_visit;
  ordered[2] cut_xray;
  vector[K_hs] hs_z;
  vector<lower=0>[K_hs] hs_local;
  vector<lower=0>[K_hs > 0] hs_global;
  vector<lower=0>[K_hs > 0] hs_slab2;
}
transformed parameters {
  vector[K] beta = append_row(beta_fixed, rep_vector(0, K_hs));
  if (K_hs > 0) {
    vector[K_hs] lambda_tilde = sqrt(hs_slab2[1] * square(hs_local)
      ./ (hs_slab2[1] + square(hs_global[1]) * square(hs_local)));
    beta[(K - K_hs + 1):K] = hs_z .* lambda_tilde * hs_global[1];
  }
  vector[N] eta = X * beta + z;
}
model {
  vector[N] shared = rep_vector(0, N);
  if (method_factor == 1) shared = delta[1] * u;

  beta_fixed ~ normal(0, sd_beta);
  hs_z ~ std_normal();
  hs_local ~ student_t(1, 0, 1);
  hs_global ~ student_t(1, 0, hs_scale_global);
  hs_slab2 ~ inv_gamma(0.5 * hs_df_slab, 0.5 * hs_df_slab * square(hs_scale_slab));
  z ~ std_normal();
  u ~ std_normal();
  alpha ~ normal(1, sd_alpha);
  phi ~ gamma(0.01, 0.01);
  lambda_extract ~ normal(0, sd_lambda_extract);
  lambda_hygiene ~ normal(0, sd_lambda_ordinal);
  lambda_visit ~ normal(0, sd_lambda_ordinal);
  lambda_xray ~ normal(0, sd_lambda_ordinal);
  delta ~ normal(0, sd_lambda_ordinal);
  cut_hygiene ~ normal(0, sd_cut);
  cut_visit ~ normal(0, sd_cut);
  cut_xray ~ normal(0, sd_cut);

  extract ~ neg_binomial_2_log(alpha + lambda_extract * eta, phi);
  if (use_hygiene == 1) hygiene ~ ordered_logistic(lambda_hygiene * eta, cut_hygiene);
  visit[idx_visit] ~ ordered_logistic(lambda_visit * eta[idx_visit] + shared[idx_visit], cut_visit);
  xray[idx_xray] ~ ordered_logistic(lambda_xray * eta[idx_xray] + shared[idx_xray], cut_xray);
}
generated quantities {
  real resid_cor_visit_xray = 0;
  vector[N] log_lik;
  array[N] int extract_rep;
  array[N] int hygiene_rep;
  array[N] int visit_rep;
  array[N] int xray_rep;
  {
    real d = 0;
    if (method_factor == 1) d = delta[1];
    resid_cor_visit_xray = square(d) / (square(d) + square(pi()) / 3);
    vector[N] mu = X * beta;
    int R = method_factor == 1 ? Q : 1;
    vector[Q * R] lp;
    for (n in 1:N) {
      // Marginal log-likelihood: z (and u) integrated out by quadrature.
      for (q in 1:Q) {
        real e = mu[n] + gh_x[q];
        real base = gh_logw[q]
          + neg_binomial_2_log_lpmf(extract[n] | alpha + lambda_extract * e, phi);
        if (use_hygiene == 1) base += ordered_logistic_lpmf(hygiene[n] | lambda_hygiene * e, cut_hygiene);
        for (r in 1:R) {
          real s = d * gh_x[r];
          real l = base + (method_factor == 1 ? gh_logw[r] : 0);
          if (visit[n] > 0) l += ordered_logistic_lpmf(visit[n] | lambda_visit * e + s, cut_visit);
          if (xray[n] > 0) l += ordered_logistic_lpmf(xray[n] | lambda_xray * e + s, cut_xray);
          lp[(q - 1) * R + r] = l;
        }
      }
      log_lik[n] = log_sum_exp(lp);

      // Replicated indicators for new draws of the latent variables.
      real e_rep = mu[n] + normal_rng(0, 1);
      real s_rep = d * normal_rng(0, 1);
      extract_rep[n] = neg_binomial_2_log_rng(alpha + lambda_extract * e_rep, phi);
      hygiene_rep[n] = ordered_logistic_rng(lambda_hygiene * e_rep, cut_hygiene);
      visit_rep[n] = ordered_logistic_rng(lambda_visit * e_rep + s_rep, cut_visit);
      xray_rep[n] = ordered_logistic_rng(lambda_xray * e_rep + s_rep, cut_xray);
    }
  }
}

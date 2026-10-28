# RQ1, measurement: do the indicators of dental neglect form one factor?
#
# Run from anywhere inside the project:  Rscript R/question_1/04_rq1_efa.R
#
# Exploratory factor analysis on the whole cohort, first on all 7 candidate
# indicators, then on the 4 core indicators used by the latent neglect model
# (R/question_1/05_rq1_neglect.R). A lavaan CFA of the core indicators checks the fit
# with and without a residual correlation between last dentist visit and
# X-ray.
#
# Correlations are polychoric (ordinal pairs), polyserial (ordinal with
# log counts) and Pearson (log counts), on pairwise-complete data. Every
# indicator is coded so that higher = more neglect.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
  library(psych)
  library(lavaan)
})

in_file <- here("data", "processed", "clean_data.rds")
out_dir <- here("results", "question_1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2026) # parallel analysis simulates random data

# Plot colours (dataviz reference palette).
surface  <- "#fcfcfb"
ink      <- "#0b0b0b"
ink_soft <- "#52514e"
ink_mute <- "#898781"
grid_col <- "#e1e0d9"
series   <- c("#2a78d6", "#eb6834")
rho_fill <- c("#2a78d6", "#f0efec", "#e34948")

theme_rq1 <- theme_minimal(base_size = 10) +
  theme(
    plot.background  = element_rect(fill = surface, colour = NA),
    panel.grid.major = element_line(colour = grid_col, linewidth = 0.3),
    panel.grid.minor = element_blank(),
    text             = element_text(colour = ink),
    axis.text        = element_text(colour = ink_soft),
    strip.text       = element_text(colour = ink, hjust = 0)
  )

# ---- Indicators ------------------------------------------------------------

labels <- c(
  hygiene_poor = "Poor hygiene",
  last_visit   = "Years since dentist",
  xray_old     = "Old / no X-ray",
  extract_log  = "Teeth to extract (log)",
  pain         = "Visit with pain",
  treat_log    = "Teeth to treat (log)",
  cbct         = "CBCT needed"
)
core <- c("hygiene_poor", "last_visit", "xray_old", "extract_log")
continuous <- c("extract_log", "treat_log")
dichotomous <- c("pain", "cbct")

d <- readRDS(in_file)
ind <- tibble(
  hygiene_poor = 4L - as.integer(d$oral_hygiene), # 1 good .. 3 poor
  last_visit   = as.integer(d$last_dentist_visit), # 1 <1 year .. 4 never
  xray_old     = as.integer(d$dental_xray),        # 1 current .. 3 never
  extract_log  = log1p(d$teeth_to_extract),
  pain         = d$visit_with_pain,
  treat_log    = log1p(d$teeth_to_treat),
  cbct         = d$cbct_needed
)

# ---- EFA -------------------------------------------------------------------

# Mixed correlation matrix and pairwise n for a set of indicators.
mixed_cor <- function(vars) {
  x <- as.data.frame(ind[vars])
  idx <- \(type) which(vars %in% type)
  r <- mixedCor(x,
    c = idx(continuous), p = idx(setdiff(vars, c(continuous, dichotomous))), d = idx(dichotomous),
    use = "pairwise", global = FALSE
  )$rho
  n <- crossprod(!is.na(as.matrix(x)))
  list(r = r[vars, vars], n = n)
}

# One EFA set: correlations, parallel analysis, 1-factor solution, fit.
# Fit statistics use the harmonic mean of the pairwise n.
run_efa <- function(vars, set) {
  m <- mixed_cor(vars)
  n_obs <- 1 / mean(1 / m$n[upper.tri(m$n)])
  pa <- fa.parallel(m$r, n.obs = n_obs, fa = "fa", plot = FALSE)
  f1 <- fa(m$r, nfactors = 1, n.obs = n_obs, fm = "minres")
  fit <- c(
    "harmonic mean pairwise n" = n_obs,
    "KMO" = KMO(m$r)$MSA,
    "Bartlett chi-square" = cortest.bartlett(m$r, n = n_obs)$chisq,
    "Bartlett p" = cortest.bartlett(m$r, n = n_obs)$p.value,
    "parallel analysis: factors" = pa$nfact,
    "chi-square" = f1$STATISTIC,
    "df" = f1$dof,
    "RMSEA" = f1$RMSEA[["RMSEA"]],
    "RMSEA 90% lower" = f1$RMSEA[["lower"]],
    "RMSEA 90% upper" = f1$RMSEA[["upper"]],
    "TLI" = f1$TLI,
    "RMSR" = f1$rms,
    "variance explained" = f1$Vaccounted["Proportion Var", 1]
  )
  list(
    set = set, vars = vars, r = m$r, n = m$n, pa = pa, f1 = f1,
    fit = tibble(set = set, measure = names(fit), value = unname(fit))
  )
}

efa_all <- run_efa(names(labels), "EFA, 7 candidate indicators")
efa_core <- run_efa(core, "EFA, 4 core indicators")

# Reliability of the core set.
omega_core <- suppressMessages(omega(efa_core$r, nfactors = 1, n.obs = efa_core$fit$value[1], plot = FALSE))
reliability <- tibble(
  set = efa_core$set,
  measure = c("omega total", "alpha"),
  value = c(omega_core$omega.tot, omega_core$alpha)
)

# ---- CFA of the core indicators --------------------------------------------

cfa_data <- ind |>
  select(all_of(core)) |>
  mutate(across(c(hygiene_poor, last_visit, xray_old), ordered))
cfa_one <- "neglect =~ hygiene_poor + last_visit + xray_old + extract_log"
cfa_models <- list(
  "CFA, 1 factor" = cfa_one,
  "CFA, 1 factor + visit-X-ray residual" = paste(cfa_one, "last_visit ~~ xray_old", sep = "\n")
)
cfa_fits <- map(cfa_models, \(m) {
  cfa(m, data = cfa_data, ordered = c("hygiene_poor", "last_visit", "xray_old"),
    estimator = "WLSMV", missing = "pairwise", std.lv = TRUE)
})
cfa_measures <- c(
  "chi-square" = "chisq.scaled", "df" = "df.scaled", "p" = "pvalue.scaled",
  "CFI" = "cfi.scaled", "TLI" = "tli.scaled", "RMSEA" = "rmsea.scaled",
  "RMSEA 90% lower" = "rmsea.ci.lower.scaled", "RMSEA 90% upper" = "rmsea.ci.upper.scaled",
  "SRMR" = "srmr"
)
cfa_fit <- imap(cfa_fits, \(f, set) {
  tibble(set = set, measure = names(cfa_measures), value = as.numeric(fitMeasures(f, cfa_measures)))
}) |>
  list_rbind()
cfa_resid <- standardizedSolution(cfa_fits[[2]]) |>
  filter(lhs == "last_visit", op == "~~", rhs == "xray_old")
cfa_fit <- bind_rows(cfa_fit, tibble(
  set = names(cfa_models)[2],
  measure = c("residual correlation visit-X-ray", "residual correlation 95% lower", "residual correlation 95% upper"),
  value = as.numeric(c(cfa_resid$est.std, cfa_resid$ci.lower, cfa_resid$ci.upper))
))

# ---- Tables ----------------------------------------------------------------

cor_long <- function(e) {
  expand_grid(var1 = e$vars, var2 = e$vars) |>
    filter(match(var1, e$vars) < match(var2, e$vars)) |>
    mutate(r = e$r[cbind(var1, var2)], n = e$n[cbind(var1, var2)])
}
write_csv(mutate(cor_long(efa_all), r = round(r, 3)), file.path(out_dir, "rq1_efa_correlations.csv"))

cfa_std <- standardizedSolution(cfa_fits[[2]]) |> filter(op == "=~")
loadings <- tibble(indicator = names(labels), label = unname(labels)) |>
  mutate(
    loading_all7 = efa_all$f1$loadings[indicator, 1],
    communality_all7 = efa_all$f1$communality[indicator],
    loading_core = efa_core$f1$loadings[match(indicator, core), 1],
    communality_core = efa_core$f1$communality[match(indicator, core)],
    loading_cfa = cfa_std$est.std[match(indicator, cfa_std$rhs)],
    core = indicator %in% core
  ) |>
  mutate(across(where(is.double), \(x) round(x, 3)))
write_csv(loadings, file.path(out_dir, "rq1_efa_loadings.csv"))

fit_table <- bind_rows(efa_all$fit, efa_core$fit, reliability, cfa_fit) |>
  mutate(value = round(value, 3))
write_csv(fit_table, file.path(out_dir, "rq1_efa_fit.csv"))

# ---- Plots -----------------------------------------------------------------

# Correlation heatmap, lower triangle, core indicators first.
heat <- cor_long(efa_all) |>
  mutate(
    var1 = factor(labels[var1], labels),
    var2 = factor(labels[var2], rev(labels))
  )
heat_plot <- ggplot(heat, aes(var1, var2, fill = r)) +
  geom_tile(colour = surface, linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f\nn=%d", r, n)), size = 2.8, colour = ink, lineheight = 0.9) +
  scale_fill_gradientn(colours = rho_fill, limits = c(-1, 1), name = "r") +
  labs(x = NULL, y = NULL, title = "Indicator correlations (polychoric / polyserial, pairwise)") +
  theme_rq1 +
  theme(panel.grid.major = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out_dir, "rq1_efa_correlations.png"), heat_plot, width = 8, height = 6.5, dpi = 300, bg = surface)

# Scree plot against parallel analysis.
scree <- map(list(efa_all, efa_core), \(e) {
  k <- seq_along(e$pa$fa.values)
  bind_rows(
    tibble(set = e$set, factor = k, eigenvalue = e$pa$fa.values, series = "Observed data"),
    tibble(set = e$set, factor = k, eigenvalue = e$pa$fa.sim, series = "Simulated random data")
  )
}) |>
  list_rbind() |>
  mutate(set = factor(set, c(efa_all$set, efa_core$set)))
write_csv(mutate(scree, eigenvalue = round(eigenvalue, 3)), file.path(out_dir, "rq1_efa_scree.csv"))
scree_plot <- ggplot(scree, aes(factor, eigenvalue, colour = series)) +
  geom_hline(yintercept = 0, colour = ink_mute, linewidth = 0.3) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~set, scales = "free_x") +
  scale_colour_manual(values = c("Observed data" = series[1], "Simulated random data" = ink_mute), name = NULL) +
  scale_x_continuous(breaks = 1:7) +
  labs(x = "Factor", y = "Eigenvalue (reduced correlation matrix)", title = "Scree plot with parallel analysis") +
  theme_rq1 +
  theme(legend.position = "top")
ggsave(file.path(out_dir, "rq1_efa_scree.png"), scree_plot, width = 9, height = 4.5, dpi = 300, bg = surface)

# One-factor loadings from both EFA sets and the CFA.
load_long <- loadings |>
  select(label, core, `EFA, 7 indicators` = loading_all7, `EFA, 4 core` = loading_core, `CFA, 4 core + residual` = loading_cfa) |>
  pivot_longer(-c(label, core), names_to = "solution", values_to = "loading") |>
  filter(!is.na(loading)) |>
  mutate(
    label = factor(label, rev(labels)),
    solution = factor(solution, c("EFA, 7 indicators", "EFA, 4 core", "CFA, 4 core + residual"))
  )
load_plot <- ggplot(load_long, aes(loading, label, colour = solution)) +
  geom_vline(xintercept = 0, colour = ink_mute, linewidth = 0.4) +
  geom_vline(xintercept = c(-0.4, 0.4), colour = ink_mute, linewidth = 0.3) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_colour_manual(values = c(series, "#1baf7a"), name = NULL) +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
  labs(x = "Standardised loading on one factor", y = NULL,
       title = "One-factor loadings", subtitle = "Thin lines at ±0.4") +
  theme_rq1 +
  theme(legend.position = "top", panel.grid.major.x = element_blank())
ggsave(file.path(out_dir, "rq1_efa_loadings.png"), load_plot, width = 8, height = 4.5, dpi = 300, bg = surface)

# ---- Console summary -------------------------------------------------------

print(as.data.frame(loadings), row.names = FALSE)
cat("\n")
print(as.data.frame(pivot_wider(fit_table, names_from = set, values_from = value)), row.names = FALSE)

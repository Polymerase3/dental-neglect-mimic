# Article 1 - Table S1: factor-analysis fit of the indicators of dental neglect.
#           Table S2: posterior summaries of the reported MIMIC model.
#
# Run from anywhere inside the project:  Rscript R/question_1/21_art1_tables_s1_s2.R
#
# Table S1 has one column per solution: EFA of the 7 candidates, EFA of the 4
# retained indicators, CFA of the 4 with one factor, and with the visit-X-ray
# residual correlation. All analyses use the 240 patients in the register,
# pairwise.
#
# Table S2 has the effects on neglect (beta, probability more neglected,
# P(beta > 0), and the probability from the model adjusted for visit year
# only), then the measurement part (RR or OR per SD, standardised loadings,
# residual correlation).
#
# Inputs (nothing is refitted):
#   results/question_1/rq1_efa_fit.csv               04_rq1_efa.R
#   results/question_1/rq1_neglect_coefficients.csv  05_rq1_neglect.R
#   results/question_1/rq1_unadjusted_comparison.csv 20_rq1_neglect_unadjusted.R
# Outputs: results/art1_MIMIC/art1_table_s1.csv and art1_table_s2.csv

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

res_dir <- here("results", "question_1")
out_dir <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sets <- c(
  "EFA, 7 candidate indicators"          = "EFA, 7 candidates",
  "EFA, 4 core indicators"               = "EFA, 4 retained",
  "CFA, 1 factor"                        = "CFA, one factor",
  "CFA, 1 factor + visit-X-ray residual" = "CFA, one factor + residual correlation"
)
n_vars <- c(7, 4, 4, 4)

fit <- read_csv(file.path(res_dir, "rq1_efa_fit.csv"), show_col_types = FALSE) |>
  pivot_wider(names_from = measure, values_from = value) |>
  mutate(set = factor(set, names(sets))) |>
  arrange(set)
stopifnot(identical(as.character(fit$set), names(sets)))

# Cell formats; an empty cell is a measure that does not apply.
f2 <- function(x) if_else(is.na(x), "", sub("-", "−", sprintf("%.2f", x)))
p_fmt <- function(p) if_else(p < 0.001, "< 0.001", sprintf("%.3f", p))
ci <- function(est, lo, hi) if_else(is.na(est), "", sprintf("%s [%s-%s]", f2(est), f2(lo), f2(hi)))

# The EFA chi-square p-values are not stored; they follow from chi-square and df.
fit <- fit |>
  mutate(
    p = coalesce(p, pchisq(`chi-square`, df, lower.tail = FALSE)),
    bartlett_df = n_vars * (n_vars - 1) / 2
  )

table_s1 <- with(fit, tribble(
  ~row, ~values,
  "Patients per pair, harmonic mean", if_else(is.na(`harmonic mean pairwise n`), "", sprintf("%.0f", `harmonic mean pairwise n`)),
  "Kaiser-Meyer-Olkin measure", f2(KMO),
  "Bartlett's test, χ² (df), p", if_else(is.na(`Bartlett chi-square`), "",
    sprintf("%.1f (%.0f), %s", `Bartlett chi-square`, bartlett_df, p_fmt(`Bartlett p`))),
  "Factors suggested by parallel analysis", if_else(is.na(`parallel analysis: factors`), "", sprintf("%.0f", `parallel analysis: factors`)),
  "One-factor model", rep("", 4),
  "  χ² (df), p", sprintf("%.1f (%.0f), %s", `chi-square`, df, p_fmt(p)),
  "  RMSEA [90% CI]", ci(RMSEA, `RMSEA 90% lower`, `RMSEA 90% upper`),
  "  TLI", f2(TLI),
  "  CFI", f2(CFI),
  "  RMSR (EFA) or SRMR (CFA)", f2(coalesce(RMSR, SRMR)),
  "  Variance explained", if_else(is.na(`variance explained`), "", sprintf("%.0f%%", 100 * `variance explained`)),
  "McDonald's ω total", f2(`omega total`),
  "Cronbach's α", f2(alpha),
  "Residual correlation, last dentist visit - X-ray [95% CI]",
    ci(`residual correlation visit-X-ray`, `residual correlation 95% lower`, `residual correlation 95% upper`)
)) |>
  unnest_wider(values, names_sep = "_")
names(table_s1) <- c("row", unname(sets))

write_csv(table_s1, file.path(out_dir, "art1_table_s1.csv"), na = "")
print(as.data.frame(table_s1), right = FALSE, row.names = FALSE)
cat("Written:", file.path(out_dir, "art1_table_s1.csv"), "\n")

# ---- Table S2 --------------------------------------------------------------

coefs <- read_csv(file.path(res_dir, "rq1_neglect_coefficients.csv"), show_col_types = FALSE) |>
  filter(model == "correlated")
unadjusted <- read_csv(file.path(res_dir, "rq1_unadjusted_comparison.csv"), show_col_types = FALSE) |>
  filter(model == "unadjusted")

# Ratios as in Figure 2; intervals with a negative bound use "to", not a hyphen.
num <- function(x) if_else(x >= 100, sprintf("%.0f", x), if_else(x >= 10, sprintf("%.1f", x), sprintf("%.2f", x)))
interval <- function(est, lo, hi) {
  if_else(lo < 0, sprintf("%s [%s to %s]", est, lo, hi), sprintf("%s [%s-%s]", est, lo, hi))
}
signed <- function(x) gsub("-", "−", sprintf("%.2f", x))
p_pos <- function(p) case_when(p > 0.995 ~ "> 0.99", p < 0.005 ~ "< 0.01", .default = sprintf("%.2f", p))

pick <- function(d, sc) {
  d |>
    filter(scale == sc) |>
    select(term, median, lower_95, upper_95, p_above_null)
}
sd_cell <- pick(coefs, "difference in latent neglect (SD)") |>
  transmute(term, beta = interval(signed(median), signed(lower_95), signed(upper_95)), p = p_pos(p_above_null))
pct_cell <- function(d) {
  pick(d, "probability more neglected (%)") |>
    transmute(term, value = sprintf("%.0f [%.0f-%.0f]", median, lower_95, upper_95))
}

# Effects: the four simple contrasts, the interaction (SD only), the one group
# contrast cited in Methods S1 (RTH seen after vs BF seen before), and the
# covariates. Age, sex and smoking are not in the visit-year-only model.
effect_rows <- c(
  "RTH vs BF, seen before treatment" = "RTH vs BF, seen before treatment",
  "RTH vs BF, seen after treatment"  = "RTH vs BF, seen after treatment",
  "After vs before, BF"              = "After vs before, BF",
  "After vs before, RTH"             = "After vs before, RTH",
  "Interaction: RTH × after"         = "Interaction, RTH × after",
  "RTH, seen after"                  = "RTH after vs BF before",
  "Age, per 10 years"                = "Age, per 10 years",
  "Male vs female"                   = "Male vs female",
  "Former vs never smoker"           = "Former vs never smoker",
  "Current vs never smoker"          = "Current vs never smoker",
  "Seen in 2025 vs 2024"             = "Seen in 2025 vs 2024",
  "Seen in 2026 vs 2024"             = "Seen in 2026 vs 2024"
)
effects <- tibble(term = names(effect_rows)) |>
  left_join(sd_cell, by = "term") |>
  left_join(pct_cell(coefs), by = "term") |>
  left_join(rename(pct_cell(unadjusted), visit_year_only = value), by = "term") |>
  mutate(term = paste0("  ", effect_rows[term]))
stopifnot(nrow(effects) == length(effect_rows), !anyNA(effects$beta))

# Measurement: per-SD ratio and standardised loading of each indicator, then
# the residual correlation.
indicators <- c(
  "Teeth to extract" = "Teeth to extract, RR", "Poor hygiene" = "Poor hygiene, OR",
  "Years since dentist" = "Years since last dentist visit, OR", "Old / no X-ray" = "Old or no dental X-ray, OR"
)
ratio <- coefs |>
  filter(scale %in% c("rate ratio per SD", "odds ratio per SD")) |>
  transmute(term, ratio = interval(num(median), num(lower_95), num(upper_95)))
std <- coefs |>
  filter(scale %in% c("standardised loading", "correlation")) |>
  transmute(term, std = interval(sprintf("%.2f", median), sprintf("%.2f", lower_95), sprintf("%.2f", upper_95)))
measurement <- tibble(term = c(names(indicators), "Residual correlation, years since dentist - X-ray")) |>
  left_join(ratio, by = "term") |>
  left_join(std, by = "term") |>
  mutate(term = paste0("  ", c(indicators, "Residual correlation, last dentist visit - X-ray")))
stopifnot(!anyNA(measurement$std), sum(is.na(measurement$ratio)) == 1)

na_blank <- function(x) coalesce(x, "")
table_s2 <- bind_rows(
  tibble(term = "Risk type × timing"),
  effects[1:6, ],
  tibble(term = "Covariates"),
  effects[7:12, ],
  tibble(term = "Measurement", beta = "RR or OR per SD of neglect [95% CrI]", p = "",
    value = "Standardised loading [95% CrI]", visit_year_only = ""),
  transmute(measurement, term, beta = ratio, value = std)
) |>
  select(term, beta, value, p, visit_year_only) |>
  mutate(across(everything(), na_blank))
names(table_s2) <- c(
  "Effect", "β, SD of neglect [95% CrI]", "Probability more neglected, % [95% CrI]", "P(β > 0)",
  "Adjusted for visit year only: probability more neglected, % [95% CrI]"
)

write_csv(table_s2, file.path(out_dir, "art1_table_s2.csv"), na = "")
print(as.data.frame(table_s2), right = FALSE, row.names = FALSE)
cat("Written:", file.path(out_dir, "art1_table_s2.csv"), "\n")

# Article 1 - Table 1: characteristics of the patients in the MIMIC model, by group.
#
# Run from anywhere inside the project:  Rscript R/question_1/16_art1_table1.R
#
# Same patients as R/question_1/05_rq1_neglect.R (exactly one of groups 1-4,
# complete covariates). One column per group plus a Total column. Percentages
# use the patients with data; "Not recorded" rows use all patients in the group.
#
# No significance tests (STROBE advises against them in descriptive tables).
# Imbalance is shown as standardised mean differences (SMD; Yang & Dalton 2012)
# for the confounders in the model and the comorbidities, with no threshold:
#   SMD, risk type - the larger of RTH vs BF before and RTH vs BF during/after;
#   SMD, timing    - the larger of during/after vs before within BF and within RTH.
# Cancer type and oncological treatment largely define the groups, and the
# indicators are the outcome (the MIMIC gives their adjusted effects), so
# neither gets an SMD.
#
# Output: results/art1_MIMIC/art1_table1.csv

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
})

in_file <- here("data", "processed", "clean_data.rds")
out_dir <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

group_levels <- c("BF before", "BF during/after", "RTH before", "RTH after")

# ---- Data ------------------------------------------------------------------

d <- readRDS(in_file) |>
  mutate(n_groups = group_bf_before + group_bf_after + group_rth_before + group_rth_after + group_other) |>
  filter(n_groups == 1, group_other == 0) |>
  drop_na(teeth_to_extract, oral_hygiene, age_at_visit, sex, smoking_status, visit_year) |>
  mutate(
    group = factor(case_when(
      group_bf_before == 1 ~ "BF before",
      group_bf_after == 1 ~ "BF during/after",
      group_rth_before == 1 ~ "RTH before",
      group_rth_after == 1 ~ "RTH after"
    ), group_levels),
    female = sex == "F",
    pack_years_ever = if_else(smoking_status %in% c("former", "current"), smoking_pack_years, NA_real_),
    visit_year = factor(visit_year),
    across(starts_with("comorb_"), \(x) x == 1)
  )

# ---- Table -----------------------------------------------------------------

# Cell formats. Percentages use the patients with data as the denominator.
n_pct <- function(x) sprintf("%d (%.0f%%)", sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
mean_sd <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
med_iqr <- function(x) {
  q <- round(quantile(x, c(0.5, 0.25, 0.75), na.rm = TRUE, names = FALSE), 1)
  sprintf("%g (%g-%g)", q[1], q[2], q[3])
}

# n (%) with value 1, one row per 0/1 column in `labels`.
flag_cells <- function(g, labels) {
  set_names(map_chr(names(labels), \(v) n_pct(g[[v]] == 1)), paste0(labels, ", n (%)"))
}

# n (%) per level of a factor, indented under its header row.
level_cells <- function(x, labels) {
  set_names(map_chr(names(labels), \(l) n_pct(x == l)), paste0("  ", labels))
}

cancer_labels <- c(
  cancer_head_neck = "Head-and-neck cancer", cancer_breast = "Breast cancer",
  cancer_other = "Other cancer"
)
prior_labels <- c(
  prior_chemotherapy = "Prior chemotherapy", prior_radiotherapy = "Prior radiotherapy (any site)",
  prior_immunotherapy = "Prior immunotherapy", prior_hormone_therapy = "Prior hormone therapy",
  prior_targeted_therapy = "Prior targeted therapy", prior_bf_denosumab = "Prior BF / denosumab"
)
comorb_labels <- c(
  comorb_diabetes = "Diabetes", comorb_hypertension = "Hypertension",
  comorb_thyroid = "Thyroid disease", comorb_cardiovascular = "Cardiovascular disease"
)

# Section and header rows are NA (blank in the CSV).
table_rows <- function(g) enframe(c(
  "Patient and visit" = NA,
  "Age, years, mean (SD)" = mean_sd(g$age_at_visit),
  "Female, n (%)" = n_pct(g$female),
  "Smoking status, n (%)" = NA,
  level_cells(g$smoking_status, c(never = "Never", former = "Former", current = "Current")),
  "Pack-years in ever smokers, median (IQR)" = med_iqr(g$pack_years_ever),
  "Visit year, n (%)" = NA,
  level_cells(g$visit_year, c("2024" = "2024", "2025" = "2025", "2026" = "2026")),
  "Cancer and oncological treatment" = NA,
  flag_cells(g, cancer_labels),
  "Two or more cancer types, n (%)" = n_pct(g$cancer_head_neck + g$cancer_breast + g$cancer_other >= 2),
  flag_cells(g, prior_labels),
  "Comorbidities" = NA,
  flag_cells(g, comorb_labels),
  "Neglect indicators" = NA,
  "Teeth to extract, median (IQR)" = med_iqr(g$teeth_to_extract),
  "Oral hygiene, n (%)" = NA,
  level_cells(g$oral_hygiene, c(good = "Good", average = "Average", poor = "Poor")),
  "Last dentist visit, n (%)" = NA,
  level_cells(g$last_dentist_visit, c("<1 year" = "< 1 year", "1-5 years" = "1-5 years", "5-10 years" = "5-10 years", never = "Never")),
  "  Not recorded" = n_pct(is.na(g$last_dentist_visit)),
  "Most recent dental X-ray, n (%)" = NA,
  level_cells(g$dental_xray, c(current = "Current", ">5 years" = "Older than 5 years", never = "Never")),
  "  Not recorded" = n_pct(is.na(g$dental_xray))
), "row", "value") |>
  mutate(id = row_number())

# Standardised mean difference between groups a and b: |difference in means| /
# pooled SD for numbers; for logicals and factors, the multivariate version over
# the category proportions, which reduces to the usual formula for a proportion.
smd <- function(x, a, b) {
  xa <- na.omit(x[d$group == a])
  xb <- na.omit(x[d$group == b])
  if (is.numeric(x)) return(abs(mean(xa) - mean(xb)) / sqrt((var(xa) + var(xb)) / 2))
  lv <- if (is.logical(x)) c(FALSE, TRUE) else levels(x)
  p <- as.vector(prop.table(table(factor(xa, lv))))[-1]
  q <- as.vector(prop.table(table(factor(xb, lv))))[-1]
  S <- (diag(p, length(p)) - outer(p, p) + diag(q, length(q)) - outer(q, q)) / 2
  sqrt(drop(t(p - q) %*% solve(S, p - q)))
}

smd_vars <- c(
  "Age, years, mean (SD)" = "age_at_visit",
  "Female, n (%)" = "female",
  "Smoking status, n (%)" = "smoking_status",
  "Pack-years in ever smokers, median (IQR)" = "pack_years_ever",
  "Visit year, n (%)" = "visit_year",
  set_names(names(comorb_labels), paste0(comorb_labels, ", n (%)"))
)
smds <- tibble(
  row = names(smd_vars),
  `SMD, risk type` = map_dbl(smd_vars, \(v) max(smd(d[[v]], "BF before", "RTH before"), smd(d[[v]], "BF during/after", "RTH after"))),
  `SMD, timing` = map_dbl(smd_vars, \(v) max(smd(d[[v]], "BF before", "BF during/after"), smd(d[[v]], "RTH before", "RTH after")))
) |>
  mutate(across(-row, \(x) sprintf("%.2f", x)))

col_names <- set_names(
  c(group_levels, "Total"),
  sprintf("%s (n = %d)", c(group_levels, "Total"), c(table(d$group), nrow(d)))
)
table1 <- d |>
  group_by(group) |>
  reframe(table_rows(pick(everything()))) |>
  mutate(group = as.character(group)) |>
  bind_rows(mutate(table_rows(d), group = "Total")) |>
  pivot_wider(names_from = group, values_from = value) |>
  rename(!!!col_names) |>
  left_join(smds, by = "row") |>
  select(-id)

write_csv(table1, file.path(out_dir, "art1_table1.csv"), na = "")

# ---- Console summary -------------------------------------------------------

cat("Patients:", nrow(d), "\n")
cat("Missing (for footnotes):\n")
print(colSums(is.na(select(d, all_of(names(comorb_labels)), last_dentist_visit, dental_xray))))
cat("Ever smokers without pack-years:", sum(d$smoking_status != "never" & is.na(d$smoking_pack_years)), "\n")
print(as.data.frame(table1), right = FALSE, row.names = FALSE)

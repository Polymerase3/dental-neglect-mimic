# dental-neglect-mimic

Code for the article on dental neglect at the first visit as a latent construct, estimated with a Bayesian MIMIC model.

- `R/` — analysis scripts
- `data/raw/` — original, untouched data (not tracked by git)
- `data/processed/` — cleaned/derived data (not tracked by git)
- `data/private/` — name → patient ID key (not tracked by git)
- `results/` — generated tables and figures (not tracked by git)

The patient data are not public, so the scripts need the raw sheet in `data/raw/` to run.

## Cleaning the data

```sh
Rscript R/01_clean_data.R
```

Reads the raw sheet from `data/raw/` and writes `data/processed/clean_data.csv` and `.rds`.

## Analysis

```sh
Rscript R/question_1/04_rq1_efa.R                  # factor analysis of the dental-neglect indicators
Rscript R/question_1/05_rq1_neglect.R              # latent dental-neglect model (Stan), checks, LOO
Rscript R/question_1/06_rq1_neglect_path.R         # path diagram of that model, two versions, SVG + 600 DPI PNG (needs Inkscape)
Rscript R/question_1/12_rq1_neglect_no_hygiene.R   # sensitivity analysis: the latent model without hygiene, compared with 05 (~1 min)
Rscript R/question_1/13_art1_flow_diagram.R        # participant flow (240 to 175): SVG, PDF, 600 DPI PNG (needs Inkscape)
Rscript R/question_1/14_rq1_candidate_effects.R    # other patient characteristics, each added on its own to the model from 05
Rscript R/question_1/15_rq1_prior_sensitivity.R    # prior sensitivity (BARG): the latent model refitted with every prior SD halved and doubled, compared with 05
Rscript R/question_1/16_art1_table1.R              # Table 1: characteristics of the patients in the model, by group
```

Run 05 before 06, 12 and 15, which read its coefficients. Outputs go to `results/question_1/` with the prefix `rq1_`; the flow diagram and Table 1 go to `results/art1_MIMIC/`.

R packages: bayesplot, dplyr, ggplot2, here, lavaan, loo, patchwork, posterior, psych, purrr, readr, rstan, stringr, tibble, tidyr.

- [CODEBOOK.md](CODEBOOK.md) — what every variable means

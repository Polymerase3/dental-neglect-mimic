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
Rscript R/question_1/17_art1_figure_s2.R           # Figure S2: posterior predictive checks, no-hygiene and prior sensitivity
Rscript R/question_1/18_art1_figure2_path.R        # Figure 2: path diagram of the reported model, SVG + PDF + 600 DPI PNG (needs Inkscape)
Rscript R/question_1/19_art1_figure1.R             # Figure 1: selection of the indicators (correlations, scree, loadings)
Rscript R/question_1/20_rq1_neglect_unadjusted.R   # the latent model adjusted for visit year only, compared with 05
Rscript R/question_1/21_art1_tables_s1_s2.R        # Tables S1 and S2: factor-analysis fit, posterior summaries of the reported model
Rscript R/question_1/22_art1_forest_figures.R      # Figure 3: effects on neglect, adjusted and adjusted for visit year only
Rscript R/question_1/23_art1_figure_s3.R           # Figure S3: effects of the other patient characteristics, one model each
```

Run 04 before 19 and 21; 05 before 06, 12, 14, 15, 17, 18 and 20, which read its fit or coefficients; 12 and 15 before 17; 20 before 21 and 22; 14 before 23. Outputs go to `results/question_1/` with the prefix `rq1_`; the flow diagram and the article tables and figures go to `results/art1_MIMIC/`.

R packages: bayesplot, cowplot, dplyr, ggplot2, here, lavaan, loo, patchwork, posterior, psych, purrr, ragg, readr, rstan, statmod, stringr, systemfonts, tibble, tidyr.

- [CODEBOOK.md](CODEBOOK.md) — what every variable means

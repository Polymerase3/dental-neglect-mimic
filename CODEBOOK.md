# Codebook — `data/processed/clean_data.csv` / `.rds`

One row per patient. Produced by `R/01_clean_data.R`; the cleaning rules are in [DECISIONS.md](DECISIONS.md).

**Conventions**
- **0/1** columns are integers: `1` = yes, `0` = no, empty/`NA` = missing.
- **Multi-choice** source columns (e.g. "2+4") are split into one 0/1 column per category, so a patient can have several 1s.
- **Factor** columns hold English labels; the `.rds` file keeps the level order shown here.
- "Source" is the column header in the raw sheet (Polish, sometimes shortened).

**Abbreviations used in the source sheet:** BF = bisphosphonates (bisfosfoniany); RTH = radiotherapy; GIS = head and neck (głowa i szyja); CHT = chemotherapy; HP = histopathology; CBCT = cone-beam CT; ORN = osteoradionecrosis; MRONJ = medication-related osteonecrosis of the jaw.

## Identification and demographics

| Variable | Type | Values | Source | Meaning |
|---|---|---|---|---|
| `patient_id` | text | `P001`, `P002`, … | Imię i nazwisko | Pseudonymous patient ID. The name → ID key is kept only in `data/private/patient_ids.csv` (not in git). |
| `age_at_visit` | numeric | years | rok urodzenia, data pierwszej wizyty | Year of first visit minus year of birth (accurate to ±1 year). |
| `visit_year` | integer | year | data pierwszej wizyty | Year of the first visit. The dental history (`last_dentist_visit`, `dental_xray`, `visit_with_pain`) became routine in July 2025 (`history_routine`), so the year explains most of the missing values there. |
| `history_routine` | 0/1 | | data pierwszej wizyty | First visit in or after July 2025, from when pain, last dentist visit and dental X-ray were recorded for every patient. Before that, pain was written down only when present (15 yes, 0 no, 79 blank). |
| `sex` | factor | `F`, `M` | płeć | Sex (K → F). |
| `referring_centre` | factor | `UCK`, `NIO`, `KCO`, `BCO`, `other` | ośrodek kierujący | Centre that referred the patient. |

## Patient group

Source: *grupa 1 - BF przed, 2- BF po, 3 - przed RTH 4 - po RTH 5 - inne*. A patient can belong to more than one group.

| Variable | Type | Meaning |
|---|---|---|
| `group_bf_before` | 0/1 | Seen before starting bisphosphonate/denosumab therapy (1) |
| `group_bf_after` | 0/1 | Seen during or after bisphosphonate/denosumab therapy (2) |
| `group_rth_before` | 0/1 | Seen before head-and-neck radiotherapy (3) |
| `group_rth_after` | 0/1 | Seen after head-and-neck radiotherapy (4) |
| `group_other` | 0/1 | Other (5) |

## Smoking

| Variable | Type | Values | Source | Meaning |
|---|---|---|---|---|
| `smoking_status` | factor | `never`, `former`, `current` | nikotynizm P1/P2/P3 | P3 = never, P2 = former, P1 = current smoker. |
| `smoking_years_since_quit` | numeric | years | nikotynizm (note after P2) | Years since a former smoker quit. Only for `former`. |
| `smoking_pack_years` | numeric | pack-years | lb paczko-lat | Pack-years. |

## Comorbidities

Source: *1 - cukrzyca, 2 - nadciśnienie 3 - tarczyca 4- choroby sercowo-naczyniowe* and *brak schorzen (1- brak)*. Other conditions listed only in the free-text column are not coded.

| Variable | Type | Meaning |
|---|---|---|
| `comorb_diabetes` | 0/1 | Diabetes (1) |
| `comorb_hypertension` | 0/1 | Hypertension (2) |
| `comorb_thyroid` | 0/1 | Thyroid disease (3) |
| `comorb_cardiovascular` | 0/1 | Cardiovascular disease (4) |
| `comorb_none` | 0/1 | No chronic disease at all (*brak schorzeń* = 1) |

## Cancer

Source for type: *1(glowa i szyja) 2) pierś 3) inne*. A patient can have more than one cancer.

| Variable | Type | Meaning |
|---|---|---|
| `cancer_head_neck` | 0/1 | Head-and-neck cancer (1) |
| `cancer_breast` | 0/1 | Breast cancer (2) |
| `cancer_other` | 0/1 | Other cancer (3) |

Head-and-neck primary site. Source: *1.gardlo srodkowe … 10) slinianka*. 0 for patients without a head-and-neck cancer.

| Variable | Type | Meaning |
|---|---|---|
| `hn_site_oropharynx` | 0/1 | Oropharynx — tonsil, base of tongue, pharynx (1) |
| `hn_site_larynx` | 0/1 | Larynx (2) |
| `hn_site_palate` | 0/1 | Palate (3) |
| `hn_site_maxilla` | 0/1 | Maxillary complex (4) |
| `hn_site_neck_nodes` | 0/1 | Cervical lymph nodes (5) |
| `hn_site_tongue` | 0/1 | Tongue (6) |
| `hn_site_cheek` | 0/1 | Cheek (7) |
| `hn_site_floor_of_mouth` | 0/1 | Floor of mouth (8) |
| `hn_site_retromolar_trigone` | 0/1 | Retromolar trigone (9) |
| `hn_site_salivary_gland` | 0/1 | Salivary gland (10) |
| `hn_site_other` | 0/1 | Other site ("inne") |

Head-and-neck region. Source: *1- jama ustna 2- gardlo 3- krtań 4 - szyja*. 0 for patients without a head-and-neck cancer.

| Variable | Type | Meaning |
|---|---|---|
| `hn_region_oral_cavity` | 0/1 | Oral cavity (1) |
| `hn_region_pharynx` | 0/1 | Pharynx (2) |
| `hn_region_larynx` | 0/1 | Larynx (3) |
| `hn_region_neck` | 0/1 | Neck (4) |
| `hn_region_maxilla` | 0/1 | Maxilla ("szczeka") |

## Oncological treatment received by the time of the visit

Source: the columns *CHT (czy było przebyte leczenie w momencie wizyty)*, *RTH*, *IMMUNO*, *HORMONO*, *TER CELOWANA*, *BISFOSFONIANY DENOSUMAB*.

| Variable | Type | Meaning |
|---|---|---|
| `prior_chemotherapy` | 0/1 | Chemotherapy |
| `prior_radiotherapy` | 0/1 | Radiotherapy (any site) |
| `prior_immunotherapy` | 0/1 | Immunotherapy |
| `prior_hormone_therapy` | 0/1 | Hormone therapy |
| `prior_targeted_therapy` | 0/1 | Targeted therapy |
| `prior_bf_denosumab` | 0/1 | Bisphosphonates or denosumab |

## Timing of the first dental consultation

Source: *kiedy odbyła się pierwsza konsultacja u mnie - 1 - przed leczeniem (…), 2 - w trakcie, 3 - po zakończonym*.

| Variable | Type | Meaning |
|---|---|---|
| `first_consult_before_tx` | 0/1 | Before oncological treatment, including after surgery or after the histopathology result (1) |
| `first_consult_during_tx` | 0/1 | During treatment (2) |
| `first_consult_after_tx` | 0/1 | After treatment had ended (3) |

## Head-and-neck radiotherapy

| Variable | Type | Values | Source | Meaning |
|---|---|---|---|---|
| `hn_rth_before` | 0/1 | | przed RTH GIS | Visit before head-and-neck radiotherapy. |
| `hn_rth_before_days` | numeric | days | ile dni | Days between the visit and the start of radiotherapy. |
| `hn_rth_after` | 0/1 | | po RTH GIS | Visit after head-and-neck radiotherapy. |
| `hn_rth_after_months` | numeric | months | ile miesięcy | Months since radiotherapy ended. |

## Bisphosphonates / denosumab

| Variable | Type | Values | Source | Meaning |
|---|---|---|---|---|
| `bf_before_prophylactic` | 0/1 | | BF przed (Profilaktycznie) | Seen before planned prophylactic BF. |
| `bf_before_therapeutic` | 0/1 | | BF przed (Leczniczo - np. przerzuty do kości) | Seen before planned therapeutic BF (e.g. bone metastases). |
| `bf_doses` | numeric | count | ile dawek | Number of BF doses received. |
| `bf_months` | numeric | months | ile m-cy BF | Months on BF. |
| `bf_osteoporosis` | 0/1 | | BF NA OSTEOPOROZĘ | BF taken for osteoporosis. |

## Dental status at the visit

| Variable | Type | Values | Source | Meaning |
|---|---|---|---|---|
| `visit_with_pain` | 0/1 | | wizyta z bólem | Patient came with pain. |
| `dental_assessment_before_tx` | factor | `none`, `adequate`, `inadequate` | ocena stom przed leczeniem - 0(brak), 1, 2 (nieprawidłowa) | Dental assessment before oncological treatment: none (0), done (1), done inadequately (2). |
| `last_dentist_visit` | factor | `<1 year`, `1-5 years`, `5-10 years`, `never` | stomatolog - (1 - do roku, 2 - 1-5l, 3 - 5-10l, 4 - nigdy) | Time since the last dentist visit. |
| `dental_xray` | factor | `current`, `>5 years`, `never` | RTG (1-aktualne, 2 - >5 lat, 3 nigdy) | Most recent dental X-ray. |
| `oral_hygiene` | factor | `poor`, `average`, `good` | H 1-zla higiena, 2 - przeciętna, 3 dobra | Oral hygiene. |
| `teeth_to_treat` | numeric | count | L - zeby do leczenia | Teeth needing treatment. |
| `teeth_to_extract` | numeric | count | E - zeby zakwalifikowane do ekstrakcji | Teeth qualified for extraction (0 also for edentulous patients). |
| `extraction_during_onco_tx` | 0/1 | | E w trakcie leczenia onko (lub po BF i RTH GIS) | Extraction needed during oncological treatment (or after BF / head-and-neck RTH). |
| `months_since_onco_tx` | numeric | months | ile m-cy od leczenia onko a ekstrakcja LUB WIZYTA U MNIE | Months from oncological treatment to the extraction or visit. 0 = after histopathology, before RTH. |
| `referred_maxillofacial` | 0/1 | | skierowano dalej do Poradni Chirurgii Szczękowo-Twarzowej | Referred on to a maxillofacial surgery clinic. |
| `cbct_needed` | 0/1 | | konieczność CBCT | CBCT was indicated (not necessarily performed). |

## Upper denture

| Variable | Type | Values | Source | Meaning |
|---|---|---|---|---|
| `upper_denture` | 0/1 | | proteza gorna | Wears an upper denture. |
| `upper_denture_years` | numeric | years | ilosc lat | Years wearing the denture. |
| `denture_night_use` | 0/1 | | użytkowanie w nocy 1- całodobowo, 0- przerwa nocna | 1 = worn day and night, 0 = removed at night. |

## Oral mucosa

| Variable | Type | Values | Source | Meaning |
|---|---|---|---|---|
| `oral_lesions` | 0/1 | | zmiany w jamie ustnej | Any mucosal lesion found (fungal infection, erythema, leukoplakia, denture stomatitis, …). |

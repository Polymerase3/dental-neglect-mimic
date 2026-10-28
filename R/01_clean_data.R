# Clean the raw patient sheet into an analysis-ready dataset.
#
# Run from anywhere inside the project:  Rscript R/01_clean_data.R
#
# Every rule applied here is explained in DECISIONS.md, and every output
# variable is described in CODEBOOK.md. Keep the three files in sync.

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(purrr)
  library(stringr)
})

raw_file <- here("data", "raw", "Kopia Kopia onkologiczny pacjenci NFZ - Arkusz1.csv")
id_file  <- here("data", "private", "patient_ids.csv")
out_dir  <- here("data", "processed")

# Raw columns kept for cleaning (short working name = original header).
# Columns not listed here (free text, "nie do analizy", unnamed) are dropped.
source_cols <- c(
  name                 = "Imię i nazwisko",
  group                = "grupa 1 - BF przed, 2- BF po, 3 - przed RTH 4 - po RTH 5 - inne",
  centre               = "ośrodek kierujący",
  visit_date           = "data pierwszej wizyty w gabinecie",
  birth_year           = "rok urodzenia",
  sex                  = "płeć",
  smoking              = "nikotynizm P1 - tak /P2 - w przeszłości (ile nie pali) /P3 - nie",
  pack_years           = "lb paczko-lat",
  comorb               = "1 - cukrzyca, 2 - nadciśnienie 3 - tarczyca 4- choroby sercowo-naczyniowe",
  comorb_none          = "brak schorzen (1- brak)",
  cancer               = "1(glowa i szyja) 2) pierś 3) inne",
  hn_site              = "1.gardlo srodkowe (migdalek, nasada jez, gardlo, ), 2 - krtań 3 - podniebienie 4 - masyw szczekowy 5) wezly chlonne szyi 6) jezyk 7) policzek 8) dno 9) trojkat 10) slinianka",
  hn_region            = "1- jama ustna 2- gardlo 3- krtań 4 - szyja",
  chemo                = "CHT (czy było przebyte leczenie w momencie wizyty)",
  radio                = "RTH",
  immuno               = "IMMUNO",
  hormone              = "HORMONO",
  targeted             = "TER CELOWANA",
  bf_denosumab         = "BISFOSFONIANY DENOSUMAB",
  first_consult        = "kiedy odbyła się pierwsza konsultacja u mnie - 1 - przed leczeniem (w tym po leczeniu chirurgicznym, wyniku badania HP), 2 - w trakcie, 3 - po zakończonym",
  rth_before           = "przed RTH GIS",
  rth_before_days      = "ile dni",
  rth_after            = "po RTH GIS",
  rth_after_months     = "ile miesięcy",
  bf_prophylactic      = "BF przed (Profilaktyktycznie)",
  bf_therapeutic       = "BF przed (Leczniczo - np. przerzuty do kości)",
  bf_doses             = "ile dawek",
  bf_months            = "ile m-cy BF",
  bf_osteoporosis      = "BF NA OSTEOPOROZĘ",
  pain                 = "wizyta z bólem",
  dental_assessment    = "ocena stom przed leczeniem - 0(brak) , 1, 2 (nieprawidłowa)",
  last_dentist         = "stomatolog - (1 - do roku,2 - 1-5l,3 - 5-10l ,4 - nigdy)",
  xray                 = "RTG (1-aktualne, 2 - >5 lat, 3 nigdy)",
  hygiene              = "H 1-zla higiena, 2 - przeciętna, 3 dobra)",
  teeth_treat          = "L - zeby do leczenia",
  teeth_extract        = "E - zeby zakwalifikowane do ekstrakcji",
  extraction_during_tx = "E w trakcie leczenia onko (lub po BF i RTH GIS)",
  months_since_tx      = "ile m-cy od leczenia onko a ekstrakcja LUB WIZYTA U MNIE ? Jak 0/1 to pacjent po histpacie - przed RTH",
  referred_mfs         = "skierowano dalej do Poradni Chirurgii Szczękowo-Twarzowej",
  cbct                 = "konieczność CBCT (ale nie u wszystkich wykonano?)",
  denture              = "proteza gorna",
  denture_years        = "ilosc lat",
  denture_night        = "użytkowanie w nocy 1- całodobowo, 0- przerwa nocna",
  oral_lesions         = "zmiany w jamie ustnej"
)

# ---- Helpers ---------------------------------------------------------------

# Leading number of a cell: "1 (ja)" -> 1, "7,5" -> 7.5, "3!" -> 3.
# Missing if the cell is empty or text, or the number is marked uncertain or
# ambiguous ("1?", "0/1", "1/2").
parse_num <- function(x) {
  m <- str_match(x, "^`?(\\d+(?:[.,]\\d+)?)\\s*([?/])?")
  as.numeric(if_else(is.na(m[, 3]), str_replace(m[, 2], ",", "."), NA_character_))
}

# Single-choice code column -> factor with the given labels; other codes -> NA.
code_factor <- function(x, labels) {
  factor(unname(labels[as.character(parse_num(x))]), levels = unname(labels))
}

# 0/1 column where "no" is written as 0: empty stays missing.
parse_binary <- function(x) {
  v <- parse_num(x)
  as.integer(if_else(v %in% c(0, 1), v, NA_real_))
}

# 1-or-blank flag column: blank -> 0; codes as written ("?" -> NA); a drug
# name written instead of the code -> 1.
parse_flag <- function(x) {
  case_when(
    x == "" ~ 0L,
    str_detect(x, "^[`\\d?]") ~ as.integer(parse_num(x) != 0),
    TRUE ~ 1L
  )
}

# Multi-choice cell ("2+3", "1,2", "2 (osteop)+4") -> one 0/1 column per code.
# `codes` maps each code (number or word) to its output column. Notes in
# brackets are ignored and repeated codes count once. A cell with no known
# code ("?", free text) -> all NA; an empty cell -> `if_empty`.
make_dummies <- function(x, codes, if_empty) {
  x <- str_to_lower(x)
  tokens <- str_extract_all(str_remove_all(x, "\\([^)]*\\)"), "\\d+|\\p{L}+")
  known <- map_lgl(tokens, \(t) any(t %in% names(codes)))
  map(names(codes), \(code) {
    v <- as.integer(map_lgl(tokens, \(t) code %in% t))
    case_when(x == "" ~ as.integer(if_empty), !known ~ NA_integer_, TRUE ~ v)
  }) |>
    set_names(unname(codes)) |>
    as_tibble()
}

# Duration in months. Bare numbers are months, "l"/"lat"/"lata" years, "dni"
# days; "0/1"-style codes take the first number. Several values joined by "+"
# or "i" (separate treatments) -> the smallest, i.e. the most recent one.
# Anything marked "?" -> NA.
parse_months <- function(x) {
  map_dbl(str_to_lower(x), \(v) {
    if (v == "" || str_detect(v, fixed("?"))) return(NA_real_)
    parts <- str_split_1(v, "\\s*(\\+|\\bi\\b)\\s*")
    n <- as.numeric(str_replace(str_extract(parts, "^\\d+(?:,\\d+)?"), ",", "."))
    unit <- case_when(
      str_detect(parts, "^\\d+\\s*l") ~ 12,
      str_detect(parts, "dni") ~ 1 / 30.44,
      TRUE ~ 1
    )
    min(n * unit)
  })
}

# Year of a visit date written in any of the sheet's formats: a 4-digit year
# if present ("2025.05.12", "maj 2024"), else the last two digits
# ("22.09.25", "20.1025"). Notes in brackets are ignored.
parse_visit_year <- function(x) {
  x <- str_remove_all(x, "\\([^)]*\\)")
  y4 <- as.integer(str_extract(x, "\\b20\\d{2}\\b"))
  y2 <- as.integer(str_match(x, "(\\d{2})\\D*$")[, 2])
  coalesce(y4, 2000L + y2)
}

# Month of a visit date, same formats: the number after a 4-digit year
# ("2025.05.12", "2025.08."), the number before one ("10.2025"), the second
# number of a day.month.year date ("12.05.25", "20.1025"), or a Polish month
# name ("maj 2024").
parse_visit_month <- function(x) {
  x <- str_to_lower(str_remove_all(x, "\\([^)]*\\)"))
  months <- c(sty = 1, lut = 2, mar = 3, kwi = 4, maj = 5, cze = 6, lip = 7, sie = 8, wrz = 9, paz = 10, lis = 11, gru = 12)
  coalesce(
    as.integer(str_match(x, "\\b20\\d{2}\\.(\\d{1,2})")[, 2]),
    as.integer(str_match(x, "^(\\d{1,2})\\.20\\d{2}\\b")[, 2]),
    as.integer(str_match(x, "^\\d{1,2}\\.(\\d{1,2})\\.?\\d{2}")[, 2]),
    as.integer(months[str_sub(str_replace(str_extract(x, "\\p{L}+"), "ź", "z"), 1, 3)])
  )
}

# Smoking: P1/P2/P3 -> current/former/never; notes dropped.
parse_smoking <- function(x) {
  code_factor(str_extract(x, "(?<=^P)[123]"), c("3" = "never", "2" = "former", "1" = "current"))
}

# Years since quitting, from the note after P2: "P2 (15l)" -> 15,
# "P2 (2m-ce)" -> 2/12. A bare number is years.
parse_years_since_quit <- function(x) {
  m <- str_match(str_to_lower(x), "^p2\\s*\\(\\s*(\\d+(?:,\\d+)?)\\s*(\\S*)")
  as.numeric(str_replace(m[, 2], ",", ".")) / if_else(str_starts(m[, 3], "m"), 12, 1)
}

# Referring centre -> UCK / NIO / KCO / BCO / other.
parse_centre <- function(x) {
  x <- str_to_upper(x)
  case_when(
    x == "" | str_starts(x, "NIE ") ~ NA_character_,
    str_starts(x, "UCK") ~ "UCK",
    x %in% c("NIO", "KCO", "BCO") ~ x,
    TRUE ~ "other"
  ) |>
    factor(levels = c("UCK", "NIO", "KCO", "BCO", "other"))
}

# Oral lesions: "0"/"brak" -> 0, a described lesion -> 1, denture-storage
# notes (not lesions) -> NA.
parse_lesions <- function(x) {
  x <- str_to_lower(x)
  case_when(
    x == "" ~ NA_integer_,
    x %in% c("0", "brak", "brak zmian") ~ 0L,
    str_detect(x, "przechowuje|corega") ~ NA_integer_,
    TRUE ~ 1L
  )
}

# Stable patient IDs: names already in the private key keep their ID, new
# names get the next number. The key is the only place names are kept.
assign_ids <- function(names, key_file) {
  key <- if (file.exists(key_file)) {
    read_csv(key_file, col_types = "cc")
  } else {
    tibble(name = character(), patient_id = character())
  }
  new <- setdiff(unique(names), key$name)
  key <- bind_rows(key, tibble(name = new, patient_id = sprintf("P%03d", nrow(key) + seq_along(new))))
  dir.create(dirname(key_file), showWarnings = FALSE)
  write_csv(key, key_file)
  key$patient_id[match(names, key$name)]
}

# ---- Read ------------------------------------------------------------------

raw <- read_csv(
  raw_file,
  col_types = cols(.default = col_character()),
  na = character(),
  name_repair = "unique_quiet"
)
names(raw) <- str_squish(names(raw))
raw <- mutate(raw, across(everything(), str_squish))

# Patients only: drop blank rows and the summary notes typed into the name
# column below the data (a "name" with no other cell filled).
raw <- filter(raw, `Imię i nazwisko` != "", if_any(-`Imię i nazwisko`, \(x) x != ""))

src <- select(raw, all_of(source_cols))

# ---- Clean -----------------------------------------------------------------

clean <- src |>
  transmute(
    patient_id = assign_ids(name, id_file),
    age_at_visit = parse_visit_year(visit_date) - parse_num(birth_year),
    visit_year = parse_visit_year(visit_date),
    # Pain, last dentist visit and X-ray are recorded for every patient from
    # July 2025; only this flag is kept, not the month.
    history_routine = as.integer(visit_year > 2025 | (visit_year == 2025 & parse_visit_month(visit_date) >= 7)),
    sex = factor(unname(c(K = "F", M = "M")[sex]), levels = c("F", "M")),
    referring_centre = parse_centre(centre),
    make_dummies(group, c(
      "1" = "group_bf_before", "2" = "group_bf_after", "3" = "group_rth_before",
      "4" = "group_rth_after", "5" = "group_other"
    ), if_empty = NA),

    smoking_status = parse_smoking(smoking),
    smoking_years_since_quit = parse_years_since_quit(smoking),
    smoking_pack_years = parse_num(pack_years),

    make_dummies(comorb, c(
      "1" = "comorb_diabetes", "2" = "comorb_hypertension",
      "3" = "comorb_thyroid", "4" = "comorb_cardiovascular"
    ), if_empty = 0),
    comorb_none = parse_flag(comorb_none),

    make_dummies(cancer, c(
      "1" = "cancer_head_neck", "2" = "cancer_breast", "3" = "cancer_other"
    ), if_empty = NA),
    make_dummies(hn_site, c(
      "1" = "hn_site_oropharynx", "2" = "hn_site_larynx", "3" = "hn_site_palate",
      "4" = "hn_site_maxilla", "5" = "hn_site_neck_nodes", "6" = "hn_site_tongue",
      "7" = "hn_site_cheek", "8" = "hn_site_floor_of_mouth",
      "9" = "hn_site_retromolar_trigone", "10" = "hn_site_salivary_gland",
      "inne" = "hn_site_other"
    ), if_empty = 0),
    make_dummies(hn_region, c(
      "1" = "hn_region_oral_cavity", "2" = "hn_region_pharynx",
      "3" = "hn_region_larynx", "4" = "hn_region_neck", "szczeka" = "hn_region_maxilla"
    ), if_empty = 0),

    prior_chemotherapy = parse_flag(chemo),
    prior_radiotherapy = parse_flag(radio),
    prior_immunotherapy = parse_flag(immuno),
    prior_hormone_therapy = parse_flag(hormone),
    prior_targeted_therapy = parse_flag(targeted),
    prior_bf_denosumab = parse_flag(bf_denosumab),

    make_dummies(first_consult, c(
      "1" = "first_consult_before_tx", "2" = "first_consult_during_tx",
      "3" = "first_consult_after_tx"
    ), if_empty = NA),

    hn_rth_before = parse_binary(rth_before),
    hn_rth_before_days = parse_num(rth_before_days),
    hn_rth_after = parse_binary(rth_after),
    hn_rth_after_months = parse_months(rth_after_months),

    bf_before_prophylactic = parse_binary(bf_prophylactic),
    bf_before_therapeutic = parse_binary(bf_therapeutic),
    bf_doses = parse_num(bf_doses),
    bf_months = parse_num(bf_months),
    bf_osteoporosis = parse_flag(bf_osteoporosis),

    visit_with_pain = parse_binary(pain),
    dental_assessment_before_tx = code_factor(dental_assessment, c(
      "0" = "none", "1" = "adequate", "2" = "inadequate"
    )),
    last_dentist_visit = code_factor(last_dentist, c(
      "1" = "<1 year", "2" = "1-5 years", "3" = "5-10 years", "4" = "never"
    )),
    dental_xray = code_factor(xray, c("1" = "current", "2" = ">5 years", "3" = "never")),
    oral_hygiene = code_factor(hygiene, c("1" = "poor", "2" = "average", "3" = "good")),
    teeth_to_treat = parse_num(teeth_treat),
    teeth_to_extract = parse_num(teeth_extract),
    extraction_during_onco_tx = parse_binary(extraction_during_tx),
    months_since_onco_tx = parse_months(months_since_tx),
    referred_maxillofacial = parse_binary(referred_mfs),
    cbct_needed = parse_binary(cbct),

    upper_denture = parse_binary(denture),
    upper_denture_years = parse_num(denture_years),
    denture_night_use = parse_binary(denture_night),
    oral_lesions = parse_lesions(oral_lesions)
  )

# ---- Write -----------------------------------------------------------------

write_csv(clean, file.path(out_dir, "clean_data.csv"), na = "")
saveRDS(clean, file.path(out_dir, "clean_data.rds"))

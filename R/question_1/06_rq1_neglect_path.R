# RQ1 - path diagram of the latent dental-neglect model (MIMIC).
#
# Run from anywhere inside the project:  Rscript R/question_1/06_rq1_neglect_path.R
#
# Reads the posterior summaries written by R/question_1/05_rq1_neglect.R (no refitting)
# and draws the reported model as a hand-built SVG in standard SEM notation:
# observed variables as rectangles, the latent construct as an ellipse,
# straight paths, and a figure legend. Inkscape renders a 600 DPI PNG and an
# SVG with the text converted to outlines, so both look the same on machines
# without the font (Arial).
#
# Two versions share the same drawing code:
#   rq1_neglect_path     title, caption and colour keys below the diagram
#   rq1_neglect_path_v2  no title or caption; colour keys in a side column
#
# Colour encodes effect size and credibility:
#   covariate -> construct  probability more neglected: red = more, blue =
#                           less; light to dark = small, moderate, large
#                           (Cohen's d 0.2 / 0.5 / 0.8 as probabilities);
#                           grey = the 95% credible interval includes 50%
#   construct -> indicator  standardised loading, teal light to dark = weak
#                           (< 0.5), moderate (0.5-0.7), strong (>= 0.7);
#                           grey = the OR / RR interval includes 1

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(purrr)
})

in_file <- here("results", "question_1", "rq1_neglect_coefficients.csv")
out_dir <- here("results", "question_1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png_dpi <- 600

# ---- Palette and type ------------------------------------------------------

# Dataviz reference palette; the three ramps pass the ordinal checks
# (monotone lightness, single hue, light end >= 2:1 on the surface).
surface      <- "#fcfcfb"
paper        <- "#ffffff"
ink          <- "#0b0b0b"
ink_soft     <- "#52514e"
rule         <- "#d9d7d0"
separator    <- "#e0ded7"
not_credible <- "#8f8d85"
more    <- c(small = "#f0918e", moderate = "#e05250", large = "#b3302f")
less    <- c(small = "#86b6ef", moderate = "#3987e5", large = "#1c5cab")
loading <- c(weak = "#5dbb96", moderate = "#1f9a6e", strong = "#0b6e4c")
serif   <- "Arial, Helvetica, sans-serif"  # one generic family throughout
sans    <- serif

# Probabilities equivalent to Cohen's d 0.2 / 0.5 / 0.8: Phi(d / sqrt(2)).
ps_cut <- 100 * pnorm(c(0.2, 0.5, 0.8) / sqrt(2))

# ---- Estimates -------------------------------------------------------------

coef <- read_csv(in_file, show_col_types = FALSE) |>
  filter(model == "correlated")
row_of <- function(term, scale) filter(coef, term == !!term, scale == !!scale)

ps_style <- function(r) {
  if (r$lower_95 <= 50 && r$upper_95 >= 50) return(list(col = not_credible, credible = FALSE))
  dist <- abs(r$median - 50)
  size <- if (dist >= ps_cut[3] - 50) "large" else if (dist >= ps_cut[2] - 50) "moderate" else "small"
  list(col = if (r$median > 50) more[[size]] else less[[size]], credible = TRUE)
}
load_style <- function(std, ratio) {
  if (ratio$lower_95 <= 1) return(list(col = not_credible, credible = FALSE, size = "not credible"))
  size <- if (std$median >= 0.7) "strong" else if (std$median >= 0.5) "moderate" else "weak"
  list(col = loading[[size]], credible = TRUE, size = size)
}

num <- function(x) ifelse(x >= 100, sprintf("%.0f", x), ifelse(x >= 10, sprintf("%.1f", x), sprintf("%.2f", x)))
pct_ci <- function(r) sprintf("%.0f%% [%.0f–%.0f]", r$median, r$lower_95, r$upper_95)
ratio_ci <- function(r, name) sprintf("%s %s [%s–%s]", name, num(r$median), num(r$lower_95), num(r$upper_95))

groups <- list(
  list(title = "Patient group", ref = "(vs BF, seen before)", rows = c(
    "BF, seen after" = "BF, seen after",
    "RTH, seen before" = "RTH, seen before",
    "RTH, seen after" = "RTH, seen after"
  )),
  list(title = "Age and sex", ref = NULL, rows = c(
    "Age, per 10 years" = "10 years older",
    "Male vs female" = "Male (vs female)"
  )),
  list(title = "Smoking", ref = "(vs never)", rows = c(
    "Former vs never smoker" = "Former smoker",
    "Current vs never smoker" = "Current smoker"
  )),
  list(title = "Visit year", ref = "(vs 2024)", rows = c(
    "Seen in 2025 vs 2024" = "Seen in 2025",
    "Seen in 2026 vs 2024" = "Seen in 2026"
  ))
)
indicators <- tibble(
  term = c("Teeth to extract", "Poor hygiene", "Years since dentist", "Old / no X-ray"),
  ratio_scale = c("rate ratio per SD", "odds ratio per SD", "odds ratio per SD", "odds ratio per SD"),
  ratio_name = c("RR", "OR", "OR", "OR")
)

# ---- SVG primitives --------------------------------------------------------

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
# `label` is markup: escape plain strings with esc() or build tspans.
text_el <- function(x, y, label, size = 14, weight = 400, fill = ink, anchor = "start",
                    family = sans, italic = FALSE, spacing = 0) {
  sprintf('<text x="%.1f" y="%.1f" font-family="%s" font-size="%.1f" font-weight="%d"%s%s fill="%s" text-anchor="%s">%s</text>',
    x, y, family, size, weight,
    if (italic) ' font-style="italic"' else "",
    if (spacing != 0) sprintf(' letter-spacing="%.2f"', spacing) else "",
    fill, anchor, label)
}
tspan <- function(label, weight = 400, fill = ink, italic = FALSE) {
  sprintf('<tspan font-weight="%d"%s fill="%s">%s</tspan>', weight, if (italic) ' font-style="italic"' else "", fill, esc(label))
}
rect_el <- function(x, y, w, h, fill, stroke = "none", sw = 1, rx = 0) {
  sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%.1f" fill="%s" stroke="%s" stroke-width="%.2f"/>',
    x, y, w, h, rx, fill, stroke, sw)
}
marker_id <- function(col) paste0("arrow-", sub("#", "", col))
path_el <- function(d, stroke, sw, start = FALSE, end = TRUE) {
  sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="%.2f"%s%s/>',
    d, stroke, sw,
    if (start) sprintf(' marker-start="url(#%s)"', marker_id(stroke)) else "",
    if (end) sprintf(' marker-end="url(#%s)"', marker_id(stroke)) else "")
}
marker_def <- function(col) {
  sprintf(paste0('<marker id="%s" viewBox="0 0 10 10" refX="9.5" refY="5" markerWidth="12" markerHeight="12" ',
    'markerUnits="userSpaceOnUse" orient="auto-start-reverse"><path d="M0,1.2 L10,5 L0,8.8 L2.4,5 Z" fill="%s"/></marker>'),
    marker_id(col), col)
}
# Approximate text width for knockouts behind labels (Arial ~0.55 em).
text_width <- function(label, size) nchar(label) * size * 0.55

# A path through `points` (rows of x, y) with its value label on a white
# knockout at `pos`. The label sits on a short horizontal stub next to its
# box, so each number stays with its own row however steep the path is.
labelled_path <- function(points, pos, label, col, credible) {
  d <- paste0("M ", paste(sprintf("%.1f %.1f", points[, 1], points[, 2]), collapse = " L "))
  w <- text_width(label, 14) + 12
  c(
    path_el(d, col, if (credible) 2.2 else 1.3),
    rect_el(pos[1] - w / 2, pos[2] - 10, w, 20, surface),
    text_el(pos[1], pos[2] + 5, esc(label), size = 14, weight = if (credible) 700 else 400,
      fill = if (credible) ink else ink_soft, anchor = "middle")
  )
}
# Length of the labelled horizontal stub next to each box; the indicator
# labels are longer, and their stub must leave room for the arrowhead.
stub_cov <- 132
stub_ind <- 205

# ---- Draw both versions ----------------------------------------------------

draw_diagram <- function(version, out_svg, out_png) {
  full <- version == 1

  # ---- Layout

  W <- if (full) 1600 else 1880
  legend_sep <- 1584  # version 2: dashed line before the key column
  legend_x <- 1614
  columns <- c(56, 600, 1000, 1544)  # covariates | construct | indicators
  title_y <- 64
  header_y <- if (full) 162 else 64

  # Covariate column: group headings, then boxes on a fixed pitch.
  box_x <- 56
  box_w <- 226
  box_h <- 38
  pitch <- 48
  cov <- list()
  heads <- list()
  y <- header_y + 84
  for (g in groups) {
    heads <- c(heads, list(list(y = y, title = g$title, ref = g$ref)))
    y <- y + 12
    for (k in seq_along(g$rows)) {
      cov <- c(cov, list(list(term = names(g$rows)[k], label = g$rows[[k]], top = y, mid = y + box_h / 2)))
      y <- y + pitch
    }
    y <- y + 28
  }
  cov_bottom <- cov[[length(cov)]]$top + box_h

  # Construct, centred on the covariate block.
  cx <- 800
  cy <- (cov[[1]]$mid + cov[[length(cov)]]$mid) / 2
  rx <- 132
  ry <- 80

  # Point on the ellipse boundary in the direction of (x, y), pulled back by `gap`.
  on_ellipse <- function(x, y, gap = 3) {
    u <- c(x - cx, y - cy)
    s <- 1 / sqrt((u[1] / rx)^2 + (u[2] / ry)^2)
    p <- c(cx, cy) + s * u
    p + gap * u / sqrt(sum(u^2))
  }

  # Indicator column.
  ind_x <- 1238
  ind_w <- 252
  ind_h <- 54
  ind_mid <- cy + c(-180, -60, 60, 180)

  footer_top <- cov_bottom + 54
  H <- if (full) footer_top + 300 else cov_bottom + 48

  # ---- Draw

  el <- character()
  add <- function(...) el <<- c(el, ...)

  add(rect_el(0, 0, W, H, surface))

  # Title, centred (version 1 only), and column headings in small capitals.
  if (full) {
    add(text_el(W / 2, title_y, esc("MIMIC model of dental neglect"), size = 32, weight = 600, family = serif, anchor = "middle"))
    add(text_el(W / 2, title_y + 30, esc(paste(
      "Multiple indicators, multiple causes: a Bayesian latent-variable model of 175 patients.",
      "Values are posterior medians with 95% credible intervals."
    )), size = 15.5, family = serif, italic = TRUE, fill = ink_soft, anchor = "middle"))
  }
  header <- list(
    list(x = mean(columns[1:2]), title = "COVARIATES", sub = "patient characteristics that shift neglect"),
    list(x = mean(columns[2:3]), title = "CONSTRUCT", sub = "dental neglect, not measured directly"),
    list(x = mean(columns[3:4]), title = "INDICATORS", sub = "clinical findings that reflect neglect"),
    if (!full) list(x = mean(c(legend_x, W - 40)), title = "KEY", sub = "colours and cut-offs")
  )
  for (h in header) {
    add(text_el(h$x, header_y, esc(h$title), size = 13, weight = 700, anchor = "middle", spacing = 2.2))
    add(text_el(h$x, header_y + 22, esc(h$sub), size = 15, family = serif, italic = TRUE, fill = ink_soft, anchor = "middle"))
  }
  for (x in c(columns[2:3], if (!full) legend_sep)) {
    add(sprintf('<line x1="%.1f" y1="%.1f" x2="%.1f" y2="%.1f" stroke="%s" stroke-width="1" stroke-dasharray="4 6"/>',
      x, header_y - 24, x, cov_bottom + 18, separator))
  }

  # Covariate group headings.
  for (h in heads) {
    label <- paste0(tspan(h$title, 700, ink_soft), if (!is.null(h$ref)) tspan(paste0(" ", h$ref), 400, ink_soft))
    add(text_el(box_x, h$y, label, size = 13.5))
  }

  markers <- c(not_credible, ink_soft)

  # Covariate paths (drawn first, so boxes sit on top) and boxes.
  for (i in seq_along(cov)) {
    r <- row_of(cov[[i]]$term, "probability more neglected (%)")
    s <- ps_style(r)
    markers <- c(markers, s$col)
    x0 <- box_x + box_w
    ym <- cov[[i]]$mid
    end <- on_ellipse(x0 + stub_cov, ym)
    add(labelled_path(rbind(c(x0, ym), c(x0 + stub_cov, ym), end), c(x0 + stub_cov / 2, ym), pct_ci(r), s$col, s$credible))
    add(rect_el(box_x, cov[[i]]$top, box_w, box_h, paper, stroke = s$col, sw = if (s$credible) 1.8 else 1.1))
    add(text_el(box_x + 13, cov[[i]]$mid + 5, esc(cov[[i]]$label), size = 15.5))
  }

  # Construct and its residual.
  z_y <- cy - ry - 80
  add(path_el(sprintf("M %.1f %.1f L %.1f %.1f", cx, z_y + 16, cx, cy - ry - 3), not_credible, 1.2))
  add(sprintf('<circle cx="%.1f" cy="%.1f" r="15" fill="%s" stroke="%s" stroke-width="1.2"/>', cx, z_y, paper, ink_soft))
  add(text_el(cx, z_y + 5.5, "ζ", size = 18, family = serif, italic = TRUE, fill = ink, anchor = "middle"))
  add(text_el(cx + 25, z_y + 5, esc("unexplained variation"), size = 14.5, family = serif, italic = TRUE, fill = ink_soft))
  add(sprintf('<ellipse cx="%.1f" cy="%.1f" rx="%.1f" ry="%.1f" fill="%s" stroke="%s" stroke-width="2.6"/>', cx, cy, rx, ry, paper, ink))
  add(text_el(cx, cy + 4, esc("Dental neglect"), size = 25, weight = 600, family = serif, anchor = "middle"))
  add(text_el(cx, cy + 27, esc("latent construct"), size = 14.5, family = serif, italic = TRUE, fill = ink_soft, anchor = "middle"))

  # Indicator paths and boxes.
  for (i in seq_len(nrow(indicators))) {
    ratio <- row_of(indicators$term[i], indicators$ratio_scale[i])
    std <- row_of(indicators$term[i], "standardised loading")
    s <- load_style(std, ratio)
    markers <- c(markers, s$col)
    ym <- ind_mid[i]
    start <- on_ellipse(ind_x - stub_ind, ym, gap = 2)
    add(labelled_path(rbind(start, c(ind_x - stub_ind, ym), c(ind_x - 1, ym)), c(ind_x - stub_ind / 2 - 6, ym),
      ratio_ci(ratio, indicators$ratio_name[i]), s$col, s$credible))
    top <- ind_mid[i] - ind_h / 2
    add(rect_el(ind_x, top, ind_w, ind_h, paper, stroke = s$col, sw = 1.8))
    add(text_el(ind_x + 13, top + 23, esc(indicators$term[i]), size = 15.5, weight = 700))
    add(text_el(ind_x + 13, top + 42, esc(sprintf("standardised loading %.2f, %s", std$median, s$size)), size = 13, fill = ink_soft))
  }

  # Residual correlation between last visit and X-ray.
  resid <- row_of("Residual correlation, years since dentist - X-ray", "correlation")
  resid_col <- if (resid$lower_95 > 0) ink_soft else not_credible
  xa <- ind_x + ind_w + 3
  add(path_el(sprintf("M %.1f %.1f C %.1f %.1f %.1f %.1f %.1f %.1f",
    xa, ind_mid[3], xa + 44, ind_mid[3], xa + 44, ind_mid[4], xa, ind_mid[4]), resid_col, 1.6, start = TRUE))
  add(text_el(ind_x + ind_w / 2 + 6, (ind_mid[3] + ind_mid[4]) / 2 + 5, paste0(
    tspan("residual ", 400, ink_soft, italic = TRUE),
    tspan(sprintf("r = %.2f [%.2f–%.2f]", resid$median, resid$lower_95, resid$upper_95), 700, ink)
  ), size = 14, anchor = "middle"))


  # Version 2: keys stacked in a side column.
  side_key <- function() {
    kx <- legend_x
    ky <- header_y + 84
    row_h <- 23
    sw_w <- 28
    add(text_el(kx, ky, esc("Probability more neglected"), size = 14, weight = 700))
    add(text_el(kx, ky + 19, esc("red more, blue less neglected"), size = 13, italic = TRUE, fill = ink_soft))
    cuts <- round(c(100 - rev(ps_cut), ps_cut))
    key <- tibble(
      col = c(more[["large"]], more[["moderate"]], more[["small"]], paper, less[["small"]], less[["moderate"]], less[["large"]]),
      name = c("large", "moderate", "small", "negligible", "small", "moderate", "large"),
      range = c(sprintf("> %d%%", cuts[6]), sprintf("%d–%d%%", cuts[5], cuts[6]), sprintf("%d–%d%%", cuts[4], cuts[5]),
        sprintf("%d–%d%%", cuts[3], cuts[4]), sprintf("%d–%d%%", cuts[2], cuts[3]), sprintf("%d–%d%%", cuts[1], cuts[2]),
        sprintf("< %d%%", cuts[1]))
    )
    y0 <- ky + 38
    for (k in seq_len(nrow(key))) {
      yk <- y0 + (k - 1) * row_h
      add(rect_el(kx, yk, sw_w, 14, key$col[k], stroke = if (key$col[k] == paper) not_credible else "none", sw = 0.8))
      add(text_el(kx + sw_w + 10, yk + 12, esc(key$name[k]), size = 13.5, fill = ink_soft))
      add(text_el(W - 40, yk + 12, esc(key$range[k]), size = 13.5, fill = ink_soft, anchor = "end"))
    }
    yn <- y0 + nrow(key) * row_h + 14
    add(text_el(kx, yn, esc("Cut-offs: Cohen's d of 0.2, 0.5"), size = 13, fill = ink_soft))
    add(text_el(kx, yn + 18, esc("and 0.8, as probabilities."), size = 13, fill = ink_soft))

    ky <- yn + 62
    add(text_el(kx, ky, esc("Standardised loading"), size = 14, weight = 700))
    add(text_el(kx, ky + 19, esc("0 unrelated, 1 tracks neglect"), size = 13, italic = TRUE, fill = ink_soft))
    load_key <- tibble(
      col = c(loading[["strong"]], loading[["moderate"]], loading[["weak"]]),
      name = c("strong", "moderate", "weak"),
      range = c("≥ 0.7", "0.5–0.7", "< 0.5")
    )
    y0 <- ky + 38
    for (k in seq_len(nrow(load_key))) {
      yk <- y0 + (k - 1) * row_h
      add(rect_el(kx, yk, sw_w, 14, load_key$col[k]))
      add(text_el(kx + sw_w + 10, yk + 12, esc(load_key$name[k]), size = 13.5, fill = ink_soft))
      add(text_el(W - 40, yk + 12, esc(load_key$range[k]), size = 13.5, fill = ink_soft, anchor = "end"))
    }

    yg <- y0 + nrow(load_key) * row_h + 34
    add(path_el(sprintf("M %.1f %.1f L %.1f %.1f", kx, yg - 4, kx + sw_w + 2, yg - 4), not_credible, 1.3))
    add(text_el(kx + sw_w + 10, yg, esc("grey: 95% interval"), size = 13.5, fill = ink_soft))
    add(text_el(kx + sw_w + 10, yg + 18, esc("includes no effect"), size = 13.5, fill = ink_soft))
  }

  # ---- Legend

  if (full) {
    add(sprintf('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="%s" stroke-width="1"/>', columns[1], footer_top, columns[4], footer_top, rule))

    # Running caption: each line is a list of (text, bold) runs.
    caption <- list(
      list(c("Arrows into the construct", 1), c(" give the probability that a patient with the characteristic is more neglected", 0)),
      list(c("than an otherwise identical patient in the reference group. 50% means no difference; 76% means the patient", 0)),
      list(c("with the characteristic is the more neglected one in about three of every four such pairs.", 0)),
      list(c("Arrows to the indicators", 1), c(" give how much a finding worsens with one standard deviation more neglect:", 0)),
      list(c("the odds of a worse category (OR) or the number of teeth to extract (RR); 1 means unrelated.", 0)),
      list(c("Values", 1), c(" are posterior medians with 95% credible intervals in brackets: given the data and the model, the true", 0)),
      list(c("value lies within the brackets with 95% probability. Coloured paths exclude no effect (50%, or 1 for OR and", 0)),
      list(c("RR); grey paths do not. ", 0), c("Reference groups:", 1), c(" BF patients seen before treatment, women, never smokers,", 0)),
      list(c("patients seen in 2024. The curved arrow is the residual correlation between last visit and X-ray, which", 0)),
      list(c("agree beyond neglect because both reflect contact with dental care.", 0))
    )
    cap_y <- footer_top + 40
    for (k in seq_along(caption)) {
      runs <- paste0(map_chr(caption[[k]], \(r) tspan(r[1], if (r[2] == "1") 600 else 400, if (r[2] == "1") ink else ink_soft)), collapse = "")
      add(text_el(columns[1], cap_y + (k - 1) * 23, runs, size = 15.5, family = serif))
    }

    # Keys.
    kx <- 1010
    key_w <- columns[4] - kx
    ky <- footer_top + 40
    add(text_el(kx, ky, esc("RULES OF THUMB"), size = 12.5, weight = 700, spacing = 2))
    ky <- ky + 30
    add(text_el(kx, ky, paste0(tspan("Probability more neglected", 700, ink), tspan("   blue less, red more", 400, ink_soft, italic = TRUE)), size = 14))
    seg_cols <- c(less[["large"]], less[["moderate"]], less[["small"]], paper, more[["small"]], more[["moderate"]], more[["large"]])
    seg_names <- c("large", "moderate", "small", "negligible", "small", "moderate", "large")
    cuts <- round(c(100 - rev(ps_cut), ps_cut))
    seg_ranges <- c(
      sprintf("< %d%%", cuts[1]), sprintf("%d–%d", cuts[1], cuts[2]), sprintf("%d–%d", cuts[2], cuts[3]),
      sprintf("%d–%d", cuts[3], cuts[4]), sprintf("%d–%d", cuts[4], cuts[5]), sprintf("%d–%d", cuts[5], cuts[6]),
      sprintf("> %d%%", cuts[6])
    )
    seg_w <- key_w / 7
    for (k in 1:7) {
      x0 <- kx + (k - 1) * seg_w
      add(rect_el(x0, ky + 10, seg_w, 12, seg_cols[k], stroke = if (k == 4) not_credible else surface, sw = if (k == 4) 0.8 else 1.5))
      add(text_el(x0 + seg_w / 2, ky + 40, esc(seg_names[k]), size = 13, fill = ink_soft, anchor = "middle"))
      add(text_el(x0 + seg_w / 2, ky + 56, esc(seg_ranges[k]), size = 13, fill = ink_soft, anchor = "middle"))
    }
    add(text_el(kx, ky + 80, esc("Cut-offs are Cohen's d of 0.2, 0.5 and 0.8, expressed as probabilities."), size = 13, fill = ink_soft))

    ky <- ky + 112
    add(text_el(kx, ky, paste0(tspan("Standardised loading", 700, ink), tspan("   0 unrelated, 1 tracks neglect exactly", 400, ink_soft, italic = TRUE)), size = 14))
    load_names <- c("weak  < 0.5", "moderate  0.5–0.7", "strong  ≥ 0.7")
    load_w <- key_w / 3
    for (k in 1:3) {
      x0 <- kx + (k - 1) * load_w
      add(rect_el(x0, ky + 10, load_w, 12, loading[[k]], stroke = surface, sw = 1.5))
      add(text_el(x0 + load_w / 2, ky + 40, esc(load_names[k]), size = 13, fill = ink_soft, anchor = "middle"))
    }

    ky <- ky + 66
    add(path_el(sprintf("M %.1f %.1f L %.1f %.1f", kx, ky - 4, kx + 34, ky - 4), not_credible, 1.3))
    add(text_el(kx + 44, ky, esc("grey: 95% interval includes no effect"), size = 13, fill = ink_soft))
    add(text_el(kx, ky + 26, esc("These are guides, not thresholds: a small effect can still matter clinically."),
      size = 14.5, family = serif, italic = TRUE, fill = ink_soft))
  } else {
    side_key()
  }

  # ---- Write

  svg <- c(
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" xml:space="preserve" width="%d" height="%.0f" viewBox="0 0 %d %.0f">', W, H, W, H),
    "<title>MIMIC model of dental neglect</title>",
    paste0("<desc>Path diagram: nine patient characteristics point to latent dental neglect, which points to four ",
      "clinical indicators; posterior medians and 95% credible intervals on every arrow.</desc>"),
    "<defs>", map_chr(unique(markers), marker_def), "</defs>",
    el,
    "</svg>"
  )
  source_svg <- tempfile(fileext = ".svg")
  writeLines(svg, source_svg, useBytes = TRUE)
  inkscape <- function(...) system2("inkscape", c(shQuote(source_svg), ...), stdout = FALSE, stderr = FALSE)
  inkscape("--export-type=png", paste0("--export-dpi=", png_dpi), paste0("--export-filename=", shQuote(out_png)))
  inkscape("--export-type=svg", "--export-plain-svg", "--export-text-to-path", paste0("--export-filename=", shQuote(out_svg)))
  cat("Written:", out_svg, "and", out_png, "\n")
}

draw_diagram(1, file.path(out_dir, "rq1_neglect_path.svg"), file.path(out_dir, "rq1_neglect_path.png"))
draw_diagram(2, file.path(out_dir, "rq1_neglect_path_v2.svg"), file.path(out_dir, "rq1_neglect_path_v2.png"))

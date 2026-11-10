# Article 1 (RQ1 MIMIC) - Figure 2: path diagram of the reported model.
#
# Run from anywhere inside the project:  Rscript R/question_1/18_art1_figure2_path.R
#
# The publication version of rq1_neglect_path_v2 (R/question_1/06_rq1_neglect_path.R):
# the same estimates and colour rules, redrawn in the plain journal style of
# R/question_1/13_art1_flow_diagram.R (Arial, white background, no title or
# explanatory text: the caption goes in the manuscript), with the colour keys
# under the diagram and the labels of the manuscript. Reads the posterior
# summaries written by R/question_1/05_rq1_neglect.R (no refitting) and has
# Inkscape write
#   art1_figure2_path.svg  text converted to outlines
#   art1_figure2_path.pdf  fonts embedded
#   art1_figure2_path.png  600 DPI
# The drawing is 680 px = 180 mm wide (full page width), so 1 px prints at
# 0.75 pt: box labels at 11 px print at about 8 pt, values and keys at 9.5-10
# px at about 7-7.5 pt.
#
# Colour encodes effect size and credibility, as in 06:
#   covariate -> construct  probability more neglected: red = more, blue =
#                           less; light to dark = small, moderate, large
#                           (Cohen's d 0.2 / 0.5 / 0.8 as probabilities);
#                           grey = the 95% credible interval includes 50%
#   construct -> indicator  standardised loading, green light to dark = weak
#                           (< 0.5), moderate (0.5-0.7), strong (>= 0.7);
#                           grey = the OR / RR interval includes 1

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(dplyr)
  library(purrr)
})

in_file <- here("results", "question_1", "rq1_neglect_coefficients.csv")
out_dir <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png_dpi <- 600

# ---- Palette and type ------------------------------------------------------

# Ramps as in 06 (dataviz reference palette, ordinal checks passed).
paper        <- "#ffffff"
ink          <- "#0b0b0b"
ink_soft     <- "#52514e"
not_credible <- "#8f8d85"
more    <- c(small = "#f0918e", moderate = "#e05250", large = "#b3302f")
less    <- c(small = "#86b6ef", moderate = "#3987e5", large = "#1c5cab")
loading <- c(weak = "#5dbb96", moderate = "#1f9a6e", strong = "#0b6e4c")
sans    <- "Arial, Helvetica, sans-serif"

# Type sizes in px (x 0.75 = pt in print).
fs_label <- 11    # box labels, column headings
fs_value <- 10    # values on the paths, group headings
fs_small <- 9.5   # secondary text and keys
fs_construct <- 13.5

# Probabilities equivalent to Cohen's d 0.2 / 0.5 / 0.8: Phi(d / sqrt(2)).
ps_cut <- 100 * pnorm(c(0.2, 0.5, 0.8) / sqrt(2))

# ---- Estimates -------------------------------------------------------------

coef <- read_csv(in_file, show_col_types = FALSE) |>
  filter(model == "correlated")
row_of <- function(term, scale) {
  r <- filter(coef, term == !!term, scale == !!scale)
  stopifnot(nrow(r) == 1)
  r
}

ps_style <- function(r) {
  if (r$lower_95 <= 50 && r$upper_95 >= 50) return(list(col = not_credible, credible = FALSE))
  dist <- abs(r$median - 50)
  size <- if (dist >= ps_cut[3] - 50) "large" else if (dist >= ps_cut[2] - 50) "moderate" else "small"
  list(col = if (r$median > 50) more[[size]] else less[[size]], credible = TRUE)
}
load_style <- function(std, ratio) {
  if (ratio$lower_95 <= 1) return(list(col = not_credible, credible = FALSE))
  size <- if (std$median >= 0.7) "strong" else if (std$median >= 0.5) "moderate" else "weak"
  list(col = loading[[size]], credible = TRUE)
}

# Ranges with hyphens, as in the manuscript.
num <- function(x) ifelse(x >= 100, sprintf("%.0f", x), ifelse(x >= 10, sprintf("%.1f", x), sprintf("%.2f", x)))
pct_ci <- function(r) sprintf("%.0f%% [%.0f-%.0f]", r$median, r$lower_95, r$upper_95)
ratio_ci <- function(r, name) sprintf("%s %s [%s-%s]", name, num(r$median), num(r$lower_95), num(r$upper_95))

# Labels as in the manuscript and Figures S1-S2; names are the terms in the CSV.
groups <- list(
  list(title = "Patient group", ref = "(vs before BF)", rows = c(
    "BF, seen after" = "During or after BF",
    "RTH, seen before" = "Before RTH",
    "RTH, seen after" = "After RTH"
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
  label = c("Teeth to extract", "Oral hygiene", "Last dentist visit", "Dental X-ray"),
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
text_el <- function(x, y, label, size = fs_label, weight = 400, fill = ink, anchor = "start", italic = FALSE) {
  sprintf('<text x="%.1f" y="%.1f" font-family="%s" font-size="%.1f" font-weight="%d"%s fill="%s" text-anchor="%s">%s</text>',
    x, y, sans, size, weight, if (italic) ' font-style="italic"' else "", fill, anchor, label)
}
tspan <- function(label, weight = 400, fill = ink) {
  sprintf('<tspan font-weight="%d" fill="%s">%s</tspan>', weight, fill, esc(label))
}
rect_el <- function(x, y, w, h, fill, stroke = "none", sw = 1) {
  sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s" stroke-width="%.2f"/>',
    x, y, w, h, fill, stroke, sw)
}
marker_id <- function(col) paste0("arrow-", sub("#", "", col))
path_el <- function(d, stroke, sw, start = FALSE, end = TRUE) {
  sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="%.2f"%s%s/>',
    d, stroke, sw,
    if (start) sprintf(' marker-start="url(#%s)"', marker_id(stroke)) else "",
    if (end) sprintf(' marker-end="url(#%s)"', marker_id(stroke)) else "")
}
marker_def <- function(col) {
  sprintf(paste0('<marker id="%s" viewBox="0 0 10 10" refX="9.5" refY="5" markerWidth="7" markerHeight="7" ',
    'markerUnits="userSpaceOnUse" orient="auto-start-reverse"><path d="M0,1.2 L10,5 L0,8.8 L2.4,5 Z" fill="%s"/></marker>'),
    marker_id(col), col)
}
# Text width in px, measured with the Arial font file (72 dpi: 1 pt = 1 px), for layout and knockouts.
text_width <- function(label, size, bold = FALSE) {
  systemfonts::string_width(label, family = "Arial", weight = if (bold) "bold" else "normal", size = size, res = 72)
}

# A path through `points` (rows of x, y) with its value label on a white
# knockout at `pos`, on a short horizontal stub next to its box.
labelled_path <- function(points, pos, label, col, credible) {
  d <- paste0("M ", paste(sprintf("%.1f %.1f", points[, 1], points[, 2]), collapse = " L "))
  w <- text_width(label, fs_value, credible) + 6
  c(
    path_el(d, col, if (credible) 1.4 else 0.9),
    rect_el(pos[1] - w / 2, pos[2] - 7, w, 14, paper),
    text_el(pos[1], pos[2] + 3.5, esc(label), size = fs_value, weight = if (credible) 700 else 400,
      fill = if (credible) ink else ink_soft, anchor = "middle")
  )
}

# ---- Layout ----------------------------------------------------------------

W <- 680
margin <- 14
heading_y <- 16

# Covariate column: group headings, then boxes on a fixed pitch.
cov_labels <- unlist(map(groups, "rows"))
box_x <- margin
box_w <- max(text_width(cov_labels, fs_label)) + 16
box_h <- 19
pitch <- 25
stub_cov <- max(text_width("76% [62-86]", fs_value, TRUE)) + 10
cov <- list()
heads <- list()
y <- heading_y + 28
for (g in groups) {
  heads <- c(heads, list(list(y = y, title = g$title, ref = g$ref)))
  y <- y + 7
  for (k in seq_along(g$rows)) {
    cov <- c(cov, list(list(term = names(g$rows)[k], label = g$rows[[k]], top = y, mid = y + box_h / 2)))
    y <- y + pitch
  }
  y <- y + 13
}
cov_bottom <- cov[[length(cov)]]$top + box_h

# Construct, centred on the covariate block.
rx <- 64
ry <- 36
cx <- box_x + box_w + stub_cov + 54 + rx
cy <- (cov[[1]]$mid + cov[[length(cov)]]$mid) / 2

# Point on the ellipse boundary in the direction of (x, y), pulled back by `gap`.
on_ellipse <- function(x, y, gap = 2) {
  u <- c(x - cx, y - cy)
  s <- 1 / sqrt((u[1] / rx)^2 + (u[2] / ry)^2)
  p <- c(cx, cy) + s * u
  p + gap * u / sqrt(sum(u^2))
}

# Indicator column: stub long enough for "OR 15.7 [4.92-160]" and the arrowhead.
stub_ind <- text_width("OR 15.7 [4.92-160]", fs_value, TRUE) + 26
ind_x <- cx + rx + 14 + stub_ind
ind_w <- max(text_width(indicators$label, fs_label, TRUE), text_width("standardised loading 0.00", fs_small)) + 14
ind_h <- 32
ind_mid <- cy + c(-120, -40, 40, 120)
stopifnot(ind_x + ind_w + 30 <= W)

key_top <- max(cov_bottom, ind_mid[4] + ind_h / 2) + 34

# ---- Draw ------------------------------------------------------------------

el <- character()
add <- function(...) el <<- c(el, ...)
markers <- c(not_credible, ink_soft)

add(rect_el(0, 0, W, 1e4, paper))  # cropped to H below

# Column headings.
add(text_el(box_x + (box_w + stub_cov) / 2, heading_y, esc("Covariates"), weight = 700, anchor = "middle"))
add(text_el(cx, heading_y, esc("Construct"), weight = 700, anchor = "middle"))
add(text_el(ind_x + ind_w / 2 - stub_ind / 2, heading_y, esc("Indicators"), weight = 700, anchor = "middle"))

# Covariate group headings.
for (h in heads) {
  label <- paste0(tspan(h$title, 700, ink_soft), if (!is.null(h$ref)) tspan(paste0(" ", h$ref), 400, ink_soft))
  add(text_el(box_x, h$y, label, size = fs_value))
}

# Covariate paths (drawn first, so boxes sit on top) and boxes.
for (i in seq_along(cov)) {
  r <- row_of(cov[[i]]$term, "probability more neglected (%)")
  s <- ps_style(r)
  markers <- c(markers, s$col)
  x0 <- box_x + box_w
  ym <- cov[[i]]$mid
  end <- on_ellipse(x0 + stub_cov, ym)
  add(labelled_path(rbind(c(x0, ym), c(x0 + stub_cov, ym), end), c(x0 + stub_cov / 2, ym), pct_ci(r), s$col, s$credible))
  add(rect_el(box_x, cov[[i]]$top, box_w, box_h, paper, stroke = s$col, sw = if (s$credible) 1.2 else 0.8))
  add(text_el(box_x + 8, ym + 4, esc(cov[[i]]$label)))
}

# Construct and its residual.
z_y <- cy - ry - 34
add(path_el(sprintf("M %.1f %.1f L %.1f %.1f", cx, z_y + 9, cx, cy - ry - 2), not_credible, 0.9))
add(sprintf('<circle cx="%.1f" cy="%.1f" r="9" fill="%s" stroke="%s" stroke-width="0.9"/>', cx, z_y, paper, ink_soft))
add(text_el(cx, z_y + 4, "ζ", size = fs_value, italic = TRUE, anchor = "middle"))
add(sprintf('<ellipse cx="%.1f" cy="%.1f" rx="%.1f" ry="%.1f" fill="%s" stroke="%s" stroke-width="1.6"/>', cx, cy, rx, ry, paper, ink))
add(text_el(cx, cy + 2, esc("Dental neglect"), size = fs_construct, weight = 700, anchor = "middle"))
add(text_el(cx, cy + 16, esc("latent construct"), size = fs_small, fill = ink_soft, anchor = "middle"))

# Indicator paths and boxes.
for (i in seq_len(nrow(indicators))) {
  ratio <- row_of(indicators$term[i], indicators$ratio_scale[i])
  std <- row_of(indicators$term[i], "standardised loading")
  s <- load_style(std, ratio)
  markers <- c(markers, s$col)
  ym <- ind_mid[i]
  start <- on_ellipse(ind_x - stub_ind, ym)
  add(labelled_path(rbind(start, c(ind_x - stub_ind, ym), c(ind_x - 1, ym)), c(ind_x - stub_ind / 2 - 4, ym),
    ratio_ci(ratio, indicators$ratio_name[i]), s$col, s$credible))
  top <- ym - ind_h / 2
  add(rect_el(ind_x, top, ind_w, ind_h, paper, stroke = s$col, sw = 1.2))
  add(text_el(ind_x + 7, top + 13.5, esc(indicators$label[i]), weight = 700))
  add(text_el(ind_x + 7, top + 26, esc(sprintf("standardised loading %.2f", std$median)), size = fs_small, fill = ink_soft))
}

# Residual correlation between last visit and X-ray.
resid <- row_of("Residual correlation, years since dentist - X-ray", "correlation")
resid_col <- if (resid$lower_95 > 0) ink_soft else not_credible
xa <- ind_x + ind_w + 2
add(path_el(sprintf("M %.1f %.1f C %.1f %.1f %.1f %.1f %.1f %.1f",
  xa, ind_mid[3], xa + 22, ind_mid[3], xa + 22, ind_mid[4], xa, ind_mid[4]), resid_col, 1, start = TRUE))
add(text_el(ind_x + ind_w / 2, (ind_mid[3] + ind_mid[4]) / 2 + 3.5, paste0(
  tspan("residual ", 400, ink_soft),
  tspan(sprintf("r = %.2f [%.2f-%.2f]", resid$median, resid$lower_95, resid$upper_95), 700, ink)
), size = fs_value, anchor = "middle"))

# ---- Keys under the diagram ------------------------------------------------

# Probability more neglected: seven bands, blue (less) to red (more).
kx <- margin
ky <- key_top
cuts <- round(c(100 - rev(ps_cut), ps_cut))
seg_cols <- c(less[["large"]], less[["moderate"]], less[["small"]], paper, more[["small"]], more[["moderate"]], more[["large"]])
seg_names <- c("large", "moderate", "small", "negligible", "small", "moderate", "large")
seg_ranges <- c(
  sprintf("< %d%%", cuts[1]), sprintf("%d-%d", cuts[1], cuts[2]), sprintf("%d-%d", cuts[2], cuts[3]),
  sprintf("%d-%d", cuts[3], cuts[4]), sprintf("%d-%d", cuts[4], cuts[5]), sprintf("%d-%d", cuts[5], cuts[6]),
  sprintf("> %d%%", cuts[6])
)
seg_w <- 54
add(text_el(kx, ky, paste0(tspan("Probability more neglected", 700, ink), tspan("  blue less, red more", 400, ink_soft)), size = fs_value))
for (k in 1:7) {
  x0 <- kx + (k - 1) * seg_w
  add(rect_el(x0, ky + 7, seg_w, 8, seg_cols[k], stroke = if (k == 4) not_credible else paper, sw = if (k == 4) 0.6 else 1))
  add(text_el(x0 + seg_w / 2, ky + 27, esc(seg_names[k]), size = fs_small, fill = ink_soft, anchor = "middle"))
  add(text_el(x0 + seg_w / 2, ky + 39, esc(seg_ranges[k]), size = fs_small, fill = ink_soft, anchor = "middle"))
}
add(text_el(kx, ky + 56, esc("Cut-offs: Cohen's d of 0.2, 0.5 and 0.8, as probabilities."), size = fs_small, fill = ink_soft))

# Standardised loading: three bands, and grey for no credible effect.
lx <- kx + 7 * seg_w + 34
lw <- (W - margin - lx) / 3
add(text_el(lx, ky, paste0(tspan("Standardised loading", 700, ink), tspan("  0 unrelated, 1 tracks neglect", 400, ink_soft)), size = fs_value))
load_names <- c("weak", "moderate", "strong")
load_ranges <- c("< 0.5", "0.5-0.7", "≥ 0.7")
for (k in 1:3) {
  x0 <- lx + (k - 1) * lw
  add(rect_el(x0, ky + 7, lw, 8, loading[[k]], stroke = paper, sw = 1))
  add(text_el(x0 + lw / 2, ky + 27, esc(load_names[k]), size = fs_small, fill = ink_soft, anchor = "middle"))
  add(text_el(x0 + lw / 2, ky + 39, esc(load_ranges[k]), size = fs_small, fill = ink_soft, anchor = "middle"))
}
add(path_el(sprintf("M %.1f %.1f L %.1f %.1f", lx, ky + 53, lx + 24, ky + 53), not_credible, 0.9))
add(text_el(lx + 30, ky + 56, esc("Grey: 95% interval includes no effect"), size = fs_small, fill = ink_soft))

H <- ky + 64

# ---- Write -----------------------------------------------------------------

el[1] <- rect_el(0, 0, W, H, paper)
svg <- c(
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" xml:space="preserve" width="%d" height="%.0f" viewBox="0 0 %d %.0f">', W, H, W, H),
  "<title>Path diagram of the dental-neglect model</title>",
  paste0("<desc>Nine patient characteristics point to latent dental neglect, which points to four clinical ",
    "indicators; posterior medians and 95% credible intervals on every arrow.</desc>"),
  "<defs>", map_chr(unique(markers), marker_def), "</defs>",
  el,
  "</svg>"
)
out <- file.path(out_dir, paste0("art1_figure2_path.", c("svg", "pdf", "png")))
source_svg <- tempfile(fileext = ".svg")
writeLines(svg, source_svg, useBytes = TRUE)
inkscape <- function(...) system2("inkscape", c(shQuote(source_svg), ...), stdout = FALSE, stderr = FALSE)
inkscape("--export-type=png", paste0("--export-dpi=", png_dpi), paste0("--export-filename=", shQuote(out[3])))
inkscape("--export-type=pdf", paste0("--export-filename=", shQuote(out[2])))
inkscape("--export-type=svg", "--export-plain-svg", "--export-text-to-path", paste0("--export-filename=", shQuote(out[1])))
cat(sprintf("Written: %s (%d x %.0f px = 180 x %.0f mm)\n", paste(basename(out), collapse = ", "), W, H, H * 180 / W))

# Article 1 (RQ1 MIMIC) - participant flow diagram: from the patients in the
# register to the patients in the model.
#
# Run from anywhere inside the project:  Rscript R/question_1/13_art1_flow_diagram.R
#
# Recounts every step from the clean data with the filters of
# R/question_1/05_rq1_neglect.R (exactly one of groups 1-4, complete
# covariates), so the diagram follows reruns of the cleaning. Draws a
# hand-built SVG in plain journal style (Arial, white boxes, straight arrows,
# no title: the caption goes in the manuscript) and has Inkscape write
#   art1_flow_diagram.svg  text converted to outlines
#   art1_flow_diagram.pdf  fonts embedded
#   art1_flow_diagram.png  600 DPI
# The drawing is 680 px = 180 mm wide (full page width), so 12-13 px type
# prints at 9-10 pt.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

in_file <- here("data", "processed", "clean_data.rds")
out_dir <- here("results", "art1_MIMIC")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
png_dpi <- 600

# ---- Palette and type ------------------------------------------------------

paper    <- "#ffffff"
ink      <- "#0b0b0b"
ink_soft <- "#52514e"
sans     <- "Arial, Helvetica, sans-serif"

# ---- Counts ----------------------------------------------------------------

covariates <- c(
  teeth_to_extract = "teeth to extract", oral_hygiene = "oral hygiene", age_at_visit = "age",
  sex = "sex", smoking_status = "smoking status", visit_year = "visit year"
)
groups <- c(
  group_bf_before = "Before BF", group_bf_after = "During or after BF",
  group_rth_before = "Before RTH", group_rth_after = "After RTH"
)

d <- readRDS(in_file) |>
  mutate(n_groups = group_bf_before + group_bf_after + group_rth_before + group_rth_after + group_other)
eligible <- filter(d, n_groups == 1, group_other == 0)
analysed <- drop_na(eligible, all_of(names(covariates)))

no_group   <- sum(is.na(d$n_groups) | d$n_groups == 0)
two_groups <- sum(d$n_groups >= 2, na.rm = TRUE)
other_only <- sum(d$n_groups == 1 & d$group_other == 1, na.rm = TRUE)
missing_by <- sort(colSums(is.na(eligible[names(covariates)])), decreasing = TRUE)
missing_by <- missing_by[missing_by > 0]
group_n    <- colSums(analysed[names(groups)])
years      <- range(d$visit_year, na.rm = TRUE)

# ---- SVG primitives --------------------------------------------------------

esc <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  gsub(">", "&gt;", x, fixed = TRUE)
}
text_el <- function(x, y, label, size = 12, weight = 400, fill = ink, anchor = "start") {
  sprintf('<text x="%.1f" y="%.1f" font-family="%s" font-size="%.1f" font-weight="%d" fill="%s" text-anchor="%s">%s</text>',
    x, y, sans, size, weight, fill, anchor, label)
}
rect_el <- function(x, y, w, h, fill, stroke = "none", sw = 1) {
  sprintf('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="%s" stroke="%s" stroke-width="%.2f"/>',
    x, y, w, h, fill, stroke, sw)
}
path_el <- function(d, end = TRUE) {
  sprintf('<path d="%s" fill="none" stroke="%s" stroke-width="1.00"%s/>',
    d, ink_soft, if (end) ' marker-end="url(#arrow)"' else "")
}
marker_def <- paste0('<marker id="arrow" viewBox="0 0 10 10" refX="9.5" refY="5" markerWidth="8" markerHeight="8" ',
  'markerUnits="userSpaceOnUse" orient="auto-start-reverse"><path d="M0,1.2 L10,5 L0,8.8 L2.4,5 Z" fill="', ink_soft, '"/></marker>')

# A box with a bold first line and plain lines under it. Its height follows
# the number of lines: first baseline 22 below the top, 17 between lines.
box_h <- function(n_lines) 35 + 17 * (n_lines - 1)
box <- function(x, y, w, lines, anchor = "middle", stroke = ink_soft, sw = 1) {
  tx <- if (anchor == "middle") x + w / 2 else x + 12
  ys <- y + 22 + 17 * (seq_along(lines) - 1)
  c(rect_el(x, y, w, box_h(length(lines)), paper, stroke, sw),
    text_el(tx, ys[1], esc(lines[1]), size = 13, weight = 700, anchor = anchor),
    map2_chr(ys[-1], lines[-1], \(y_line, line) text_el(tx, y_line, esc(line), anchor = anchor)))
}

# ---- Layout ----------------------------------------------------------------

register_lines <- c(sprintf("%d patients in the register", nrow(d)), sprintf("First visits %d–%d", years[1], years[2]))
ex1_lines <- c(
  sprintf("%d excluded", no_group + two_groups + other_only),
  sprintf("%d with no group recorded", no_group),
  sprintf("%d recorded in more than one group", two_groups),
  sprintf("%d in the “other” group only", other_only)
)
eligible_lines <- c(sprintf("%d eligible", nrow(eligible)), "Exactly one of groups 1–4")
ex2_lines <- c(
  sprintf("%d excluded: missing covariate", nrow(eligible) - nrow(analysed)),
  sprintf("%d %s", missing_by, covariates[names(missing_by)])
)
analysed_lines <- sprintf("%d analysed", nrow(analysed))

W <- 680
m <- 16       # outer margin and gap between group boxes
gap <- 20     # vertical gap between rows
main_x <- m
main_w <- 280
cx <- main_x + main_w / 2
side_x <- 360
side_w <- W - m - side_x
grp_w <- (W - 2 * m - 3 * m) / 4
grp_x <- m + (0:3) * (grp_w + m)
grp_cx <- grp_x + grp_w / 2

y_register <- m
y_ex1 <- y_register + box_h(length(register_lines)) + gap
y_eligible <- y_ex1 + box_h(length(ex1_lines)) + gap
y_ex2 <- y_eligible + box_h(length(eligible_lines)) + gap
y_analysed <- y_ex2 + box_h(length(ex2_lines)) + gap
y_bus <- y_analysed + box_h(1) + gap
y_groups <- y_bus + gap
H <- y_groups + box_h(2) + m

el <- c(
  rect_el(0, 0, W, H, paper),
  path_el(sprintf("M %.1f %.1f L %.1f %.1f", cx, y_register + box_h(length(register_lines)), cx, y_eligible)),
  path_el(sprintf("M %.1f %.1f L %.1f %.1f", cx, y_ex1 + box_h(length(ex1_lines)) / 2, side_x, y_ex1 + box_h(length(ex1_lines)) / 2)),
  path_el(sprintf("M %.1f %.1f L %.1f %.1f", cx, y_eligible + box_h(length(eligible_lines)), cx, y_analysed)),
  path_el(sprintf("M %.1f %.1f L %.1f %.1f", cx, y_ex2 + box_h(length(ex2_lines)) / 2, side_x, y_ex2 + box_h(length(ex2_lines)) / 2)),
  path_el(sprintf("M %.1f %.1f L %.1f %.1f", cx, y_analysed + box_h(1), cx, y_bus), end = FALSE),
  path_el(sprintf("M %.1f %.1f L %.1f %.1f", grp_cx[1], y_bus, grp_cx[4], y_bus), end = FALSE),
  map_chr(grp_cx, \(x) path_el(sprintf("M %.1f %.1f L %.1f %.1f", x, y_bus, x, y_groups))),
  box(main_x, y_register, main_w, register_lines),
  box(side_x, y_ex1, side_w, ex1_lines, anchor = "start"),
  box(main_x, y_eligible, main_w, eligible_lines),
  box(side_x, y_ex2, side_w, ex2_lines, anchor = "start"),
  box(main_x, y_analysed, main_w, analysed_lines, stroke = ink, sw = 2),
  unlist(map(1:4, \(i) box(grp_x[i], y_groups, grp_w, c(groups[[i]], sprintf("n = %d", group_n[[i]])))))
)

# ---- Write -----------------------------------------------------------------

svg <- c(
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>",
  sprintf('<svg xmlns="http://www.w3.org/2000/svg" xml:space="preserve" width="%d" height="%.0f" viewBox="0 0 %d %.0f">', W, H, W, H),
  "<title>Participant flow</title>",
  sprintf("<desc>Participant flow: %d patients in the register, %d eligible, %d analysed in four groups.</desc>",
    nrow(d), nrow(eligible), nrow(analysed)),
  "<defs>", marker_def, "</defs>",
  el,
  "</svg>"
)
source_svg <- tempfile(fileext = ".svg")
writeLines(svg, source_svg, useBytes = TRUE)
out <- file.path(out_dir, paste0("art1_flow_diagram.", c(svg = "svg", pdf = "pdf", png = "png")))
inkscape <- function(...) system2("inkscape", c(shQuote(source_svg), ...), stdout = FALSE, stderr = FALSE)
inkscape("--export-type=png", paste0("--export-dpi=", png_dpi), paste0("--export-filename=", shQuote(out[3])))
inkscape("--export-type=pdf", paste0("--export-filename=", shQuote(out[2])))
inkscape("--export-type=svg", "--export-plain-svg", "--export-text-to-path", paste0("--export-filename=", shQuote(out[1])))
cat("Written:", paste(out, collapse = ", "), "\n")

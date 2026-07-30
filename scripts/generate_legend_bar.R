# ==============================================================================
# SCRIPT: SPLIT-SCALE COLOUR BAR FOR GIS LEGEND (QGIS-READY)
# Project:  Accumulated Thermal Stress in Posidonia oceanica (Balearic Islands)
# Author:   Yeray Gonzalez-Marrero
# Purpose:  Generate a publication-quality horizontal colour bar with a split
#           scale for use as a legend in QGIS map layouts. The split scale
#           allocates 60% of the bar width to the ecologically critical DHD
#           range (0–7), where both management thresholds lie, and 40% to the
#           extreme stress range (7–55.6), corresponding to the 2022 heatwave.
#           This design ensures that the threshold transitions are visually
#           prominent without compressing the high-stress range.
# Outputs:  outputs/maps/legend_DHD_split.png  (300 dpi, for QGIS import)
#           outputs/maps/legend_DHD_split.svg  (vector, for Inkscape editing)
# ==============================================================================

# Delete workspace
rm(list = ls())


library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

# ==============================================================================
# SECTION 1: INDICATOR PARAMETERS
# ==============================================================================
# These values must match the thresholds in outputs/tables/management_thresholds.xlsx
# and the colour ramp used in Script 04 (Concept 5).
alert_threshold    <- 4.9      # Green-to-Yellow transition (75th percentile of sub-lethal DHD)
critical_threshold <- 5.6      # Yellow-to-Red transition (zero-crossing of GAMM prediction curve)
split_at           <- 7.0      # Scale break point: separates the threshold zone from the extreme zone
max_dhd            <- 55.6195  # Maximum observed DHD (2022 heatwave year)

# Proportional width allocation for each bar segment
frac_left  <- 0.60   # 60% of total width covers DHD 0 – split_at (threshold zone)
frac_right <- 0.40   # 40% of total width covers DHD split_at – max_dhd (extreme zone)

# ==============================================================================
# SECTION 2: COLOUR GRADIENT DEFINITION
# ==============================================================================
# The gradient uses stepped transitions at the management thresholds to produce
# visually distinct zones (Green / Yellow / Red) whilst retaining within-zone
# variability through a continuous gradient. Duplicate colour values at each
# threshold create an abrupt colour shift rather than a gradual blend.
#
# Colour anchors (calibrated to the Xarxa de Monitoratge field data):
#   #2e7d32 — Dark green  (DHD = 0, no thermal stress)
#   #fbc02d — Yellow      (DHD = alert threshold, sub-lethal stress onset)
#   #d3402f — Red         (DHD = critical threshold, net regression onset)
#   #5d0000 — Dark red    (DHD = max, extreme heatwave conditions)

colour_min      <- "#2e7d32"   # Dark green  — DHD = 0
colour_alert    <- "#fbc02d"   # Yellow      — DHD = alert threshold
colour_critical <- "#d3402f"   # Red         — DHD = critical threshold
colour_max      <- "#5d0000"   # Dark red    — DHD = max_dhd

gradient_colours <- c(
  "#2e7d32",          # Dark green — DHD = 0
  "#29dc2f",          # Light green — just below alert threshold
  "#fbc02d",          # Yellow — alert threshold (abrupt transition)
  "#fbc02d",          # Yellow — just above alert threshold (holds yellow zone)
  "#d3402f",          # Red — critical threshold (abrupt transition)
  "#d3402f",          # Red — just above critical threshold (holds red zone)
  "#5d0000"           # Dark red — maximum DHD
)

# Normalised anchor positions (0–1) mapped to the full DHD range (0–max_dhd).
# The near-duplicate values at each threshold create the step-like transitions.
gradient_values <- rescale(
  c(0,
    alert_threshold - 0.05,   # just below alert: still green
    alert_threshold,           # alert: switches to yellow
    alert_threshold + 0.4,    # solidly yellow zone
    critical_threshold,        # critical: switches to red
    critical_threshold + 0.05, # solidly red zone
    max_dhd),
  to   = c(0, 1),
  from = c(0, max_dhd)
)

# Helper function: return the mapped colour for a given DHD value
colour_for_dhd <- function(val) {
  cmap <- scales::col_numeric(palette = gradient_colours, domain = c(0, max_dhd))
  cmap(val)
}

# ==============================================================================
# SECTION 3: LEFT SEGMENT (DHD 0 → split_at, 60% of bar width)
# ==============================================================================
# Rendered as a dense sequence of geom_tile() rectangles to approximate a
# smooth gradient. The two management thresholds are marked with white dashed
# vertical lines and annotated above the bar.
n_steps <- 500   # Number of gradient steps (higher = smoother gradient)

df_left <- data.frame(
  x   = seq(0, frac_left, length.out = n_steps),
  val = seq(0, split_at,  length.out = n_steps),
  y   = 0.5
)

p_left <- ggplot(df_left, aes(x = x, y = y, fill = val)) +
  geom_tile(height = 1, width = frac_left / n_steps) +
  scale_fill_gradientn(
    colours = gradient_colours,
    values  = gradient_values,
    limits  = c(0, max_dhd),
    guide   = "none"
  ) +
  # Vertical threshold markers
  geom_vline(xintercept = (alert_threshold    / split_at) * frac_left,
             colour = "white", linewidth = 0.9, linetype = "dashed") +
  geom_vline(xintercept = (critical_threshold / split_at) * frac_left,
             colour = "white", linewidth = 0.9, linetype = "dashed") +
  # Alert threshold annotation (positioned to the left of the line)
  annotate("text",
           x     = (alert_threshold / split_at) * frac_left - 0.008,
           y     = 1.55,
           label = sprintf("%.4f DHD", alert_threshold),
           hjust = 1, vjust = 0, size = 4.2, fontface = "bold",
           colour = "#222222") +
  annotate("text",
           x     = (alert_threshold / split_at) * frac_left - 0.008,
           y     = 1.18,
           label = "Sub-lethal stress onset",
           hjust = 1, vjust = 0, size = 3.8, colour = "#555555") +
  # Critical threshold annotation (positioned to the left to avoid clipping)
  annotate("text",
           x     = (critical_threshold / split_at) * frac_left - 0.008,
           y     = 1.55,
           label = sprintf("%.4f DHD", critical_threshold),
           hjust = 0, vjust = 0, size = 4.2, fontface = "bold",
           colour = "#222222") +
  annotate("text",
           x     = (critical_threshold / split_at) * frac_left - 0.008,
           y     = 1.18,
           label = "Net regression onset",
           hjust = 0, vjust = 0, size = 3.8, colour = "#555555") +
  # Numeric axis ticks at threshold values and segment boundary
  scale_x_continuous(
    limits = c(0, frac_left),
    breaks = c(0,
               (alert_threshold    / split_at) * frac_left,
               (critical_threshold / split_at) * frac_left,
               frac_left),
    labels = c("0",
               sprintf("%.1f", alert_threshold),
               sprintf("%.1f", critical_threshold),
               sprintf("%.0f", split_at)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(-0.5, 2.2), expand = c(0, 0)) +
  theme_void() +
  theme(
    axis.text.x          = element_text(size = 9.5, colour = "#333333",
                                        margin = margin(t = 3)),
    axis.ticks.x         = element_line(colour = "#666666", linewidth = 0.4),
    axis.ticks.length.x  = unit(2.5, "pt"),
    plot.margin          = margin(t = 36, r = 0, b = 14, l = 8)
  )

# ==============================================================================
# SECTION 4: RIGHT SEGMENT (DHD split_at → max_dhd, 40% of bar width)
# ==============================================================================
# Covers the extreme stress range observed during the 2022 heatwave.
# No threshold annotations are needed here; only numeric axis ticks at
# representative DHD values (20, 40, and the maximum).
df_right <- data.frame(
  x   = seq(0, frac_right, length.out = n_steps),
  val = seq(split_at, max_dhd, length.out = n_steps),
  y   = 0.5
)

p_right <- ggplot(df_right, aes(x = x, y = y, fill = val)) +
  geom_tile(height = 1, width = frac_right / n_steps) +
  scale_fill_gradientn(
    colours = gradient_colours,
    values  = gradient_values,
    limits  = c(0, max_dhd),
    guide   = "none"
  ) +
  scale_x_continuous(
    limits = c(0, frac_right),
    breaks = c(0,
               (20 - split_at) / (max_dhd - split_at) * frac_right,
               (40 - split_at) / (max_dhd - split_at) * frac_right,
               frac_right),
    labels = c("", "20", "40", sprintf("%.1f", max_dhd)),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(-0.5, 2.2), expand = c(0, 0)) +
  theme_void() +
  theme(
    axis.text.x          = element_text(size = 9.5, colour = "#333333",
                                        margin = margin(t = 3)),
    axis.ticks.x         = element_line(colour = "#666666", linewidth = 0.4),
    axis.ticks.length.x  = unit(2.5, "pt"),
    plot.margin          = margin(t = 36, r = 8, b = 14, l = 0)
  )

# ==============================================================================
# SECTION 5: SCALE BREAK SYMBOL (//)
# ==============================================================================
# A pair of diagonal line segments between the two bar segments indicates the
# discontinuity in the scale. This is the standard cartographic convention for
# split-scale legends.
p_break <- ggplot() +
  annotate("segment", x = 0.2, xend = 0.8, y = 0.2, yend = 0.8,
           colour = "#888888", linewidth = 0.7) +
  annotate("segment", x = 0.4, xend = 1.0, y = 0.2, yend = 0.8,
           colour = "#888888", linewidth = 0.7) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(t = 28, r = 0, b = 14, l = 0))

# ==============================================================================
# SECTION 6: LEGEND TITLE
# ==============================================================================
legend_title <- textGrob(
  "Accumulated heat (°C \u00b7 days)",
  gp = gpar(fontsize = 12, fontface = "bold", col = "#111111")
)

# ==============================================================================
# SECTION 7: COMPOSE AND EXPORT
# ==============================================================================
# The three bar segments and the title are assembled using gridExtra::arrangeGrob().
# Layout: [left_segment | break_symbol | right_segment] on row 1,
#         [legend_title (spanning all columns)]          on row 2.
layout_matrix <- rbind(c(1, 2, 3),
                       c(4, 4, 4))

legend_grob <- arrangeGrob(
  p_left, p_break, p_right,
  legend_title,
  layout_matrix = layout_matrix,
  widths         = c(frac_left, 0.04, frac_right),
  heights        = c(0.88, 0.12)
)

# Export as PNG (for direct import into QGIS Print Layout)
ggsave("outputs/maps/legend_DHD_split.png",
       plot   = legend_grob,
       width  = 15,
       height = 4.5,
       units  = "cm",
       dpi    = 300,
       bg     = "white")

# Export as SVG (for vector editing in Inkscape or Illustrator)
ggsave("outputs/maps/legend_DHD_split.svg",
       plot   = legend_grob,
       width  = 18,    # slightly wider to prevent label clipping
       height = 5.5,
       units  = "cm",
       bg     = "white")

cat("[OK] Legend saved to:\n")
cat("     outputs/maps/legend_DHD_split.png  (import into QGIS Print Layout via Add Item > Add Image)\n")
cat("     outputs/maps/legend_DHD_split.svg  (open in Inkscape for further editing)\n")

# ============================================================
# ERS/NDSU plotting framework for the gwkit examples
# ============================================================
# Self-contained copy of the plotting style used in the
# warming-impacts-alfalfa study, so the example .Rmd files render with the
# same look without depending on that sibling repository. Source this file from
# each example's setup chunk: source("_ers_framework.R").
#
#   ers_theme()         - the USDA-REE-ERS MTED ggplot2 theme.
#   gw_diverging_map()  - signed continuous surface on a state basemap
#                         (red = negative, blue = positive), symmetric 98% clip.
#   state_class_map()   - categorical state choropleth (NDSU palette).
#
# ers_theme() is copied from https://github.com/USDA-REE-ERS/MTED-Theme.
# ============================================================

# NDSU / ERS categorical palette (Rust -> Night)
ndsu_palette <- c(
  "#BE5E27", "#FFC425", "#FEF389", "#BED73B", "#A0BD78", "#00583D",
  "#003524", "#D7E5C8", "#9DD9F7", "#51ABA0", "#0F374B")

# --- ers_theme(): USDA-REE-ERS MTED theme (copied 08/01/2025) -----------------
ers_theme <- function() {
  ggplot2::theme(
    line   = ggplot2::element_line(colour = "black", linewidth = 0.5, linetype = 1,
                                   lineend = "butt"),
    rect   = ggplot2::element_rect(fill = "white", colour = NA, linewidth = 0.5,
                                   linetype = 1),
    text   = ggplot2::element_text(family = "sans", face = "plain", colour = "black",
                                   size = 9, hjust = 0.5, vjust = 0.5, angle = 0,
                                   lineheight = 0.9,
                                   margin = ggplot2::margin(0, 0, 0, 0, "pt")),
    axis.title.y          = ggplot2::element_blank(),
    axis.ticks            = ggplot2::element_blank(),
    axis.line             = ggplot2::element_blank(),
    legend.key            = ggplot2::element_rect(fill = "white", colour = NA),
    legend.key.size       = ggplot2::unit(1.2, "lines"),
    legend.title          = ggplot2::element_text(hjust = 0, size = 9),
    legend.text           = ggplot2::element_text(size = 9),
    legend.position       = "bottom",
    legend.justification  = "center",
    panel.border          = ggplot2::element_blank(),
    panel.grid.major.x    = ggplot2::element_blank(),
    panel.grid.minor.x    = ggplot2::element_blank(),
    panel.grid.major.y    = ggplot2::element_line(colour = "grey92"),
    panel.grid.minor.y    = ggplot2::element_line(linewidth = ggplot2::rel(0.5)),
    plot.title            = ggplot2::element_text(face = "bold", size = 10.5, hjust = 0,
                                                  vjust = 1,
                                                  margin = ggplot2::margin(5.5, 0, 5.5, 0, "pt")),
    plot.title.position   = "plot",
    plot.subtitle         = ggplot2::element_text(hjust = 0, vjust = 1,
                                                  margin = ggplot2::margin(0, 0, 5, 0, "pt")),
    plot.caption          = ggplot2::element_text(size = 8, hjust = 0, vjust = 1,
                                                  margin = ggplot2::margin(5.5, 0, 0, 0, "pt")),
    plot.caption.position = "plot",
    strip.background      = ggplot2::element_rect(fill = "grey85", colour = "grey20"),
    strip.text            = ggplot2::element_text(colour = "grey10", size = ggplot2::rel(0.8),
                                                  margin = ggplot2::margin(4.4, 4.4, 4.4, 4.4, "pt"))
  )
}

# --- gw_diverging_map(): signed continuous surface on a state basemap ---------
# `sf_states` is an sf state layer carrying `fill_col`. Mirrors the alfalfa
# 101_gw_scalar_consensus_map idiom: grey basemap, red-white-blue diverging fill
# centred at 0, symmetric clip at the 98th percentile of |value|.
gw_diverging_map <- function(sf_states, fill_col, title, subtitle = NULL,
                             legend = "Local\ncoefficient", caption = NULL,
                             limits = NULL) {
  v   <- as.data.frame(sf_states)[, fill_col]
  lim <- if (is.null(limits)) c(-1, 1) * stats::quantile(abs(v), 0.98, na.rm = TRUE) else limits
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf_states, colour = "black", fill = "grey95", linewidth = 0.2) +
    ggplot2::geom_sf(data = sf_states, ggplot2::aes(fill = .data[[fill_col]]),
                     colour = "grey30", linewidth = 0.1) +
    ggplot2::scale_fill_gradient2(low = "#B2182B", mid = "#F7F7F7", high = "#2166AC",
                                  midpoint = 0, limits = lim, oob = scales::squish,
                                  name = legend) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "", y = "", caption = caption) +
    ers_theme() + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9)) +
    ggplot2::coord_sf()
}

# --- gw_sequential_map(): non-diverging level surface (NDSU greens) -----------
gw_sequential_map <- function(sf_states, fill_col, title, subtitle = NULL,
                              legend = "Level", caption = NULL) {
  ggplot2::ggplot() +
    ggplot2::geom_sf(data = sf_states, colour = "black", fill = "grey95", linewidth = 0.2) +
    ggplot2::geom_sf(data = sf_states, ggplot2::aes(fill = .data[[fill_col]]),
                     colour = "grey30", linewidth = 0.1) +
    ggplot2::scale_fill_gradient(low = "#D7E5C8", high = "#003524", name = legend) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "", y = "", caption = caption) +
    ers_theme() + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9)) +
    ggplot2::coord_sf()
}

# --- state_class_map(): categorical state choropleth (NDSU palette) -----------
# Trimmed plot_us_states_choropleth: fills states by a categorical column and
# labels each with its two-letter abbreviation. `palette` is a named vector
# mapping class label -> hex.
state_class_map <- function(sf_states, class_col, title, subtitle = NULL,
                            legend = "Class", palette = NULL, caption = NULL) {
  if (is.null(palette))
    palette <- stats::setNames(ndsu_palette[seq_along(unique(sf_states[[class_col]]))],
                               sort(unique(as.character(sf_states[[class_col]]))))
  lab <- if ("state_abbv" %in% names(sf_states)) "state_abbv" else NULL
  g <- ggplot2::ggplot(sf_states) +
    ggplot2::geom_sf(ggplot2::aes(fill = .data[[class_col]]), colour = "grey30",
                     linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = palette, na.value = "grey95", name = legend)
  if (!is.null(lab))
    g <- g + ggplot2::geom_sf_text(ggplot2::aes(label = .data[[lab]]), size = 2.3,
                                   colour = "black")
  g +
    ggplot2::labs(title = title, subtitle = subtitle, x = "", y = "", caption = caption) +
    ers_theme() + ggplot2::theme_bw() +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   legend.position = "right",
                   plot.title = ggplot2::element_text(size = 11, face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 9)) +
    ggplot2::coord_sf()
}

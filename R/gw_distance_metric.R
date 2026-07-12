#' GW distance metric presets for GWmodel
#'
#' Provides a curated set of distance metric presets (Minkowski family and great-circle)
#' for \pkg{GWmodel}. Each preset specifies \code{(p, theta, longlat)} for
#' \code{GWmodel::gw.dist()}.
#'
#' @return A named list of presets, each entry a \code{list(p, theta, longlat)}.
#' @export
gw_distance_metric_presets <- function() {
  list(
    # Euclidean / L2
    "Euclidean"                       = list(p = 2.0,  theta = 0.0, longlat = FALSE),
    "Euclidean (rotated theta=0.8)"       = list(p = 2.0,  theta = 0.8, longlat = FALSE), # rotation no-op for p=2

    # Manhattan / L1
    "Manhattan"                       = list(p = 1.0,  theta = 0.0, longlat = FALSE),
    "Manhattan (rotated theta=0.5)"       = list(p = 1.0,  theta = 0.5, longlat = FALSE),

    # Minkowski (general Lp)
    "Minkowski p=1.5"                 = list(p = 1.5,  theta = 0.0, longlat = FALSE),
    "Minkowski p=1.5 (rotated theta=0.8)" = list(p = 1.5,  theta = 0.8, longlat = FALSE),
    "Minkowski p=3"                   = list(p = 3.0,  theta = 0.0, longlat = FALSE),
    "Minkowski p=3 (rotated theta=0.8)"   = list(p = 3.0,  theta = 0.8, longlat = FALSE),

    # Chebyshev / L_inf (approx via large p)
    "Chebyshev (approx L_inf, p = 10)"     = list(p = 10.0, theta = 0.0, longlat = FALSE),

    # Geodesic
    "Great Circle"                    = list(p = 2.0,  theta = 0.0, longlat = TRUE)
  )
}

#' Resolve a GW distance metric preset
#' @param name Character scalar. One of \code{gw_distance_metric_names()}.
#' @param stop_on_error Logical. If \code{TRUE}, throw for unknown names; else \code{NULL}.
#' @return \code{list(p, theta, longlat)} or \code{NULL}.
#' @export
resolve_distance_metric <- function(name, stop_on_error = TRUE) {
  presets <- gw_distance_metric_presets()
  dm <- presets[[name]]
  if (is.null(dm) && isTRUE(stop_on_error)) {
    stop("Unknown `distance_metric`: ", name,
         "\nAvailable: ", paste(names(presets), collapse = ", "))
  }
  dm
}

#' List valid GW distance metric names
#' @return Character vector of valid preset names.
#' @export
gw_distance_metric_names <- function() {
  names(gw_distance_metric_presets())
}

#' gwkit: convenience wrappers for geographically weighted analysis
#'
#' @description
#' `gwkit` is a thin toolkit layer over `GWmodel` and related spatial packages
#' (`sf`, `terra`, `sp`). It bundles reproducible geographically weighted
#' workflows so that common tasks - fitting geographically weighted summary
#' statistics at points, polygons, or counties, sweeping distance-metric and
#' kernel settings, and reconciling the results into a single per-unit consensus
#' - can be run with a couple of calls. It wraps, rather than reimplements, the
#' underlying estimators.
#'
#' @section Main function families:
#' \describe{
#'   \item{Geographically weighted summary statistics}{`estimate_gwss()`.}
#'   \item{Geographically weighted spatial lag}{`estimate_gwlag()`.}
#'   \item{Geographically weighted regression (plain, fixed-effects, mean/variance)}{
#'     `estimate_gwr()`.}
#'   \item{Distance-metric presets}{`gw_distance_metric_presets()`,
#'     `gw_distance_metric_names()`, `resolve_distance_metric()`.}
#'   \item{Consensus across the kernel x distance-metric domain}{
#'     `gw_consensus_scalar()`, `gw_consensus_class()`.}
#' }
#'
#' @keywords internal
#' @import data.table
"_PACKAGE"
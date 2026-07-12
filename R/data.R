#' U.S. state agricultural land and crop-insurance panel
#'
#' A state x Census-of-Agriculture-year panel of agricultural land value/area,
#' crop-insurance coverage, and Federal Crop Insurance Program (FCIP) risk and
#' participation measures, assembled from the USDA NASS Census of Agriculture and
#' RMA extracts published in the `USFarmSafetyNetLab` data release. Bundled as
#' example data for the geographically weighted tools in gwkit (a real, gap-prone
#' state-level surface for `estimate_gwr()`, `estimate_gwlag()`, `estimate_gwss()`,
#' and the consensus tools). Join to a state polygon layer (e.g.
#' `urbnmapr::get_urbn_map("states", sf = TRUE)`) by state FIPS.
#'
#' The example analyses use it to ask where farmland trades at a *risk discount*
#' and whether FCIP offsets it: regress `log(ag_land_value_per_acre)` on the FCIP
#' base rate (`fcip_base_rate`, a risk proxy), participation (`fcip_adoption`),
#' and their interaction.
#'
#' @format A `data.table` with one row per (`state_code` x `census_year`), keyed
#'   on both, and the columns:
#' \describe{
#'   \item{state_code}{State FIPS code (integer).}
#'   \item{census_year}{Census of Agriculture year (integer).}
#'   \item{crop_insurance_acres}{Agricultural land enrolled in crop insurance (acres).}
#'   \item{cropland_acres}{Cropland (acres).}
#'   \item{ag_land_value}{Total value of agricultural land and buildings (USD).}
#'   \item{ag_land}{Total agricultural land (acres).}
#'   \item{fcip_base_rate}{Acreage-weighted FCIP base premium rate near the census
#'     year (a proxy for assessed production risk).}
#'   \item{fcip_acres}{FCIP net reported acres near the census year.}
#'   \item{fcip_adoption}{FCIP participation: `pmin(fcip_acres / cropland_acres, 1)`.}
#'   \item{ag_land_value_per_acre}{`ag_land_value / ag_land` (USD/acre).}
#' }
#'
#' @source USDA NASS Census of Agriculture and USDA RMA, via the piggyback release
#'   \href{https://github.com/ftsiboe/USFarmSafetyNetLab}{ftsiboe/USFarmSafetyNetLab}:
#'   `agCensusAcres.rds` and `agCensusInsurance.rds` (tag `nass_extracts`),
#'   `fcip_aph_base_rate.rds` (tag `adm_extracts`), and `sobscc_all.rds` (tag
#'   `sob`). Built by `data-raw/scripts/internal_datasets.R`.
#' @examples
#' data(us_state_ag_census)
#' us_state_ag_census[, participation := crop_insurance_acres / cropland_acres]
#' head(us_state_ag_census)
"us_state_ag_census"

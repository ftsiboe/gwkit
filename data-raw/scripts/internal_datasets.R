# ============================================================
# Internal dataset builder: us_state_ag_census
# ============================================================
# Assembles a state x Census-of-Agriculture-year panel of agricultural land and
# crop-insurance coverage, used as bundled example / reference data for the
# geographically weighted tools in gwkit (e.g. a real, gap-prone state-level
# surface to demonstrate estimate_gwss() smoothing / imputation).
#
# Source (USDA NASS Census of Agriculture + USDA RMA, via USFarmSafetyNetLab):
#   piggyback repo `ftsiboe/USFarmSafetyNetLab`
#     - agCensusAcres.rds      (tag nass_extracts) -> ag land value + ag land acres
#     - agCensusInsurance.rds   (tag nass_extracts) -> crop-insured + cropland acres (STATE)
#     - fcip_aph_base_rate.rds  (tag adm_extracts)  -> FCIP APH base rates
#     - sobscc_all.rds          (tag sob)           -> FCIP net reported acres (SOB)
#
# FORMAT: a data.table, one row per (state_code x census_year), keyed on both:
#   state_code             integer  State FIPS code.
#   census_year            integer  Census of Agriculture year.
#   crop_insurance_acres   numeric  Agricultural land enrolled in crop insurance (acres).
#   cropland_acres         numeric  Cropland (acres).
#   ag_land_value          numeric  Total value of agricultural land and buildings (USD).
#   ag_land                numeric  Total agricultural land (acres).
#   fcip_base_rate         numeric  Acreage-weighted FCIP base premium rate (risk proxy).
#   fcip_acres             numeric  FCIP net reported acres near the census year.
#   fcip_adoption          numeric  FCIP participation = pmin(fcip_acres/cropland_acres, 1).
#   ag_land_value_per_acre numeric  ag_land_value / ag_land (USD/acre).
#
# REBUILD: source("data-raw/scripts/internal_datasets.R")  # writes data/us_state_ag_census.rda
# (documented in R/data.R). See data-raw/examples/ for worked analyses.
# ============================================================

rm(list = ls(all = TRUE)); gc()
library(data.table)

# Large piggyback downloads can exceed the default 60s timeout on slow links.
# (gwkit itself no longer sets this globally on load - see R/zzz.R.)
options(timeout = max(600L, getOption("timeout")))

temp_dir <- tempdir()

# --- Download the raw NASS Census extracts -----------------------------------
for (ff in list(
  c("fcip_aph_base_rate.rds", "adm_extracts"),
  c("sobscc_all.rds", "sob"),
  c("agCensusInsurance.rds", "nass_extracts"),
  c("agCensusAcres.rds",     "nass_extracts"))) {
  piggyback::pb_download(
    file = ff[1], dest = temp_dir,
    repo = "ftsiboe/USFarmSafetyNetLab", tag = ff[2], overwrite = TRUE)
}

# --- Agricultural land value + acres, aggregated to state x census_year ------
agCensusAcres <- readRDS(file.path(temp_dir, "agCensusAcres.rds"))
agCensusAcres <- agCensusAcres[
  , .(ag_land_value = sum(ag_land_value, na.rm = TRUE),
      ag_land       = sum(ag_land,       na.rm = TRUE)),
  by = c("state_code", "census_year")]

# --- State-level crop-insured + cropland acres -------------------------------
agCensusInsurance <- readRDS(file.path(temp_dir, "agCensusInsurance.rds"))
agCensusInsurance <- agCensusInsurance[agg_level_desc %in% "STATE"]
agCensusInsurance <- agCensusInsurance[short_desc %in% c("AG LAND, CROP INSURANCE - ACRES", "AG LAND, CROPLAND - ACRES")]
agCensusInsurance <- agCensusInsurance[domain_desc %in% c("TOTAL")]
agCensusInsurance <- agCensusInsurance[domaincat_desc %in% c("NOT SPECIFIED")]
agCensusInsurance[
  , short_desc := factor(short_desc, c("AG LAND, CROP INSURANCE - ACRES", "AG LAND, CROPLAND - ACRES"),
                         c("crop_insurance_acres", "cropland_acres"))]
agCensusInsurance[, c("data_source","domaincat_desc","domain_desc","agg_level_desc","asd_code", "county_code") := NULL]
agCensusInsurance <- agCensusInsurance |> tidyr::spread(short_desc, value)
agCensusInsurance$census_year <- agCensusInsurance$year

# --- State-level fcip variables -------------------------------
fcip_rates <- readRDS(file.path(temp_dir, "fcip_aph_base_rate.rds"))
fcip_rates[, census_year := NA_real_]

sobscc_all <- readRDS(file.path(temp_dir, "sobscc_all.rds"))

sobscc_all <- sobscc_all[
  , .(fcip_acres = sum(net_reported_quantity, na.rm = TRUE)),
  by = c("state_code","county_code", "commodity_year","commodity_code")]

fcip_data <- as.data.table(dplyr::inner_join(sobscc_all, fcip_rates))

# Map each FCIP commodity_year to the nearest census wave (+/- 2 years).
for (yr in unique(agCensusInsurance$census_year)) {
  fcip_data[commodity_year %in% (yr-2):(yr+2), census_year := yr]
}

fcip_data <- fcip_data[
  , .(fcip_base_rate = weighted.mean(x=tau_adm, w=fcip_acres, na.rm = TRUE),
      fcip_acres = sum(fcip_acres, na.rm = TRUE)),
  by = c("state_code", "census_year")]

fcip_data <- fcip_data[!census_year %in% NA]

# --- Join and keep the documented analytic columns ---------------------------
us_state_ag_census <- as.data.table(dplyr::inner_join(
  agCensusInsurance, agCensusAcres, by = c("state_code", "census_year")))
us_state_ag_census[, year := NULL]
us_state_ag_census <- us_state_ag_census[
  , .(state_code, census_year, crop_insurance_acres, cropland_acres,
      ag_land_value, ag_land)]

# Augment with the FCIP risk/participation measures. LEFT join so the full
# Census panel is retained (state x year rows with no FCIP match carry NA fcip_*).
n_rows <- nrow(us_state_ag_census)
us_state_ag_census <- as.data.table(dplyr::left_join(
  us_state_ag_census, fcip_data, by = c("state_code", "census_year")))
n_na <- sum(is.na(us_state_ag_census$fcip_base_rate))
if (n_na > 0L)
  message(sprintf(
    "internal_datasets: %d of %d state x year rows have no FCIP match (fcip_* = NA).",
    n_na, n_rows))

us_state_ag_census[, fcip_adoption          := pmin(fcip_acres / cropland_acres, 1)]
us_state_ag_census[, ag_land_value_per_acre := ag_land_value / ag_land]

# Key the shipped object on the documented panel identifiers (kept through save).
data.table::setkey(us_state_ag_census, state_code, census_year)

# --- Persist as exported package data (data/us_state_ag_census.rda) ----------
usethis::use_data(us_state_ag_census, overwrite = TRUE)

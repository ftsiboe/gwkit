#' @name PACKAGE_GLOBALVARIABLES
#' @title global_names
#' @description A combined dataset for global_names
#' @format A list of global names.
#' @source Internal innovation
PACKAGE_GLOBALVARIABLES <-  strsplit(
  " . ..class_levels ..coords county_fips latitude_flag longitude_flag
    model_estimator polygon_fips value queen_agreement queen_order queen_value",
  "\\s+"
)[[1]]


#' @name PACKAGE_GLOBALVARIABLES
#' @title global_names
#' @description A combined dataset for global_names
#' @format A list of global names.
#' @source Internal innovation
PACKAGE_GLOBALVARIABLES <-  strsplit(
  ". ..class_levels ..final_cols ..keep_cols ALL complete.cases obs",
  "\\s+"
)[[1]]


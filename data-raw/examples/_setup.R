# ============================================================
# Shared setup for the gwkit example .Rmd files
# ============================================================
# Loads gwkit and the bundled ERS plotting framework. When knitted from inside
# the source tree it uses devtools::load_all() so the freshly built
# data/us_state_ag_census.rda is visible WITHOUT reinstalling the package;
# otherwise it falls back to the installed package.
#
# Prerequisite: build the dataset once (needs network + piggyback):
#   source("data-raw/scripts/internal_datasets.R")   # writes data/us_state_ag_census.rda
# ============================================================

.desc <- "../../DESCRIPTION"
if (file.exists(.desc) &&
    grepl("^Package: gwkit", readLines(.desc, n = 1L)) &&
    requireNamespace("devtools", quietly = TRUE)) {
  suppressMessages(devtools::load_all("../..", quiet = TRUE))   # dev tree: exposes data/
} else {
  library(gwkit)                                                # installed package
}

if (!exists("us_state_ag_census")) {
  utils::data("us_state_ag_census", package = "gwkit")
  if (!exists("us_state_ag_census"))
    stop("Dataset 'us_state_ag_census' not found. Build it first:\n",
         "  source(\"data-raw/scripts/internal_datasets.R\")")
}

source("_ers_framework.R")     # ers_theme(), gw_diverging_map(), state_class_map(), ...

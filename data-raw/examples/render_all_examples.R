#!/usr/bin/env Rscript
# ============================================================
# Render every gwkit example to github_document (.md)
# ============================================================
# One entry point that knits all six example .Rmd files in order. Each example
# sources _setup.R (which load_all()s the dev tree and the ERS framework), so the
# only prerequisite is that the dataset has been built once:
#
#   source("data-raw/scripts/internal_datasets.R")   # writes data/us_state_ag_census.rda
#
# Run it any of these ways:
#   Rscript data-raw/examples/render_all_examples.R
#   source("data-raw/examples/render_all_examples.R")           # from anywhere
#   (or setwd("data-raw/examples"); source("render_all_examples.R"))
# ============================================================

examples <- c(
  "01_estimate_gwr.Rmd",
  "02_estimate_gwlag.Rmd",
  "03_estimate_gwss.Rmd",
  "04_gw_consensus_scalar.Rmd",
  "05_gw_consensus_class.Rmd",
  "06_building_blocks.Rmd")

# --- locate the examples directory --------------------------------------------
# Works via Rscript, source(), OR pasting into the console: tries the detected
# script dir first, then getwd() and the usual sub-paths, and picks whichever
# actually contains the example files.
.script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- sub("^--file=", "", a[grepl("^--file=", a)])
  if (length(f)) return(dirname(normalizePath(f)))            # Rscript
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of)))        # source()
  NULL
}

.find_examples_dir <- function(files) {
  cands <- unique(Filter(Negate(is.null), c(
    .script_dir(),
    getwd(),
    file.path(getwd(), "data-raw", "examples"),   # run from package root
    file.path(getwd(), "examples"))))
  hit <- Filter(function(d) all(file.exists(file.path(d, files))), cands)
  if (length(hit)) normalizePath(hit[[1]]) else NULL
}

# --- preflight checks ---------------------------------------------------------
if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("Package 'rmarkdown' is required to render the examples.")

ex_dir <- .find_examples_dir(examples)
if (is.null(ex_dir))
  stop("Could not locate the example .Rmd files.\n",
       "Set the working directory to the package root (or to data-raw/examples) ",
       "and re-run, e.g.:\n",
       "  setwd(\"data-raw/examples\"); source(\"render_all_examples.R\")")

pkg_root <- normalizePath(file.path(ex_dir, "..", ".."), mustWork = FALSE)
if (!file.exists(file.path(pkg_root, "data", "us_state_ag_census.rda")))
  stop("Dataset not built yet. From the package root run:\n",
       "  source(\"data-raw/scripts/internal_datasets.R\")\n",
       "then re-run this script.")

# --- render each example (in the examples dir so _setup.R paths resolve) ------
ok <- character(0); failed <- character(0)
for (f in examples) {
  message("\n=== Rendering ", f, " ===")
  res <- tryCatch({
    rmarkdown::render(
      input       = file.path(ex_dir, f),
      output_format = "github_document",
      knit_root_dir = ex_dir,     # chunks (and source("_setup.R")) run here
      quiet       = TRUE)
    TRUE
  }, error = function(e) { message("  ! failed: ", conditionMessage(e)); FALSE })
  if (isTRUE(res)) ok <- c(ok, f) else failed <- c(failed, f)
}

# --- summary ------------------------------------------------------------------
message("\n============================================================")
message("Rendered ", length(ok), "/", length(examples), " examples.")
if (length(failed)) {
  msg <- paste0("Failed:\n  ", paste(failed, collapse = "\n  "))
  if (interactive()) stop(msg) else { message(msg); quit(status = 1L, save = "no") }
}
message("All examples rendered to .md in ", ex_dir)

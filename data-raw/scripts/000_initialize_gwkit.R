rm(list = ls(all = TRUE)); gc()

# Clean generated artifacts
unlink(c(
  "NAMESPACE",
  list.files("./data", full.names = TRUE),
  list.files("./man",  full.names = TRUE)
))

# if(toupper(as.character(Sys.info()[["sysname"]])) %in% "WINDOWS"){
#   source( file.path(dirname(dirname(getwd())),"code-library.R"))
#   list_function <- c(
#     file.path(codeLibrary,"plot/ers_theme.R"),
#     file.path(codeLibrary,"plot/plot_helpers.R")
#   )
#   file.copy(from= list_function, to = "R/", overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
# }

# Sanity pass through R/ sources: shows any non-ASCII characters per file
for (i in list.files("R", full.names = TRUE)) {
  print(paste0("********************", i, "********************"))
  tools::showNonASCIIfile(i)
}

# Rebuild documentation from roxygen comments
devtools::document()

# Check man pages only (faster than full devtools::check)
devtools::check_man()

# Build PDF manual into the current working directory
devtools::build_manual(path = getwd())

# Optional: run tests / full package check (uncomment when needed)
# devtools::test()
devtools::check()


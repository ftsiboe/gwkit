# Register non-standard-evaluation symbols (data.table `:=`/`.SD`, dplyr verbs)
# so that R CMD check does not flag "no visible binding for global variable".
# This runs at namespace load; PACKAGE_GLOBALVARIABLES is defined in
# PACKAGE_GLOBALVARIABLES.R (collated first).
#
# NOTE: gwkit deliberately does NOT mutate global session options (scipen,
# timeout, future.globals.maxSize, ...) on load. A package should not change a
# user's global state; set such options in your own script if you need them.
if (getRversion() >= "2.15.1") utils::globalVariables(PACKAGE_GLOBALVARIABLES)

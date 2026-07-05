# Install the packages required to develop and run CalciumInsights locally.
cran <- getOption("repos")
if (is.null(cran) || identical(unname(cran["CRAN"]), "@CRAN@")) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

required <- c(
  "config", "golem", "shiny", "shinyjs", "DT", "ggplot2", "dplyr",
  "vroom", "jsonlite", "pracma", "prospectr", "latex2exp",
  "gridExtra", "htmltools", "pkgload", "testthat", "remotes", "httpuv"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing)) {
  install.packages(missing, dependencies = TRUE)
} else {
  message("All required packages are already installed.")
}

message(
  "\nDependencies are ready. Open CalciumInsights.Rproj and run ",
  "shiny::runApp(), or source('dev/run_dev.R')."
)

#' Access files installed with the CalciumInsights package
#'
#' @param ... Character vectors identifying a path inside the package.
#' @noRd
app_sys <- function(...) {
  system.file(..., package = "CalciumInsights")
}

#' Read the golem configuration
#'
#' @param value Configuration key.
#' @param config Active configuration.
#' @param use_parent Whether config should search parent directories.
#' @param file Configuration file.
#' @noRd
get_golem_config <- function(
  value,
  config = Sys.getenv(
    "GOLEM_CONFIG_ACTIVE",
    Sys.getenv("R_CONFIG_ACTIVE", "default")
  ),
  use_parent = TRUE,
  file = app_sys("golem-config.yml")
) {
  config::get(
    value = value,
    config = config,
    file = file,
    use_parent = use_parent
  )
}

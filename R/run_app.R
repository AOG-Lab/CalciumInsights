#' Run the CalciumInsights Shiny application
#'
#' @param onStart Optional function called before the application starts.
#' @param options Named list of options passed to [shiny::shinyApp()].
#' @param enableBookmarking Bookmarking mode.
#' @param uiPattern URL pattern for the application.
#' @param ... Values stored as golem options.
#'
#' @return A Shiny application object, normally run interactively.
#' @export
run_app <- function(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
) {
  current_limit <- getOption("shiny.maxRequestSize", 0)
  requested_limit <- 200 * 1024 * 1024

  if (!is.numeric(current_limit) || length(current_limit) != 1L ||
      !is.finite(current_limit) || current_limit < requested_limit) {
    base::options(shiny.maxRequestSize = requested_limit)
  }

  golem::with_golem_options(
    app = shiny::shinyApp(
      ui = app_ui,
      server = app_server,
      onStart = onStart,
      options = options,
      enableBookmarking = enableBookmarking,
      uiPattern = uiPattern
    ),
    golem_opts = list(...)
  )
}

#' Add external resources to the application
#'
#' @import shiny
#' @importFrom golem favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {

  addResourcePath(
    prefix = "www",
    directoryPath = app_sys("app/www")
  )

  tags$head(
    favicon(),

    bundle_resources(
      path = app_sys("app/www"),
      app_title = "CalciumInsights"
    ),

    tags$style(HTML("
      body {
        background-color: #ffffff;
      }

      .navbar-brand {
        font-weight: 700;
        letter-spacing: 0.3px;
      }

      iframe {
        display: block;
      }

      .tab-content {
        padding-top: 15px;
      }
    "))
  )
}


#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {

  options(
    spinner.color = "#337ab7",
    spinner.color.background = "#ffffff",
    spinner.size = 2
  )

  tagList(
    golem_add_external_resources(),

    fluidPage(
      navbarPage(
        title = "CalciumInsights",
        id = "main_navbar",
        collapsible = TRUE,

        tabPanel(
          title = "Home",
          icon = icon("home"),
          tags$iframe(
            src = "www/index.html",
            height = "900px",
            width = "100%",
            style = "border: none;"
          )
        ),

        tabPanel(
          title = "FFT Denoising Analysis",
          icon = icon("chart-line"),
          mod_Denoising_data_ui("Denoising_data_1")
        )
      )
    )
  )
}

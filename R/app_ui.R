#' Add external resources used by the application
#'
#' @noRd
golem_add_external_resources <- function() {
  www_path <- app_sys("app", "www")

  if (!nzchar(www_path) || !dir.exists(www_path)) {
    stop(
      "CalciumInsights web resources were not found. ",
      "Reinstall the package from a complete source checkout.",
      call. = FALSE
    )
  }

  shiny::addResourcePath(prefix = "www", directoryPath = www_path)

  shiny::tags$head(
    shiny::tags$style(shiny::HTML("
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

#' CalciumInsights user interface
#'
#' @param request Internal parameter supplied by Shiny.
#' @noRd
app_ui <- function(request) {
  shiny::tagList(
    golem_add_external_resources(),

    shiny::navbarPage(
      title = "CalciumInsights",
      id = "main_navbar",
      theme = NULL,
      selected = "home",

      header = shiny::tagList(
        shinyjs::useShinyjs(),

        shiny::tags$script(shiny::HTML("
          $(document).on('shiny:connected', function() {

            // Hide the Home tab from the visible navigation menu
            $('a[data-value=\"home\"]').parent().hide();

            // Make the CalciumInsights brand behave as the Home button
            $('.navbar-brand').on('click', function(e) {
              e.preventDefault();

              // Activate the hidden Home tab
              $('a[data-value=\"home\"]').tab('show');

              // Update Shiny navbar input value
              Shiny.setInputValue('main_navbar', 'home', {priority: 'event'});
            });

          });
        "))
      ),

      shiny::tabPanel(
        title = "Home",
        value = "home",
        shiny::tags$iframe(
          src = "www/index.html?v=20260521",
          width = "100%",
          height = "950px",
          style = "border: none;"
        )
      ),

      shiny::tabPanel(
        title = "FFT + Baseline Analysis",
        value = "fft",
        mod_Denoising_data_fft_ui("fft_module")
      ),

      shiny::tabPanel(
        title = "Wavelet Ridgewalking",
        value = "wavelet",
        mod_wavelet_ridgewalking_ui("wavelet_module")
      )

      # The Method Comparison module remains intentionally inactive,
      # matching the current standalone Shiny app.
    )
  )
}

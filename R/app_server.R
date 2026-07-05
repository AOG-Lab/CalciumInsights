#' CalciumInsights server
#'
#' @param input,output,session Internal parameters supplied by Shiny.
#' @noRd
app_server <- function(input, output, session) {
  mod_Denoising_data_fft_server("fft_module")
  mod_wavelet_ridgewalking_server("wavelet_module")

  # The Method Comparison module remains intentionally inactive,
  # matching the current standalone Shiny app.
}

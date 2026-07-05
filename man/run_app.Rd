\name{run_app}
\alias{run_app}
\title{Run the CalciumInsights Shiny Application}
\usage{
run_app(
  onStart = NULL,
  options = list(),
  enableBookmarking = NULL,
  uiPattern = "/",
  ...
)
}
\arguments{
  \item{onStart}{Optional function called before the application starts.}
  \item{options}{Named list of options passed to \code{shiny::shinyApp()}.}
  \item{enableBookmarking}{Bookmarking mode.}
  \item{uiPattern}{URL pattern for the application.}
  \item{...}{Values stored as golem options.}
}
\value{
A Shiny application object, normally run interactively.
}
\description{
Starts the FFT and Wavelet Ridgewalking version of CalciumInsights.
}

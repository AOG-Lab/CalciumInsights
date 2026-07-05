options(golem.app.prod = FALSE)

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Install 'pkgload' before running the development app.", call. = FALSE)
}

pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)
run_app()

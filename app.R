# Run the package directly from a cloned or downloaded repository.
if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop(
    "Package 'pkgload' is required for repository mode. ",
    "Run source('install_dependencies.R') first.",
    call. = FALSE
  )
}

pkgload::load_all(export_all = FALSE, helpers = FALSE, quiet = TRUE)
run_app()

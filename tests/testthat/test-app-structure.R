test_that("core application functions are available", {
  expect_true(is.function(CalciumInsights::run_app))
  expect_true(is.function(getFromNamespace("app_ui", "CalciumInsights")))
  expect_true(is.function(getFromNamespace("app_server", "CalciumInsights")))
  expect_true(is.function(getFromNamespace(
    "mod_Denoising_data_fft_server",
    "CalciumInsights"
  )))
  expect_true(is.function(getFromNamespace(
    "mod_wavelet_ridgewalking_server",
    "CalciumInsights"
  )))
})

test_that("installed web resources are present", {
  www <- system.file("app", "www", package = "CalciumInsights")
  expect_true(dir.exists(www))
  expect_true(file.exists(file.path(www, "index.html")))
  expect_true(file.exists(file.path(www, "images", "logo.png")))
})

# test_log_normalization.R

context("log_normalization")

test_that("Test log_normalization", {
  dat.path <- system.file("extdata", "test_dat.rda", package = "SAVERg")
  out.path <- system.file("extdata", "lognorm_out.rda", package = "SAVERg")
  load(dat.path)
  load(out.path)
  expect_that(log_normalization(test_dat), equals(lognorm_out))
})

# [END]

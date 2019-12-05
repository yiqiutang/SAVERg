# test_preprocessing.R

context("preprocessing")

test_that("Test preprocessing", {
  dat.path <- system.file("extdata", "test_dat.rda", package = "SAVERg")
  out.path <- system.file("extdata", "prepro_out.rda", package = "SAVERg")
  load(dat.path)
  load(out.path)
  expect_that(preprocessing(test_dat), throws_error())
  expect_that(preprocessing(test_dat, percent = 0.1), equals(prepro_out))
})

# [END]

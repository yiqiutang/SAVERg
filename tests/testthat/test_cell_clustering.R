# test_cell_clustering.R

context("cell_clustering")

test_that("Test cell clustering", {
  expect_that(cell_clustering(ipsc_saver, -0.1), throws_error())
  expect_that(cell_clustering(ipsc_saver, 1.1), throws_error())
})

# [END]

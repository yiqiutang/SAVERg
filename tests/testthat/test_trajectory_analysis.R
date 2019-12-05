# test_trajectory_analysis.R

context("trajectory_analysis")

test_that("Test trajectory analysis", {
  deng_cellLabels <- factor(colnames(deng_saver),
                            levels=c('zygote', 'early 2-cell', 'mid 2-cell',
                                     'late 2-cell', '4-cell', '8-cell',
                                     '16-cell', 'early blastocyst',
                                     'mid blastocyst', 'late blastocyst'))
  expect_that(trajectory_analysis(deng_saver, deng_cellLabels),
              equals(list(cor.kendall = 0.8113, POS = 0.9434)))
})

# [END]

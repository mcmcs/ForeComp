test_that("dm.test.bt.fb returns empty output for unsupported cl", {
  set.seed(11)
  d = rnorm(20)

  out = dm.test.bt.fb(d, M = 3, cl = 0.01)
  expect_length(out, 0)
})

test_that("dm.test.im validates q and handles edge cases", {
  set.seed(12)
  d = rnorm(20)

  expect_error(dm.test.im(d, q = 1), "greater than or equal to 2")
  expect_error(dm.test.im(d, q = 2.5), "integer")
  expect_error(dm.test.im(d, q = NA), "integer")

  expect_warning(
    out_large <- dm.test.im(d, q = length(d) + 1),
    "exceeds sample size"
  )
  expect_true(is.na(out_large$rej))
  expect_true(is.na(out_large$stat))
  expect_true(is.na(out_large$pval))

  expect_warning(
    out_uneven <- dm.test.im(rnorm(10), q = 3),
    "not divisible"
  )
  expect_named(out_uneven, c("rej", "stat", "pval"), ignore.order = FALSE)
  expect_true(is.logical(out_uneven$rej))
  expect_true(is.numeric(out_uneven$stat))
  expect_true(is.numeric(out_uneven$pval))
})

test_that("Plot_Tradeoff validates required arguments", {
  expect_error(
    Plot_Tradeoff(data = mikedata, f1 = "f1", f2 = "f2"),
    "supplied altogether"
  )

  expect_error(
    Plot_Tradeoff(
      data = mikedata,
      f1 = "f1",
      f2 = "f2",
      y = "y",
      n_sim = 1,
      m_set = 1L,
      cl = 0.07,
      verbose = FALSE,
      no_m_label = TRUE
    ),
    "0.05 or 0.10"
  )

  expect_error(
    Plot_Tradeoff(
      data = mikedata,
      f1 = "f1",
      f2 = "f2",
      y = "y",
      n_sim = 1,
      m_set = c(1, 2.5),
      cl = 0.05,
      verbose = FALSE,
      no_m_label = TRUE
    ),
    "positive integers"
  )

  expect_error(
    Plot_Tradeoff(
      data = mikedata,
      f1 = "f1",
      f2 = "f2",
      y = "y",
      n_sim = 1,
      m_set = c(1, Inf),
      cl = 0.05,
      verbose = FALSE,
      no_m_label = TRUE
    ),
    "positive integers"
  )

  expect_error(
    Plot_Tradeoff(
      data = mikedata,
      f1 = "f1",
      f2 = "f2",
      y = "y",
      n_sim = 1,
      m_set = c(0, 1),
      cl = 0.05,
      verbose = FALSE,
      no_m_label = TRUE
    ),
    "positive integers"
  )
})

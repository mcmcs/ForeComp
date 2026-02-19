test_that("loss.diff.p computes expected values across loss types", {
  y = c(1, 2, 3)
  f1 = c(1, 1, 4)
  f2 = c(0, 2, 2)

  expect_equal(loss.diff.p(y, f1, f2, type = "quad"), c(-1, 1, 0))
  expect_equal(loss.diff.p(y, f1, f2, type = "abs"), c(-1, 1, 0))
  expect_equal(
    loss.diff.p(y, f1, f2, type = "check", tau = 0.25),
    c(-0.25, 0.25, 0.50)
  )

  e1 = y - f1
  e2 = y - f2
  expected_linex = (exp(0.5 * e1) - 0.5 * e1 - 1) - (exp(0.5 * e2) - 0.5 * e2 - 1)
  expect_equal(loss.diff.p(y, f1, f2, type = "linex", c = 0.5), expected_linex)
})

test_that("loss.diff.p returns NA for invalid loss type", {
  captured = capture.output(out <- loss.diff.p(c(1, 2), c(1, 2), c(1, 2), type = "invalid"))
  expect_true(length(captured) >= 1)
  expect_true(is.na(out))
  expect_length(out, 1)
})

test_that("default bandwidth rules match explicit bandwidth choices", {
  set.seed(7)
  d = rnorm(73)
  n = length(d)

  m_llsw = ceiling(1.3 * sqrt(n))
  m_nw = ceiling(4 * (n / 100)^(2 / 9))
  m_ci = floor(sqrt(n))
  b_ewc = max(1, floor(0.4 * n^(2 / 3)))
  m_wpe = floor(n^(1 / 3))

  out_bt_default = dm.test.bt(d, M = NA, Mopt = 2, cl = 0.05)
  out_bt_explicit = dm.test.bt(d, M = m_nw, cl = 0.05)
  expect_equal(out_bt_default$stat, out_bt_explicit$stat, tolerance = 1e-12)
  expect_equal(out_bt_default$pval, out_bt_explicit$pval, tolerance = 1e-12)
  expect_identical(out_bt_default$rej, out_bt_explicit$rej)

  out_bt_ci_default = dm.test.bt(d, M = NA, Mopt = 4, cl = 0.05)
  out_bt_ci_explicit = dm.test.bt(d, M = m_ci, cl = 0.05)
  expect_equal(out_bt_ci_default$stat, out_bt_ci_explicit$stat, tolerance = 1e-12)
  expect_equal(out_bt_ci_default$pval, out_bt_ci_explicit$pval, tolerance = 1e-12)
  expect_identical(out_bt_ci_default$rej, out_bt_ci_explicit$rej)

  out_fb_default = dm.test.bt.fb(d, M = NA, Mopt = 1, cl = 0.05)
  out_fb_explicit = dm.test.bt.fb(d, M = m_llsw, cl = 0.05)
  expect_equal(out_fb_default$stat, out_fb_explicit$stat, tolerance = 1e-12)
  expect_identical(out_fb_default$rej, out_fb_explicit$rej)

  out_ewc_default = dm.test.ewc.fb(d, B = NA, Bopt = 1, cl = 0.05)
  out_ewc_explicit = dm.test.ewc.fb(d, B = b_ewc, cl = 0.05)
  expect_equal(out_ewc_default$stat, out_ewc_explicit$stat, tolerance = 1e-12)
  expect_equal(out_ewc_default$pval, out_ewc_explicit$pval, tolerance = 1e-12)
  expect_identical(out_ewc_default$rej, out_ewc_explicit$rej)

  out_wpe_default = dm.test.wpe.fb(d, M = NA, Mopt = 1, cl = 0.05)
  out_wpe_explicit = dm.test.wpe.fb(d, M = m_wpe, cl = 0.05)
  expect_equal(out_wpe_default$stat, out_wpe_explicit$stat, tolerance = 1e-12)
  expect_equal(out_wpe_default$pval, out_wpe_explicit$pval, tolerance = 1e-12)
  expect_identical(out_wpe_default$rej, out_wpe_explicit$rej)
})

test_that("Bartlett DM tests handle large M without recycling warnings", {
  set.seed(77)
  d = rnorm(40)
  n = length(d)
  m_big = n
  lag_eff = n - 1

  covv = stats::acf(d, lag.max = lag_eff, type = "covariance", plot = FALSE)$acf[, , 1]
  expected_var = (covv[1] + 2 * sum((1 - ((1:lag_eff) / m_big)) * covv[-1])) / n
  expected_stat = mean(d) / sqrt(expected_var)

  out_bt_big = expect_warning(dm.test.bt(d, M = m_big, cl = 0.05), NA)
  expect_equal(out_bt_big$stat, expected_stat, tolerance = 1e-12)

  out_fb_big = expect_warning(dm.test.bt.fb(d, M = m_big, cl = 0.05), NA)
  expect_equal(out_fb_big$stat, expected_stat, tolerance = 1e-12)

  b = m_big / n
  crit = 1.9600 + 2.9694 * b + 0.4160 * b^2 - 0.5324 * b^3
  expect_identical(out_fb_big$rej, abs(out_fb_big$stat) > crit)
})

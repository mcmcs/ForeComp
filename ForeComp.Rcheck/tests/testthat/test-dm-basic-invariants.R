test_that("DM-style tests satisfy basic invariants and return contracts", {
  set.seed(42)
  d = rnorm(80)

  specs = list(
    list(
      name = "dm.test.r",
      fn = function(x) dm.test.r(x, h = 1, cl = 0.05),
      has_pval = TRUE,
      signed_stat = TRUE
    ),
    list(
      name = "dm.test.r.m",
      fn = function(x) dm.test.r.m(x, h = 1, cl = 0.05),
      has_pval = TRUE,
      signed_stat = TRUE
    ),
    list(
      name = "dm.test.bt",
      fn = function(x) dm.test.bt(x, M = 5, cl = 0.05),
      has_pval = TRUE,
      signed_stat = TRUE
    ),
    list(
      name = "dm.test.bt.fb",
      fn = function(x) dm.test.bt.fb(x, M = 5, cl = 0.05),
      has_pval = FALSE,
      signed_stat = TRUE
    ),
    list(
      name = "dm.test.ewc.fb",
      fn = function(x) dm.test.ewc.fb(x, B = 4, cl = 0.05),
      has_pval = TRUE,
      signed_stat = TRUE
    ),
    list(
      name = "dm.test.wpe.fb",
      fn = function(x) dm.test.wpe.fb(x, M = 4, cl = 0.05),
      has_pval = TRUE,
      signed_stat = TRUE
    ),
    list(
      name = "dm.test.im",
      fn = function(x) dm.test.im(x, q = 5, cl = 0.05),
      has_pval = TRUE,
      signed_stat = TRUE
    ),
    list(
      name = "dm.test.cnr.t",
      fn = function(x) dm.test.cnr.t(x, q = 5, cl = 0.05),
      has_pval = TRUE,
      signed_stat = FALSE
    ),
    list(
      name = "dm.test.cnr.w",
      fn = function(x) dm.test.cnr.w(x, q = 5, cl = 0.05),
      has_pval = TRUE,
      signed_stat = FALSE
    )
  )

  for (spec in specs) {
    out = spec$fn(d)
    out_scaled = spec$fn(10 * d)
    out_neg = spec$fn(-d)

    expected_names = c("rej", "stat")
    if (isTRUE(spec$has_pval)) {
      expected_names = c(expected_names, "pval")
    }

    expect_named(out, expected_names, ignore.order = FALSE, info = spec$name)
    expect_true(is.logical(out$rej) || is.na(out$rej), info = spec$name)
    expect_true(is.numeric(out$stat) || is.na(out$stat), info = spec$name)

    if (isTRUE(spec$has_pval)) {
      expect_true(is.numeric(out$pval) || is.na(out$pval), info = spec$name)
      expect_true(
        (is.na(out$pval) || (out$pval >= 0 && out$pval <= 1)),
        info = spec$name
      )
    }

    expect_equal(out$stat, out_scaled$stat, tolerance = 1e-10, info = spec$name)
    if (isTRUE(spec$has_pval)) {
      expect_equal(out$pval, out_scaled$pval, tolerance = 1e-12, info = spec$name)
    }
    expect_identical(out$rej, out_scaled$rej, info = spec$name)

    if (isTRUE(spec$signed_stat)) {
      expect_equal(out$stat, -out_neg$stat, tolerance = 1e-10, info = spec$name)
    } else {
      expect_equal(out$stat, out_neg$stat, tolerance = 1e-10, info = spec$name)
      expect_true(out$stat >= 0 || is.na(out$stat), info = spec$name)
    }
    if (isTRUE(spec$has_pval)) {
      expect_equal(out$pval, out_neg$pval, tolerance = 1e-12, info = spec$name)
    }
    expect_identical(out$rej, out_neg$rej, info = spec$name)
  }
})

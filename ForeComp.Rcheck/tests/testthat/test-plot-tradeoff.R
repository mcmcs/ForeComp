run_plot_tradeoff = function(...) {
  plot_file = tempfile(fileext = ".pdf")
  grDevices::pdf(plot_file)
  on.exit(grDevices::dev.off(), add = TRUE)
  suppressWarnings(Plot_Tradeoff(...))
}

test_that("Plot_Tradeoff returns consistent output and test statistics", {
  set.seed(123)
  out = run_plot_tradeoff(
    data = mikedata,
    f1 = "f1",
    f2 = "f2",
    y = "y",
    n_sim = 10,
    m_set = c(1L, 2L, 3L),
    cl = 0.05,
    verbose = FALSE,
    no_m_label = TRUE
  )

  expect_length(out, 2)
  expect_s3_class(out[[1]], "ggplot")
  expect_s3_class(out[[2]], "data.frame")

  plotting_data = out[[2]]
  required_columns = c(
    "M",
    "b_size_distortion",
    "b_power_loss",
    "v_hypothesis_test_b",
    "v_test_statistic_b",
    "v_hypothesis_test_dm",
    "v_test_statistic_dm"
  )
  expect_true(all(required_columns %in% names(plotting_data)))

  d = (mikedata$y - mikedata$f1)^2 - (mikedata$y - mikedata$f2)^2
  series_length = length(d)
  m_default = min(ceiling(1.3 * sqrt(series_length)), series_length - 1)
  expect_true(m_default %in% plotting_data$M)

  for (i in seq_len(nrow(plotting_data))) {
    m_choice = plotting_data$M[i]
    expected = dm.test.bt.fb(d, M = m_choice, cl = 0.05)

    expect_equal(
      plotting_data$v_test_statistic_b[i],
      as.numeric(expected$stat),
      tolerance = 1e-10
    )

    expected_label = if (is.na(expected$rej)) {
      NA_character_
    } else if (isTRUE(expected$rej)) {
      "cross"
    } else {
      "circle"
    }

    if (is.na(expected_label)) {
      expect_true(is.na(plotting_data$v_hypothesis_test_b[i]))
    } else {
      expect_identical(plotting_data$v_hypothesis_test_b[i], expected_label)
    }
  }

  expect_true(all(is.finite(plotting_data$b_size_distortion)))
  expect_true(all(plotting_data$b_size_distortion >= -1 & plotting_data$b_size_distortion <= 1))
  expect_true(all(is.finite(plotting_data$b_power_loss)))
  expect_true(all(plotting_data$b_power_loss >= -1 & plotting_data$b_power_loss <= 1))
})

test_that("Plot_Tradeoff applies na_handling options correctly", {
  tmp = mikedata[1:20, c("f1", "f2", "y")]
  tmp$f1[c(3, 7)] = NA_real_

  expect_error(
    run_plot_tradeoff(
      data = tmp,
      f1 = "f1",
      f2 = "f2",
      y = "y",
      n_sim = 3,
      m_set = c(1L, 2L),
      cl = 0.05,
      verbose = FALSE,
      no_m_label = TRUE,
      na_handling = "error"
    ),
    "missing/non-finite values"
  )

  set.seed(321)
  out_zero = run_plot_tradeoff(
    data = tmp,
    f1 = "f1",
    f2 = "f2",
    y = "y",
    n_sim = 3,
    m_set = c(1L, 2L),
    cl = 0.05,
    verbose = FALSE,
    no_m_label = TRUE,
    na_handling = "zero"
  )
  expect_s3_class(out_zero[[1]], "ggplot")

  set.seed(321)
  out_drop = run_plot_tradeoff(
    data = tmp,
    f1 = "f1",
    f2 = "f2",
    y = "y",
    n_sim = 3,
    m_set = c(1L, 2L),
    cl = 0.05,
    verbose = FALSE,
    no_m_label = TRUE,
    na_handling = "drop"
  )
  expect_s3_class(out_drop[[1]], "ggplot")
})

test_that("Plot_Tradeoff local regression check supports n_sim = 500", {
  skip_on_cran()

  set.seed(2026)
  out = run_plot_tradeoff(
    data = mikedata,
    f1 = "f1",
    f2 = "f2",
    y = "y",
    n_sim = 500,
    m_set = c(1L, 3L, 5L),
    cl = 0.05,
    verbose = FALSE,
    no_m_label = TRUE
  )
  plotting_data = out[[2]]

  expect_true(all(is.finite(plotting_data$b_size_distortion)))
  expect_true(all(plotting_data$b_size_distortion >= -1 & plotting_data$b_size_distortion <= 1))
  expect_true(all(is.finite(plotting_data$b_power_loss)))
  expect_true(all(plotting_data$b_power_loss >= -1 & plotting_data$b_power_loss <= 1))
})

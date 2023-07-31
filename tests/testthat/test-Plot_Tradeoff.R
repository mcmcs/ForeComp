

# =====================================================
# Test 1

test_that("Size distortion cannot be smaller than -5%, which is the pre-specified confidence level for this package", {

  spf_step1 <- TBILL$SPFfor_Step1
  na_spf_step1 <-is.na(spf_step1)
  spf_step1_imputed <- ifelse(na_spf_step1, 0, spf_step1)

  TBILL$SPFfor_Step1 <- spf_step1_imputed

  nc_step1 <- TBILL$NCfor_Step1
  na_nc_step1 <-is.na(nc_step1)
  nc_step1_imputed <- ifelse(na_nc_step1, 0, nc_step1)

  TBILL$NCfor_Step1 <- nc_step1_imputed

  realiz_step1 <- TBILL$Realiz1
  na_realiz_step1 <-is.na(realiz_step1)
  realiz_step1_imputed <- ifelse(na_realiz_step1, 0, realiz_step1)

  TBILL$Realiz1 <- realiz_step1_imputed

  expect_true(min(Plot_Tradeoff(data = TBILL,
                                f1   = "SPFfor_Step1",
                                f2   = "NCfor_Step1",
                                y    = "Realiz1")[[2]]$b_size_distortion) >= -0.05)

})

# =====================================================
# Test 2
test_that("Check the dimension of the output.", {

  spf_step1 <- TBILL$SPFfor_Step1
  na_spf_step1 <-is.na(spf_step1)
  spf_step1_imputed <- ifelse(na_spf_step1, 0, spf_step1)

  TBILL$SPFfor_Step1 <- spf_step1_imputed

  nc_step1 <- TBILL$NCfor_Step1
  na_nc_step1 <-is.na(nc_step1)
  nc_step1_imputed <- ifelse(na_nc_step1, 0, nc_step1)

  TBILL$NCfor_Step1 <- nc_step1_imputed

  realiz_step1 <- TBILL$Realiz1
  na_realiz_step1 <-is.na(realiz_step1)
  realiz_step1_imputed <- ifelse(na_realiz_step1, 0, realiz_step1)

  TBILL$Realiz1 <- realiz_step1_imputed

  m_set = c(1,2,3,4,5);

  output = Plot_Tradeoff(data = TBILL,
                f1   = "SPFfor_Step1",
                f2   = "NCfor_Step1",
                y    = "Realiz1",
                m_set = m_set);

  expect_true(nrow(output[[2]]) == length(m_set));

})

# =====================================================
# Test 3
test_that("Check that the M labels are plotted by default.", {

  spf_step1 <- TBILL$SPFfor_Step1
  na_spf_step1 <-is.na(spf_step1)
  spf_step1_imputed <- ifelse(na_spf_step1, 0, spf_step1)

  TBILL$SPFfor_Step1 <- spf_step1_imputed

  nc_step1 <- TBILL$NCfor_Step1
  na_nc_step1 <-is.na(nc_step1)
  nc_step1_imputed <- ifelse(na_nc_step1, 0, nc_step1)

  TBILL$NCfor_Step1 <- nc_step1_imputed

  realiz_step1 <- TBILL$Realiz1
  na_realiz_step1 <-is.na(realiz_step1)
  realiz_step1_imputed <- ifelse(na_realiz_step1, 0, realiz_step1)

  TBILL$Realiz1 <- realiz_step1_imputed

  output = Plot_Tradeoff(data = TBILL,
                f1   = "SPFfor_Step1",
                f2   = "NCfor_Step1",
                y    = "Realiz1"
                );

  expect_false(is.null(output[[1]]$layers[[3]][["constructor"]][[2]][["label"]]))
})

# =====================================================
# Test 4
test_that("Check points have no labels when no_m_label = TRUE.", {

  spf_step1 <- TBILL$SPFfor_Step1
  na_spf_step1 <-is.na(spf_step1)
  spf_step1_imputed <- ifelse(na_spf_step1, 0, spf_step1)

  TBILL$SPFfor_Step1 <- spf_step1_imputed

  nc_step1 <- TBILL$NCfor_Step1
  na_nc_step1 <-is.na(nc_step1)
  nc_step1_imputed <- ifelse(na_nc_step1, 0, nc_step1)

  TBILL$NCfor_Step1 <- nc_step1_imputed

  realiz_step1 <- TBILL$Realiz1
  na_realiz_step1 <-is.na(realiz_step1)
  realiz_step1_imputed <- ifelse(na_realiz_step1, 0, realiz_step1)

  TBILL$Realiz1 <- realiz_step1_imputed

  output = Plot_Tradeoff(data = TBILL,
                f1   = "SPFfor_Step1",
                f2   = "NCfor_Step1",
                y    = "Realiz1",
                no_m_label = TRUE
                );

  expect_error(is.null(output[[1]]$layers[[3]]$computed_mapping$label))

})

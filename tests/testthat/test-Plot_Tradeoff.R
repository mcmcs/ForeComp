
# =====================================================
# Test 1
test_that("Size distortion cannot be smaller than 5% confidence level", {
  expect_true(min(Plot_Tradeoff(data = TBILL,
                                f1   = "SPFfor_Step1",
                                f2   = "NCfor_Step1",
                                y    = "Realiz1")[[2]]$b_size_distortion) >= -0.05)

})

# =====================================================
# Test 2
test_that("Check the dimension of the output.", {

  m_set = c(1,2,3,4,5);

  output = Plot_Tradeoff(data = TBILL,
                f1   = "SPFfor_Step1",
                f2   = "NCfor_Step1",
                y    = "Realiz1",
                m_set = m_set);

  expect_true(nrow(output[[2]]) == length(m_set));

})

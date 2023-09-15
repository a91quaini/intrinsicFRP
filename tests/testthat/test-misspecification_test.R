test_that("Test HJMisspecificationTest", {

  factors = factors[,-1]
  returns = returns[,-1]

  expect_no_error(HJMisspecificationTest(returns, factors))

})

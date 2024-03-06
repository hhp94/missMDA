test_that("compare modded_imputeMCA and imputeMCA", {
  test_env <- new.env()
  data(vnf, envir = test_env)
  vnf <- get("vnf", envir = test_env)
  vnf0 <- vnf[-c(841,1073), ]
  # For these tests we expect number of different values to be less than 1 value
  # test 1: ncp = 0
  t0 <- fit_compare_fns(vnf0, ncp = 0)
  diff <- sum(t0$test != t0$modded_test)
  expect_lt(diff, 2L)
  
  # test 2: ncp = 1
  t1 <- fit_compare_fns(vnf0, ncp = 1)
  diff <- sum(t1$test != t1$modded_test)
  expect_lt(diff, 2L)
  
  # test 3: ncp = 2
  t2 <- fit_compare_fns(vnf0, ncp = 2)
  diff <- sum(t2$test != t2$modded_test)
  expect_lt(diff, 2L)
})

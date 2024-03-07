test_that("modded_estim_ncpMCA works", {
  vnf <- get_data_clean("vnf")
  vnf0 <- vnf[-c(841, 1073),]
  set.seed(1234)
  result <- estim_ncpMCA(
    vnf0,
    ncp.min = 2,
    ncp.max = 4,
    nbsim = 25,
    pNA = 0.01
  )

  set.seed(1234)
  result1 <- modded_estim_ncpMCA(
    vnf0,
    ncp.min = 2,
    ncp.max = 4,
    nbsim = 25,
    pNA = 0.01
  )

  diff_criterion <- all(abs(result$criterion - result1$criterion) < 0.001)
  expect_true(diff_criterion)
})

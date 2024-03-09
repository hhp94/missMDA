test_that("modded_estim_ncpMCA works", {
  # Make sure `modded_estim_ncpMCA` is the same as `estim_ncpMCA`
  vnf <- get_data_clean("vnf")
  vnf0 <- vnf[-c(841, 1073),]
  set.seed(1234)
  result <- estim_ncpMCA(
    vnf0,
    ncp.min = 2,
    ncp.max = 4,
    nbsim = 15,
    pNA = 0.005
  )

  set.seed(1234)
  result1 <- modded_estim_ncpMCA(
    vnf0,
    ncp.min = 2,
    ncp.max = 4,
    nbsim = 15,
    pNA = 0.005
  )

  expect_equal(result$ncp, result1$ncp)
  expect_true(all(near(result$criterion, result1$criterion)))
})

test_that("prodna1 and generate_k_fold", {
  vnf <- get_data_clean("vnf")
  vnf0 <- vnf[-c(841, 1073),]
  set.seed(1234)
  expect_error(prodna1(100, 100, 0), regex = "Increase")
  expect_error(prodna1(100, 100, 0.00000001), regex = "Increase")

  set.seed(1234)
  expect_error(generate_k_fold(vnf0, nrow(vnf0), ncol(vnf0), pNA = 1))
  expect_error(generate_k_fold(vnf0, nrow(vnf0), ncol(vnf0), pNA = 0.99, max_iterations = 1))
})

test_that("compare modded_imputeMCA and imputeMCA", {
  vnf <- get_data_clean("vnf")
  vnf0 <- vnf[-c(841,1073), ]
  # For these tests we expect number of different values to be less than 1 value
  # test 1: ncp = 0
  t0 <- fit_compare_fns(vnf0, ncp = 0)
  diff <- sum(t0$test != t0$modded_test)
  expect_equal(diff, 0)

  # test 2: ncp = 1
  t1 <- fit_compare_fns(vnf0, ncp = 1)
  diff <- sum(t1$test != t1$modded_test)
  expect_equal(diff, 0)

  # test 3: ncp = 2
  t2 <- fit_compare_fns(vnf0, ncp = 2)
  diff <- sum(t2$test != t2$modded_test)
  expect_equal(diff, 0)
})

test_that("test different modded_imputeMCA svd functions tall", {
  vnf <- get_data_clean("vnf")
  vnf0 <- vnf[-c(841,1073), ]

  for(i in seq_len(3)) {
    # Tall
    expect_equal(
      imputeMCA(vnf0, ncp = i)$completeObs,
      modded_imputeMCA(vnf0, ncp = i, svd_fns = "svd")$completeObs
    )
    expect_equal(
      modded_imputeMCA(vnf0, ncp = i, svd_fns = "svd")$completeObs,
      modded_imputeMCA(vnf0, ncp = i, svd_fns = "bootSVD")$completeObs
    )
  }
})

test_that("test different modded_imputeMCA svd functions wide", {
  sim_wide <- get_data_clean("sim_wide")
  # Wide
  expect_equal(
    imputeMCA(sim_wide$results, ncp = 2)$completeObs,
    modded_imputeMCA(sim_wide$results, ncp = 2, svd_fns = "svd")$completeObs
  )
  expect_equal(
    modded_imputeMCA(sim_wide$results, ncp = 2, svd_fns = "svd")$completeObs,
    modded_imputeMCA(sim_wide$results, ncp = 2, svd_fns = "bootSVD")$completeObs
  )
})

test_that("underflow caused by too many ncp", {
  vnf <- get_data_clean("vnf")
  vnf0 <- vnf[-c(841,1073), ]

  expect_error(imputeMCA(vnf0, ncp = 21))
  expect_no_error(modded_imputeMCA(vnf0, ncp = 21))
})

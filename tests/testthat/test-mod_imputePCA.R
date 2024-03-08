test_that("methyl_imputePCA", {
  orange <- get_data_clean("orange")
  dig <- 7
  for(i in seq_len(3)) {
    res.comp <- lapply(imputePCA(orange, ncp = i), round, digits = dig)
    expect_warning(methyl_imputePCA(as.matrix(orange), ncp = i))
    
    # Warning suppressed because this function is designed for fat methyldata
    res.comp1 <- suppressWarnings(lapply(methyl_imputePCA(as.matrix(orange), ncp = i), round, digits = dig))
    expect_identical(res.comp, res.comp1)
  }
})

test_that("methyl_imputePCA", {
  orange <- get_data_clean("orange")
  dig <- 7
  
  nb <- lapply(estim_ncpPCA(orange, ncp.min = 0, ncp.max = 4), round, digits = dig)
  # Warning suppressed because this function is designed for fat methyldata
  nb1 <- suppressWarnings(lapply(
    methyl_estim_ncpPCA(as.matrix(orange), ncp.min = 0, ncp.max = 4),
    round,
    digits = dig
  ))
  
  expect_identical(nb, nb1)
})

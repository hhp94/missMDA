test_that("methyl_imputePCA", {
  x <- new.env()
  data("orange", envir = x)
  orange <- get("orange", envir = x)
  dig <- 7
  res.comp <-
    imputePCA(orange, ncp = 2) |> lapply(round, digits = dig)
  # Warning suppressed because this function is designed for fat methyldata
  res.comp1 <- suppressWarnings(lapply(methyl_imputePCA(as.matrix(orange), ncp = 2), round, digits = dig))
  
  expect_identical(res.comp, res.comp1)
})

test_that("methyl_imputePCA", {
  x <- new.env()
  data("orange", envir = x)
  orange <- get("orange", envir = x)
  dig <- 7
  
  nb <-
    lapply(estim_ncpPCA(orange, ncp.min = 0, ncp.max = 4), round, digits = dig)
  # Warning suppressed because this function is designed for fat methyldata
  nb1 <- suppressWarnings(lapply(
    methyl_estim_ncpPCA(as.matrix(orange), ncp.min = 0, ncp.max = 4),
    round,
    digits = dig
  ))
  
  expect_identical(nb, nb1)
})

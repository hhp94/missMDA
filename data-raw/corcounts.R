## code to prepare `corcounts` dataset goes here
# sim_tall ####################
n_var.sim_tall <- 15
n_row.sim_tall <- 700
rate.sim_tall <- rep(1.25, times = n_var.sim_tall)

sim_fns <- function(n_row, n_var, rate, prop, cut = 4) {
  d <- simstudy::genCorGen(
    n = n_row,
    nvars = n_var,
    params1 = rate,
    dist = "poisson",
    rho = 0.7,
    corstr = "ar1",
    wide = TRUE
  )
  d$id <- NULL
  d <- as.matrix(d)
  d[d > cut] <- cut
  n_entry <- ncol(d) * nrow(d)
  index <- sample.int(n_entry, size = floor(n_entry * prop))
  truth <- d[index]
  d[index] <- NA
  results <- as.data.frame(d) |>
    dplyr::mutate(dplyr::across(.cols = dplyr::everything(), factor))
  list(results = results, truth = truth, index = index)
}

set.seed(1234)
sim_tall <-
  sim_fns(
    n_row = n_row.sim_tall,
    n_var = n_var.sim_tall,
    rate = rate.sim_tall,
    prop = 0.025
  )

sim_tall
# sim_wide ####################
n_var.sim_wide <- 300
n_row.sim_wide <- 100
rate.sim_wide <- rep(1.25, times = n_var.sim_wide)

set.seed(1234)
sim_wide <-
  sim_fns(
    n_row = n_row.sim_wide,
    n_var = n_var.sim_wide,
    rate = rate.sim_wide,
    prop = 0.01
  )

sim_wide
# apply(sim_tall$results, 2, as.numeric) |> PCA(scale.unit = F) |> summary()
# sim_tall$results |> as.numeric()|> princomp() |> summary()
usethis::use_data(sim_tall, overwrite = TRUE)
usethis::use_data(sim_wide, overwrite = TRUE)

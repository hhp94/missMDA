modded_imputeMCA_p <- function(don, ncp, row.w, verbose, printm, ...) {
  if (verbose) {
    cat(paste(printm, "...", sep = ""))
  }
  res <- modded_imputeMCA(don = don, ncp = ncp, row.w = row.w)
  return(res)
}

# Define the normalization function separately
normalize_cols <- function(x, col.suppr) {
  if (sum(x[1:col.suppr[1]]) != 1) {
    x[1:col.suppr[1]] <- x[1:col.suppr[1]] / sum(x[1:col.suppr[1]])
  }
  for (i in 2:length(col.suppr)) {
    x[(col.suppr[i - 1] + 1):(col.suppr[i])] <- x[(col.suppr[i - 1] + 1):(col.suppr[i])] /
      sum(x[(col.suppr[i - 1] + 1):col.suppr[i]])
  }
  return(x)
}

# Main function
normtdc <- function(tab.disj, data.na, ...) {
  # scale a fuzzy table
  tab.disj[tab.disj < 0] <- 0
  tab.disj[tab.disj > 1] <- 1

  col.suppr <- cumsum(vapply(data.na, function(x) {
    nlevels(x)
  }, numeric(1)))

  tab.disj <- t(apply(tab.disj, 1, normalize_cols, col.suppr = col.suppr))

  return(tab.disj)
}

draw <- function(tabdisj, Don, ...) {
  # draw from a scaled fuzzy table
  nbdummy <- vapply(Don, nlevels, numeric(1))
  vec <- c(0, cumsum(nbdummy))

  Donres <- Don
  for (i in seq_along(Don)) {
    Donres[, i] <- as.factor(levels(Don[, i])[
      apply(
        tabdisj[, (vec[i] + 1):vec[i + 1]],
        1,
        function(x) {
          sample(1:length(x), size = 1, prob = x)
        }
      )
    ])
    # Keep the order of the levels from original data
    Donres[, i] <- factor(Donres[, i], levels(Don[, i]))
  }

  return(don.imp = Donres)
}

modded_MIMCA <- function(X, nboot = 100, ncp, coeff.ridge = 1, threshold = 1e-06, maxiter = 1000, verbose = FALSE) {
  if (verbose) {
    temp <- if (coeff.ridge == 1) {
      "regularized"
    } else if (coeff.ridge == 0) {
      "EM"
    } else {
      paste("coeff.ridge=", coeff.ridge)
    }
    cat("Multiple Imputation using", temp, "MCA using", nboot, "imputed arrays", "\n")
  }
  validate_MCA(X)

  n <- nrow(X)
  Boot <- matrix(sample(1:n, size = nboot * n, replace = T), n, nboot)
  Weight <- matrix(1 / (n * 1000), n, nboot, dimnames = list(1:n, paste("nboot=", 1:nboot, sep = "")))
  Boot.table <- apply(Boot, 2, table)
  for (i in 1:nboot) Weight[names(Boot.table[[i]]), i] <- Boot.table[[i]]
  Weight <- sweep(Weight, 2, STATS = colSums(Weight), FUN = "/")
  Weight <- as.data.frame(Weight)
  if (future::nbrOfWorkers() == 1) {
    f_mapply <- mapply
  } else {
    f_mapply <- future.apply::future_mapply
  }
  res.imp <- f_mapply(
    # Each column of weight is implicitly passed to row.w. SMH
    Weight,
    FUN = modded_imputeMCA_p,
    MoreArgs = list(
      don = X,
      ncp = ncp,
      verbose = verbose
    ),
    printm = as.character(1:nboot),
    SIMPLIFY = FALSE,
    future.seed = TRUE
  )
  tdc.imp <- lapply(res.imp, \(x) x[["tab.disj"]])
  res.comp <- lapply(res.imp, \(x) x[["completeObs"]])
  tdc.norm <- f_mapply(
    FUN = normtdc,
    tab.disj = tdc.imp,
    data.na = res.comp,
    SIMPLIFY = F,
    future.seed = TRUE
  )
  X.imp <- f_mapply(
    FUN = draw,
    tabdisj = tdc.norm,
    Don = res.comp,
    SIMPLIFY = F,
    future.seed = TRUE
  )
  if (verbose) {
    cat("\ndone!\n")
  }
  return(list(res.MI = X.imp))
}

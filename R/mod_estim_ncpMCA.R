prodna1 <- function(x, noNA) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- rep(FALSE, n * p)
  n_miss <- floor(n * p * noNA)
  stopifnot("Increase pNA st floor(n * p * noNA) > 0" = n_miss > 0)
  NAloc[sample(n * p, n_miss)] <- TRUE
  x[matrix(NAloc, nrow = n, ncol = p)] <- NA
  return(x)
}

generate_k_fold <- function(don, pNA, max_iterations = 50) {
  stopifnot(0 < pNA, pNA < 1)
  compteur <- 1
  # Calculate the levels of each column before amputation
  levels_before <- sum(sapply(don, nlevels))

  while (compteur <= max_iterations) {
    donNA <- prodna1(don, pNA) # Amputate the data
    donNA <- droplevels(donNA) # Drop unused level
    # If we drop a level of any column by prodna, levels_after will be less than
    # levels_before
    levels_after <- sum(sapply(donNA, nlevels))

    if (levels_before == levels_after) {
      return(donNA)
    }

    compteur <- compteur + 1
  }

  stop(
    paste(
      "It is too difficult to suppress some cells.",
      "Maybe several categories are taken by only 1 individual.",
      "You should suppress these variables, try method.cv='loo', or lower pNA.",
      sep = "\n"
    )
  )
}

MCA_kfold_crit <- function(nbaxes, vrai.tab, const, don, donNA, threshold) {
  tab.disj.comp <- modded_imputeMCA(
    donNA,
    ncp = nbaxes, threshold = threshold
  )$tab.disj
  crit <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE) / const
  return(crit)
}

MCA_kfold_sim <- function(
    sim, vrai.tab, n_NA_init, don, pNA, ncp.min, ncp.max, threshold, pb, nbsim, verbose) {
  donNA <- generate_k_fold(don = don, pNA = pNA)
  const <- (sum(is.na(FactoMineR::tab.disjonctif(donNA))) - n_NA_init)
  res <- sapply(
    ncp.min:ncp.max,
    \(ncp) {
      MCA_kfold_crit(
        ncp,
        vrai.tab = vrai.tab, const = const, don = don, donNA = donNA, threshold = threshold
      )
    }
  )

  if (verbose) {
    setTxtProgressBar(pb, sim / nbsim * 100)
  }
  return(res)
}

#' Perform `modded_estim_ncpMCA` with K-fold CV
#'
#' @description
#' Estimates the number of dimensions to retain in Multiple Correspondence Analysis (MCA) using K-fold cross-validation.
#'
#' @details
#' This function performs K-fold cross-validation to estimate the number of dimensions to retain in MCA.
#' It introduces missing values into the input data.frame using the `prodna1` function.
#' For each simulation and each number of dimensions (ncp), it computes the imputed data.frame using the `modded_imputeMCA` function.
#' The sum of squared differences between the imputed and true tables is calculated and stored in the res matrix.
#' The criterion values are computed by taking the mean of the sum of squared differences across simulations for each ncp.
#' The estimated number of dimensions (ncp) is determined by selecting the ncp value with the minimum criterion value.
#'
#' @inheritParams modded_estim_ncpMCA
#'
#' @return a list containing the estimated ncp and the corresponding criterion values.
estim_ncpMCA_kfold <- function(don, ncp.min, ncp.max, nbsim, pNA, threshold, verbose) {
  stopifnot(
    "ncp.min should be integer" = ncp.min == as.integer(ncp.min),
    "ncp.max should be integer" = ncp.max == as.integer(ncp.max),
    "ncp.max should be >= ncp.min" = ncp.max >= ncp.min
  )
  vrai.tab <- tab.disjonctif(don)
  n_NA_init <- sum(is.na(vrai.tab))

  if (verbose) {
    pb <- txtProgressBar(min = 1 / nbsim * 100, max = 100, style = 3)
  }

  res <- vapply(
    seq_len(nbsim),
    function(x) {
      MCA_kfold_sim(
        sim = x,
        vrai.tab = vrai.tab, n_NA_init = n_NA_init, don = don, pNA = pNA,
        ncp.min = ncp.min, ncp.max = ncp.max, threshold = threshold, pb = pb,
        nbsim = nbsim, verbose = verbose
      )
    },
    FUN.VALUE = numeric(length(ncp.min:ncp.max))
  )

  if (verbose) {
    close(pb)
  }

  crit <- apply(res, 1, mean, na.rm = TRUE)
  names(crit) <- c(ncp.min:ncp.max)
  result <- list(ncp = as.integer(which.min(crit) + ncp.min - 1), criterion = crit)
  return(result)
}

#' Perform `modded_estim_ncpMCA` with LOO-CV
#'
#' @description
#' Estimates the number of dimensions to retain in Multiple Correspondence Analysis (MCA) using LOO-CV.
#'
#' @details
#' This function performs leave-one-out cross-validation to estimate the number of dimensions to retain in MCA.
#' It iterates over each element of the input data.frame and introduces a missing value at that position.
#' For each missing value and each number of dimensions (ncp), it computes the imputed value using the `modded_imputeMCA` function.
#' The imputed values are stored in the tab.disj.hat matrix.
#' The criterion values are computed by taking the mean of the squared differences between the imputed and true tables.
#' The estimated number of dimensions (ncp) is determined by selecting the ncp value with the minimum criterion value.
#'
#' @inheritParams modded_estim_ncpMCA
#'
#' @return a list containing the estimated ncp and the corresponding criterion values.
estim_ncpMCA_loo <- function(don, ncp.min, ncp.max, threshold, verbose) {
  vrai.tab <- tab.disjonctif(don)
  if (verbose) pb <- txtProgressBar(min = 0, max = 100, style = 3)
  crit <- NULL
  tab.disj.hat <- vrai.tab
  col.in.indicator <- c(0, sapply(don, nlevels))

  for (nbaxes in ncp.min:ncp.max) {
    for (i in 1:nrow(don)) {
      for (j in seq_len(ncol(don))) {
        if (!is.na(don[i, j])) {
          donNA <- as.matrix(don)
          donNA[i, j] <- NA
          if (!any(unlist(sapply(donNA, summary)) == 0)) {
            for (k in 1:ncol(donNA)) donNA[, k] <- as.factor(as.character(donNA[, k]))
            tab.disj.hat[i, (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])] <- modded_imputeMCA(donNA, ncp = nbaxes, threshold = threshold)$tab.disj[i, (cumsum(col.in.indicator)[j] + 1):(cumsum(col.in.indicator)[j + 1])]
          }
        }
      }
      if (verbose) setTxtProgressBar(pb, round((((1:length(ncp.min:ncp.max))[which(nbaxes == (ncp.min:ncp.max))] - 1) * nrow(don) + i) / (length(ncp.min:ncp.max) * nrow(don)) * 100))
    }
    crit <- c(crit, mean((tab.disj.hat - vrai.tab)^2, na.rm = TRUE))
  }
  if (verbose) close(pb)

  names(crit) <- c(ncp.min:ncp.max)
  return(list(ncp = as.integer(which.min(crit) + ncp.min - 1), criterion = crit))
}

#' Estimate the Best `ncp` for `modded_imputeMCA` with cross-validation
#'
#' @description
#' Estimates the number of dimensions to retain in Multiple Correspondence Analysis
#'
#' @details
#' This function estimates the number of dimensions to retain in MCA using the specified cross-validation method.
#' It calls either the `estim_ncpMCA_kfold` or `estim_ncpMCA_loo` function based on the method.cv argument.
#' This function first validates the input data.frame using the validate_MCA function.
#' It computes the true table of the input data.frame using the `tab.disjonctif` function.
#' This function then passes the necessary arguments to the selected cross-validation function and returns the result.
#'
#' @param don: The input data.frame as a matrix or data frame.
#' @param ncp.min: The minimum number of dimensions to consider (default = 0).
#' @param ncp.max: The maximum number of dimensions to consider (default = 5).
#' @param method.cv: The cross-validation method to use, either "Kfold" or "loo" (default = "Kfold").
#' @param nbsim: The number of simulations to perform (default = 100).
#' @param pNA: The proportion of missing values to introduce in each simulation (default = 0.05).
#' @param threshold: The threshold value for convergence in the imputeMCA function (default = 1e-4).
#' @param verbose: Logical value indicating whether to display progress information (default = TRUE).
#'
#' @return a list containing the estimated ncp and the corresponding criterion values.
#' @export
modded_estim_ncpMCA <- function(don, ncp.min = 0, ncp.max = 5,
                                method.cv = c("Kfold", "loo"), nbsim = 100, pNA = 0.05,
                                threshold = 1e-4, verbose = TRUE) {
  method.cv <- match.arg(method.cv, c("loo", "Kfold", "kfold", "LOO"), several.ok = T)[1]
  method.cv <- tolower(method.cv)

  validate_MCA(don)

  if (method.cv == "kfold") {
    return(estim_ncpMCA_kfold(don, ncp.min, ncp.max, nbsim, pNA, threshold, verbose))
  }

  if (method.cv == "loo") {
    return(estim_ncpMCA_loo(don, ncp.min, ncp.max, threshold, verbose))
  }
}

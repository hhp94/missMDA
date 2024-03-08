#' Insert NA into data.frame
#'
#' @description
#' Randomly introduces missing values (NA) into a given data.frame
#'
#' @details
#' This function creates a matrix of logical values (NAloc) with the same dimensions as the input data.frame.
#' For each column in the data.frame, it randomly selects a subset of row indices where the corresponding element is not an empty string.
#' The selected indices are used to assign TRUE values in the corresponding positions of the NAloc matrix.
#' The input data.frame is then updated by setting the elements corresponding to the TRUE values in NAloc as missing (NA).
#' The modified data.frame with missing values is returned.
#'
#' @param x The input data.frame
#' @param noNA The proportion of missing values to introduce, as a numeric value between 0 and 1.
#'
#' @return The input data.frame with randomly introduced missing values
prodna1 <- function(x, noNA) {
  n <- nrow(x)
  p <- ncol(x)
  NAloc <- matrix(FALSE, nrow = n, ncol = p)
  for (j in seq_len(p)) {
    na_indices <- sample(seq_len(n)[x[, j] != ""], floor(n * noNA))
    NAloc[na_indices, j] <- TRUE
  }
  x[NAloc] <- NA
  return(x)
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
  vrai.tab <- tab.disjonctif(don)
  res <- matrix(NA, ncp.max - ncp.min + 1, nbsim)
  if (verbose) { 
    pb <- txtProgressBar(min = 1 / nbsim * 100, max = 100, style = 3)
  }

  for (sim in seq_len(nbsim)) {
    donNA <- prodna1(don, pNA)
    for (i in seq_len(ncol(don))) {
      donNA[, i] <- as.factor(as.character(donNA[, i]))
    }
    
    for (nbaxes in ncp.min:ncp.max) {
      tab.disj.comp <- imputeMCA(donNA, ncp = nbaxes, threshold = threshold)$tab.disj
      res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE) / (sum(is.na(tab.disjonctif(donNA))) - sum(is.na(tab.disjonctif(don))))
    }
    
    if (verbose) { 
      setTxtProgressBar(pb, sim / nbsim * 100) 
    }
  }
  if (verbose) { 
    close(pb)
  }

  crit <- apply(res, 1, mean, na.rm = TRUE)
  names(crit) <- c(ncp.min:ncp.max)
  result <- list(ncp = as.integer(which.min(crit) + ncp.min - 1), criterion = crit)
  return(result)
}

# Sequential over both
# for (sim in seq_len(nbsim)) {
#   donNA <- prodna1(don, pNA)
#   for (i in seq_len(ncol(don))) {
#     donNA[, i] <- as.factor(as.character(donNA[, i]))
#   }
#   
#   for (nbaxes in ncp.min:ncp.max) {
#     tab.disj.comp <- imputeMCA(donNA, ncp = nbaxes, threshold = threshold)$tab.disj
#     res[nbaxes - ncp.min + 1, sim] <- sum((tab.disj.comp - vrai.tab)^2, na.rm = TRUE) / (sum(is.na(tab.disjonctif(donNA))) - sum(is.na(tab.disjonctif(don))))
#   }
#   
#   if (verbose) { 
#     setTxtProgressBar(pb, sim / nbsim * 100) 
#   }
# }

# Parallel over chunks of nbsim

# Parallel over the ncp
# res[, sim] <- future.apply::future_sapply(
#   ncp.min:ncp.max,
#   function(nbaxes) {
#     tab.disj.comp <- modded_imputeMCA(donNA, ncp = nbaxes, threshold = threshold)$tab.disj
#     crit <- sum((tab.disj.comp - vrai.tab) ^ 2, na.rm = TRUE) /
#       (sum(is.na(tab.disjonctif(donNA))) - sum(is.na(tab.disjonctif(don))))
#     return(crit)
#   },
#   USE.NAMES = FALSE,
#   future.packages = c("FactoMineR"),
#   future.seed = TRUE
# )



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

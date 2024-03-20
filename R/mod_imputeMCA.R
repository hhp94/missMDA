#' Modded imputeMCA by Refractoring
#'
#' Mod the imputeMCA function by omitting weights and sup variables
#'
#' @param don a data.frame with categorical variables containing missing values
#' @param ncp integer corresponding to the number of dimensions used to predict the missing entries
#' @param coeff.ridge 1 by default to perform the regularized imputeMCA algorithm;
#' useful only if method="Regularized". Other regularization terms can be
#' implemented by setting the value to less than 1 in order to regularized less
#' (to get closer to the results of the EM method) or more than 1 to regularized
#' more (to get closer to the results of the proportion imputation)
#' @param threshold the threshold for assessing convergence
#' @param seed integer, by default seed = NULL implies that missing values are
#' initially imputed by the proportion of the category for the categorical
#' variables coded with indicator matrices of dummy variables. Other values
#' leads to a random initialization
#' @param svd_fns which svd function to use "bootSVD"(default), "corpcor", or "svd"
#' @param maxiter integer, maximum number of iterations for the regularized iterative
#' MCA algorithm
#'
#' @return list of tab.disj and completeObs. See ?imputeMCA
#' @export
modded_imputeMCA <-
  function(don,
           ncp = 2,
           coeff.ridge = 1,
           threshold = 1e-6,
           seed = NULL,
           svd_fns = c("bootSVD", "corpcor", "svd"),
           maxiter = 1000) {
    svd_fns <- match.arg(svd_fns)
    svd_fns <- switch(svd_fns,
      bootSVD = bootSVD_wrap,
      corpcor = corpcor_wrap,
      svd = svd,
      stop("Not Supported")
    )

    ########## Debut programme principal
    validate_MCA(don)

    don <- droplevels(don)
    ncol_don <- ncol(don)
    nrow_don <- nrow(don)

    row.w <- rep(1 / nrow_don, nrow_don)

    if (ncp == 0) {
      tab.disj <- tab.disjonctif.prop(don, NULL)
      compObs <- find.category.1(don, tab.disj)
      return(list(tab.disj = tab.disj, completeObs = compObs))
    }

    # Convert to dummy matrix and get the coordinate of the missing values
    tab.disj.NA <- tab.disjonctif(don)
    tab.disj.comp <- tab.disjonctif.prop(don, seed, row.w = row.w)
    stopifnot("tab.disjonctif.prop should not return any NA" = !anyNA(tab.disj.comp))
    # Repeatedly calculated values
    ncol_tab <- ncol(tab.disj.comp)

    # Cannot have too many ncp for regularization to work
    k_j <- ncol_tab - ncol_don
    regularize <- TRUE
    if (ncp >= k_j) {
      regularize <- FALSE
    }

    ncp_vec <- seq_len(ncp)

    # Initialize
    hidden <- which(is.na(tab.disj.NA))
    tab.disj.rec.old <- tab.disj.comp

    continue <- TRUE
    nbiter <- 0

    while (continue) {
      nbiter <- nbiter + 1
      stopifnot("maxiter reached" = nbiter <= maxiter)
      # weights are always > 0, ncol_don always > 0. So if value is smaller than
      # zero then it's not because of ncol_don

      sum_weighted <- colSums(tab.disj.comp * row.w[1])
      M <- sum_weighted / ncol_don
      if (any(M < 0)) {
        stop(
          paste(
            "The algorithm fails to converge. Choose a number of components (ncp) less or equal than ",
            ncp - 1,
            " or a number of iterations (maxiter) less or equal than ",
            maxiter - 1,
            sep = ""
          )
        )
      }

      Z <- t(t(tab.disj.comp) / sum_weighted)
      Z <- t(t(Z) - colSums(Z * row.w[1]))
      Zscale <- t(t(Z) * sqrt(M))

      # Run svd based on the Zscale matrix
      # svd.Zscale <- FactoMineR::svd.triplet(Zscale, row.w = row.w, ncp = ncp)
      svd.Zscale <- modded_svd.triplet(Zscale, row.w = row.w, ncp = ncp, svd_fns = svd_fns)

      # Regularization is meaningless if ncp is too high. In which case moyeig to zero
      if (regularize == FALSE) {
        moyeig <- 0
      } else {
        # Regularizing
        if (nrow(don) > (k_j)) {
          moyeig <- mean(svd.Zscale$vs[-c(ncp_vec, (k_j + 1):ncol_tab)]^2)
        } else {
          moyeig <-
            mean(svd.Zscale$vs[-c(ncp_vec, ncol_tab:length(svd.Zscale$vs))]^2)
        }
        moyeig <- min(moyeig * coeff.ridge, svd.Zscale$vs[ncp + 1]^2)
      }
      eig.shrunk <- ((svd.Zscale$vs[ncp_vec]^2 - moyeig) / svd.Zscale$vs[ncp_vec])
      rec <- tcrossprod(
        t(t(svd.Zscale$U[, ncp_vec, drop = FALSE]) * eig.shrunk),
        svd.Zscale$V[, ncp_vec, drop = FALSE]
      )

      tab.disj.rec <- t(t(rec) / sqrt(M)) + matrix(1, nrow(rec), ncol(rec))
      tab.disj.rec <- t(t(tab.disj.rec) * sum_weighted)

      diff <- tab.disj.rec - tab.disj.rec.old
      diff[hidden] <- 0
      relch <- sum(diff^2 * row.w)
      tab.disj.rec.old <- tab.disj.rec
      tab.disj.comp[hidden] <- tab.disj.rec[hidden]
      continue <- (relch > threshold) && (nbiter < maxiter)
      # End of while loop
    }

    compObs <- find.category.1(don, tab.disj.comp)
    return(list(tab.disj = tab.disj.comp, completeObs = compObs))
  }

#' Assign Category Based on Most Likely tabdisj
#'
#' @param X original data.frame
#' @param tabdisj tabdisj object
#'
#' @return imputed data.frame with values based on most likely category
find.category.1 <- function(X, tabdisj) {
  nbdummy <- rep(1, ncol(X))
  is.quali <- which(!unlist(lapply(X, is.numeric)))
  nbdummy[is.quali] <-
    unlist(lapply(X[, is.quali, drop = FALSE], nlevels))
  vec <- c(0, cumsum(nbdummy))
  Xres <- X
  for (i in is.quali) {
    temp <-
      as.factor(levels(X[, i])[apply(tabdisj[, (vec[i] + 1):vec[i + 1]], 1, which.max)])
    Xres[, i] <- factor(temp, levels(X[, is.quali][, i]))
  }
  return(Xres)
}

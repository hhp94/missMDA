#' Modded imputeMCA using sparsesvd
#'
#' Mod the imputeMCA function by omitting weights and sup variables
#'
#' @param don a data.frame with categorical variables containing missing values
#' @param ncp integer corresponding to the number of dimensions used to predict the missing entries
#' @param method "Regularized" by default or "EM"
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
#' @param maxiter integer, maximum number of iterations for the regularized iterative
#' MCA algorithm
#'
#' @return list of tab.disj and completeObs. See ?imputeMCA
#' @export
modded_imputeMCA <-
  function(don,
           ncp = 2,
           method = c("Regularized", "EM"),
           coeff.ridge = 1,
           threshold = 1e-6,
           seed = NULL,
           maxiter = 1000) {

    moy.p <- function(V, poids) {
      res <- sum(V * poids, na.rm = TRUE) / sum(poids[!is.na(V)])
      res
    }
    
    ########## Debut programme principal
    stopifnot(is.data.frame(don), all(!c("tbl_df", "tbl") %in% class(don)))
    don_class <- unique(lapply(don, class))
    stopifnot(length(don_class) == 1, don_class[[1]] == "factor")
    miss_matrix <- is.na(don)
    stopifnot(
      "All missing column(s) detected" = all(colSums(miss_matrix) < nrow(miss_matrix)),
      "All missing row(s) detected" = all(rowSums(miss_matrix) < ncol(miss_matrix))
    )
    
    method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
    method <- tolower(method)
    don <- droplevels(don)
    
    # row.w <- rep(1 / nrow(don), nrow(don))
    row.w <- rep(1, nrow(don))
    
    if (ncp == 0) {
      tab.disj <- tab.disjonctif.prop(don, NULL, row.w = row.w)
      compObs <- find.category.1(don, tab.disj)
      return(list(tab.disj = tab.disj, completeObs = compObs))
    }
    
    # Convert to dummy matrix and get the coordinate of the missing values
    tab.disj.NA <- tab.disjonctif(don)
    tab.disj.comp <- tab.disjonctif.prop(don, seed, row.w = row.w)

    # Initialize
    hidden <- which(is.na(tab.disj.NA))
    tab.disj.rec.old <- tab.disj.comp

    continue <- TRUE
    nbiter <- 0

    while (continue) {
      nbiter <- nbiter + 1
      # weights are always > 0, ncol(don) always > 0. So if value is smaller than
      # zero then it's not because of ncol(don)
      
      # sum_weighted
      M <- apply(tab.disj.comp, 2, moy.p, row.w) / ncol(don)
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

      Z <- t(t(tab.disj.comp) / apply(tab.disj.comp, 2, moy.p, row.w))
      Z <- t(t(Z) - apply(Z, 2, moy.p, row.w))
      Zscale <- t(t(Z) * sqrt(M))

      # Run svd based on the Zscale matrix
      svd.Zscale <- FactoMineR::svd.triplet(Zscale, row.w = row.w, ncp = ncp)
      moyeig <- 0
      
      NcolZscale <- ncol(Zscale)
      
      # Regularizing
      if (nrow(don) > (NcolZscale - ncol(don))) {
        moyeig <- mean(svd.Zscale$vs[-c(seq_len(ncp), (NcolZscale - ncol(don) + 1):NcolZscale)] ^ 2)
      } else {
        moyeig <- mean(svd.Zscale$vs[-c(seq_len(ncp), NcolZscale:length(svd.Zscale$vs))]^2)
      }
      moyeig <- min(moyeig * coeff.ridge, svd.Zscale$vs[ncp + 1]^2)
      
      if (method == "em") {
        moyeig <- 0
      }

      eig.shrunk <- ((svd.Zscale$vs[seq_len(ncp)]^2 - moyeig) / svd.Zscale$vs[seq_len(ncp)])

      rec <- tcrossprod(
        t(t(svd.Zscale$U[, seq_len(ncp), drop = FALSE]) * eig.shrunk),
        svd.Zscale$V[, seq_len(ncp), drop = FALSE]
      )

      tab.disj.rec <- t(t(rec) / sqrt(M)) + matrix(1, nrow(rec), ncol(rec))
      tab.disj.rec <- t(t(tab.disj.rec) * apply(tab.disj.comp, 2, moy.p, row.w))

      diff <- tab.disj.rec - tab.disj.rec.old
      diff[hidden] <- 0
      relch <- sum(diff^2 * row.w)
      tab.disj.rec.old <- tab.disj.rec
      tab.disj.comp[hidden] <- tab.disj.rec[hidden]
      continue <- (relch > threshold) & (nbiter < maxiter)
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

#' Compare Results of imputeMCA and modded_imputeMCA
#'
#' @param df data frame
#' @param ... arguments passed to imputeMCA and modded_imputeMCA 
#'
#' @return test and modded_test
fit_compare_fns <- function(df, ...) {
    list(
      test = imputeMCA(df, ...)$completeObs,
      modded_test = modded_imputeMCA(df, ...)$completeObs
    )
  }
mod_impute <-
  function(X,
           ncp,
           scale = TRUE,
           method = NULL,
           threshold = 1e-6,
           init = 1,
           maxiter = 1000,
           coeff.ridge = 1,
           nrX,
           ncX,
           ...) {
    # Initialization
    stopifnot("X should be a matrix" = is.matrix(X))    
    nb.iter <- 1
    old <- Inf
    objective <- 0
    ncp <- min(ncp, ncol(X), nrow(X) - 1)
    ncp_seq <- seq_len(ncp)
    missing <- which(is.na(X))

    # mean.p <- apply(X, 2, mod_moy.p, row.w)
    mean.p <- colMeans(X, na.rm = TRUE)
    Xhat <- t(t(X) - mean.p) # Centering

    # et <- apply(Xhat, 2, mod_ec, row.w)
    et <- matrixStats::colSds(Xhat, na.rm = TRUE)

    if (scale) {
      Xhat <- t(t(Xhat) / et)
    }
    if (any(is.na(X))) {
      Xhat[missing] <- 0
    }
    ## random initialization
    if (init > 1) {
      Xhat[missing] <- rnorm(length(missing))
    }

    fittedX <- Xhat

    if (ncp == 0) {
      nb.iter <- 0
    }

    if (method == "em") {
      stop("Do not use method em to prevent overfit")
      sigma2 <- 0
    }

    # Refit PCA
    while (nb.iter > 0) {
      Xhat[missing] <- fittedX[missing]

      if (scale) {
        Xhat <- t(t(Xhat) * et)
      }
      Xhat <- t(t(Xhat) + mean.p)

      # mean.p <- apply(Xhat, 2, mod_moy.p, row.w)
      mean.p <- colMeans(Xhat, na.rm = TRUE)
      Xhat <- t(t(Xhat) - mean.p)

      # et <- apply(Xhat, 2, mod_ec, row.w)
      et <- matrixStats::colSds(Xhat, na.rm = TRUE)
      if (scale) {
        Xhat <- t(t(Xhat) / et)
      }

      # svd.res <- FactoMineR::svd.triplet(Xhat, ncp = ncp)
      svd.res <- mod_svd(Xhat, ncp = ncp, svd_fns = corpcor_wrap, ...)

      sigma2 <- nrX * ncX /
        min(ncX, nrX - 1) *
        sum(
          (svd.res$vs[-c(ncp_seq)]^2) /
            ((nrX - 1) * ncX - (nrX - 1) * ncp - ncX * ncp + ncp^2)
        )

      # sigma2 is the minimum of the eigen value of the ncp + 1 value.
      sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp + 1]^2)
      lambda.shrinked <- (svd.res$vs[ncp_seq]^2 - sigma2) / svd.res$vs[ncp_seq]

      fittedX <- tcrossprod(
        t(t(svd.res$U[, ncp_seq, drop = FALSE]) * lambda.shrinked),
        svd.res$V[, ncp_seq, drop = FALSE]
      )

      diff <- Xhat - fittedX
      diff[missing] <- 0
      objective <- sum(diff^2)
      criterion <- abs(1 - objective / old)
      old <- objective
      nb.iter <- nb.iter + 1
      if (!is.nan(criterion)) {
        if ((criterion < threshold) && (nb.iter > 5)) {
          nb.iter <- 0
        }
        if ((objective < threshold) && (nb.iter > 5)) {
          nb.iter <- 0
        }
      }
      if (nb.iter > maxiter) {
        nb.iter <- 0
        warning(paste("Stopped after ", maxiter, " iterations"))
      }
    }

    # End of While Loop

    if (scale) {
      Xhat <- t(t(Xhat) * et)
    }
    Xhat <- t(t(Xhat) + mean.p)
    completeObs <- X
    completeObs[missing] <- Xhat[missing]

    if (scale) {
      fittedX <- t(t(fittedX) * et)
    }
    fittedX <- t(t(fittedX) + mean.p)

    result <- list()
    result$completeObs <- completeObs
    result$fittedX <- fittedX

    return(result)
  }

#' [imputePCA()] for methylation data
#'
#' @return imputed matrix
#' @export methyl_imputePCA
methyl_imputePCA <- function(X,
                             ncp = 2,
                             scale = TRUE,
                             method = c("Regularized", "EM"),
                             coeff.ridge = 1,
                             threshold = 1e-6,
                             nb.init = 1,
                             maxiter = 1000,
                             ...) {
  #### Main program
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"),
    several.ok = T
  )[1]
  obj <- Inf
  method <- tolower(method)
  nrX <- nrow(X)
  ncX <- ncol(X)
  # Weight X for the SVD
  if (ncX < nrX) {
    warning("Methyldata should be fat")
  }
  
  if (ncp > min(nrX - 2, ncX - 1)) {
    stop("ncp is too large")
  }

  for (i in seq_len(nb.init)) {
    if (!any(is.na(X))) {
      message("X is completely observed")
      return(X)
    }
    # Fitting --------------------------------------------------------
    res.impute <-
      mod_impute(
        X,
        ncp = ncp,
        scale = scale,
        method = method,
        threshold = threshold,
        init = i,
        maxiter = maxiter,
        coeff.ridge = coeff.ridge,
        nrX = nrX,
        ncX = ncX,
        ...
      )
    cur_obj <- mean((res.impute$fittedX[!is.na(X)] - X[!is.na(X)])^2)
    if (cur_obj < obj) {
      res <- res.impute
      obj <- cur_obj
    }
  }
  return(res)
}

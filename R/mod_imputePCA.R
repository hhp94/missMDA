mod.tryCatch.W.E1 <- function(expr) {
  W <- NULL

  w.handler <- function(w) {
    W <<- w
    invokeRestart("muffleWarning")
  }

  return(
    list(
      val = withCallingHandlers(
        tryCatch(expr, error = function(e) {
          e
        }),
        warning = w.handler
      ),
      warning = W
    )
  )
}

# From dplyr
near <- function (x, y, tol = .Machine$double.eps^0.5)  {
  abs(x - y) < tol
}

mod_svd <- function(X, ncp = Inf, ...) {
  ncp <- min(ncp, nrow(X) - 1, ncol(X))

  # Weight X for the SVD
  if (ncol(X) < nrow(X)) {
    warning("Methyldata should be fat")
  }

  # Try performing SVD
  # Returned object is a list of 3 matrices from svd decomposition
  ## Note that svd return all d. irlba only calculate the required ones
  ## and leave the rest as NA

  # svd.usuelle <- mod.tryCatch.W.E1(svd(t(X), nu = ncp, nv = ncp))$val
  svd.usuelle <- mod.tryCatch.W.E1(fast.svd.wrap(X = t(X), ncp = ncp))$val

  if (names(svd.usuelle)[[1]] == "message") {
    warning("svd returned message")
    # If a message is returned, then rotate the matrix and do SVD again
    # svd.usuelle <- mod.tryCatch.W.E1(svd(X, nu = ncp, nv = ncp))$val
    svd.usuelle <- mod.tryCatch.W.E1(fast.svd.wrap(X = X, ncp = ncp))$val

    if (names(svd.usuelle)[[1]] == "d") {
      warning("rotated svd returned message")
      # Because the matrix is rotated, we have to swap the u and v
      aux <- svd.usuelle$u
      svd.usuelle$u <- svd.usuelle$v
      svd.usuelle$v <- aux
    } else {
      warning("Manually computing SVD. Should not be necessary")
      bb <- eigen(crossprod(t(X), t(X)), symmetric = TRUE)
      svd.usuelle <- vector(mode = "list", length = 3)
      svd.usuelle$d[svd.usuelle$d < 0] <- 0
      svd.usuelle$d <- sqrt(svd.usuelle$d)
      svd.usuelle$v <- bb$vec[, seq_len(ncp)]
      svd.usuelle$u <-
        t(t(crossprod(X, svd.usuelle$v)) / svd.usuelle$d[seq_len(ncp)])
    }
  }

  U <- svd.usuelle$v
  V <- svd.usuelle$u

  # Reverses the sign of the columns with sign = -1
  mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
  mult[near(mult, 0)] <- 1

  # Re multiply the U V matrix by the weight because it was weighted earlier
  U <- t(t(U) * mult)
  V <- t(t(V) * mult)

  # In real useage, the last couple of vs can be extremely small. In which case
  # the author is re-multiplying that with the values in the U and V matrix.
  vs <- svd.usuelle$d[seq_len(min(ncol(X), nrow(X) - 1))]
  num <- which(vs[seq_len(ncp)] < 1e-15)

  if (length(num) == 1) {
    U[, num] <- U[, num, drop = FALSE] * vs[num]
    V[, num] <- V[, num, drop = FALSE] * vs[num]
  }

  if (length(num) > 1) {
    U[, num] <- t(t(U[, num]) * vs[num])
    V[, num] <- t(t(V[, num]) * vs[num])
  }

  return(list(vs = vs, U = U, V = V))
}

fast.svd.wrap <- function(X, ncp) {
  s <- seq_len(ncp)
  d <- corpcor::fast.svd(m = X)
  return(
    list(
      d = d$d,
      u = d$u[, s, drop = FALSE],
      v = d$v[, s, drop = FALSE]
    )
  )
}

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
      svd.res <- mod_svd(Xhat, ncp = ncp, ...)

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

  if (ncp > min(nrow(X) - 2, ncol(X) - 1)) {
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

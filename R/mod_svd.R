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

mod_svd <- function(X, ncp = Inf, svd_fns, ...) {
  nrX <- nrow(X)
  ncX <- ncol(X)
  ncp <- min(ncp, nrX - 1, ncX)
  
  # Try performing SVD
  # Returned object is a list of 3 matrices from svd decomposition
  ## Note that svd return all d. irlba only calculate the required ones
  ## and leave the rest as NA
  
  # svd.usuelle <- mod.tryCatch.W.E1(svd(t(X), nu = ncp, nv = ncp))$val
  svd.usuelle <- mod.tryCatch.W.E1(svd_fns(X = t(X), ncp = ncp))$val
  
  if (names(svd.usuelle)[[1]] == "message") {
    warning("svd returned message")
    # If a message is returned, then rotate the matrix and do SVD again
    # svd.usuelle <- mod.tryCatch.W.E1(svd(X, nu = ncp, nv = ncp))$val
    svd.usuelle <- mod.tryCatch.W.E1(svd_fns(X = X, ncp = ncp))$val
    
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
  vs <- svd.usuelle$d[seq_len(min(ncX, nrX - 1))]
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



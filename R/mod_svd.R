tryCatch.W.E.1 <- function(expr) {
  W <- NULL
  w.handler <- function(w) {
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(
    value = withCallingHandlers(tryCatch(expr, error = function(e) e), warning = w.handler),
    warning = W
  )
}

perform_svd <- function(X, svd_fns, ncp) {
  svd.usuelle <- tryCatch.W.E.1(svd_fns(X, nu = ncp, nv = ncp))$val
  if (names(svd.usuelle)[[1]] == "message") {
    svd.usuelle <- tryCatch.W.E.1(svd_fns(t(X), nu = ncp, nv = ncp))$val
    if (names(svd.usuelle)[[1]] == "d") {
      aux <- svd.usuelle$u
      svd.usuelle$u <- svd.usuelle$v
      svd.usuelle$v <- aux
    } else {
      bb <- eigen(crossprod(X, X), symmetric = TRUE)
      svd.usuelle <- vector(mode = "list", length = 3)
      svd.usuelle$d[svd.usuelle$d < 0] <- 0
      svd.usuelle$d <- sqrt(svd.usuelle$d)
      svd.usuelle$v <- bb$vec[, seq_len(ncp)]
      svd.usuelle$u <- t(t(crossprod(t(X), svd.usuelle$v)) / svd.usuelle$d[seq_len(ncp)])
    }
  }
  return(svd.usuelle)
}

modded_svd.triplet <- function(X, row.w = NULL, svd_fns = NULL, ncp = Inf) {
  if (is.null(row.w)) row.w <- rep(1 / nrow(X), nrow(X))
  ncp <- min(ncp, nrow(X) - 1, ncol(X))
  row.w <- row.w / sum(row.w)
  X <- X * sqrt(row.w)

  if (ncol(X) < nrow(X)) {
    svd.usuelle <- perform_svd(X, svd_fns, ncp)
    U <- svd.usuelle$u
    V <- svd.usuelle$v

    mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
    mult[mult == 0] <- 1
    V <- t(t(V) * mult)
    U <- t(t(U) * mult)

    U <- U / sqrt(row.w)
  } else {
    svd.usuelle <- perform_svd(t(X), svd_fns, ncp)
    U <- svd.usuelle$v
    V <- svd.usuelle$u

    mult <- sign(as.vector(crossprod(rep(1, nrow(V)), as.matrix(V))))
    mult[mult == 0] <- 1
    V <- t(t(V) * mult)
    U <- t(t(U) * mult) / sqrt(row.w)
  }

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

  res <- list(vs = vs, U = U, V = V)
  return(res)
}

corpcor_wrap <- function(X, nu, nv, ncp = NULL) {
  if(is.null(ncp)) { ncp <- nu }
  s <- seq_len(ncp)
  d <- corpcor::fast.svd(X)
  return(
    list(
      d = d$d,
      u = d$u[, s, drop = FALSE],
      v = d$v[, s, drop = FALSE]
    )
  )
}

bootSVD_wrap <- function(X, nu, nv, ncp = NULL){
  if(is.null(ncp)) { ncp <- nu }
  s <- seq_len(ncp)
  d <- bootSVD::fastSVD(X)
  return(
    list(
      d = d$d,
      u = d$u[, s, drop = FALSE],
      v = d$v[, s, drop = FALSE]
    )
  )
}
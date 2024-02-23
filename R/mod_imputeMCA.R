#' Impute categorical dataset
#'
#' Impute the missing values of a categorical dataset using Multiple
#' Correspondence Analysis (MCA). Can be used as a preliminary step before
#' performing MCA on an incomplete dataset.
#'
#' Impute the missing entries of a categorical data using the iterative MCA
#' algorithm (method="EM") or the regularised iterative MCA algorithm
#' (method="Regularized"). The (regularized) iterative MCA algorithm first
#' consists in coding the categorical variables using the indicator matrix of
#' dummy variables. Then, in the initialization step, missing values are
#' imputed with initial values such as the proportion of the category for each
#' category using the non-missing entries. This imputation corresponds also to
#' using the algorithm with ncp=0 and is sometimes called in the literature the
#' "missing fuzzy average method". If the argument seed is set to a specific
#' value, a random initialization is performed: random values are drawn in such
#' a way that the constraint that the sum of the entries corresponding to one
#' individual and one variable is equal to one in the indicator matrix of dummy
#' variables.  The second step of the (regularized) iterative MCA algorithm
#' consists in performing MCA on the completed dataset. Then, it imputes the
#' missing values with the (regularized) reconstruction formulae of order ncp
#' (the fitted matrix computed with ncp components for the (regularized) scores
#' and loadings). These steps of estimation of the parameters via MCA and
#' imputation of the missing values using the (regularized) fitted matrix are
#' iterate until convergence. \cr We advice to use the regularized version of
#' the algorithm to avoid the overfitting problems which are very frequent when
#' there are many missing values. In the regularized algorithm, the singular
#' values of the MCA are shrinked. \cr The number of components ncp used in the
#' algorithm can be selected using the function ncpMCA. A small number of
#' components can also be seen as a way to regularize more and consequently may
#' be advices to get more stable predictions. \cr The output of the algorithm
#' can be used as an input of the MCA function of the FactoMineR package in
#' order to perform MCA on an incomplete dataset.
#'
#' @param don a data.frame with categorical variables containing missing values
#' @param ncp integer corresponding to the number of dimensions used to predict
#' the missing entries
#' @param method "Regularized" by default or "EM"
#' @param row.w row weights (by default, a vector of 1 for uniform row weights)
#' @param coeff.ridge 1 by default to perform the regularized imputeMCA
#' algorithm; useful only if method="Regularized". Other regularization terms
#' can be implemented by setting the value to less than 1 in order to
#' regularized less (to get closer to the results of the EM method) or more
#' than 1 to regularized more (to get closer to the results of the proportion
#' imputation)
#' @param threshold the threshold for assessing convergence
#' @param seed integer, by default seed = NULL implies that missing values are
#' initially imputed by the proportion of the category for the categorical
#' variables coded with indicator matrices of dummy variables. Other values
#' leads to a random initialization
#' @param maxiter integer, maximum number of iterations for the regularized
#' iterative MCA algorithm
#' @return \item{tab.disj}{The imputed indicator matrix; the observed values
#' are kept for the non-missing entries and the missing values are replaced by
#' the predicted ones.  The imputed values are real numbers but they but they
#' met the constraint that the sum of the entries corresponding to one
#' individual and one variable is equal to one. Consequently they can be seen
#' as degree of membership to the corresponding category. }
#' \item{completeObs}{The categorical imputed dataset; the observed values are
#' kept for the non-missing entries and the missing values are replaced by the
#' predicted ones. Missing values are imputed with the most plausible
#' categories according to the values in the tab.disj output}
#' @author Francois Husson \email{husson@@agrocampus-ouest.fr} and Julie Josse
#' \email{julie.josse@@polytechnique.edu}
#' @seealso \code{\link{estim_ncpMCA}},\cr
#' \href{https://www.youtube.com/watch?v=_Wa6R4PM9dY&list=PLnZgp6epRBbQzxFnQrcxg09kRt-PA66T_&index=1}{Video
#' showing how to perform MCA on an incomplete dataset}
#' @references Josse, J., Chavent, M., Liquet, B. and Husson, F. (2010).
#' Handling missing values with Regularized Iterative Multiple Correspondence
#' Analysis, Journal of Clcassification, 29 (1), pp. 91-116.\cr Josse, J. and
#' Husson, F. missMDA (2016). A Package for Handling Missing Values in
#' Multivariate Data Analysis. Journal of Statistical Software, 70 (1), pp 1-31
#' <doi:10.18637/jss.v070.i01>
#' @keywords models multivariate imputation
#' @examples
#' \dontrun{
#' data(vnf)
#' ## First the number of components has to be chosen
#' ##   (for the reconstruction step)
#' ## nb <- estim_ncpMCA(vnf,ncp.max=5) ## Time-consuming, nb = 4
#'
#' ## Impute the indicator matrix and perform a MCA
#' res.impute <- imputeMCA(vnf, ncp = 4)
#'
#' ## The imputed indicator matrix can be used as an input of the MCA function of the
#' ## FactoMineR package to perform the MCA on the incomplete data vnf
#' res.mca <- MCA(vnf, tab.disj = res.impute$tab.disj)
#'
#' ## With supplementary variables (var 11 to 14), impute the active ones
#' res.impute <- imputeMCA(vnf[, 1:10], ncp = 4)
#' res.mca <- MCA(vnf, tab.disj = res.impute$tab.disj, quali.sup = 11:14)
#' }
#'
#' @export
sum_score_imputeMCA <- function(
    don, ncp = 2, method = c("Regularized", "EM"), row.w = NULL, coeff.ridge = 1, 
    threshold = 1e-6, seed = NULL, maxiter = 1000
  ) {
  
  ind.sup <- NULL
  quanti.sup <- NULL
  quali.sup <- NULL
  
  moy.p <- function(V, poids) {
    res <- sum(V * poids, na.rm = TRUE) / sum(poids[!is.na(V)])
  }
  
  find.category <- function(X, tabdisj) {
    nbdummy <- rep(1, ncol(X))
    is.quali <- which(!unlist(lapply(X, is.numeric)))
    nbdummy[is.quali] <- unlist(lapply(X[, is.quali, drop = FALSE], nlevels))
    vec <- c(0, cumsum(nbdummy))
    Xres <- X
    for (i in is.quali) {
      temp <- as.factor(levels(X[, i])[apply(tabdisj[, (vec[i] + 1):vec[i + 1]], 1, which.max)])
      Xres[, i] <- factor(temp, levels(X[, is.quali][, i]))
    }
    return(Xres)
  }

  # Convert data frame of factors to matrix of dummies with reference level and rename them 
  tab.disjonctif.NA <- function(tab) {
    modalite.disjonctif <- function(i) {
      moda <- tab[, i]
      nom <- names(tab)[i]
      n <- length(moda)
      moda <- as.factor(moda)
      x <- matrix(0, n, length(levels(moda)))
      ind <- (1:n) + n * (unclass(moda) - 1)
      indNA <- which(is.na(ind))
      x[(1:n) + n * (unclass(moda) - 1)] <- 1
      x[indNA, ] <- NA
      if ((ncol(tab) != 1) & (levels(moda)[1] %in% c(1:nlevels(moda), "n", "N", "y", "Y"))) {
        dimnames(x) <- list(row.names(tab), paste(nom, levels(moda), sep = "."))
      } else {
        dimnames(x) <- list(row.names(tab), levels(moda))
      }
      return(x)
    }
    
    res <- lapply(1:ncol(tab), modalite.disjonctif)
    res <- as.matrix(data.frame(res, check.names = FALSE))
    
    return(res)
  }
  
  ########## Debut programme principal
  stopifnot("This function doesn't support tibbles" = all(!class(don) %in% c("tbl_df", "tbl")))
  stopifnot("This function only support data.frame with only factors"=all(apply(don, 2, is.factor)))
  stopifnot("There are rows with all missing in the data"=all(rowSums(is.na(don)) < ncol(don)))
  stopifnot("There should be more than 1 columns" = ncol(don) > 1) # because no sup variables
  stopifnot("There should be more than 1 rows" = nrow(don) > 1)
  
  method <- match.arg(method, c("Regularized", "regularized", "EM", "em"), several.ok = T)[1]
  method <- tolower(method)
  don <- droplevels(don)
  
  # For weighted svds. Remove for now.
  # if (is.null(row.w)) row.w <- rep(1 / nrow(don), nrow(don))
  # if (!is.null(ind.sup)) row.w[ind.sup] <- 1e-08 / length(row.w) # divide by the number of individuals to be sure that the weight will be small
  Vec <- rep("var", ncol(don))
  NLevels <- sapply(don, nlevels)
  
  # For supplementary variables. Remove for now
  # if (!is.null(quali.sup)) Vec[quali.sup] <- "quali.sup"
  # if (!is.null(quanti.sup)) {
  #   Vec[quanti.sup] <- "quanti.sup"
  #   NLevels[NLevels == 0] <- 1
  # }
  TabDisjMod <- rep(Vec, NLevels)

  if (ncp == 0) {
    tab.disj <- FactoMineR::tab.disjonctif.prop(don, NULL, row.w = row.w)
    compObs <- find.category(don, tab.disj)
    return(list(tab.disj = tab.disj, completeObs = compObs))
  }

  # Convert the factor data frame into one hot encoded matrix
  tab.disj.NA <- tab.disjonctif.NA(don)
  # Get index of missing value based on new data
  hidden <- which(is.na(tab.disj.NA))
  # Initiate the missing where the missing is calculated as proportion of the 
  # category for the categorical variables coded with indicator matrices of dummy variables
  tab.disj.comp <- FactoMineR::tab.disjonctif.prop(don, seed, row.w = row.w)
  
  # Set up for the while loop
  tab.disj.rec.old <- tab.disj.comp

  continue <- TRUE
  nbiter <- 0

  while (continue) {
    nbiter <- nbiter + 1
    if (length(quali.sup) > 0) tab.disj.comp[, TabDisjMod[TabDisjMod != "quanti.sup"] == "quali.sup"] <- tab.disj.comp[, TabDisjMod[TabDisjMod != "quanti.sup"] == "quali.sup"] * 1e-8
    M <- apply(tab.disj.comp, 2, moy.p, row.w) / ncol(don)
    if (any(M < 0)) stop(paste("The algorithm fails to converge. Choose a number of components (ncp) less or equal than ", ncp - 1, " or a number of iterations (maxiter) less or equal than ", maxiter - 1, sep = ""))

    Z <- t(t(tab.disj.comp) / apply(tab.disj.comp, 2, moy.p, row.w))
    Z <- t(t(Z) - apply(Z, 2, moy.p, row.w))
    Zscale <- t(t(Z) * sqrt(M))

    svd.Zscale <- FactoMineR::svd.triplet(Zscale, row.w = row.w, ncp = ncp)
    moyeig <- 0
    if (length(quanti.sup) + length(quali.sup) > 0) {
      NcolZscale <- sum(TabDisjMod == "var")
    } else {
      NcolZscale <- ncol(Zscale)
    }
    #  if (nrow(don)>(NcolZscale-ncol(don))) moyeig <- mean(svd.Zscale$vs[-c(1:ncp,(NcolZscale-ncol(don)-length(quali.sup)- length(quanti.sup)+1):NcolZscale)]^2)
    if (nrow(don) > (NcolZscale - ncol(don))) {
      moyeig <- mean(svd.Zscale$vs[-c(1:ncp, (NcolZscale - (ncol(don) - length(quali.sup) - length(quanti.sup)) + 1):NcolZscale)]^2)
    } else {
      moyeig <- mean(svd.Zscale$vs[-c(1:ncp, NcolZscale:length(svd.Zscale$vs))]^2)
    }
    moyeig <- min(moyeig * coeff.ridge, svd.Zscale$vs[ncp + 1]^2)
    if (method == "em") moyeig <- 0
    eig.shrunk <- ((svd.Zscale$vs[1:ncp]^2 - moyeig) / svd.Zscale$vs[1:ncp])

    if (ncp == 1) {
      rec <- tcrossprod(svd.Zscale$U[, 1] * eig.shrunk, svd.Zscale$V[, 1])
    } else {
      rec <- tcrossprod(t(t(svd.Zscale$U[, 1:ncp, drop = FALSE]) * eig.shrunk), svd.Zscale$V[, 1:ncp, drop = FALSE])
    }

    tab.disj.rec <- t(t(rec) / sqrt(M)) + matrix(1, nrow(rec), ncol(rec))
    tab.disj.rec <- t(t(tab.disj.rec) * apply(tab.disj.comp, 2, moy.p, row.w))

    diff <- tab.disj.rec - tab.disj.rec.old
    diff[hidden] <- 0
    relch <- sum(diff^2 * row.w)
    tab.disj.rec.old <- tab.disj.rec
    tab.disj.comp[hidden] <- tab.disj.rec[hidden]
    if (length(quali.sup) > 0) tab.disj.comp[, TabDisjMod[TabDisjMod != "quanti.sup"] == "quali.sup"] <- tab.disj.comp[, TabDisjMod[TabDisjMod != "quanti.sup"] == "quali.sup"] * 1e+08
    continue <- (relch > threshold) & (nbiter < maxiter)
  }
  if (is.null(quanti.sup)) {
    compObs <- find.category(don, tab.disj.comp)
  } else {
    compObs <- find.category(don[, -quanti.sup], tab.disj.comp)
    aux <- don
    aux[, -quanti.sup] <- compObs
    compObs <- aux

    Tabaux <- cbind.data.frame(tab.disj.comp, don[, quanti.sup, drop = FALSE])
    ordre <- c(which(TabDisjMod != "quanti.sup"), which(TabDisjMod == "quanti.sup"))
    tab.disj.comp <- Tabaux[, order(ordre)]
  }

  return(list(tab.disj = tab.disj.comp, completeObs = compObs))
}

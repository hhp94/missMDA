validate_MCA <- function(don) {
  stopifnot("X must be a data.frame"=is.data.frame(don) && all(!c("tbl_df", "tbl") %in% class(don)))
  don <- droplevels(don)
  pass_class <- vapply(lapply(don, class), \(x) {"factor" %in% x | "ordered" %in% x}, logical(1))
  stopifnot("Only factors are allowed in X"=all(pass_class))
  miss_matrix <- is.na(don)
  stopifnot("X must have more than 1 rows and cols"=ncol(don) > 1 && nrow(don) > 1)
  stopifnot(
    "All missing column(s) detected" = all(colSums(miss_matrix) < nrow(miss_matrix)),
    "All missing row(s) detected" = all(rowSums(miss_matrix) < ncol(miss_matrix))
  )
  single_level_cols <- vapply(don, \(x) { length(levels(x)) == 1 }, logical(1))
  stopifnot("Factors with only 1 level detected"=sum(single_level_cols) == 0)
}

get_data_clean <- function(name) {
  e <- new.env()
  utils::data(list = name, envir = e)
  return(get(name, envir = e))
}

#' Compare Results of imputeMCA and modded_imputeMCA
#'
#' @param df data frame
#' @param ... arguments passed to imputeMCA and modded_imputeMCA
#'
#' @return test and modded_test
fit_compare_fns <- function(df, ...) {
  args <- c(list(don = df), list(...))
  test_args <- args
  test_args[["svd_fns"]] <- NULL
  test <- do.call("imputeMCA", args = test_args)$completeObs

  modded_test_args <- args
  modded_test_args[["row.w"]] <- NULL
  modded_test <- do.call("modded_imputeMCA", args = modded_test_args)$completeObs
  return(list(test = test, modded_test = modded_test))
}

# From dplyr
near <- function (x, y, tol = .Machine$double.eps^0.5) {
  abs(x - y) < tol
}

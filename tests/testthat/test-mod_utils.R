test_that("validate_MCA handles valid inputs correctly", {
  # Create valid test data
  valid_df <- data.frame(
    col1 = factor(c("a", "b", "c")),
    col2 = factor(c("x", "y", "x")),
    col3 = factor(c("m", "n", "m"))
  )

  # Should not throw any errors
  expect_no_error(validate_MCA(valid_df))
})

test_that("validate_MCA rejects non-data.frame inputs", {
  # Test matrix
  mat <- matrix(1:4, nrow = 2)
  expect_error(validate_MCA(mat), "X must be a data.frame")

  # Test tibble
  tbl <- data.frame(a = factor(c("x", "y")), b = factor(c("m", "n")))
  class(tbl) <- c("tbl", class(tbl))
  expect_error(validate_MCA(tbl), "X must be a data.frame")
})

test_that("validate_MCA rejects non-factor columns", {
  # Mixed column types
  mixed_df <- data.frame(
    a = factor(c("x", "y")),
    b = c(1, 2),
    stringsAsFactors = FALSE
  )
  expect_error(validate_MCA(mixed_df), "Only factors are allowed in X")

  # Character columns
  char_df <- data.frame(
    a = c("x", "y"),
    b = c("m", "n"),
    stringsAsFactors = FALSE
  )
  expect_error(validate_MCA(char_df), "Only factors are allowed in X")
})

test_that("validate_MCA checks dimensions", {
  # Single row
  one_row <- data.frame(
    a = factor("x"),
    b = factor("y")
  )
  expect_error(validate_MCA(one_row), "X must have more than 1 rows and cols")

  # Single column
  one_col <- data.frame(
    a = factor(c("x", "y"))
  )
  expect_error(validate_MCA(one_col), "X must have more than 1 rows and cols")
})

test_that("validate_MCA handles missing values correctly", {
  # Missing column
  na_col_df <- data.frame(
    a = factor(c("x", "y", NA, NA)),
    b = factor(c(NA, NA, NA, NA))
  )
  expect_error(validate_MCA(na_col_df), "All missing column\\(s\\) detected")

  # Missing row
  na_row_df <- data.frame(
    a = factor(c(NA, "y", "z")),
    b = factor(c(NA, "m", "n")),
    c = factor(c(NA, "p", "q"))
  )
  expect_error(validate_MCA(na_row_df), "All missing row\\(s\\) detected")

  # Partial missing values (should pass)
  partial_na_df <- data.frame(
    a = factor(c(NA, "y", "z")),
    b = factor(c("x", NA, "n")),
    c = factor(c("p", "q", NA))
  )
  expect_no_error(validate_MCA(partial_na_df))
})

test_that("validate_MCA checks for single-level factors", {
  # Single level factor
  single_level_df <- data.frame(
    a = factor(c("x", "x", "x")),
    b = factor(c("y", "z", "y"))
  )
  expect_error(validate_MCA(single_level_df), "Factors with only 1 level detected")

  # Multiple single level factors
  multi_single_level_df <- data.frame(
    a = factor(c("x", "x", "x")),
    b = factor(c("y", "y", "y")),
    c = factor(c("z", "z", "z"))
  )
  expect_error(validate_MCA(multi_single_level_df), "Factors with only 1 level detected")
})

test_that("validate_MCA properly handles droplevels", {
  # Create data frame with unused levels
  df_with_unused_levels <- data.frame(
    a = factor(c("x", "y"), levels = c("x", "y", "z")),
    b = factor(c("m", "n"), levels = c("m", "n", "o"))
  )

  # Should not throw any errors
  expect_no_error(validate_MCA(df_with_unused_levels))

  # Verify that levels are dropped
  result <- droplevels(df_with_unused_levels)
  expect_equal(levels(result$a), c("x", "y"))
  expect_equal(levels(result$b), c("m", "n"))
})

test_that("validate_MCA works with different factor encodings", {
  # Ordered factors
  ordered_df <- data.frame(
    a = ordered(c("low", "med", "high")),
    b = ordered(c("small", "large", "medium"))
  )
  expect_no_error(validate_MCA(ordered_df))

  # Mixed ordered and unordered factors
  mixed_ordered_df <- data.frame(
    a = ordered(c("low", "med", "high")),
    b = factor(c("x", "y", "z"))
  )
  expect_no_error(validate_MCA(mixed_ordered_df))
})

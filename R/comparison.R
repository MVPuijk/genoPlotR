################################################################################
# comparison class and methods
################################################################################

#' comparison class and methods
#' 
#' `comparison` objects are collections of similarities between two DNA 
#' segments. These functions are class functions to create, convert, and test
#' `comparison` objects. 
#' 
#' @details
#' A `comparison` object requires at least the `start1`, `end1`, 
#' `start2`, and `end2` columns. A color can be provided as well using the `col`
#' column. Additional numeric columns can be used for color-coding using
#' `[gradient_color_scheme]`, or the `"gradient"` color scheme in the
#' `global_color_scheme` argument in `[plot_gene_map]`.
#' 
#' These objects inherit from `data.table`, and
#' can therefore be manipulated using `data.table` syntax.
#' 
#' `is.comparison` returns `TRUE` if the object tested is a
#' `comparison` object.
#' 
#' @name comparison
#' @returns `as.comparison` and `comparison` both return a 
#' `comparison` object.
#' @returns `is.comparison` returns a logical.
#' @export
#' 
#' @param x An object to be coerced. Can be a `list`, `data.frame` or 
#' `data.table` object. See details for necessary columns.
#' @param df Same as `x`.
#' 
#' @author Lionel Guy
#' 
#' @seealso [read_comparison], [print_comparison], [trim.comparison],
#' [reverse.comparison], [update_comparisons]
#' 
#' @examples
#' ## Get some values
#' starts1 <- c(2, 1000, 1050)
#' ends1 <- c(600, 800, 1345)
#' starts2 <- c(50, 800, 1200)
#' ends2 <- c(900, 1100, 1322)
#' 
#' ## From a data.frame
#' comparison1 <- as.comparison(data.frame(start1 = starts1, end1 = ends1,
#'                                         start2 = starts2, end2 = ends2))
#' print_comparison(comparison1)
#' is.comparison(comparison1)
#' is.data.frame(comparison1)
#' print_comparison(data.frame(start1 = starts1, end1 = ends1,
#'                             start2 = starts2, end2 = ends2))
#' 
#' ## From a list
#' print_comparison(list(start1 = starts1, end1 = ends1,
#'                       start2 = starts2, end2 = ends2))
#' 
#' 
#' ## Printing out a comparison like this can occasionally throw an error:
#' \dontrun{
#' comparison1
#' print(comparison1)
#' }
#' ## This can happen when print.comparison is loaded from the testthat package
#' ## To avoid this, use print_comparison as shown above
#' 
comparison <- function(x) {
  # if data.frame
  if (is.data.table(x)) {
    return(as.comparison(x))
  } else if (is.data.frame(x)) {
    return(as.comparison(as.data.table(x)))
  } else if (is.list(x)) {
    return(as.comparison(as.data.table(x)))
  } else {
    stop("Cannot coerce class ", class(x), " to comparison")
  }
}

#' @name comparison
#' @export
#' 
as.comparison <- function(df) {
  # check for class comparison, list, df, dt
  if (is.comparison(df)) return(df)
  if (is.data.frame(df) && !is.data.table(df)) df <- as.data.table(df)
  if (is.list(df) && !is.data.table(df)) df <- as.data.table(df)
  if (is.data.table(df)) {
    # check for columns
    names <- c("start1", "end1", "start2", "end2")
    if (!identical(names(df)[1:4], names)) {
      stop("Column names must start with 'start1', 'end1', 'start2', 'end2'")
    }
    if (!all(sapply(df, function(i) all(!is.null(i))))) {
      stop("NULL values not allowed")
    }
    if (!is.numeric(df$start1) | 
        !is.numeric(df$end1) |
        !is.numeric(df$start2) | 
        !is.numeric(df$end2)
        ) {
      stop("'start1', 'end1', 'start2', and 'end2' must be numeric")
    }
    # add direction column
    df$direction <- ifelse(
      sign(df$start1 - df$end1) * sign(df$start2 - df$end2) > 0,
      1,
      -1
    )
    # check color (character and default)
    if (is.factor(df$col)) df$col <- as.character(df$col)
  } else {
    stop("Unable to handle this format")
  }
  class(df) <- c("comparison", class(df))
  df
}

#' @name comparison
#' @export
#' 
#' @param comparison An object to test.
#'
is.comparison <- function(comparison) {
  inherits(comparison, "comparison")
}

#' @rdname range
#' @export
#' 
range.comparison <- function(x, overall = TRUE, ...) {
  comparison <- x
  if (overall) {
    range <- range(comparison$start1, comparison$end1,
                   comparison$start2, comparison$end2,
                   na.rm = FALSE
                   )
  } else {
    xlim1 <- range(comparison$start1, comparison$end1, na.rm = FALSE)
    xlim2 <- range(comparison$start2, comparison$end2, na.rm = FALSE)
    range <- data.frame(xlim1 = xlim1, xlim2 = xlim2)
  }
  range
}

#' Print out a comparison
#' 
#' Converts objects to a `data.table` object and prints it out
#' 
#' @details
#' The object provided is printed out as if it were `data.table` object without
#' explicitly converting it. As such, this function returns a `data.table`
#' object, not a `comparison`. 
#' 
#' This function was written to avoid potential errors, as
#' using the generic print function on a `comparison` object can cause an error
#' if the `print.comparison` method from the `testthat` package is loaded.
#' Outside of this situation, printing out a `comparison` object as normal
#' should function properly.
#' 
#' @returns A `data.table` object.
#' @export
#' 
#' @param comparison A `comparison` object or any other object that can be
#' printed out as a `data.table` object without further conversions.
#' 
#' @author Mike Puijk
#' 
#' @seealso [comparison]
#' 
#' @examples
#' ## Get some values
#' starts1 <- c(2, 1000, 1050)
#' ends1 <- c(600, 800, 1345)
#' starts2 <- c(50, 800, 1200)
#' ends2 <- c(900, 1100, 1322)
#' 
#' ## From a data.frame
#' comparison1 <- as.comparison(data.frame(start1 = starts1, end1 = ends1,
#'                                         start2 = starts2, end2 = ends2))
#' print_comparison(comparison1)
#' 
#' ## Printing out a comparison like this can occasionally throw an error:
#' \dontrun{
#' comparison1
#' print(comparison1)}
#' ## This can happen when print.comparison is loaded from the testthat package
#' ## To avoid this, use print_comparison as shown above
#' 
print_comparison <- function(comparison) {
  class(comparison) <- c("data.table", "data.frame")
  print(comparison)
}

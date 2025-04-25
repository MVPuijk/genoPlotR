################################
# dna_seg class and methods
################################

#' dna_seg class and methods
#' 
#' `dna_seg` objects are collections of genes or other elements along a genome,
#' to be represented on a map. These functions are class functions to create,
#' convert, and test `dna_seg` objects.
#' 
#' @details
#' A `dna_seg` object is a `data.table` that must start with the columns `name`,
#' `start`, `end`, and `strand`. The arguments listed above add more
#' columns, but only if those were not already present in the data provided by
#' `x`. Any number of additional columns can be added to `dna_seg` objects, to
#' act as metadata or as ways to group or identify features.
#' 
#' `dna_seg` and `as.dna_seg` can both be used to create `dna_seg` objects from
#' list, `data.frame`, or `data.table` objects. Alternatively, 
#' [read_dna_seg_from_file] can be used to create them from a file.
#' 
#' These objects inherit from `data.table`, and
#' can therefore be manipulated using `data.table` syntax.
#' 
#' `is.dna_seg` returns `TRUE` if the object tested is a `dna_seg` object.
#' 
#' @name dna_seg
#' @returns `as.dna_seg` and `dna_seg` both return a `dna_seg` object.
#' @returns `is.dna_seg` returns a logical.
#' @export
#' 
#' @param x An object to be coerced. Can be a `data.frame`, `data.table`, 
#' or `list` object. See details for necessary columns.
#' @param df Same as `x`.
#' @param ... Arguments to pass to [dna_seg], see arguments below.
#' 
#' @author Lionel Guy, Mike Puijk
#' 
#' @seealso [read_dna_seg], [trim.dna_seg], [reverse.dna_seg], [gene_types],
#' [update_dna_segs]
#' 
#' @examples
#' ## Generate data
#' names1 <- c("feat1", "feat2", "feat3")
#' starts1 <- c(2, 1000, 1050)
#' ends1 <- c(600, 800, 1345)
#' strands1 <- c("-", -1, 1)
#' cols1 <- c("blue", "grey", "red")
#' 
#' ## Create data.frame
#' df1 <- data.frame(name = names1, start = starts1, end = ends1,
#'                   strand = strands1, col = cols1)
#' 
#' ## Create dna_seg
#' dna_seg1 <- dna_seg(df1)
#' dna_seg1
#' as.dna_seg(df1)
#' 
#' ## Test dna_seg
#' is.dna_seg(dna_seg1)
#' 
dna_seg <- function(x, ...) {
  # support for list ?
  if (is.data.table(x)) {
    return(as.dna_seg(x, ...))
  } else if (is.data.frame(x)) {
    return(as.dna_seg(as.data.table(x)))
  } else if (is.list(x)) {
    return(as.dna_seg(as.data.table(x)))
  } else {
    stop("Cannot coerce class ", class(x), " to dna_seg")
  }
}

#' @name dna_seg
#' @export
#' 
#' @param col A character vector of colors, of either length one or the same
#' length as `x`. Determines the outline of features, or their overall
#' color for shapes (see [gene_types]) that have no fill color (e.g. text,
#' points).
#' @param fill Same as `col`, but determines the fill color of features.
#' @param lty A vector of either length 1 or the same length as `x`. Determines
#' the line type (see [par]) of the features.
#' @param lwd A numeric vector of either length 1 or the same length as `x`. 
#' Determines the line width (see [par]) of the features.
#' @param pch A vector of either length 1 or the same length as `x`. Determines
#' the shape (see [points]) of the features when they are represented by points.
#' @param cex A numeric vector of either length 1 or same the length as `x`.
#' Determines the size multiplier of features when they are represented by 
#' text or points.
#' @param region_plot A character vector of either length 1 or the same length
#' as `x`. Determines if features are going to be plotted when going for
#' regional plotting. See [plot_gene_map] for more details.
#' @param gene_type A character vector of either length 1 or same the length as 
#' `x`. Determines the gene type (i.e. shape) of the features. See [gene_types].
#' @param ordered Logical. If `TRUE`, orders the `dna_seg` on the `start` and
#' `end` columns.
#' 
# convert to dna_seg format. 
as.dna_seg <- function(
  df,
  col = "grey20",
  fill = "grey80",
  lty = 1,
  lwd = 1,
  pch = 8,
  cex = 1,
  gene_type = "arrows",
  region_plot = "NA",
  ordered = TRUE
) {
  
  # check for class dna_seg, list, df, dt
  if (is.dna_seg(dna_seg)) return(df)
  if (is.data.frame(df) && !is.data.table(df)) df <- as.data.table(df)
  if (is.list(df) && !is.data.table(df)) df <- as.data.table(df)
  if (is.data.table(df)) {
    # check that it has rows
    if (nrow(df) < 1) {
      stop("Number of rows is 0, check data input")
    }
    # check for columns
    names <- c("name", "start", "end", "strand")
    if (!identical(names(df)[1:4], names)) {
      stop("Column names must start with 'name', 'start', 'end', 'strand'")
    }
    if (!all(sapply(df, function(x) all(!is.null(x))))) {
      stop("NULL values not allowed")
    }
    if (!is.numeric(df$start) | !is.numeric(df$end)) {
      stop("'start' and 'end' columns must be numeric")
    }
    if (is.factor(df$name)) df$name <- as.character(df$name)
    if (is.factor(df$strand)) df$strand <- as.character(df$strand)
    if (is.factor(df$col)) df$col <- as.character(df$col)
    if (is.factor(df$fill)) df$fill <- as.character(df$fill)
    if (is.factor(df$gene_type)) df$gene_type <- as.character(df$gene_type)
    if (is.factor(df$region_plot)) df$region_plot <- as.character(df$region_plot)
    # care for strand
    if (is.character(df$strand)) {
      df$strand[df$strand == "+"] <- 1
      df$strand[df$strand == "-"] <- -1
      df$strand <- as.numeric(df$strand)
    }
    if (!all(df$strand %in% c(-1, 1))) {
      stop("'strand' column must be one of 1, -1, -, +")
    }
    # gene_type: not given in gene
    if (is.null(df$gene_type)) {
      if (gene_type == "auto") gene_type <- auto_gene_type(nrow(df))
      df$gene_type <- gene_type
    }
    # region_plot check
    if (is.null(df$region_plot)) {
      df$region_plot <- region_plot
    } else {
      df$region_plot <- as.character(df$region_plot)
      df$region_plot[df$region_plot %ilike%
                       "^true$|^plot$|^y$|^t$|^yes$"] <- "TRUE"
      df$region_plot[!df$region_plot %ilike%
                       "^true$|^plot$|^y$|^t$|^yes$|^start$|^end$"] <- "NA"
      df$gene_type[df$region_plot %ilike% "^start$|^end$"] <- "boundaries"
    }
    # col
    if (is.null(df$col)) df$col <- col
    if (is.null(df$fill)) df$fill <- fill
    # lwd & lty
    if (is.null(df$lty)) df$lty <- lty
    if (is.null(df$lwd)) {
      df$lwd <- lwd
      df$lwd[df$gene_type == "boundaries"] <- lwd * 3
    } 
    if (is.null(df$pch)) df$pch <- pch
    if (is.null(df$cex)) df$cex <- cex
    # check for correct argument types
    if (!is.character(df$name)) stop("'name' column must be character")
    if (!(is.numeric(df$start) && is.numeric(df$end))) {
      stop("Non-numeric 'start' or 'end'")
    }
    if (!all(is.numeric(c(df$lwd, df$lty, df$pch, df$cex)))) {
      stop("'lwd', 'lty', 'pch', and 'cex' must be numeric")
    }
    if (!is.character(df$gene_type)) {
      stop("'gene_type' must be a character vector, made of: ",
           paste(gene_types(), collapse = ", "), ", or a function name")
    }
    
    # helper function that checks if all colors are recognizable
    colorcheck <- function(df_colors, column) {
      
      col_is_rec <- sapply(df_colors, function(x) tryCatch({
        is.matrix(col2rgb(x))},
        error = function(e) FALSE
      ))
      
      if (!all(col_is_rec)) {
        # try to parse unrecognizable colors
        not_rec <- col_is_rec[col_is_rec == FALSE]
        not_rec_names <- names(not_rec)
        comma_grep <- grep(",", not_rec_names, fixed = TRUE)
        both_grep <- grep("[, ]", not_rec_names)
        if (length(both_grep) > 0) {
          for (i in both_grep) {
            if (any(i == comma_grep)) {
              splitted <- unlist(strsplit(not_rec_names[i], ","))
            } else {
              splitted <- unlist(strsplit(not_rec_names[i], " "))
            }
            num_splitted <- suppressWarnings(as.numeric(splitted))
            if (!length(num_splitted) %in% c(3,4) || any(is.na(num_splitted))) {
              stop("All values in '", column, "' must be valid colors")
            } else {
              num_list <- as.list(num_splitted)
              num_list$maxColorValue <- 255
              tryCatch({
                rgb_colors <- do.call(rgb, num_list)
              }, error = function(e) {
                stop("All values in '", column, "' must be valid colors")
              })
              not_rec_names[i] <- rgb_colors
            }
          }
        } else {
          stop("All values in '", column, "' must be valid colors")
        }
        df_colors[col_is_rec == FALSE] <- not_rec_names
      }
      
      col_is_rec <- sapply(df_colors, function(x) tryCatch({
        is.matrix(col2rgb(x))},
        error = function(e) FALSE
      ))
      if (!all(col_is_rec)) {
        stop("All values in '", column, "' must be valid colors")
      }
      
      df_colors
    }
    
    df$col <- colorcheck(df$col, "col")
    df$fill <- colorcheck(df$fill, "fill")
    
  } else {
    stop("Unable to handle this format")
  }
  if (ordered) df <- df[order(start, end)]
  
  class(df) <- c("dna_seg", class(df))
  return(df)
}

#' @name dna_seg
#' @export
#' 
#' @param dna_seg An object to test.
#' 
is.dna_seg <- function(dna_seg) {
  inherits(dna_seg, "dna_seg")
}

#' @rdname range
#' @export
#' 
range.dna_seg <- function(x, ...) {
  dna_seg <- x
  range(dna_seg$start, dna_seg$end, na.rm = FALSE)
}


# # merges two dna_segs. So far returns the minimal common set of columns, maybe
# # changed in future
# c.dna_seg <- function(...) {
#   # a little helper function to grab colnames
#   #fold <- function(f, x, L) (for (e in L) x <- f(x, e))
#   fold <- function(x, fun) {
#     if (length(x) == 1) return(fun(x))
#     accumulator <- fun(x[[1]], x[[2]])
#     if (length(x) == 2) return(accumulator)
#     for (i in 3:length(x)) {
#       accumulator <- fun(accumulator, x[[i]])
#     }
#     accumulator
#   }
#   # parse args
#   x <- list(...)
#   n <- length(x)
#   if (!all(sapply(x, is.dna_seg)))
#     stop("All elements must be of class dna_seg")
#   if (n == 1) return(x)
#   cols <- lapply(x, names)
#   com_cols <- fold(cols, intersect)
#   dna_seg <- x[[1]][, ..com_cols]
#   for (i in 2:n) {
#     dna_seg <- rbindlist(list(dna_seg, x[[i]][, ..com_cols]))
#   }
#   # rbind from data.table makes its arguments lose their original class
#   class(dna_seg) <- c("dna_seg", class(dna_seg))
#   dna_seg
# }

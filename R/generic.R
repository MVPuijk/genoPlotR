################################################################################
# S3 method: generic functions
################################################################################

#' Trimming data frames using a numeric interval
#' 
#' Trims data frames with 2 or more `numeric` columns using a `numeric` interval. 
#' Returns rows with `numeric` values that fall within the interval.
#' 
#' @details
#' By default it uses the first 2 `numeric` columns in `x`.
#' If `x` is a `dna_seg` object, it uses the `start` and `end` columns.
#' If `x` is a `comparison` object, it filters using the `start1` and `end1`
#' columns with the `xlim1` argument, and `start2` and `end2` with the `xlim2`
#' argument.
#' If `x` is `annotation` object, it uses the `x1` and `x2` columns.
#' If `x` is a `seg_plot` object, the function uses the `xargs` variable from
#' `x` to define which vectors determine the x position (they should be the
#' same length). Then, all the variables (including those inside the `gp` 
#' variable) that are the same length as the x vectors are trimmed, so that only
#' the rows for which the x values are inside the `xlim` argument are kept.
#' 
#' @rdname trim
#' @returns An object with the same type as `x`, with the rows (or subset) 
#' corresponding to the given interval
#' @export
#' 
#' @param x An object to trim, generally a `data.frame`, `data.table`, or 
#' `seg_plot` object
#' @param xlim A `numeric` of length 2. In a general case, the rows whose values
#' are included in this interval are returned
#' @param xlim1 A `numeric` of length 2, used only for `comparison` objects, 
#' the interval to filter the first side.
#' @param xlim2 A `numeric` of length 2, used only for `comparison` objects, 
#' the interval to filter the second side.
#' @param ... Unused.
#' 
#' @author Lionel Guy
#' 
#' @seealso [dna_seg], [comparison], [annotation], [seg_plot]
#' 
#' @examples
#' ## Load
#' data(barto)
#' xlim_ref <- c(10000, 45000)
#' ## Seg 2 (ref)
#' barto$dna_segs[[2]] <- trim(barto$dna_segs[[2]], xlim = xlim_ref)
#' ## Seg 1
#' barto$comparisons[[1]] <- trim(barto$comparisons[[1]], xlim2 = xlim_ref)
#' xlim1 <- range(barto$comparisons[[1]], overall = FALSE)$xlim1
#' barto$dna_segs[[1]] <- trim(barto$dna_segs[[1]], xlim = xlim1)
#' ## Seg 3
#' barto$comparisons[[2]] <- trim(barto$comparisons[[2]], xlim1 = xlim_ref)
#' xlim3 <- range(barto$comparisons[[2]], overall = FALSE)$xlim2
#' barto$dna_segs[[3]] <- trim(barto$dna_segs[[3]], xlim = xlim3)
#' ## Seg 4
#' barto$comparisons[[3]] <- trim(barto$comparisons[[3]], xlim1 = xlim3)
#' xlim4 <- range(barto$comparisons[[3]], overall = FALSE)$xlim2
#' barto$dna_segs[[4]] <- trim(barto$dna_segs[[4]], xlim = xlim4)
#' ## Plot
#' plot_gene_map(barto$dna_segs, barto$comparisons)
#' 
#' ## With seg_plot
#' x <- 1:20
#' y <- rnorm(20)
#' sp <- seg_plot(func = pointsGrob,
#'                args = list(x = x, y = y, gp = gpar(col = 1:20, cex = 1:3)))
#' ## Trim 
#' sp_trim <- trim(sp, c(3, 10))
#' str(sp_trim)
#' range(sp_trim$arg$x)
#' 
trim <- function(x, ...) {
  UseMethod("trim")
}

#' @rdname trim
#' @export
trim.default <- function(x, xlim = NULL, ...) {
  if (!is.numeric(xlim)) stop("'xlim' must be numeric")
  if (length(xlim) != 2) stop("'xlim' must be length 2")
  num_col <- which(sapply(x, is.numeric))
  if (length(num_col) < 2) stop ("'x' must have at least 2 numeric columns")
  x <- x[x[,num_col[1]] >= xlim[1] & x[,num_col[2]] <= xlim[2]]
  x
}

#' @rdname trim
#' @export
# trim dna_seg given x limits
trim.dna_seg <- function(x, xlim = NULL, ...) {
  dna_seg <- x
  if (!is.null(xlim)) {
    if (!is.numeric(xlim)) stop("'xlim' must be numeric")
    if (length(xlim) != 2) stop("'xlim' must be length 2")
    dna_seg <- dna_seg[dna_seg$start >= xlim[1] & dna_seg$end <= xlim[2],]
  }
  dna_seg
}

#' @rdname trim
#' @export
# trim comparison given x limits
trim.comparison <- function(x,
                            xlim1 = c(-Inf, Inf),
                            xlim2 = c(-Inf, Inf),
                            ...) {
  comparison <- x
  if (!is.null(xlim1) && !is.null(xlim2)) {
    if (!is.numeric(xlim1) || !is.numeric(xlim2)) {
      stop("'xlims' must be numeric")
    }
    if (length(xlim1) != 2 || length(xlim2) != 2) {
      stop("'xlims' must be length 2")
    }
    # testing to include overlapping comps
    # direction 1
    comparison$start1[comparison$start1 < xlim1[1] &
                        comparison$end1 > xlim1[1]
                      ] <- xlim1[1]
    comparison$end1[comparison$start1 < xlim1[2] &
                      comparison$end1 > xlim1[2]
                    ] <- xlim1[2]
    comparison$start2[comparison$start2 < xlim2[1] &
                        comparison$end2 > xlim2[1]
                      ] <- xlim2[1]
    comparison$end2[comparison$start2 < xlim2[2] &
                      comparison$end2 > xlim2[2]
                    ] <- xlim2[2]
    # direction -1
    comparison$start1[comparison$start1 > xlim1[2] &
                        comparison$end1 < xlim1[2]
                      ] <- xlim1[2]
    comparison$end1[comparison$start1 > xlim1[1] &
                      comparison$end1 < xlim1[1]
                    ] <- xlim1[1]
    comparison$start2[comparison$start2 > xlim2[2] &
                        comparison$end2 < xlim2[2]
                      ] <- xlim2[2]
    comparison$end2[comparison$start2 > xlim2[1] &
                      comparison$end2 < xlim2[1]
                    ] <- xlim2[1]
    comparison <- comparison[comparison$start1 >= xlim1[1] &
                               comparison$end1 <= xlim1[2] &
                               comparison$start2 >= xlim2[1] &
                               comparison$end2 <= xlim2[2],
                             ]
  }
  comparison
}

#' @rdname trim
#' @export
#' 
trim.annotation <- function(x, xlim = NULL, ...) {
  annotation <- x
  xlim <- as.numeric(xlim)
  if (!is.null(xlim)) {
    if (!is.numeric(xlim)) stop("'xlim' must be numeric")
    if (length(xlim) != 2) stop("'xlim' must be length 2")
    # to be accepted, x1 > xlim1 and, if x2 = NA, xlim1 also < xlim1 or,
    # x2 < xlim2
    annotation <- 
      annotation[annotation$x1 >= xlim[1] & 
                   ((is.na(annotation$x2) & annotation$x1 <= xlim[2]) |
                      (!is.na(annotation$x2) & annotation$x2 <= xlim[2])) ,
                 ]
  }
  annotation
}

#' @rdname trim
#' @export
#' 
# Tries to trim correctly the seg_plot objects...
# Made difficult because of the complexity of the objects
trim.seg_plot <- function(x, xlim = NULL, ...) {
  # Check that we are in the right class
  if (!is.seg_plot(x)) stop("'x' must be a seg_plot object")
  xlim <- as.numeric(xlim)
  # If xlim is null, return whole object. Else, trim
  if (!is.null(xlim)) {
    if (!is.numeric(xlim)) stop("'xlim' must be numeric")
    if (length(xlim) != 2) stop("'xlim' must be length 2")
    mainargs <- names(x$args)
    gpargs <- NULL
    # Define arguments to filter
    if ("gp" %in% mainargs && inherits(x$args$gp, "gpar")) {
      mainargs <- mainargs[mainargs != "gp"]
      gpargs <- names(x$args$gp)
    }
    ##std_xargs <- c("x", "x0", "x1", "x2", "v")
    xargs <- mainargs[mainargs %in% x$xargs]
    if (length(xargs) > 0) {
      # Check that all xargs have the same length
      nrows <- length(x$args[[xargs[1]]])
      if (length(xargs) > 1) {
        sapply(xargs[2:length(xargs)], function(sp) {
          if (length(x$args[[sp]]) != nrows) {
            stop(paste("Argument", sp, "has a different length as argument",
                       xargs[1]))
          }
        })
      }
      # Get the correct indexes
      idx <- rep(TRUE, nrows)
      for (i in 1:length(xargs)) {
        idx <- idx & x$args[[xargs[i]]] >=
          xlim[1] & x$args[[xargs[i]]] <= xlim[2]
      }
      # Filter everything...
      for (arg in mainargs) {
        if (length(x$args[[arg]]) == nrows) {
          x$args[[arg]] <- x$args[[arg]][idx]
        }
      }
      if (!is.null(gpargs)) {
        for (arg in gpargs) {
          if (length(x$args$gp[[arg]]) == nrows) {
            x$args$gp[[arg]] <- x$args$gp[[arg]][idx]
          }
        }
      }
    }
  }
  x
}

#' Reverse objects
#' 
#' Reverses objects, mainly meant for `dna_seg` and `comparison` objects.
#' 
#' @name reverse
#' @returns An object with the same type as `x`.
#' @export
#' 
#' @param x An object to reverse.
#' @param side Used only when `x` is a `comparison` object, the side that 
#' should be reversed. If `side = 1`, the first side will be reversed. If 
#' `side = 2`, the second side will be reversed. If `side < 1`, no side is 
#' reversed. If `side > 2`, both sides are reversed.
#' @param ... Unused.
#' 
#' @author Lionel Guy
#' 
#' @seealso [dna_seg], [comparison]
#' 
#' @examples
#' ## load data
#' data(three_genes)
#' dna_segs <- three_genes$dna_segs
#' comparisons <- three_genes$comparisons
#' 
#' ## on dna_seg
#' dna_segs[[1]]
#' reverse(dna_segs[[1]])
#' ## on comparison
#' reverse(comparisons[[2]], side = 1)
#' reverse(comparisons[[2]], side = 3)
#' 
#' ## With mauve backbone
#' data(mauve_bbone)
#' ## Plot
#' plot_gene_map(dna_segs = mauve_bbone$dna_segs,
#'               comparisons = mauve_bbone$comparisons,
#'               alpha_comparisons = 0.4)
#' 
#' ## Reverse B_bacilliformis, and the corresponding comparison (first "side")
#' mauve_bbone$dna_segs[[1]] <- reverse(mauve_bbone$dna_segs[[1]])
#' mauve_bbone$comparisons[[1]] <- reverse(mauve_bbone$comparisons[[1]], 1)
#' plot_gene_map(dna_segs = mauve_bbone$dna_segs,
#'               comparisons = mauve_bbone$comparisons,
#'               alpha_comparisons = 0.4)
#' 
reverse <- function(x, ...) {
  UseMethod("reverse")
}

#' @rdname reverse
#' @export
#' 
# default method
reverse.default <- function(x, ...) {
  num_col <- which(sapply(x, is.numeric))
  if (length(num_col) < 2) stop ("'x' must have at least 2 numeric columns")
  tmp <- -x[,num_col[2]]
  x[,num_col[2]] <- -x[,num_col[1]]
  x[,num_col[1]] <- tmp
  x
}

#' @rdname reverse
#' @export
#' 
# reverse dna_seg
reverse.dna_seg <- function(x, ...) {
  dna_seg <- x
  start <- -dna_seg$end
  dna_seg$end <- -dna_seg$start
  dna_seg$start <- start
  dna_seg$strand <- -dna_seg$strand
  dna_seg
}

#' @rdname reverse
#' @export
#' 
# reverses a comparison. side <1 for no sides, 1 for first,
# 2 for second, >2 for both
reverse.comparison <- function(x, side = 0, ...) {
  comparison <- x
  if (side > 1) {
    comparison$start2 <- -comparison$start2
    comparison$end2 <- -comparison$end2
  }
  if (side == 1 || side > 2) {
    comparison$start1 <- -comparison$start1
    comparison$end1 <- -comparison$end1    
  }
  # if only one side changed, change direction
  if (side == 1 || side == 2) {
    comparison$direction <- -comparison$direction
  }
  comparison
}

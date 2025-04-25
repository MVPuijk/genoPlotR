################################
# Annotation class and methods
################################

#' annotation class and methods
#' 
#' `annotation` objects are used to describe `dna_seg` objects. It consists
#' of labels that are attached to a single position or a range.
#' 
#' @details
#' An `annotation` object is a `data.table` with the columns `x1`, `x2`, 
#' `text`, `col`, and `rot`. These provide the first (or only) position, an
#' optional second position, the text, the color, and the rotation of the 
#' annotation, respectively. 
#' 
#' When plotted using [plot_gene_map], it will add an annotation row on top
#' of the corresponding DNA segment. For any given row of an `annotation` 
#' object, if `x2` is `NA`, the label will be attached to the position given by
#' `x1`. If `x2` is a number instead, a range will be plotted between these 
#' two values and the label itself will be attached in the middle of this range.
#' 
#' `is.annotation` returns `TRUE` if the object tested is an `annotation`
#' object.
#' 
#' @name annotation
#' @returns `as.annotation` and `annotation` both return an
#' `annotation` object.
#' @returns `is.annotation` returns a logical.
#' @export
#' 
#' @param x1 A numeric vector giving the first or only position of the label.
#' @param x2 A numeric vector of the same length as `x1` or `NA`, providing
#' an optional secondary position for the label. 
#' @param text A character vector of the same length as `x1`, providing the text
#' of the labels.
#' @param rot A numeric vector of the same length as `x1`, providing the
#' rotation of the labels in degrees.
#' @param col A character vector of the same length as `x1`, providing the color
#' of the labels.
#' 
#' @author Lionel Guy
#' 
#' @seealso [plot_gene_map], [middle]
#' 
#' @examples
#' ## Loading data
#' data(three_genes)
#' dna_segs <- three_genes$dna_segs
#' comparisons <- three_genes$comparisons
#' 
#' ## Calculating middle positions
#' mid_pos <- middle(dna_segs[[1]])
#' 
#' ## Create first annotation
#' annot1 <- annotation(x1 = mid_pos, text = dna_segs[[1]]$name)
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               annotations = annot1)
#' 
#' ## Exploring options
#' annot2 <- annotation(x1 = c(mid_pos[1], dna_segs[[1]]$end[2]),
#'                      x2 = c(NA, dna_segs[[1]]$end[3]),
#'                      text = c(dna_segs[[1]]$name[1], "region1"),
#'                      rot = c(30, 0),
#'                      col = c("grey", "black"))
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               annotations = annot2, annotation_height = 1.3)
#'               
#' ## Using a bigger dataset from a 4-genome comparison
#' data(barto)
#' 
#' ## Showing several subsegments
#' xlims2 <- list(c(1445000, 1415000, 1380000, 1412000),
#'                c(  10000,   45000,   50000,   83000, 90000, 120000),
#'                c(  15000,   36000,   90000,  120000, 74000,  98000),
#'                c(   5000,    82000))
#' ## Adding annotations for all genomes, allow segments to be placed out
#' ## of the longest segment
#' annots <- lapply(barto$dna_segs, function(x) {
#'   mid <- middle(x)
#'   annot <- annotation(x1 = mid, text = x$name, rot = 30)
#'   # removing gene names starting with "B" and keeping 1 in 4
#'   idx <- grep("^[^B]", annot$text, perl = TRUE)
#'   annot[idx[idx %% 4 == 0], ]
#' })
#' plot_gene_map(dna_segs = barto$dna_segs,
#'               comparisons = barto$comparisons,
#'               annotations = annots,
#'               xlims = xlims2,
#'               limit_to_longest_dna_seg = FALSE,
#'               dna_seg_scale = TRUE)
#' 
#' ## Annotations on all the segments
#' annots <- lapply(dna_segs, function(x) {
#'   mid <- middle(x)
#'   annot <- annotation(x1 = mid, text = x$name, rot = 30)
#' })
#' plot_gene_map(dna_segs = dna_segs,
#'               comparisons = comparisons,
#'               annotations = annots,
#'               annotation_height = 1.8,
#'               annotation_cex = 1)
#' 
annotation <- function(x1, x2 = NA, text, rot = 0, col = "black") {
  if (missing(x1) || missing(text)) stop("'x1' and 'text' must be provided")
  if (!is.numeric(x1)) stop("'x1' must be numeric")
  if (!(all(is.na(x2)) || is.numeric(x2))) stop("'x2' must be numeric or NA")
  if (!is.character(text)) stop("'text' must be character")
  as.annotation(
    data.table(x1 = x1, x2 = x2, text = text, stringsAsFactors = FALSE),
    rot = rot,
    col = col
  )
}

#' @name annotation
#' @export
#' 
#' @param df A `data.frame` to convert to an `annotation` object. Must have at
#' least the columns `x1` and `text`.
#' 
as.annotation <- function(df, x2 = NA, rot = 0, col = "black") {
  if (is.annotation(df)) return(df)
  if (!all(c("x1", "text") %in% names(df))) {
    stop("'df' must have at least an 'x1' and 'text' columns")
  }
  # attributes x2, col and arg to all rows if not defined
  if (is.null(df$x2)) df$x2 <- x2
  if (is.null(df$color)) df$color <- col
  if (is.null(df$rot)) df$rot <- rot
  df <- as.data.table(df)
  class(df) <- c("annotation", class(df))
  df
}

#' @name annotation
#' @export
#' 
#' @param annotation An object to test.
#' 
is.annotation <- function(annotation) {
  inherits(annotation, "annotation")  
}

#' Range calculation
#' 
#' Calculate the range of `dna_seg`, `comparison` or `annotation` objects.
#' 
#' @name range
#' @rdname range
#' @returns A numeric of length 2, unless `x` is a `comparison` object
#' and `overall` is `FALSE`, in which case it will be a `data.frame` object with
#' two rows and two named columns, `xlim1` and `xlim2`
#' @export
#' 
#' @param x An object to calculate the range from.
#' @param overall Used only when `x` is a `comparison` object, a logical. If
#' `FALSE`, will calculate the range of each side separately instead of
#' calculating the overall range of `x`.
#' @param ... Unused.
#' 
#' @author Lionel Guy
#' 
#' @seealso [dna_seg], [comparison], [annotation], [trim]
#' 
#' @examples
#' ## Create dna_segs and comparison
#' df1 <- data.frame(name = c("feat1", "feat2", "feat3"),
#'                   start = c(2, 1000, 1050),
#'                   end = c(600, 800, 1345),
#'                   strand = c(-1, -1, 1))
#' df2 <- data.frame(name = c("feat1", "feat2", "feat3"),
#'                   start = c(50, 800, 1200),
#'                   end = c(900, 1100, 1322),
#'                   strand = c(-1, 1, 1))
#' dna_seg1 <- dna_seg(df1)
#' dna_seg2 <- dna_seg(df2)
#' df3 <- data.frame(start1 = dna_seg1$start,
#'                   end1 = dna_seg1$end,
#'                   start2 = dna_seg2$start,
#'                   end2 = dna_seg2$end)

#' comparison1 <- comparison(df3)
#' 
#' ## Range calculation
#' range(dna_seg1)
#' range(comparison1)
#' range(comparison1, overall = FALSE)
#' 
# calculates dna_seg range
range.annotation <- function(x, ...) {
  annotation <- x
  range(annotation$x1, annotation$x2, na.rm = TRUE)
}



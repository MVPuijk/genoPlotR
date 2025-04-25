################################################################################
# Color schemes
################################################################################

#' Generate a gradient color scheme
#' 
#' A color scale is generated based on a numeric vector, with the darker colors 
#' corresponding to the higher values, or the lower values if 
#' `decreasing = TRUE`.
#' 
#' `apply_color_scheme` is synonymous with `gradient_color_scheme`, and is
#' included for backwards compatibility.
#' 
#' @details
#' Made for use alongside `comparison` objects, where it is especially 
#' useful to visualize representations of sequence similarities, like BLAST 
#' percentage identity values or e-values. 
#' 
#' For the red-blue color scheme, red colors correspond to comparisons that go
#' in the same direction (1), where the blue colors correspond to comparisons
#' that go in the opposite direction (-1). 
#' 
#' @returns A character vector of colors, with the same length as `x`.
#' @export
#' @aliases apply_color_scheme
#' 
#' @param x A numeric vector to generate colors for.
#' @param direction A numeric vector composed of `-1` and `1` values, denoting
#' the direction of the comparison. Ignored unless the color scheme is a
#' red-blue color scheme.
#' @param color_scheme A character string, one of: `"red_blue"`, `"blue_red"`,
#' `"gray"`, or `"grey"`.
#' @param decreasing Logical. Are the values of `x` representing a relationship
#' that gets stronger as the number goes down (e.g. e-values, gaps, mismatches)? 
#' @param rng Numeric of length 2. Gives the lower and upper limits to apply
#' the gradient to.
#' @param alpha A single numeric value between 0 and 1, or `FALSE`. Determines
#' the transparency applied to the color scheme, 0 being fully transparent, and
#' 1 being fully opaque.
#' @param transparency Deprecated, included for backwards compatibility.
#' When provided, replaces `alpha`.
#' 
#' @author Lionel Guy
#' 
#' @seealso [sequential_color_scheme], [uniform_color_scheme], [plot_gene_map], 
#' [comparison]
#' 
#' @examples
#' ## Load data
#' data(three_genes)
#' dna_segs <- three_genes$dna_segs
#' comparisons <- three_genes$comparisons
#' 
#' ## Color schemes
#' ## Greys
#' comparisons[[1]]$values <- c(70, 80, 90)
#' comparisons[[1]]$col <- gradient_color_scheme(
#'   comparisons[[1]]$values,
#'   color_scheme = "grey")
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons)
#' ## Red-blue
#' comparisons[[1]]$col <- gradient_color_scheme(
#'   comparisons[[1]]$values,
#'   direction = comparisons[[1]]$direction,
#'   color_scheme = "red_blue")
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons)
#' ## Decreasing
#' comparisons[[1]]$col <- gradient_color_scheme(
#'   comparisons[[1]]$values,
#'   direction = comparisons[[1]]$direction,
#'   color_scheme = "red_blue",
#'   decreasing = TRUE)
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons)
#' ## Range
#' comparisons[[1]]$col <- gradient_color_scheme(
#'   comparisons[[1]]$values,
#'   direction = comparisons[[1]]$direction,
#'   color_scheme = "red_blue",
#'   rng = c(30,100))
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons)
#' ## Transparency
#' x1 <- seq(100, 600, by = 50)
#' x2 <- seq(1100, 700, by = -50)
#' comparisons[[2]] <- as.comparison(data.frame(
#'   start1 = c(x1, x2),
#'   end1 = c(x1 + 250, x2 + 300),
#'   start2 = c(x1 + 150, x2 - 300) + 2000,
#'   end2 = c(x1 + 250, x2 - 500) + 2000))
#' comparisons[[1]]$col <- gradient_color_scheme(
#'   comparisons[[1]]$values,
#'   color_scheme = "grey",
#'   alpha = 0.8)
#' comparisons[[2]]$col <- gradient_color_scheme(
#'   1:nrow(comparisons[[2]]),
#'   direction = comparisons[[2]]$direction,
#'   color_scheme = "blue_red")
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons)
#' 
#' comparisons[[1]]$col <- gradient_color_scheme(
#'   comparisons[[1]]$values,
#'   color_scheme = "grey",
#'   alpha = 1)
#' comparisons[[2]]$col <- gradient_color_scheme(
#'   1:nrow(comparisons[[2]]),
#'   direction = comparisons[[2]]$direction,
#'   color_scheme = "blue_red",
#'   alpha = 0.2)
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons)
#' 
# apply color scheme to a numeric vector
# decreasing: relationship gets stronger with decreasing values?
# TRUE for, typically, e-values, gaps, mism, ... 190
# FALSE for bit scores, per_id, aln_length
gradient_color_scheme <- function(
  x,
  direction = NULL,
  color_scheme = "grey",
  decreasing = FALSE,
  rng = NULL,
  alpha = 0.5,
  transparency # backwards compatibility
) {
                                  
  # check arguments
  # if length is 0, return length 0
  if (!is.null(direction) && length(direction) == 0) return (character(0))
  # if x is null and direction is not, get x to 1s (mainly for blue/red)
  if (is.null(x) && !is.null(direction)) {
    x <- rep(1, length(direction))
  }
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (is.null(rng)) rng <- range(x)
  if (!missing(transparency)) alpha <- transparency
  if (is.null(alpha)) alpha <- FALSE
  if (!(is.logical(alpha) && !(alpha)) && !is.numeric(alpha)) {
    stop ("'alpha' must be FALSE or numeric")
  }
  col <- rep(grey(0.5), length(x))
  # red blue
  if (any(color_scheme %in% c("red_blue", "blue_red"))) {
    if (is.null(direction)) direction <- rep(1, length(x))
    blues <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6",
               "#4292C6", "#2171B5", "#08519C", "#08306B")
    reds  <- c("#FFF5F0", "#FEE0D2", "#FCBBA1", "#FC9272", "#FB6A4A",
               "#EF3B2C", "#CB181D", "#A50F15", "#67000D")
    if (diff(rng) == 0) { # case: only one value:
      level <- rep(9, length(x))
    } else { # case: several values
      level <- round(((x-rng[1])/diff(rng))*8+1)
    }
    if (decreasing) level <- -level + 10
    col[direction == 1] <- reds[level[direction == 1]] 
    col[direction == -1] <- blues[level[direction == -1]]
  } else if (any(color_scheme %in% c("grey", "gray", "grays", "greys"))) {
    # grey: between 0.25 and 0.75
    if (diff(rng) == 0) { # case: only one value:
      col <- rep(grey(0.5), length(x))
    } else { # case: several values
      level <- 0.75-((x-rng[1])/diff(rng))*0.5
      if (decreasing) level <- -level+1
      col <- grey(level)
    }
  } else {
    stop('\'color_scheme\' must be either "red_blue" or "grey"')
  }
  if (alpha && alpha != 1) {
    # Convert ratio into hexadecimal
    if (alpha > 1 || alpha < 0) stop("'alpha' must be between 0 and 1")
    tpc <- floor(alpha * 256)
    tpc <- sprintf("%X", tpc)
    if (nchar(tpc) == 1) tpc <- paste("0", tpc, sep = "")
    col <- paste(col, tpc, sep = "")
  }
  
  col
}

#' @name gradient_color_scheme
#' 
# Included for backwards compatibility
apply_color_scheme <- function(
  x,
  direction = NULL,
  color_scheme = "grey",
  decreasing = FALSE,
  rng = NULL,
  alpha = 0.5,
  transparency # backwards compatibility
) {
                               
  if (!missing(transparency)) alpha <- transparency
  col <- gradient_color_scheme(x = x,
                               direction = direction,
                               color_scheme = color_scheme,
                               decreasing = decreasing,
                               rng = rng,
                               alpha = alpha
                               )
  
  col
}
  
#' Apply a sequential color scheme to dna_segs and comparisons
#' 
#' Generates and applies a sequential color scheme to a list of `dna_seg` and
#' `comparison` objects. It does this by taking the colors that are already
#' there and transferring those over to any features connected to it through
#' the comparisons. For example, if a feature from a single `dna_seg` has red
#' as its `fill` attribute, the comparisons that can be linked to this feature
#' will become red as well. This is then followed up by updating any `dna_seg`
#' features linked to those comparisons, and so on.
#' 
#' @details
#' The existing colors from the input `dna_seg` and 
#' `comparison` objects are transferred over to the next object in the
#' plotting order, with the exception of their default colors, provided by
#' `default_color`. As `comparison` objects only have a single color
#' attribute `col`, those will be updated using the column provided by
#' `color_var` from the `dna_segs`, while the `dna_segs` themselves will be
#' updated using the `col` column from the comparisons regardless of 
#' `color_var`.
#' 
#' The objects are linked together through shared
#' values. The columns for these shared values are determined by the `seg_id`
#' and `comparison_id` arguments, for the `dna_segs` and `comparisons`, 
#' respectively. `comparison_id` refers to 2 columns, and defaults to `"auto"`, 
#' which will attempt to determine which columns to use automatically.
#' If for example, `comparison_id` is set as `"name"`, it will look for the
#' `"name1"` and `"name2"` columns to match to the `seg_id` in the `dna_segs`
#' above, and under it, respectively.
#' 
#' @returns A list with 2 named elements: `dna_segs` and `comparisons`, which
#' are both lists containing the `dna_seg` and `comparison` objects, 
#' respectively.
#' @export
#' 
#' @param dna_segs A list of `dna_seg` objects.
#' @param comparisons A list of `comparison` objects.
#' @param seg_id The name of a `dna_seg` column, whose values will be used to
#' make the links to the `comparisons`.
#' @param comparison_id The shared name of the `comparison` columns, whose 
#' values will be used to make the links to the `dna_segs`. See details.
#' @param color_var A character string denoting which color attribute to update
#' for the `dna_segs`, one of: `"fill"`, `"col"`.
#' @param default_color A character string providing a default color, must be 
#' either `NULL` or a valid color. The color given by this argument will be
#' ignored when updating, never overwriting any other color.
#' @param both_directions Logical. If `FALSE`, the color scheme will be applied
#' sequentially in plotting order, starting from the first `dna_seg`. When 
#' `both_directions` is `TRUE`, it will then additionally update each `dna_seg`
#' and `comparison` in reverse plotting order.
#' 
#' @author Mike Puijk
#' 
#' @seealso [gradient_color_scheme], [uniform_color_scheme], [plot_gene_map],
#' [dna_seg], [comparison]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("1A", "1B", "1C")
#' names2 <- c("2A", "2C", "2B")
#' names3 <- c("3B", "3A", "3C")
#' 
#' ## Make dna_segs with some alternate colors
#' dna_seg1 <- dna_seg(data.frame(name = names1,
#'                                start = (1:3) * 3,
#'                                end = (1:3) * 3 + 2,
#'                                strand = rep(1, 3),
#'                                fill = c("darkred", "grey80", "darkblue")))
#' dna_seg2 <- dna_seg(data.frame(name = names2,
#'                                start = (1:3) * 3,
#'                                end = (1:3) * 3 + 2,
#'                                strand = rep(1, 3),
#'                                fill = c("grey80", "grey80", "darkgreen")))
#' dna_seg3 <- dna_seg(data.frame(name = names3,
#'                                start = (1:3) * 3,
#'                                end = (1:3) * 3 + 2,
#'                                strand = rep(1, 3)))
#' 
#' ## Make comparisons
#' comp1 <- comparison(data.frame(start1 = c(3, 6, 9), end1 = c(5, 8, 11),
#'                                start2 = c(3, 9, 6), end2 = c(5, 11, 8),
#'                                name1 = c("1A", "1B", "1C"), 
#'                                name2 = c("2A", "2B", "2C"),
#'                                direction = c(1, 1, 1)))
#' comp2 <- comparison(data.frame(start1 = c(3, 9, 6), end1 = c(5, 11, 8),
#'                                start2 = c(6, 3, 9), end2 = c(8, 5, 11),
#'                                name1 = c("2A", "2B", "2C"), 
#'                                name2 = c("3A", "3B", "3C"),
#'                                direction = c(1, 1, 1)))
#' 
#' ## Before adding a color scheme
#' plot_gene_map(dna_segs = list(dna_seg1, dna_seg2, dna_seg3), 
#'               comparisons = list(comp1, comp2),
#'               alpha_comparisons = 0.6)
#' 
#' ## Sequential color scheme without going both directions
#' full_data <- sequential_color_scheme(list(dna_seg1, dna_seg2, dna_seg3),
#'                                      comparisons = list(comp1, comp2),
#'                                      seg_id = "name",
#'                                      both_directions = FALSE)
#' plot_gene_map(dna_segs = full_data$dna_segs, 
#'               comparisons = full_data$comparisons,
#'               alpha_comparisons = 0.6)
#' 
#' ## Sequential color scheme with both directions
#' full_data <- sequential_color_scheme(list(dna_seg1, dna_seg2, dna_seg3),
#'                                      comparisons = list(comp1, comp2),
#'                                      seg_id = "name")
#' plot_gene_map(dna_segs = full_data$dna_segs, 
#'               comparisons = full_data$comparisons,
#'               alpha_comparisons = 0.6)
#' 
sequential_color_scheme <- function(
  dna_segs,
  comparisons, 
  seg_id = "locus_id",
  comparison_id = "auto",
  color_var = "fill",
  default_color = "grey80",
  both_directions = TRUE
) {
  
  sequential_updates(dna_segs = dna_segs,
                     comparisons = comparisons, 
                     seg_id = seg_id,
                     comparison_id = comparison_id,
                     color_var = color_var,
                     default_color = default_color,
                     both_directions = both_directions,
                     update_region_plot = FALSE,
                     update_positions = FALSE
                     )
  
}

#' Apply a uniform color scheme to dna_segs and comparisons
#' 
#' Generates and applies a uniform color scheme to a list of `dna_seg` and/ or
#' `comparison` objects. This is done by generating a set of colors based
#' on a set of unique values taken from a given (shared) column.
#' 
#' @details
#' The amount of colors generated is equal to the amount of unique values
#' present in the `id_column` in the data set(s). 
#' 
#' The `colors` argument can be
#' left as `NULL`, in which case one of 3 different palettes will be used, based
#' on the amount of colors that need to be generated. If not `NULL`, then 
#' `colors` must be either a character string representing a known palette 
#' (see [hcl.pals] and [palette.pals]), or a character vector with values
#' that can be recognized as colors (i.e. `"red"` or `"#88CCEE"`). If there are
#' not enough colors in this vector, then they will be duplicated and a warning
#' will be given. 
#' 
#' @returns If only `dna_segs` is provided, a list of `dna_seg` objects. 
#' @returns If only `comparisons` is provided, a list of `comparison`
#' objects.
#' @returns If both `dna_segs` and `comparisons` are provided, a list with 2
#' named elements: `dna_segs` and `comparisons`, which are both lists containing
#' the `dna_seg` and `comparison` objects, respectively.
#' 
#' @export
#' 
#' @param dna_segs A list of `dna_seg` objects.
#' @param comparisons A list of `comparison` objects.
#' @param id_column The name of a column, whose unique values will be used to 
#' generate a color scheme. Must exist in all provided data sets.
#' @param ids A character vector of values from the `id_column` to generate 
#' colors for. Any values found in the `id_column` that are not in this set will
#' be ignored.
#' @param cluster_ids Logical. If `TRUE`, numbered values may be clustered 
#' together, sharing a single color. Specifically, it will look for values that
#' end in a number, and when these values are the same once the number is 
#' removed, it will cluster them together. For example: `"lac_1"`, `"lac-2"`, 
#' and `"lac"` will all be given the same color. Finally, the value `"-"` 
#' (which is usually the value in the `gene` column for hypothetical proteins)
#' will also be given its own dedicated color.
#' @param colors Choice of colors to use, can be a palette, or a vector of
#' colors. See details. 
#' @param color_var A character string denoting which color attribute to update
#' for the `dna_segs`, one of: `"fill"`, `"col"`.
#' @param alpha_dna_segs A single numeric value between 0 and 1, or `NULL`. 
#' Determines the transparency applied to the colors used for the `dna_segs`, 0
#' being fully transparent, and 1 being fully opaque.
#' @param alpha_comparisons A single numeric value between 0 and 1, or `NULL`. 
#' Determines the transparency applied to the colors used for the `comparisons`,
#' 0 being fully transparent, and 1 being fully opaque.
#' 
#' @author Mike Puijk
#' 
#' @seealso [gradient_color_scheme], [sequential_color_scheme], [plot_gene_map],
#' [dna_seg], [comparison]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("A", "B", "C")
#' names2 <- c("A", "C", "B")
#' names3 <- c("B", "A", "C")
#' 
#' ## Make dna_segs with some alternate colors
#' dna_seg1 <- dna_seg(data.frame(name = names1,
#'                                start = (1:3) * 3,
#'                                end = (1:3) * 3 + 2,
#'                                strand = rep(1, 3)))
#' dna_seg2 <- dna_seg(data.frame(name = names2,
#'                                start = (1:3) * 3,
#'                                end = (1:3) * 3 + 2,
#'                                strand = rep(1, 3)))
#' dna_seg3 <- dna_seg(data.frame(name = names3,
#'                                start = (1:3) * 3,
#'                                end = (1:3) * 3 + 2,
#'                                strand = rep(1, 3)))
#' 
#' ## Make comparisons
#' comp1 <- comparison(data.frame(start1 = c(3, 6, 9), end1 = c(5, 8, 11),
#'                                start2 = c(3, 9, 6), end2 = c(5, 11, 8),
#'                                name = c("A", "B", "C"),
#'                                direction = c(1, 1, 1)))
#' comp2 <- comparison(data.frame(start1 = c(3, 9, 6), end1 = c(5, 11, 8),
#'                                start2 = c(6, 3, 9), end2 = c(8, 5, 11),
#'                                name = c("A", "B", "C"),
#'                                direction = c(1, 1, 1)))
#' 
#' ## Before adding a color scheme
#' plot_gene_map(dna_segs = list(dna_seg1, dna_seg2, dna_seg3), 
#'               comparisons = list(comp1, comp2),
#'               alpha_comparisons = 0.6)
#' 
#' ## Apply uniform color scheme
#' full_data <- uniform_color_scheme(list(dna_seg1, dna_seg2, dna_seg3),
#'                                   comparisons = list(comp1, comp2),
#'                                   id_column = "name")
#' plot_gene_map(dna_segs = full_data$dna_segs, 
#'               comparisons = full_data$comparisons,
#'               alpha_comparisons = 0.6)
#' 
#' ## Use ids and colors to specify the values and colors to use
#' full_data <- uniform_color_scheme(list(dna_seg1, dna_seg2, dna_seg3),
#'                                   comparisons = list(comp1, comp2),
#'                                   id_column = "name",
#'                                   ids = c("A", "C"),
#'                                   colors = c("purple", "#F8A19F"))
#' plot_gene_map(dna_segs = full_data$dna_segs, 
#'               comparisons = full_data$comparisons,
#'               alpha_comparisons = 0.6)
#' 
#' 
uniform_color_scheme <- function(
  dna_segs = NULL,
  comparisons = NULL,
  id_column,
  ids = NULL,
  cluster_ids = TRUE,
  colors = NULL,
  color_var = "fill",
  alpha_dna_segs = NULL,
  alpha_comparisons = NULL
) {
  
  # Check mandatory arguments
  if (is.null(dna_segs) & is.null(comparisons)) {
    stop("Either 'dna_segs' or 'comparisons' must be provided")
  }
  # Check dna_segs and set flags for later
  if (!is.null(dna_segs)) {
    if (!is.list(dna_segs) || !all(sapply(dna_segs, is.dna_seg))) {
      stop("'dna_segs' must be a list of dna_seg objects")
    }
    if (!all(sapply(dna_segs, function(x) id_column %in% names(x)))) {
      stop("'id_column' must refer to a column that is present in all ", 
           "provided dna_seg objects")
    }
    dna_segs <- copy(dna_segs)
  }
  # Check comparisons and set flags for later
  if (!is.null(comparisons)) {
    if (!is.list(comparisons) || !all(sapply(comparisons, is.comparison))) {
      stop("'comparisons' must be a list of comparison objects")
    }
    if (!all(sapply(comparisons, function(x) id_column %in% names(x)))) {
      stop("'id_column' must refer to a column that is present in all ", 
           "provided comparison objects")
    }
    comparisons <- copy(comparisons)
  }
  
  # Prepare ids and determine amount of colors
  if (!is.null(ids)) {
    if (!is.character(ids)) stop("'ids' must be a character vector")
    id_vector <- unique(ids)
  } else {
    if (!is.null(dna_segs)) {
      seg_id_list <- sapply(dna_segs,
                            function(x) x[, c(id_column), with = FALSE])
    } else {
      seg_id_list <- NULL
    }
    
    if (!is.null(comparisons)) {
      comp_id_list <- sapply(comparisons,
                             function(x) x[, c(id_column), with = FALSE])
    } else {
      comp_id_list <- NULL
    }
    id_vector <- unique(unlist(c(seg_id_list, comp_id_list)))
    id_vector <- id_vector[id_vector != "NA"]
    id_vector <- sort(id_vector)
  }
  
  if (cluster_ids) {
    # Remove trailing number ranges, then check for duplicates
    trimmed_ids <- sub("[_-]?[0-9]+$", "", id_vector)
    any_dups <- any(duplicated(trimmed_ids))
    if (any_dups) {
      # If there are duplicates, only make colors for the unique ones
      original_ids <- id_vector
      id_vector <- unique(trimmed_ids)
    }
  }
  
  n_colors <- length(id_vector)
  if (cluster_ids && any("-" == id_vector)) n_colors <- n_colors - 1

  # Check colors and make vector of colors
  if (is.null(colors)) {
    # Is colors not provided? Use a default palette
    if (n_colors <= 11) {
      # A colorblind safe palette
      color_vec <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", 
                     "#AA4499", "#44AA99", "#999933", "#882255", "#661100", 
                     "#6699CC")[1:n_colors]
    } else if (n_colors <= 34) {
      # A very divergent palette, first 2 colors are omitted (both grey)
      color_vec <- palette.colors(n_colors + 2,
                                  palette = "polychrome 36")[3:(n_colors + 2)]
    } else {
      max_n <- ifelse(n_colors > 266, 266, n_colors)
      color_vec <- rep_len(hcl.colors(max_n, palette = "dark 3"), n_colors)
    }
  } else if (length(colors) == 1 &
             all(tolower(colors) %in% c(tolower(hcl.pals()),
                                        tolower(palette.pals())
                                        ))
  ) {
    # Is colors a single character string and a known palette, use that
    colors <- tolower(colors)
    if (any(colors == tolower(palette.pals()))) {
      max_n <- length(unique(palette.colors(n_colors, palette = colors)))
      
      color_vec <- palette.colors(n_colors, palette = colors, recycle = TRUE)
    } else {
      max_n <- 0
      # Find where the palette stops producing unique colors consistently
      while (n_colors > max_n) {
        if (length(unique(hcl.colors(max_n, palette = colors))) < max_n) {
          max_n <- max_n - 1
          break
        }
        max_n <- max_n + 1
      }
      color_vec <- rep_len(hcl.colors(max_n, palette = colors), n_colors)
    }
    if (n_colors > max_n) {
      warning('The maximum amount of colors for "', colors, '" is ', max_n,
              ' but ', n_colors, ' were needed, so colors were duplicated ', 
              'as needed')
    }
  } else if (all(sapply(colors,
                        function(x) tryCatch(is.matrix(col2rgb(x)),
                                             error = function(e) FALSE
                                             )
                        ))
             ) {
    # Does colors consist of elements that can be recognized as colors by R?
    if (length(colors) < n_colors) {
      color_vec <- rep_len(colors, length.out = n_colors)
      warning(colors, " colors were provided but ", n_colors, " are needed, ",
              "so the provided colors were duplicated")
    } else {
      color_vec <- colors[1:n_colors]
    }
  } else {
    stop("If 'colors' is provided, it must be a single character consisting ", 
         "of a known palette (one of palette.pals() or hcl.pals()), or a ", 
         "vector of valid colors")
  }
  
  if (cluster_ids) {
    if (any("-" == id_vector)) {
      # "-" isn't explicitly NA, in the gene column it often represents
      # hypothetical proteins, so a special case is made here to color those
      id_vector <- id_vector[-which("-" == id_vector)]
      id_vector <- c(id_vector, "-")
      color_vec <- c(color_vec, "#D2D0B6")
    }
    
    if (any_dups) {
      # If there were duplicates among the trimmed IDs, extend color_vec to 
      # account for the grouped IDs
      color_vec <- sapply(
        1:length(original_ids),
        function (x) color_vec[which(trimmed_ids[x] == id_vector)]
      )
      id_vector <- original_ids
    }
  }
  
  # Apply alpha to color vector as needed for dna_segs
  if (!is.null(dna_segs)) {
    if (!is.null(alpha_dna_segs)) {
      # Convert color_vec into hexadecimal
      rgb_matrix <- col2rgb(color_vec)
      hex_color_vec <- rgb(red = rgb_matrix[1, ], green = rgb_matrix[2, ], 
                           blue = rgb_matrix[3, ], maxColorValue = 255)
      # Convert ratio into hexadecimal
      if (alpha_dna_segs > 1 | alpha_dna_segs < 0) {
        stop("'alpha_dna_segs' must be between 0 and 1")
      }
      tpc <- floor(alpha_dna_segs * 256)
      tpc <- sprintf("%X", tpc)
      if (nchar(tpc) == 1) tpc <- paste("0", tpc, sep = "")
      seg_color_table <- data.table(V1 = id_vector,
                                    V2 = paste(hex_color_vec, tpc, sep = ""))
    } else {
      seg_color_table <- data.table(V1 = id_vector, V2 = color_vec)
    }
    setnames(seg_color_table,
             old = c("V1", "V2"),
             new = c(id_column, color_var)
             )
  }
  # Do the same for comparisons
  if (!is.null(comparisons)) {
    if (!is.null(alpha_comparisons)) {
      # Convert color_vec into hexadecimal
      rgb_matrix <- col2rgb(color_vec)
      hex_color_vec <- rgb(red = rgb_matrix[1, ], green = rgb_matrix[2, ], 
                           blue = rgb_matrix[3, ], maxColorValue = 255)
      # Convert ratio into hexadecimal
      if (alpha_comparisons > 1 || alpha_comparisons < 0) {
        stop("'alpha_comparisons' must be between 0 and 1")
      }
      tpc <- floor(alpha_comparisons * 256)
      tpc <- sprintf("%X", tpc)
      if (nchar(tpc) == 1) tpc <- paste("0", tpc, sep = "")
      comp_color_table <- data.table(V1 = id_vector,
                                     V2 = paste(hex_color_vec, tpc, sep = ""))
    } else {
      comp_color_table <- data.table(V1 = id_vector, V2 = color_vec)
    }
    setnames(comp_color_table,
             old = c("V1", "V2"),
             new = c(id_column, color_var)
             )
  }

  # Loop through and update dna_segs if provided
  if (!is.null(dna_segs)) {
    for (i in 1:length(dna_segs)) {
      if (color_var == "col") {
        seg_merge <- seg_color_table[dna_segs[[i]],
                                     on = id_column
                                     ][is.na(col),
                                       col := i.col
                                       ][, i.col := NULL]
      } else if (color_var == "fill") {
        seg_merge <- seg_color_table[dna_segs[[i]],
                                     on = id_column
                                     ][is.na(fill),
                                       fill := i.fill
                                       ][, i.fill := NULL]
      }
      setcolorder(seg_merge, neworder = names(dna_segs[[i]]))
      dna_segs[[i]] <- as.dna_seg(seg_merge)
    }
  }
  
  # Loop through and update comparisons if provided
  if (!is.null(comparisons)) {
    for (i in 1:length(comparisons)) {
      if (color_var == "col") {
        comp_merge <- comp_color_table[comparisons[[i]],
                                       on = id_column
                                       ][is.na(col),
                                         col := i.col
                                         ][, i.col := NULL]
      } else if (color_var == "fill") {
        comp_merge <- comp_color_table[comparisons[[i]],
                                       on = id_column
                                       ][!is.na(fill),
                                         col := fill
                                         ][,
                                           fill := NULL
                                           ][is.na(col), col := "NA"]
      }
      setcolorder(comp_merge, neworder = names(comparisons[[i]]))
      comparisons[[i]] <- as.comparison(comp_merge)
    }
  }
  
  # Return dna_segs and comparisons as appropriate
  if (is.null(dna_segs)) return(comparisons)
  if (is.null(comparisons)) return(dna_segs)
  list(dna_segs = dna_segs, comparisons = comparisons)
}

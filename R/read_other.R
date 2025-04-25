################################################################################
# Other read functions
################################################################################

#' Parse `dna_seg` offset values for plotting
#' 
#' This function takes a character string and separates it to parse offset
#' values for plotting in [plot_gene_map]. The offset values determine the
#' amount of space before each subsegment of each plotted `dna_seg`.
#' 
#' @details
#' The format of `offsets` corresponds to the output of the `print_offsets`
#' argument of [plot_gene_map].
#' `offsets` must be a character string containing 1 or more offset values for
#' each `dna_seg`. `dna_segs` are separated by commas, in plotting order,
#' with their subsegments separated by spaces. For each `dna_seg`, if a single
#' value is provided, then that value will become the offset for its first 
#' subsegment. If more values are provided for a `dna_seg`, then there must 
#' be as many values as there are subsegments for that dna_seg.',
#' 
#' @returns A list of offset values that can be passed directly to the `offsets`
#' argument of [plot_gene_map]. Each element of the list is a numeric vector
#' that corresponds to the offset values of one `dna_seg`.
#' @export
#' 
#' @param offsets A character string with offset values to parse. See details.
#' 
#' @author Mike Puijk
#' 
#' @seealso [read_offsets_file], [plot_gene_map]
#' 
#' @examples
#' ## Parse offset values for 2 dna_segs
#' ## dna_seg 1: its 2 subsegments will be offset by 150 and 1300, respectively
#' ## dna_seg 2: only its first subsegment  will be offset by 300
#' read_offsets("150 1300,300")
#' 
read_offsets <- function(offsets) {
  lapply(strsplit(trimws(unlist(strsplit(offsets, ","))), " "),
         function(x) as.numeric(x))
}

#' Parse `dna_seg` offset values for plotting from a file
#' 
#' This function parses a file containing offset
#' values for plotting in [plot_gene_map]. The offset values determine the
#' amount of space before each subsegment of each plotted `dna_seg`.
#' 
#' @details
#' The file provided by `offsets` must contain 1 or more offset values for
#' each `dna_seg`. Each line must contain offset values for 1 `dna_seg`, in 
#' plotting order, with their subsegments separated by spaces. For each 
#' `dna_seg`, if a single value is provided, then that value will become the
#' offset for its first subsegment. If more values are provided for a `dna_seg`,
#' then there must be as many values as there are subsegments for that 
#' `dna_seg`.
#' 
#' @returns A list of offset values that can be passed directly to the `offsets`
#' argument of [plot_gene_map]. Each element of the list is a numeric vector
#' that corresponds to the offset values of one `dna_seg`.
#' @export
#' 
#' @param offsets A character string containing a file path, or a file 
#' connection. See details for file format.
#' 
#' @author Mike Puijk
#' 
#' @seealso [read_offsets], [plot_gene_map]
#' 
#' @examples
#' ## Parse offset values for 2 dna_segs from a file
#' offsets_file <- system.file('extdata/offsets.txt', package = 'genoPlotR')
#' 
#' ## When this file is parsed:
#' ## dna_seg 1: its 2 subsegments will be offset by 150 and 1300, respectively
#' ## dna_seg 2: only its first subsegment  will be offset by 300
#' 
#' ## Example of an offsets file:
#' cat(readLines(offsets_file), sep = "\n")
#' 
#' ## Parse offsets file
#' read_offsets_file(offsets_file)
#' 
#' ## This is equivalent to:
#' read_offsets("150 1300,300")
#' 
read_offsets_file <- function(offsets) {
  lapply(strsplit(trimws(readLines(offsets)), " "), function (x) as.numeric(x))
}

#' Parse `xlims`, a set of `dna_seg` positions for plotting
#' 
#' This function takes a character string containing a set of `xlims`, positions
#' to print for `dna_segs`. Multiple start and end positions can be provided for
#' each `dna_seg`, forming subsegments to plot.
#' 
#' @details
#' The format of `xlims` corresponds to the output of the `print_xlims`
#' argument of [plot_gene_map].
#' `xlims` must be a character string containing a set of `xlims` for each 
#' `dna_seg`. For each `dna_seg`, start and end positions must be provided for
#' each subsegment, all separated by spaces. The `dna_segs` themselves are 
#' separated by commas, in plotting order. An example of the format is provided
#' below.
#' 
#' @returns A list of `xlims` that can be passed directly to the `xlims`
#' argument of [plot_gene_map]. Each element of the list is a numeric vector
#' that corresponds to the `xlims` of one `dna_seg`.
#' @export
#' 
#' @param xlims A character string with `xlims` to parse. See details.
#' 
#' @author Mike Puijk
#' 
#' @seealso [read_xlims_file], [plot_gene_map]
#' 
#' @examples
#' ## Parse xlims for 2 dna_segs
#' ## For dna_seg 1: 2 subsegments, 1st starts at 1 and ends at 10000
#' ## while the 2nd subsegment starts at 155000 and ends at 154000
#' ## For dna_seg 2: 1 subsegment, starting at 500 and ending at 10500
#' read_xlims("1 10000 155000 154000,500 10500")
#' 
read_xlims <- function(xlims) {
  lapply(strsplit(trimws(unlist(strsplit(xlims, ","))), " "),
         function(x) as.numeric(x))
}

#' Parse `xlims` from a file, a set of `dna_seg` positions for plotting
#' 
#' This function parses a file containing a set of `xlims`, positions
#' to print for `dna_segs`. Multiple start and end positions can be provided for
#' each `dna_seg`, forming subsegments to plot.
#' 
#' @details
#' The file format of `xlims` corresponds to the output of the `outfile_xlims`
#' argument of [plot_gene_map]. The format in question is a 
#' tab-separated file containing `xlims` for each `dna_seg`. Each row
#' represents one subsegment and has 3 columns, containing the start position,
#' end position, and `dna_seg` label, respectively. An example of the format:
#' 
#' \tabular{lllll}{
#' x0    \tab x1    \tab seg_label\cr
#' 2000  \tab 12000 \tab Genome_2 \cr
#' 1000  \tab 10000 \tab Genome_1 \cr
#' 6000  \tab 16000 \tab Genome_3 \cr
#' 30000 \tab 40000 \tab Genome_2 \cr
#' }
#' 
#' @returns A list of `xlims` that can be passed directly to the `xlims`
#' argument of [plot_gene_map]. Each element of the list is a numeric vector
#' that contains the positions to plot of one `dna_seg`.
#' @returns If `reorder_dna_segs = TRUE`, a list with 2 named elements will 
#' be returned instead (`xlims` 
#' and either `dna_segs` or `seg_labels`, depending on the `dna_segs` argument).
#' @export
#' 
#' @param dna_segs Either a character vector containing `dna_seg` labels, or a
#' list of `dna_seg` objects.
#' @param xlims A character string containing a file path, or a file 
#' connection. See details for file format.
#' @param reorder_dna_segs Logical. If `TRUE`, then the order of `dna_segs` will
#' be changed to match the order in which they first appear in the `xlims` file.
#' Additionally, any `dna_seg` that was not present in the `xlims` file will be 
#' removed. Then, this function will return a list with the reordered 
#' `dna_segs` and parsed `xlims` as elements, instead of just the `xlims`.
#' 
#' @author Mike Puijk
#' 
#' @seealso [read_xlims], [plot_gene_map]
#' 
#' @examples
#' ## Parse xlims for 3 dna_segs from a file
#' dna_seg_labels <- c("Genome_1", "Genome_2", "Genome_3")
#' xlims_file <- system.file('extdata/xlims.tab', package = 'genoPlotR')
#' 
#' ## Example of an xlims file:
#' cat(readLines(xlims_file), sep = "\n")
#' 
#' ## Parse xlims file, keep order provided by dna_seg_labels
#' read_xlims_file(dna_seg_labels, xlims_file)
#' 
#' ## Parse xlims file, reorder dna_segs based on the xlims file
#' xlims_read <- read_xlims_file(dna_seg_labels, xlims_file, 
#'                               reorder_dna_segs = TRUE)
#' xlims_read$xlims
#' xlims_read$seg_labels
#' 
read_xlims_file <- function(dna_segs, xlims, reorder_dna_segs = FALSE) {
  # Check mandatory arguments
  if (missing(dna_segs)) stop("'dna_segs' must be provided")
  if (is.character(dna_segs)) {
    seg_labels <- dna_segs
    return_list <- FALSE
  } else if (is.list(dna_segs) && all(sapply(dna_segs, is.dna_seg))) {
    seg_labels <- get_seg_labels(dna_segs, seg_labels = NULL)
    return_list <- TRUE
    dna_segs <- copy(dna_segs)
    if (is.null(seg_labels)) {
      stop("dna_seg labels could not be determined from this dna_seg list")
    }
  } else {
    stop("'dna_segs' must be a list of dna_seg objects or a character ",
         "vector containing the seg labels (i.e. names(dna_segs))")
  }
  
  # If xlims argument is a character string, it is assumed to be a file path
  if (missing(xlims)) stop("'xlims' must be provided")
  if (is.character(xlims) && file.exists(xlims)) {
    xlims <- fread(file = xlims) 
  } else if (is.data.frame(xlims)) {
    xlims <- xlims
  } else {
    stop("'xlims' must be either a character string (specifying a file path) ",
         "or a data.frame")
  }
  
  xlim_labels <- unique(xlims$seg_label)
  if (reorder_dna_segs) {
    if (!all(xlim_labels %in% seg_labels)) {
      ### Add more information to this error message. Which are missing?
      stop("xlim seg_labels do not match dna_seg seg_labels")
    }
    # Loops through xlim labels, takes xlims for each, formats them into list 
    xlims <- lapply(
      xlim_labels,
      function(x) as.numeric(strsplit(
        paste(
          xlims[seg_label == x, c("x0"), with = FALSE][,x0],
          xlims[seg_label == x, c("x1"), with = FALSE][,x1],
          sep = ",", collapse = ","
        ), split = ",")[[1]])
    )
    if (return_list) {
      dna_return <- lapply(xlim_labels, function(x) dna_segs[[x]])
      names(dna_return) <- xlim_labels
      return(list(dna_segs = dna_return, xlims = xlims))
    } else {
      return(list(seg_labels = xlim_labels, xlims = xlims))
    }
    
  } else {
    if (length(xlim_labels) != length(seg_labels) | 
        !all(seg_labels %in% xlim_labels)) {
      stop("xlim seg_labels do not match dna_seg seg labels")
    }
    # Loops through seg_labels, takes xlims for each seg, formats them into list 
    return(lapply(
      seg_labels,
      function(x) as.numeric(strsplit(
        paste(
          xlims[seg_label == x, c("x0"), with = FALSE][,x0],
          xlims[seg_label == x, c("x1"), with = FALSE][,x1],
          sep = ",", collapse = ","
        ), split = ",")[[1]])
    ))
  }
}

################################################################################
# File reading functions: read mauve backbone
################################################################################

#' Read mauve backbone files to create dna_segs and comparisons
#' 
#' Parses Progressive Mauve backbone files to create a list of `dna_segs`, and a
#' matching list of `comparisons` between them. 
#' 
#' @details
#' Mauve Backbone files are tabular files that summarize similarities between
#' genomes in blocks. Each genome has 2 columns containing the start and end
#' coordinates of each block, respectively. The header, if present, uses
#' sequence numbers instead of genome names. See 
#' https://darlinglab.org/mauve/user-guide/files.html for more info on the file
#' format. This function should be able to read both `progressiveMauve` and 
#' `mauveAligner` outputs.
#' 
#' @returns A list with 2 named elements: `dna_segs` and `comparisons`, which
#' are both lists containing the `dna_seg` and `comparison` objects, 
#' respectively.
#' @export
#' 
#' @param file A character string containing a file path, or a file connection.
#' @param ref A numeric indicating which DNA segment to use as the reference, 
#' i.e. which one will have its blocks in order.
#' @param gene_type A character string, determines how genes are visualized.
#' Must be a valid gene type (see [gene_types]).
#' @param header Logical. If `TRUE`, parses the first line of the file as a
#' header containing column names. 
#' @param filter_low A numeric, if larger than 0, all blocks smaller than this
#' number will be filtered out.
#' @param common_blocks_only Logical. If TRUE, reads only common blocks (core 
#' blocks).
#' @param ... Further arguments to pass to [as.dna_seg].
#' 
#' @author Lionel Guy, Jens Roat Kultima
#' 
#' @references Mauve: <https://darlinglab.org/mauve/mauve.html>
#' 
#' @seealso [dna_seg], [comparison]
#' 
#' @examples
#' ## Mauve backbone
#' bbone_file <- system.file('extdata/barto.backbone', package = 'genoPlotR')
#' bbone <- read_mauve_backbone(bbone_file)
#' names <- c("B_bacilliformis", "B_grahamii", "B_henselae", "B_quintana")
#' names(bbone$dna_segs) <- names
#' 
#' ## Plot
#' plot_gene_map(dna_segs = bbone$dna_segs, comparisons = bbone$comparisons)
#' 
#' ## Using filter_low & changing reference sequence
#' bbone <- read_mauve_backbone(bbone_file, ref = 2, filter_low = 2000) 
#' names(bbone$dna_segs) <- names
#' plot_gene_map(dna_segs = bbone$dna_segs, comparisons = bbone$comparisons)
#' 
read_mauve_backbone <- function(
  file,
  ref = 1,
  gene_type = "side_blocks",
  header = TRUE,
  filter_low = 0,
  common_blocks_only = TRUE,
  ...
) {
  
  blocks <- read.table(file, stringsAsFactors = FALSE, header = header)
  n_orgs <- ncol(blocks) / 2
  n_blocks <- nrow(blocks)
  # check ref parameter
  if (!(ref %in% 1:n_orgs)) {
    stop("'ref' must be a positive integer =< row number")
  }
  # ordering from ref:
  for (i in 1:n_blocks) {
    if (blocks[i, ref * 2] < 0) blocks[i,] <- -blocks[i,]
  }
  if (any(blocks[, ref * 2] < 0)) {
    stop("Not all rows in ref columns are positive. Contact author.")
  }
  blocks <- blocks[do.call(order, list(blocks[,ref*2])),]
  # filter, if needed
  sizes <- matrix(NA, nrow = n_blocks, ncol = n_orgs)
  for (i in 1:n_orgs) {
    sizes[, i] <- abs(blocks[, 2 * i]) - abs(blocks[, 2 * i - 1])
  }
  min_sizes <- apply(sizes, 1, min, na.rm = TRUE)
  rows_to_keep <- rep(TRUE, n_blocks)
  if (common_blocks_only) {
    rows_to_keep <- (min_sizes > 0)
  }
  if (is.numeric(filter_low) && filter_low > 1) {
    rows_to_keep <- rows_to_keep & (min_sizes >= filter_low)
  }
  blocks <- blocks[rows_to_keep,]
  n_blocks <- nrow(blocks)

  # check that there are rows remaining
  if (nrow(blocks) < 1) {
    stop("No row left in data. Check data and/or lower filter_low argument")
  }
  # prepare objects
  dna_segs <- list()
  comparisons <- list()
  # colors: rainbow starting from the beginning
  col <- rainbow(n = n_blocks)
  # run through organisms
  for (i in 1:n_orgs) {
    # prepare dna_seg
    s <- blocks[,2*i-1]
    e <- blocks[,2*i]
    strand <- sign(s)
    df <- data.frame(name = paste("block", 1:n_blocks, sep = "_"),
                     start = abs(s), end = abs(e), strand = strand,
                     stringsAsFactors = FALSE)
    df$col <- col
    if (all(df$strand == 0)) {
      stop("One or several of the columns is composed only of 0. Check data")
    }
    dna_segs[[i]] <- as.dna_seg(df[df$strand != 0,], gene_type = gene_type, ...)
    # prepare comparison (not with i=1)
    if (i > 1) {
      s0 <- blocks[, 2 * (i - 1) - 1]
      e0 <- blocks[, 2 * (i - 1)]
      strand0 <- sign(s0)
      df <- data.frame(start1 = ifelse(strand0 == 1, abs(s0), abs(e0)),
                       end1 = ifelse(strand0 == 1, abs(e0), abs(s0)),
                       start2 = ifelse(strand == 1, abs(s), abs(e)),
                       end2 = ifelse(strand == 1, abs(e), abs(s)))
      comparison <- as.comparison(df[df$start1 != 0 & df$start2 != 0,])
      # apply red_blue color scheme
      #if (i == n_orgs) browser()
      comparison$col <- gradient_color_scheme(x = NULL,
                                              direction = comparison$direction,
                                              color_scheme = "red_blue")
      comparisons[[i-1]] <- comparison
    }
  }
  list(dna_segs = dna_segs, comparisons = comparisons)
}

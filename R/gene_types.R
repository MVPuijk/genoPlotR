################################################################################
# Plot helpers: list and determine gene "types"
################################################################################

#' Lists gene_types for dna_seg objects
#' 
#' Returns a character vector containing the available gene types when plotting
#' `dna_seg` features. 
#' 
#' @details
#' `dna_segs` may contain the `gene_type` column, which determines the shape
#' of `dna_seg` features when they are plotted using [plot_gene_map]. 
#' Elements in this column should either be one of the predefined gene types 
#' returned by this function, or they must refer to a graphical function with 
#' exactly the same name that returns a `grob` or a `gList` object.
#' 
#' @returns A character vector.
#' @export
#' 
#' @param auto Logical. If `TRUE`, includes the `"auto"` gene type in the
#' output, which has [plot_gene_map] determine what gene type to use
#' automatically.
#' 
#' @author Lionel Guy
#' 
#' @seealso [plot_gene_map], [dna_seg]
#' 
#' @examples
#' ## To view pre-coded gene types:
#' gene_types()
#' 
#' ## Load data
#' data(barto)
#' n <- length(gene_types(auto = FALSE))
#' 
#' ## Get a small subset from the barto dataset
#' dna_seg <- barto$dna_segs[[3]][1:n,]
#' plot_gene_map(list(dna_seg))
#' 
#' ## Change gene_types and plot again
#' dna_seg$gene_type <- gene_types(auto = FALSE)
#' dna_seg$fill <- rainbow(n)
#' dna_seg_r <- dna_seg
#' dna_seg_r$strand <- -dna_seg$strand
#' 
#' ## Add an annotation
#' annot <- annotation(middle(dna_seg), text = dna_seg$gene_type, rot = 45,
#'                     col = dna_seg$col)
#' 
#' ## Plot
#' plot_gene_map(list(dna_seg, dna_seg_r), annotations = list(annot, annot),
#'               annotation_height = 5, dna_seg_line = grey(0.7))
#' 
gene_types <- function(auto = TRUE) {
  types <- c("arrows", "headless_arrows",
             "blocks", "bars", "points", "text", "lines", 
             "side_blocks", "side_bars", "side_points", "side_text",
             "side_lines",
             "introns", "exons", "side_exons",
             "boundaries"
             )
  if (auto) types <- c("auto", types)
  types
}

#' Automatically determine gene type
#' 
#' Determines the gene type to use based on how many genes are plotted.
#' 
#' @returns A character string containing a gene type.
#' @noRd
#' 
#' @param n_genes A numeric vector containing the amount of genes. The maximum
#' is taken when multiple values are provided.
#' 
auto_gene_type <- function(n_genes) {
  if (max(n_genes) > 1000) {
    gene_type <- "side_bars"
  } else if (max(n_genes) > 100) {
    gene_type <- "side_blocks"
  } else {
    gene_type <- "side_bars"
  }
  gene_type
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genoPlotR.R
\docType{data}
\name{barto}
\alias{barto}
\title{Comparison of 4 Bartonella genomes}
\format{
\code{barto}, a list of 3 data frame lists, representing the
four genomes and their pairwise comparisons:
\describe{
\item{\code{dna_segs}}{A list of 4 \code{dna_seg} objects, containing all the
protein-coding genes for each genome. Obtained by reading ptt files
downloaded from NCBI with \link{read_dna_seg_from_ptt}.}
\item{\code{comparisons}}{A list of 3 \code{comparison} objects, obtained by doing
genome-to-genome (fasta files) BLASTS, and then reading the resulting tab
files with \link{read_comparison_from_blast}.}
\item{\code{rnt_segs}}{A list of 4 \code{dna_seg} objects, containing all the RNA
genes of the four genomes. Obtained by reading rnt files downloaded from
NCBI with \link{read_dna_seg_from_ptt}.}
}
}
\usage{
data(barto)
}
\description{
A comparison of 4 Bartonella genomes by BLAST.
}
\examples{
data(barto)
plot_gene_map(barto$rnt_segs, barto$comparisons, gene_type = "blocks")

}
\keyword{datasets}

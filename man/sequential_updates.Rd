% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edit_data.R
\name{sequential_updates}
\alias{sequential_updates}
\title{Update dna_segs and comparisons sequentially in plotting order}
\usage{
sequential_updates(
  dna_segs,
  comparisons,
  seg_id = "locus_id",
  comparison_id = "auto",
  color_var = NULL,
  default_color = "grey80",
  update_region_plot = TRUE,
  update_positions = FALSE,
  both_directions = TRUE
)
}
\arguments{
\item{dna_segs}{A list of \code{dna_seg} objects.}

\item{comparisons}{A list of \code{comparison} objects.}

\item{seg_id}{The name of a \code{dna_seg} column, whose values will be used to
make the links to the \code{comparisons}.}

\item{comparison_id}{The shared name of the \code{comparison} columns, whose
values will be used to make the links to the \code{dna_segs}. See details.}

\item{color_var}{A character string denoting which color attribute to update,
must be one of: \code{"fill"}, \code{"col"}, or left as \code{NULL}, which means no color
will be updated.}

\item{default_color}{A character string providing a default color, must be
either \code{NULL} or a valid color. The color given by this argument will be
ignored when updating, never overwriting any other color.}

\item{update_region_plot}{Logical. If \code{TRUE}, updates the \code{region_plot}
attribute of \code{dna_segs} (and of \code{comparisons} but only to pass this
information to the \code{dna_segs} that follow), which determines whether the
neighborhood of these features is plotted in a regional plot.}

\item{update_positions}{Logical. If \code{TRUE}, updates the plotted positions
of the comparisons to match the \code{dna_segs}. \code{start1} and \code{end1} will be
updated using the \code{dna_seg} above the \code{comparison} in plotting order,
while \code{start2} and \code{end2} will be updated using the \code{dna_seg} under the
\code{comparison} in plotting order.}

\item{both_directions}{Logical. If \code{FALSE}, updates are applied
sequentially in plotting order, starting from the first \code{dna_seg}. When
\code{both_directions} is \code{TRUE}, it will then additionally update each \code{dna_seg}
and \code{comparison} in reverse plotting order.}
}
\value{
A list with 2 named elements: \code{dna_segs} and \code{comparisons}, which
are both lists containing the \code{dna_seg} and \code{comparison} objects,
respectively.
}
\description{
Takes a list of \code{dna_seg} and a list of \code{comparison} objects, and then
updates both sequentially, in plotting order. This can be used to update the
positions of the comparisons, as well as the color and \code{region_plot}
attributes of both the \code{dna_segs} and comparisons. It does this by taking the
values of features for these attributes and transferring those over to
comparisons directly connected to it, and then to the features connected to
those comparisons. For example, if a feature from a single \code{dna_seg} has red
as its \code{fill} attribute, the comparisons that can be linked to this feature
will become red as well. This is then followed up by updating any \code{dna_seg}
features linked to those comparisons, and so on.
}
\details{
When updating colors, the existing colors from the input \code{dna_seg} and
\code{comparison} objects are transferred over to the next object in the
plotting order, with the exception of their default colors, provided by
\code{default_color}. As \code{comparison} objects only have a single color
attribute \code{col}, those will be updated using the column provided by
\code{color_var} from the \code{dna_segs}, while the \code{dna_segs} themselves will be
updated using the \code{col} column from the comparisons regardless of \code{color_var}
(unless it is left as \code{NULL} to avoid updating colors entirely).

The objects are linked together through shared
values. The columns for these shared values are determined by the \code{seg_id}
and \code{comparison_id} arguments, for the \code{dna_segs} and \code{comparisons},
respectively. \code{comparison_id} refers to 2 columns, and defaults to \code{"auto"},
which will attempt to determine which columns to use automatically.
If for example, \code{comparison_id} is set as \code{"name"}, it will look for the
\code{"name1"} and \code{"name2"} columns to match to the \code{seg_id} in the \code{dna_segs}
above, and under it, respectively.

There are some circumstances with very interconnected comparisons where
you might want to set \code{both_directions} to \code{FALSE} to avoid transferring
over the \code{region_plot} attribute to too many features.
}
\examples{
## Prepare dna_seg
names1 <- c("1A", "1B", "1C")
names2 <- c("2A", "2C", "2B")
names3 <- c("3B", "3A", "3C")

## Make dna_segs with some alternate colors
dna_seg1 <- dna_seg(data.frame(name=names1, start=(1:3)*3, end=(1:3)*3+2,
                               strand=rep(1, 3),
                               fill=c("darkred", "grey80", "darkblue")))
dna_seg2 <- dna_seg(data.frame(name=names2, start=(1:3)*3, end=(1:3)*3+2,
                               strand=rep(1, 3),
                               fill=c("grey80", "grey80", "darkgreen")))
dna_seg3 <- dna_seg(data.frame(name=names3, start=(1:3)*3, end=(1:3)*3+2,
                               strand=rep(1, 3)))

## Make comparisons
comp1 <- comparison(data.frame(start1=c(3,6,9), end1=c(5,8,11),
                               start2=c(3,9,6), end2=c(5,11,8),
                               name1=c("1A", "1B", "1C"), 
                               name2=c("2A", "2B", "2C"),
                               direction=c(1,1,1)))
comp2 <- comparison(data.frame(start1=c(3,9,6), end1=c(5,11,8),
                               start2=c(6,3,9), end2=c(8,5,11),
                               name1=c("2A", "2B", "2C"), 
                               name2=c("3A", "3B", "3C"),
                               direction=c(1,1,1)))

## Before applying sequential_updates
plot_gene_map(dna_segs=list(dna_seg1, dna_seg2, dna_seg3), 
              comparisons=list(comp1, comp2),
              alpha_comparisons=0.6)

## Applying sequential_updates without updating in both directions
full_data <- sequential_updates(list(dna_seg1, dna_seg2, dna_seg3),
                                comparisons=list(comp1, comp2),
                                seg_id = "name",
                                color_var = "fill",
                                both_directions = FALSE)
plot_gene_map(dna_segs=full_data$dna_segs,
              comparisons=full_data$comparisons,
              alpha_comparisons=0.6)

## Applying sequential_updates using both directions
full_data <- sequential_updates(list(dna_seg1, dna_seg2, dna_seg3),
                                comparisons=list(comp1, comp2),
                                seg_id = "name",
                                color_var = "fill")
plot_gene_map(dna_segs=full_data$dna_segs, 
              comparisons=full_data$comparisons,
              alpha_comparisons=0.6)

}
\seealso{
\link{update_comparisons}, \link{update_dna_segs}, \link{dna_seg}, \link{comparison}
}
\author{
Mike Puijk
}

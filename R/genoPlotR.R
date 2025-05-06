################################################################################
# PLOT GENE MAPS                                                               #
################################################################################
# A R framework to plot comparison of gene stretches, a la Artemis, but
# with production-like graphics, and a static interface
################################################################################
# Plotting function
################################################################################


#' Plot gene and genome maps
#' 
#' This plotting function represents linear DNA segments and comparisons
#' between them. It will plot one line per DNA segment, eventually separated
#' by the comparisons. In addition, a tree can be plotted on the left of the
#' plot, as well annotations on each DNA segment. This function is based on the
#' grid package, and as such the plots can be placed into
#' other graphics or modified.
#' 
#' @details
#' One line is plotted per `dna_seg`. The shape of the `dna_seg` elements is
#' determined by the `gene_type` column, unless they are explicitly overwritten
#' using the `gene_type` argument. See [gene_types] for more details.  When 
#' provided, `comparisons` are placed
#' in between the `dna_segs`, annotations are displayed above each `dna_seg`,
#' and accompanying data for each `dna_seg` can be plotted using `seg_plots` 
#' (see [seg_plot] for more details).
#' 
#' A phylogenetic tree (a `phylog` object from package `ade4`) can be drawn at
#' the left of the plot (see [newick2phylog]). The tree does not need to be 
#' ordered in the same way as the `dna_seg_labels`, but a permutation of the
#' tree with that order should exist. If the tree is large, the number of 
#' permutations can become too large, causing the function to stop
#' (>100000 permutations). The solution is then to provide `dna_segs` that are
#' ordered in the same manner as the tree labels, or vice-versa (see 
#' [permute_dna_segs]). There is (experimental) support for branch annotations.
#' These are given in the Newick tree, directly after the parenthesis closing a
#' node. They can be characters or integers, but so far [newick2phylog] does
#' not support the `.` character in any labels. Tags will be ignored if they
#' start with "I", and trimmed if they start with "X".
#' 
#' There are two ways to draw a legend using this function:
#' 
#' The first is a 
#' customizable legend, which can be drawn by providing a vector of labels to 
#' the `legend_labels` argument and a vector of colors to link to these labels 
#' to the `legend_colors` argument. 
#' 
#' The other way is through the `legend_column`
#' argument, which must be `"auto"`, or the name of a column present in each
#' `dna_seg`. If `"auto"`, the column to use will be the same as the column
#' that was used for the `"uniform"` color scheme if this was applied, the 
#' `gene` column if this was present, or finally the `name` column if neither
#' of the previous options are possible. The legend will take all possible 
#' values for the chosen column, alongside the color of the first instance it
#' can find for each possible value. If `legend_labels` was provided, only
#' values found in that character vector are considered for the legend.
#' 
#' The `xlims` allow the user to plot subsegments of a `dna_seg`. `xlims` 
#' consists of a list composed of as many numeric vectors as there are 
#' `dna_segs`. Each of these numeric vectors give pairs of left and right 
#' borders, and gives the direction. For example, `c(1,2,6,4)` will plot two 
#' subsegments for a single `dna_seg`, the first subsegment will go from 1 to 2
#' which is plotted left to right and the second subsegment will go from 4 to 6,
#' plotted right to left. `-Inf` and `Inf` values are accepted. `NULL` values
#' will result in plotting the whole segment.
#' 
#' Alternatively, the `region_size` argument can be used to make regional plots
#' in a more automated fashion. Any `dna_seg` features with the value `"TRUE"`
#' in their `region_plot` column will be plotted and become the center of a
#' subsegment. Portions of the `dna_seg` to the left and right of these features
#' will be included and the size of these portions is determined by
#' `region_size`. Any overlapping subsegments will be merged, and features
#' marked as `boundaries` will end subsegments early (e.g. chromosome or contig
#' boundaries). Functions like [edit_dna_segs] can be used to alter which 
#' features to focus on by altering the `region_plot` column, but `dna_segs`
#' can be altered in the same manner as `data.table` objects as well, see the
#' `data.table` package for more details.
#' 
#' `offsets` allows the user to define the placement of the subsegments
#' (as defined by `xlims` or `region_size). If 
#' `offsets` is a list, each element represents 1 `dna_seg` and must be a
#' numeric vector giving gap sizes between the subsegments of that `dna_seg`,
#' including the first one (which will be the distance between the left
#' border of the plot and the first subsegment).
#' Each element of this list must be the same length as the number
#' of subsegments (see `xlims` and details). If `offsets` is a numeric
#' vector, then those numbers will give the distance before the first subsegment
#' for each `dna_seg`, while the others gaps remain static.  If `offsets` is 
#' `NULL`, then the gaps are optimized to minimize `comparison` length.
#' 
#' `dna_seg_line` determines whether a line should be drawn through each
#' `dna_seg`, and if so, what color. If only value is provided, then this
#' value will be repeated for each `dna_seg`. If `dna_seg_line` is a logical 
#' vector, `TRUE` will default to drawing a black line for that `dna_seg`,
#' and `FALSE` will result in no line. If `dna_seg_line` is a character vector,
#' `"FALSE"` will still result in no line, but any other value will be
#' interpreted as a color choice for the line.
#' 
#' `global_color_scheme` applies a color scheme to the `dna_segs` and/or
#' `comparisons`, overriding the colors present in those objects. There are 3
#' options: `"uniform"`,`"gradient"`, and `"sequential"`. These color schemes
#' can be further customized using the other color scheme arguments.
#' 
#' * `"uniform"`: Applies a different color for each possible value of a
#'   given column (specified by `color_scheme_column`) from the `dna_segs` and/
#'   or `comparisons` (specified by `color_scheme_dataset`). Only the values
#'   for features shown in the plot are included. The `color_scheme_column` must
#'   be a column present in all `dna_segs` (`comparisons` are skipped if it is
#'   not present there). If `color_scheme_column = "auto"`, it will determine
#'   which column to use, prioritizing column names related to orthology and 
#'   groups, followed by the `gene` column and finally the `gene_type` column. 
#'   The choice of colors is dependent on the `color_scheme_colors` argument. 
#'   This can be a color palette, or a character vector of colors recognizable 
#'   by R. If `color_scheme_colors = "auto"`, a palette will be chosen based on 
#'   the amount of distinct colors that are required. See [uniform_color_scheme]
#'   for the function that is used to apply this color scheme.
#' * `"gradient"`: Applies a gradient color scheme to the `comparisons` 
#'   based on a numerical column (specified by `color_scheme_column`) that is 
#'   present in all the `comparisons`. The gradient is dependent on 
#'   `color_scheme_colors`, which can be `"red_blue"`, `"blue_red"`, or `"gray"` 
#'   (`"auto"` will default to `"red_blue"`). If `color_scheme_column = "auto"`, 
#'   it will determine which column to use, prioritizing column names that would
#'   be present if the `comparisons` were parsed from BLAST results. The
#'   direction of the gradient
#'   is dependent on `gradient_scheme_direction`, which should be `"increasing"`
#'   for variables that represent a relationship that increases as the numbers
#'   go up (e.g. bit score, alignment length), `"decreasing"` for variables that
#'   represent a relationship that decreases as the numbers go up (e.g. e-value,
#'   gaps, mismatches), or `"auto"`, which will attempt to determine this 
#'   automatically depending on the chosen column. See [gradient_color_scheme]
#'   for the function that is used to apply this color scheme.
#' * `"sequential"`: Transfers over any colors already present in the
#'   `dna_segs` and `comparisons` and copies them over to linked features.
#'   Features are linked through shared identifiers (specified by 
#'   `color_scheme_column`) and `comparisons` that connect them.
#'   See [sequential_color_scheme] for
#'   the function that is used to apply this color scheme.
#' 
#' 
#' @returns A lattice graphic is plotted on the current device. The function
#' itself returns nothing (invisible `NULL`).
#' @export
#' 
#' @param dna_segs A list of `dna_seg` objects. Mandatory.
#' @param comparisons A list of `comparison` objects. If provided, they will
#' plotted between the `dna_segs`. The number of `comparisons` should be 1 less
#' than the number of `dna_segs`.
#' @param tree A tree, in the form of a [phylog] object. If provided, will be
#' plotted on the left. See details.
#' @param tree_width A single numeric, giving the width of the tree area in
#' the plot, in inches. By default it will take 20% of the total plotting area.
#' @param tree_branch_labels_cex A single numeric, giving a size
#' multiplier for possible node annotations of the provided tree (if present).
#' @param tree_scale Logical. If `TRUE`, plots a scale for the tree. 
#' @param legend_column A character string, must be either `"auto"`, or refer to
#' the name of a column present in each `dna_seg`. When provided, will attempt
#' to create a plot legend based on the given column. See details.
#' @param legend_labels A character vector of labels to display in the legend.
#' See details.
#' @param legend_colors A character vector of colors to use for the legend, but
#' only when `legend_labels` is provided. See details.
#' @param annotations An `annotation` or a list of `annotation` objects with 
#' the same length as `dna_segs`. If provided, plots annotations above the 
#' `dna_seg`(s).
#' @param annotation_height A single numeric, giving the height reserved
#' for plotting `annotations`. For comparison, the height of a `dna_seg` is 1.
#' @param annotation_cex A single numeric, giving a size multiplier for the
#' annotations.
#' @param seg_plots A `seg_plot` or a list of `seg_plot` objects with the same
#' length as `dna_segs`. If provided, plots additional data above the 
#' `dna_segs`. See [seg_plot] for more information and some examples.
#' @param seg_plot_height A single numeric, giving the height of the `seg_plot`
#' regions, measured in the unit provided by `seg_plot_height_unit`.
#' @param seg_plot_height_unit The unit of the height of the `seg_plot` regions.
#' Must be a valid unit, see the `grid` documentation for more details. If
#' this argument is set to `"null"`, then the height will be calculated as a
#' proportion of the `comparison` region (i.e. 0.5 means the `seg_plot` region
#' will be half the size of a `comparison`).
#' @param seg_plot_yaxis Can be `NULL`, `FALSE` or a numeric. In the first two
#' cases, no y-axis is drawn for the `seg_plots`. If numeric, an axis is drawn
#' with approximately that number of ticks.
#' @param seg_plot_yaxis_cex A single numeric, giving a size multiplier for the
#' `seg_plot` y-axis. 
#' @param region_size A single numeric or numeric vector with the same length
#' as `dna_segs`, providing the neighbourhood size to use for creating regional
#' plots. Ignored if `xlims` are provided. See details.
#' @param xlims A list with as many elements as there are `dna_segs`, or
#' `NULL`. If `NULL`, the whole DNA segment will be represented. If a list is
#' provided, each element of the list must be a numeric vector, representing 
#' pairs of left and right limits for each subsegment. See details.
#' @param print_xlims Logical. If `TRUE`, prints out the `xlims` (start and
#' end coordinates of each subsegment of each `dna_seg`).
#' @param outfile_xlims A file path. If provided, the `xlims` (start and
#' end coordinates of each subsegment of each `dna_seg`) are written to this
#' file.
#' @param offsets A list or numeric vector with the same length as `dna_segs`,
#' or `NULL`, giving the distance before and between subsegments. Each element
#' of this list must be the same length as the number
#' of subsegments (see `xlims` and details). If `offsets` is `NULL`, then the
#' gaps are optimized to minimize `comparison` length. See details.
#' @param print_offsets Logical. If `TRUE`, prints out the `offsets`, the gap
#' lengths between subsegments for each `dna_seg`.
#' @param minimum_gap_size A single numeric, giving the minimum gap size
#' between subsegments, proportional to the plot region (e.g. `0.03` means the
#' width of the gaps will be at least 3% of the overall plot width).
#' @param fixed_gap_length Logical. If `TRUE`, then the gaps between
#' subsegments will all have the same fixed length instead of optimizing the
#' gap size to minimize `comparison` length.
#' @param limit_to_longest_dna_seg Logical. If `TRUE`, restricts the plot width
#' to the length of the longest `dna_seg`. If `FALSE`, the sizes of the shorter
#' `dna_segs` can be extended to better fit the `comparisons`, but this can
#' lead to extremely wide plots.
#' @param main A character string that gives the main title of the plot.
#' @param main_pos A character string that gives the position of the plot title.
#' Must be one of `"centre"`, `"left"`, or `"right"`.
#' @param dna_seg_labels A character vector with the same length as `dna_segs`,
#' or `NULL`. If `NULL`, the names of the `dna_segs` will be determined
#' automatically where possible (e.g. if `dna_segs` is a named list). Labels
#' are optional, but must be provided or findable if a tree is provided.
#' @param dna_seg_label_cex A single numeric, giving a size multiplier for the
#' `dna_seg` labels.
#' @param dna_seg_label_col A character vector, providing the color(s) for the
#' `dna_seg` labels. Must provide either 1 color or as many as there are 
#' `dna_segs`.
#' @param gene_type A character string, determines the gene type (i.e. shape) 
#' of the features, overriding the `gene_type` column of the `dna_segs`.
#' See [gene_types].
#' @param arrow_head_len A single numeric, giving the length of the arrow heads
#' for the `"arrows"` and `"headless_arrows"` gene types. At maximum, the arrow
#' heads extend to half of the total length. Can be set to `Inf` to force this
#' behavior.
#' @param dna_seg_line A vector, either logical or character, with a length of
#' either 1 or the amount of `dna_segs`. Determines whether a line should be
#' drawn through the middle of each `dna_seg`, and if so, what color. See 
#' details.
#' @param scale Logical. If `TRUE`, a scale will be displayed on the plot.
#' @param dna_seg_scale A single logical that determines whether a scale
#' should be drawn under each `dna_seg`, or a logical vector of the same
#' length as `dna_segs`, for making this choice for each `dna_seg` separately.
#' @param n_scale_ticks A single numeric, giving the approximate number of
#' ticks to display on the longest segment.
#' @param scale_cex A single numeric, giving a size multiplier for the scale
#' labels.
#' @param global_color_scheme A character string, adding a color scheme to
#' the `dna_segs` and/ or `comparisons`. Must be one of: `"uniform"`,
#' `"gradient"`, or `"sequential"`. See details.
#' @param color_scheme_column A character string, must be either `"auto"`, or
#' refer to the name of a column. Depending on `global_color_scheme`, the colors
#' will be determined based on the values found in this column. See details.
#' @param color_scheme_colors A color scheme. If
#' `global_color_scheme = "gradient"`, then it must be one of: `"red_blue"`,
#' `"blue_red"`, or `"gray"`. If `global_color_scheme = "uniform"`, it can be
#' a palette or vector of colors to use. See details.
#' @param color_scheme_dataset A character string, a choice of which data to
#' apply the color scheme to. Must be one of: `"auto"`, `"dna_segs"`, or 
#' `"comparisons"`. If `"auto"`, the color scheme will be applied to both,
#' unless it is not possible to apply it to the `comparisons`. Only applies
#' when `global_color_scheme` is `"uniform"` or `"sequential"`.
#' @param gradient_scheme_direction A character string, indicating the
#' direction of the scale used in the gradient color scheme. Must be one of:
#' `"increasing"`, `"decreasing"`, or `"auto"`. See details.
#' @param alpha_dna_segs A single numeric value between 0 and 1, or `NULL`. 
#' Determines the transparency applied to the `dna_segs`, 0 being fully 
#' transparent, and 1 being fully opaque. This overrides any existing alpha
#' values. If `NULL`, no change is made.
#' @param alpha_comparisons A single numeric value between 0 and 1, or `NULL`. 
#' Determines the transparency applied to the `comparisons`, 0 being fully 
#' transparent, and 1 being fully opaque. This overrides any existing alpha
#' values. If `NULL`, no change is made.
#' @param plot_new Logical. If `TRUE`, uses `grid.newpage()` to produce a new
#' plot. If `FALSE`, integrates it on the current plot.
#' @param outfile A file path. If provided, the plot will be saved to this
#' file instead of the regular output.
#' @param outfile_format A character string, giving the file format for the
#' saved plot `outfile`. Must be one of: `"pdf"`, `"png"`, or `"bmp"`.
#' @param outfile_height A single numeric, giving the height of the plot 
#' `outfile` in inches, or `"auto"`. If `"auto"`, an appropriate height will be
#' approximated automatically.
#' @param outfile_width A single numeric, giving the width of the plot 
#' `outfile` in inches, or `"auto"`. If `"auto"`, an appropriate width will be
#' approximated automatically.
#' @param debug A numeric. If larger than `0`, only that number of elements
#'  will be plotted for each `dna_seg` and `comparison`.
#' @param verbose Logical. If `TRUE`, reports the timing of various steps.
#' @param ... Further arguments to be passed to user-defined graphical
#' functions.
#' 
#' @author Lionel Guy `<lionel.guy@ebc.uu.se>`, Jens Roat Kultima. Mike Puijk
#' 
#' @seealso [dna_seg] and [comparison] for the base objects;
#' [read_dna_seg_from_file], [read_comparison_from_file], and 
#' [read_orthogroup_from_file] to read from files; [annotation] to annotate 
#' `dna_segs`; [seg_plot] to draw plots above `dna_segs`; [gene_types] for 
#' `gene_type` argument; [uniform_color_scheme], [gradient_color_scheme], and
#' [uniform_color_scheme] for color schemes
#' 
#' @examples
#' old.par <- par(no.readonly = TRUE)
#' data("three_genes")
#' dna_segs <- three_genes$dna_segs
#' comparisons <- three_genes$comparisons
#' 
#' ## Segments only
#' plot_gene_map(dna_segs = dna_segs) 
#' 
#' ## With comparisons
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons) 
#' 
#' ## Tree
#' names <- c("A_aaa", "B_bbb", "C_ccc")
#' names(dna_segs) <- names
#' tree <- ade4::newick2phylog("(((A_aaa:4.2,B_bbb:3.9):3.1,C_ccc:7.3):1);")
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               tree = tree)
#' ## Increasing tree width
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               tree = tree, tree_width = 3)
#' ## Annotations on the tree
#' tree2 <- ade4::newick2phylog("(((A_aaa:4.2,B_bbb:3.9)97:3.1,C_ccc:7.3)78:1);")
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               tree = tree2, tree_width = 3)
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               tree = tree2, tree_width = 3, tree_branch_labels_cex = 0.6)
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               tree = tree2, tree_width = 3, tree_branch_labels_cex = 0)
#' 
#' ## Annotation
#' ## Calculating middle positions
#' mid_pos <- middle(dna_segs[[1]])
#' 
#' # Create first annotation
#' annot1 <- annotation(x1 = mid_pos, text = dna_segs[[1]]$name)
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               annotations = annot1)
#' 
#' ## Exploring options
#' annot2 <- annotation(x1 = c(mid_pos[1], dna_segs[[1]]$end[2]),
#'                      x2 = c(NA, dna_segs[[1]]$end[3]),
#'                      text = c(dna_segs[[1]]$name[1], "region1"),
#'                      rot = c(30, 0), col = c("grey", "black"))
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               annotations = annot2, annotation_height = 1.3)
#' 
#' ## xlims
#' ## Just reversing 1 segment
#' plot_gene_map(dna_segs, comparisons,
#'               xlims = list(NULL, NULL, c(Inf,-Inf)),
#'               dna_seg_scale = TRUE)
#' ## Removing one gene
#' plot_gene_map(dna_segs, comparisons,
#'               xlims = list(NULL, NULL, c(-Inf,2800)),
#'               dna_seg_scale = TRUE)
#' 
#' ## offsets
#' offsets <- c(0, 0, 0)  
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               offsets = offsets)
#' offsets <- c(200, 400, 0)  
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               offsets = offsets)
#' 
#' ## main
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               main = "Comparison of A, B and C")
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               main = "Comparison of A, B and C", main_pos = "left")
#' 
#' ## dna_seg_labels
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               dna_seg_labels = c("Huey", "Dewey", "Louie"))
#' 
#' ## dna_seg_labels size
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               dna_seg_labels = c("Huey", "Dewey", "Louie"),
#'               dna_seg_label_cex = 2)
#' 
#' ## dna_seg_line
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               dna_seg_line = c("FALSE", "red", grey(0.6)))
#' 
#' ## gene_type
#' plot_gene_map(dna_segs = dna_segs, comparisons = comparisons,
#'               gene_type = "side_blocks")
#' 
#' ##
#' ## From here on, using a bigger dataset from a 4-genome comparison
#' ##
#' data("barto")
#' ## Adding a tree
#' tree <- ade4::newick2phylog("(BB:2.5,(BG:1.8,(BH:1,BQ:0.8):1.9):3);")
#' ## Showing only subsegments
#' xlims1 <- list(c(1380000, 1445000),
#'                c(10000, 83000),
#'                c(15000, 98000),
#'                c(5000, 82000))
#' ## Reducing dataset size for speed purpose
#' for (i in 1:length(barto$dna_segs)) {
#'   barto$dna_segs[[i]] <- trim(barto$dna_segs[[i]], xlim = xlims1[[i]])
#'   if (i < length(barto$dna_segs)) {
#'     barto$comparisons[[i]] <- trim(barto$comparisons[[i]],
#'                                    xlim1 = xlims1[[i]], xlims1[[i+1]])
#'   }
#' }
#' plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree,
#'               xlims = xlims1,
#'               dna_seg_scale = TRUE)
#' ## Showing several subsegments per genome
#' xlims2 <- list(c(1445000, 1415000, 1380000, 1412000),
#'                c(  10000,   45000,   50000,   83000, 90000, 120000),
#'                c(  15000,   36000,   90000,  120000, 74000,  98000),
#'                c(   5000,   82000))
#' plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree,
#'               xlims = xlims2,
#'               dna_seg_scale = TRUE)
#' ## Hand-made offsets: size of all gaps
#' offsets2 <- list(c(10000, 10000),
#'                  c(2000, 2000, 2000),
#'                  c(10000, 5000, 2000),
#'                  c(10000))
#' plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree,
#'               xlims = xlims2,
#'               offsets = offsets2,
#'               dna_seg_scale = TRUE)
#' 
#' ## dna_seg_scale, global_color_scheme, size, number, color of dna_seg_scale,
#' ## size of dna_seg_scale labels
#' plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree,
#'               xlims = xlims2,
#'               dna_seg_scale = c(TRUE, FALSE, FALSE, TRUE),
#'               scale = FALSE,
#'               dna_seg_label_cex = 1.4,
#'               dna_seg_label_col = c("black", "grey", "blue", "red"),
#'               global_color_scheme = "gradient",
#'               alpha_comparisons = 0.5,
#'               n_scale_ticks = 3, scale_cex = 1)
#' 
#' ##
#' ## Exploring and modifying a previously plotted gene map plot
#' ##
#' plot_gene_map(barto$dna_segs, barto$comparisons, tree = tree,
#'               xlims = xlims2, offsets = offsets2, dna_seg_scale = TRUE)
#' ## View viewports
#' current.vpTree()
#' ## Go down to one of the viewports, add an xaxis, go back up to root viewport
#' downViewport("dna_seg_scale.3.2")
#' grid.rect()
#' upViewport(0)
#' ## Get all the names of the objects
#' grobNames <- getNames()
#' grobNames
#' ## Change the color of the scale line
#' grid.edit("scale.lines", gp = gpar(col = "grey"))
#' ## Remove first dna_seg_lines
#' grid.remove("dna_seg_line.1.1")
#' 
#' ##
#' ## Plot genoPlotR logo
#' ##
#' col_vec <- c("#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
#'              "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC")
#' cex <- 2.3
#' ## First segment 
#' start1 <- c(150, 390, 570, 270, 530)
#' end1   <- c(  1, 490, 690, 270, 530)
#' ## Second segment
#' start2 <- c(100, 520, 550, 330)
#' end2   <- c(240, 420, 650, 330)
#' ## dna_segs
#' ds1 <- as.dna_seg(data.frame(name = c("", "", "", "geno", "R"),
#'                              start = start1, end = end1, strand = rep(1, 5),
#'                              fill = col_vec[c(2, 6, 1, 8, 9)]
#'                              ))
#' ds_genoR <- edit_dna_segs(ds1, ids = data.frame(id = c("geno", "R"), 
#'                                                 cex = c(2.3, 2.3),
#'                                                 gene_type = c("text", "text")
#'                                                 ))
#' ds2 <- as.dna_seg(data.frame(name = c("", "", "", "Plot"),
#'                              start = start2, end = end2,
#'                              strand = rep(1, 4),
#'                              fill = col_vec[c(5, 3, 7, 1)]
#'                              ))
#' ds_Plot <- edit_dna_segs(ds2, ids = data.frame(id = "Plot", 
#'                                                cex = 2.3,
#'                                                gene_type = "text"
#'                                                ))
#' ## comparison
#' c1 <- as.comparison(data.frame(start1 = start1[1:3], end1 = end1[1:3],
#'                                start2 = start2[1:3], end2 = end2[1:3],
#'                                fill = grey(c(0.6, 0.8, 0.5))))
#' ## Generate genoPlotR logo
#' \dontrun{
#'   pdf("logo.pdf", h = 0.7, w = 3)
#' }
#' par(fin = c(0.7, 3))
#' plot_gene_map(dna_segs = list(ds_genoR, ds_Plot),
#'               comparisons = list(c1), scale = FALSE, dna_seg_scale = FALSE,
#'               dna_seg_line = grey(0.7), offsets = c(-20,160))
#' \dontrun{
#'   dev.off()
#' }
#' par(old.par)
#' 
plot_gene_map <- function(
  dna_segs,           # list of dna_segs, mandatory
  comparisons = NULL, # list of comparisons between dna_segs, optional
  tree = NULL,        # tree to plot, phylog object, optional
  tree_width = NULL,  # tree width in inches
  tree_branch_labels_cex = NULL, # relative size of tree node annotations
  tree_scale = FALSE, # plot tree scale?
  legend_column = NULL, # dna_seg/ comparison column to base legend on
  legend_labels = NULL, # vector of labels to include in legend
  legend_colors = NULL, # a vector of colors to be used for legend_labels
  annotations = NULL,      # list of dna_seg annotations, optional
  annotation_height = 0.1, # height of annotation line
  annotation_cex = 0.8,    # relative size of annotations
  seg_plots = NULL,               # user-defined plots
  seg_plot_height = 3,            # height of seg_plots
  seg_plot_height_unit = "lines", # unit of preceding height variable
  seg_plot_yaxis = 3,             # if non-NULL or non false, number of ticks
  seg_plot_yaxis_cex = scale_cex, # relative size of seg_plot yaxis
  region_size = NULL,   # region size for regional plots, sets xlims
  xlims = NULL,         # plot parts of dna_segs, ignores region_size argument
  print_xlims = FALSE,  # prints the xlims out to terminal
  outfile_xlims = NULL, # file path to write xlims out to for later use
  offsets = NULL,           # regulates manual gap alignment between dna_segs
  print_offsets = FALSE,    # prints the offsets to terminal
  minimum_gap_size = 0.03,  # minimum gap between dna_segs (% of entire plot)
  fixed_gap_length = FALSE, # make gap length uniform?
  limit_to_longest_dna_seg = TRUE, # limit gap sizes of longest dna_seg?
  main = NULL,           # main title
  main_pos = "centre",   # centre, left, right
  dna_seg_labels = NULL,       # labels on dna_segs
  dna_seg_label_cex = 0.9,     # relative size of dna_seg labels
  dna_seg_label_col = "black", # color of dna_seg labels
  gene_type = NULL,      # if not NULL, resets gene_type (feature shape)
  arrow_head_len = 200,  # length of arrow heads for arrow gene types
  dna_seg_line = TRUE,   # draw a line on each dna_seg? or its color
  scale = TRUE,          # draw a scale in the bottom right?
  dna_seg_scale = FALSE, # draw a scale on each dna_seg?
  n_scale_ticks = 7,     # approximate number of tick marks for dna_seg_scale
  scale_cex = 0.6,       # relative size of text on scale
  global_color_scheme = NULL,         # choice of color scheme to apply
  color_scheme_column = "auto",       # column to use as basis for color scheme
  color_scheme_colors = "auto",       # colors for color scheme
  color_scheme_dataset = "auto",      # data to apply color scheme to 
  gradient_scheme_direction = "auto", # direction of the gradient color scheme
  alpha_dna_segs = NULL,    # dna_segs transparency, number between 0 and 1
  alpha_comparisons = NULL, # comparisons transparency, number between 0 and 1
  plot_new = TRUE, # FALSE to integrate on a bigger plot
  outfile = NULL,          # file path to save plot to 
  outfile_format = "pdf",  # file type for saved plots (pdf, png, bmp)
  outfile_height = "auto", # height of saved plot file in inches
  outfile_width = "auto",  # width of saved plot file in inches
  debug = 0,
  verbose = FALSE, # reports timings of various steps
  ...
) {
  
  #----------------------------------------------------------------------------#
  # check arguments
  #----------------------------------------------------------------------------#
  # objects #
  # dna_segs
  if (verbose) cat_time("Starting plot_gene_map.")
  if (missing(dna_segs)) stop("'dna_segs' must be provided")
  if (!is.list(dna_segs) || !all(sapply(dna_segs, is.dna_seg))) {
    stop("'dna_segs' must be a list of dna_seg objects")
  }
  n_dna_segs <- length(dna_segs)
  n_rows <- 3*n_dna_segs-1
  
  # comparisons
  n_comparisons <- length(comparisons)
  if (n_comparisons > 0) {
    if (!is.list(comparisons) || !all(sapply(comparisons, is.comparison))) {
      stop("'comparisons' must be a list of comparison objects")
    }
    # check that there are enough comparisons compared to dna segments
    if (!(n_dna_segs - n_comparisons == 1)) {
      stop("The amount of comparison objects must be 1 less than ", 
           "the amount of dna_seg objects")
    }
  }
  
  # check seg labels and assign them if they were left as NULL
  seg_labels <- get_seg_labels(dna_segs, dna_seg_labels)
  
  # check seg_label colors
  if (length(dna_seg_label_col) == 1) {
    dna_seg_label_col <- rep(dna_seg_label_col, n_dna_segs)
  }
  else if (!length(dna_seg_label_col) == n_dna_segs) {
    stop("'dna_seg_label_col' must be either a single color or ",
         "a vector of colors with the same length as 'dna_segs'")
  }
  
  # check tree
  if (!is.null(tree)) {
    if (!inherits(tree, "phylog")) stop("'tree' must be a phylog object")
    # check correspondence between labels (given by names(dna_segs) or
    # seg_labels) and tree leaves
    if (is.null(seg_labels)) {
      stop("If 'tree' is provided, label names must be provided via named ",
           "list dna_segs or 'dna_seg_labels'")
    }
    # check that number of leaves corresponds to number of segs
    if (length(tree$leaves) != n_dna_segs) {
      stop("Number of leaves in 'tree' not equal to the number of dna_segs")
    }
    if (!all(gsub(".", "_", seg_labels, fixed = TRUE) %in% names(tree$leaves))) {
      stop("The tree leaf labels must correspond to the dna_seg labels")
    }
    # check whether nodes have added labels
    if (is.null(tree_branch_labels_cex)) {
      if (length(grep("^[^I]", names(tree$nodes)))) {
        tree_branch_labels_cex <- 0.8
      } else {
        tree_branch_labels_cex <- 0
      }
    }
  }

  # check seg_plots (user-defined plots)
  if (!is.null(seg_plots)) {
    # if there is only one seg_plot, put it on the top
    if (is.seg_plot(seg_plots)) {
      s_plot <- seg_plots
      seg_plots <- c(list(s_plot), rep(list(NULL), n_dna_segs - 1))
    } else if (length(seg_plots) == n_dna_segs) {
      if (!all(sapply(seg_plots, function(x) is.seg_plot(x) || is.null(x)))) {
        stop("All elements of 'seg_plots' must be NULL or seg_plot objects")
      }
    } else {
      stop ("'seg_plots' must be of same length as 'dna_segs'")
    }
    seg_plot_h <- ifelse(sapply(seg_plots, is.null), 0, seg_plot_height)
  } else {
    seg_plot_h <- rep(0, n_dna_segs)
  }
  # check seg_plot_yaxis
  if (is.null(seg_plot_yaxis) || !is.numeric(seg_plot_yaxis) ||
      is.null(seg_plots)
      ) {
    seg_plot_yaxis <- 0
  }
  
  # check annotation
  if (!is.null(annotations)) {
    # if there is only one annotation, put it on the top
    if (is.annotation(annotations)) {
      annot <- annotations
      annotations <- c(list(annot), rep(list(NULL), n_dna_segs-1))
    } else if (length(annotations) == n_dna_segs) {
      if (!all(sapply(annotations,
                      function(x) is.annotation(x) || is.null(x)
                      ))
          ) {
        stop("All elements of 'annotations' must be NULL or annotation objects")
      }
    } else {
      stop("'annotation' must be of same length as 'dna_segs'")
    }
    annot_h <- ifelse(sapply(annotations, is.null), 0, annotation_height)
  } else {
    annot_h <- rep(0, n_dna_segs)
  }
  
  if (verbose) cat_time("Determining xlims.")

  # graphical arguments #
  # check xlims
  if (!is.null(xlims)) {
    # User-supplied xlims 
    if (length(xlims) != n_dna_segs && length(xlims) != 1) {
      stop("'xlims' must be of length 1 or the same length as 'dna_segs'")
    }
    if (!is.list(xlims) || 
        !all(sapply(xlims,function(x) (length(x) %% 2) == 0))
        ) {
      stop("'xlims' must be a list whose elements all contain ",
           "an even number of elements")
    }
    if (length(xlims) == 1) {
      # Only 1 set of xlims was provided, so duplicate it for all dna_segs
      xlims <- lapply(seq(1, n_dna_segs), function(x) xlims[[x]] <- xlims[[1]])
    }
    for (i in 1:length(xlims)) {
      xlim <- xlims[[i]]
      # check numeric, replace null and inf values
      rng <- range(dna_segs[[i]])
      seg_min <- rng[1] - 0.01*diff(rng)
      seg_max <- rng[2] + 0.01*diff(rng)
      # replace NULL value by min-max
      if (is.null(xlim)) xlim <- c(seg_min, seg_max)
      # replace -Inf and Inf by min max
      xlim[xlim == -Inf] <- seg_min
      xlim[xlim == Inf] <- seg_max
      # check that all are numeric
      if (!is.numeric(xlim)) stop("All elements of 'xlims' must be numeric")
      # transform to data.frame
      xlim <- data.frame(matrix(xlim, ncol = 2, byrow = TRUE))
      names(xlim) <- c("x0", "x1")
      # check the strand
      xlim$strand <- ifelse(xlim$x0 < xlim$x1, 1, -1)
      # sort x0 x1
      for (j in 1:nrow(xlim)) xlim[j,1:2] <- sort(as.numeric(xlim[j,1:2]))
      xlims[[i]] <- xlim
    }
  } else if (!is.null(region_size)) {
    region_size <- unlist(region_size)
    if (!is.numeric(region_size) |
        (length(region_size) != 1) & length(region_size) != n_dna_segs
        ) {
      stop("'region_size' must be a numeric vector or list of length 1 or ",
           "the same length as 'dna_segs'")
    }
    if (length(region_size) == 1) {
      # A single region size was provided, so duplicate it for all dna_segs
      region_size <- rep(unlist(region_size), n_dna_segs)
    }
    # calculate xlims of gene regions
    if (!all(sapply(dna_segs, function(x) ("region_plot" %in% names(x) &&
                                           any(x[, region_plot == "TRUE"]))))
        ) {
      stop("If 'xlims' specifies a neighbourhood size, then each dna_seg ",
           "must have at least 1 feature to plot, specified by the ",
           "'region_plot' column"
           )
    }
    xlims <- list()
    for (i in 1:n_dna_segs) {
      to_plot <- dna_segs[[i]][region_plot %ilike% "^TRUE$",
                               .(start, end)]
      bounds <- dcast(dna_segs[[i]][region_plot %like% "start|end",
                                    .(seq_origin, start, region_plot)
                                    ],
                      seq_origin ~ region_plot,
                      value.var = "start"
                      )
      if (nrow(bounds) > 0) {
        # Find and then add sequence boundaries for each feature to plot
        to_plot <- cbind(
          to_plot,
          bounds[
            sapply(1:nrow(to_plot), function(x) which(
            to_plot$start[x] >= bounds$start & to_plot$end[x] <= bounds$end
          )), .(seqStart = start, seqEnd = end)]
        )
        # Using sequence boundaries, determine xlims
        to_plot[, let(
          xlimstart = (max((.SD$start - region_size[i]), .SD$seqStart ) - 5), 
          xlimend = (min((.SD$end + region_size[[i]]), .SD$seqEnd) + 5)
          ),
          by = seq_len(nrow(to_plot))
        ]
      } else { # No sequence boundaries present
        rng <- range(dna_segs[[i]])
        seg_min <- rng[1] - 0.02*diff(rng)
        seg_max <- rng[2] + 0.02*diff(rng)
        to_plot[, let(
          xlimstart = max((.SD$start - region_size[i] - 5), seg_min), 
          xlimend = min((.SD$end + region_size[i] + 5), seg_max)), 
          by = seq_len(nrow(to_plot))
        ]
      }
      # Merge overlapping neighbourhoods, add strand column
      xlim <- data.frame(
        to_plot[
          ,
          .(x0 = min(xlimstart), x1 = max(xlimend)),
          by = .(group = cumsum(c(1, tail(xlimstart, - 1) > head(xlimend, -1))))
        ][
          ,
          .(x0, x1, strand = fifelse(x0<x1, 1, -1))
        ]
      )
      for (j in 1:nrow(xlim)) xlim[j, 1:2] <- sort(as.numeric(xlim[j, 1:2]))
      xlims[[i]] <- xlim
    }
  } else {
    xlims <- list()
    for (i in 1:n_dna_segs) {
      rng <- range(dna_segs[[i]])
      xlims[[i]] <- data.frame(x0 = rng[1] - 0.01 * diff(rng),
                               x1 = rng[2] + 0.01 * diff(rng),
                               strand = 1
                               )
    }
  }
  
  if (!is.null(outfile_xlims)) {
    xlim_labels <- get_seg_labels(dna_segs, seg_labels, required = TRUE)
    named_xlims <- copy(xlims)
    for (i in 1:n_dna_segs) {
      named_xlims[[i]] <- as.data.table(named_xlims[[i]])
      named_xlims[[i]][strand == -1, c("x0", "x1") := .(x1, x0)]
      named_xlims[[i]]$strand <- NULL
      named_xlims[[i]]$seg_label <- xlim_labels[i]
    }
    named_xlims <- rbindlist(named_xlims)
    fwrite(named_xlims, file = outfile_xlims, sep = "\t", eol = "\n")
  }
  
  if (print_xlims) {
    formatted_xlims <- character(n_dna_segs)
    for (i in 1:n_dna_segs) {
      temp_xlims <- as.data.table(copy(xlims[[i]]))
      temp_xlims[strand == -1, c("x0", "x1") := .(x1, x0)]
      formatted_xlims[[i]] <- paste(unlist(temp_xlims[, 1], use.names = FALSE),
                                    unlist(temp_xlims[, 2], use.names = FALSE),
                                    collapse = " "
                                    )
    }
    cat_time("xlims:\n",
             paste0(formatted_xlims, collapse = ","),
             "\n",
             sep = ""
             )
  }
  
  if (verbose) cat_time("Finished determining xlims.")
  
  # check offsets
  if (!is.null(offsets) && length(offsets) != n_dna_segs) {
    stop("'offsets' must be of the same length as 'dna_segs'")
  }

  # check main_pos
  if (main_pos == "centre") {
    main_x <- 0.5
    main_just <- "centre"
  } else if (main_pos == "left") {
    main_x <- 0
    main_just <- "left"
  } else if (main_pos == "right") {
    main_x <- 1
    main_just <- "right"
  } else {
    stop("'main_pos' must be one of: centre, left, right")
  }

  # check dna_seg_line
  if (is.logical(dna_seg_line)) {
    dna_seg_line <- as.character(dna_seg_line)
    dna_seg_line[dna_seg_line == "TRUE"] <- "black"
  }
  if (!is.character(dna_seg_line)) {
    stop("'dna_seg_line' must be either logical, or ",
         "a character string providing a color")
  }
  if (length(dna_seg_line) == 1) {
    dna_seg_line <- rep(dna_seg_line, n_dna_segs)
  } else if (length(dna_seg_line) != n_dna_segs) {
    stop("'dna_seg_line' must be of length 1 or the same length as 'dna_segs'")
  }
  
  # check gene_type
  if (!is.null(gene_type) && !(gene_type %in% gene_types())) {
    stop("'gene_type' must be one of: ", paste(gene_types(), collapse = ", "))
  }
  
  # check dna_seg_scale. Must be logical
  if (is.logical(dna_seg_scale)) {
    # if length 1, make dna_seg_scale the same length as dna_segs
    if (length(dna_seg_scale) == 1) {
      dna_seg_scale <- rep(dna_seg_scale, n_dna_segs)
    }
    # if a different length, must be of the same length as n_dna_segs
    else if (length(dna_seg_scale) != n_dna_segs) {
      stop("'dna_seg_scale' must be of length 1 or",
           "the same length as 'dna_segs'")
    }
  } else {
    stop("'dna_seg_scale' must be logical")
  }
  
  # check alpha_dna_segs
  if (!is.null(alpha_dna_segs)) {
    if (suppressWarnings(is.na(as.double(alpha_dna_segs)))) {
      stop("'alpha_dna_segs' must be a number between 0 and 1")
    }
    alpha_dna_segs <- as.double(alpha_dna_segs)
    if (alpha_dna_segs > 1 | alpha_dna_segs < 0) {
      stop("'alpha_dna_segs' must be a number between 0 and 1")
    }
  }
  
  
  # check alpha_comparisons
  if (!is.null(alpha_comparisons)) {
    if (suppressWarnings(is.na(as.double(alpha_comparisons)))) {
      stop("'alpha_comparisons' must be a number between 0 and 1")
    }
    alpha_comparisons <- as.double(alpha_comparisons)
    if (alpha_comparisons > 1 | alpha_comparisons < 0) {
      stop("'alpha_comparisons' must be a number between 0 and 1")
    }
  }
  
  #----------------------------------------------------------------------------#
  # plotting options
  #----------------------------------------------------------------------------#
  # dna_seg lines height. 1 line for the dna_seg, 0.5 in addition for
  # the dna_seg_scale, if needed. 1 null in between for comparisons
  # to rewrite sometimes
  h <- rep(1, n_rows)
  h[seq(2, n_rows, by = 3)] <- 1 + scale_cex*dna_seg_scale + annot_h
  h[seq(1, n_rows, by = 3)] <- seg_plot_h
  dna_seg_heights <- unit(
    h,
    c(rep(c(seg_plot_height_unit, "lines", "null"), n_dna_segs), "lines")
  )
  

  
  #----------------------------------------------------------------------------#
  # prepare plotting frame & filter objects 
  #----------------------------------------------------------------------------#
  # further calculations on xlims & gaps
  for (i in 1:n_dna_segs) {
    xlims[[i]]$length <- xlims[[i]]$x1 - xlims[[i]]$x0
  }
  # default gap_length is a 20th of the max length
  def_gap_length <- max(sapply(xlims, function(x) sum(x$length))
                        ) * minimum_gap_size
  unpadded_lengths <- sapply(xlims, function(x)
                             sum(x$length) + (nrow(x) - 1) * def_gap_length)
  longest_seg <- which.max(unpadded_lengths)
  max_length <- unpadded_lengths[longest_seg]
  scale_unit <- diff(pretty(c(0, max_length), n = n_scale_ticks + 2)[1:2])

  # trim comparisons #
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      # concatenate all possible combinations
      comp1 <- comparisons[[i]][0,]
      for (j in 1:nrow(xlims[[i]])) {
        for (k in 1:nrow(xlims[[i+1]])) {
          comp1 <- rbindlist(list(comp1, trim.comparison(
            comparisons[[i]],
            xlim1 = as.numeric(xlims[[i]][j, c("x0", "x1")]),
            xlim2 = as.numeric(xlims[[i + 1]][k, c("x0", "x1")])
          )))
        }
      }
      # rbind from data.table makes its arguments lose their original class
      class(comp1) <- c("comparison", class(comp1))
      comparisons[[i]] <- comp1
    }
  }
  
  # check global_color_scheme
  if (!is.null(global_color_scheme)) {
    if (length(global_color_scheme) == 4) {
      # For backwards compatibility with the old version of this argument
      if (verbose) {
        warning("'global_color_scheme' has a length of 4, where 1 was ",
                'expected, so it is treated as "gradient", using these 4',
                'elements as the other arguments for this color scheme.')
      }
      color_scheme_column <- global_color_scheme[1]
      gradient_scheme_direction <- global_color_scheme[2]
      color_scheme_colors <- global_color_scheme[3]
      alpha_comparisons <- as.numeric(global_color_scheme[4])
      global_color_scheme <- "gradient"
    }
    if (!any(global_color_scheme == c("gradient", "uniform", "sequential"))) {
      stop("'global_color_scheme' must be one of: ",
           "gradient, uniform, sequential")
    }
    
    if (global_color_scheme == "gradient" & !n_comparisons > 0) {
      stop("When 'global_color_scheme' is ",
           '"gradient", comparisons must be provided')
    }
    
    # color_scheme_column: choice of column
    if (global_color_scheme == "gradient") {
      if (color_scheme_column == "auto") {
        names_1 <- names(comparisons[[1]])
        # collect numerical columns
        num_cols <- lapply(comparisons,
                           function(x) names(x)[sapply(x, is.numeric)]
                           )
        shared_num_cols <- names(
          which(table(unlist(num_cols)) == length(num_cols))
        )
        shared_num_cols <- shared_num_cols[
          !shared_num_cols %in% c("start1", "start2", "end1", "end2")
        ]
        # take if present, per_id, bit_score, or e_value, in that order
        # else take first shared numerical column from first comparison
        relevant_colnames <- c("per_id", "bit_score", "e_value")
        relevant_cols <- relevant_colnames %in% tolower(shared_num_cols)
        if (any(relevant_cols)) {
          color_scheme_column <- relevant_colnames[which(relevant_cols)[1]]
          color_scheme_column <- shared_num_cols[
            which(tolower(shared_num_cols) == color_scheme_column)
          ]
        } else {
          color_scheme_column <- names_1[names_1 %in% shared_num_cols][1]
        }
      } else if (!all(sapply(comparisons,
                             function(x) color_scheme_column %in% names(x)))
                 ) {
        stop("'color_scheme_column' must be either \"auto\" or a column name. ",
             "A column name was given, but not all comparison objects ",
             "contained the given column name")
      }
    } else if (global_color_scheme == "uniform") {
      # deal with color_scheme_dataset first: which data are we modifying
      if (!any(color_scheme_dataset == c("auto", "dna_segs", "comparisons"))) {
        stop("'color_scheme_dataset' must be one of: auto, dna_segs, ",
             "comparisons")
      }
      if (color_scheme_dataset == "comparisons" & !n_comparisons > 0) {
        stop("Argument 'color_scheme_dataset' was \"comparisons\", but no ",
             "comparisons were found")
      }
      if (color_scheme_column == "auto") {
        if (color_scheme_dataset == "comparisons") {
          names_1 <- names(comparisons[[1]])
          cols <- lapply(comparisons, function(x) names(x))
        } else {
          names_1 <- names(dna_segs[[1]])
          cols <- lapply(dna_segs, function(x) names(x))
        }
        shared_cols <- names(which(table(unlist(cols)) == length(cols)))
        shared_cols <- shared_cols[!shared_cols %in% c("start", "end", 
                                                       "strand", "length")]
        # take if present, anything related to groups, clusters, or homology
        relevant_colnames <- c("orthogroup", "cog", "kog", "og", "orthology", 
          "group", "homolog", "ortholog", "cluster", "gene", "gene_type")
        relevant_cols <- relevant_colnames %in% tolower(shared_cols)
        if (any(relevant_cols)) {
          color_scheme_column <- relevant_colnames[which(relevant_cols)[1]]
          color_scheme_column <- shared_cols[
            which(tolower(shared_cols) == color_scheme_column)
          ]
        } else {
          color_scheme_column <- names_1[names_1 %in% shared_cols][1]
        }
        if (color_scheme_dataset == "comparisons") {
          # Currently does not matter, here to make sure these objects exist
          column_in_comps <- TRUE 
          column_in_segs <- FALSE
        } else if (color_scheme_dataset == "dna_segs") {
          # Currently does not matter, here to make sure these objects exist
          column_in_segs <- TRUE
          column_in_comps <- FALSE
        } else {
          column_in_segs <- TRUE
          column_in_comps <- all(sapply(
            comparisons,
            function(x) color_scheme_column %in% names(x)
          ))
        }
      } else {
        # Check to see if the chosen column is present
        if (color_scheme_dataset == "comparisons" & 
            !all(sapply(comparisons,
                        function(x) color_scheme_column %in% names(x)
                        ))
            ) {
          stop("'color_scheme_column' must be either \"auto\" or a column ",
               "name. A column name was given, but not all comparison ",
               "objects contained the given column name"
               )
        } else if (color_scheme_dataset == "dna_segs" &
                   !all(sapply(dna_segs,
                               function(x) color_scheme_column %in% names(x)
                               ))
                   ) {
          stop("'color_scheme_column' must be either \"auto\" or a column ",
               "name. A column name was given, but not all dna_seg ",
               "objects contained the given column name"
               )
        } else if (color_scheme_dataset == "auto") {
          if (n_comparisons > 0) {
            column_in_comps <- all(sapply(
              comparisons,
              function(x) color_scheme_column %in% names(x)
            ))
          } else {
            column_in_comps <- FALSE
            color_scheme_dataset <- "dna_segs"
          }
          column_in_segs <- all(sapply(
            dna_segs,
            function(x) color_scheme_column %in% names(x)
          ))
          if (!column_in_comps & !column_in_segs) {
            stop("'color_scheme_column' must be either \"auto\" or a column ",
                 "name. A column name was given, but not all comparison ",
                 "and/or dna_seg objects contained the given column name"
                 )
          }
        }
      }
    } else if (global_color_scheme == "sequential") {
      # deal with color_scheme_dataset first: which data are we modifying
      if (!any(color_scheme_dataset == c("auto", "dna_segs", "comparisons"))) {
        stop("'color_scheme_dataset' must be one of: auto, dna_segs, ",
             "comparisons")
      }
      if (color_scheme_dataset == "comparisons" & !n_comparisons > 0) {
        stop("Argument 'color_scheme_dataset' was \"comparisons\", but no ",
             "comparisons were found")
      }
      if (color_scheme_dataset == "auto" & !n_comparisons > 0) {
        color_scheme_dataset <- "dna_segs"
      }
      if (color_scheme_column == "auto") {
        # collect shared dna_seg columns 
        names_1 <- names(dna_segs[[1]])
        cols <- lapply(dna_segs, function(x) names(x))
        shared_cols <- names(which(table(unlist(cols)) == length(cols)))
        shared_cols <- shared_cols[!shared_cols %in% c("start", "end", 
                                                       "strand", "length")]
        relevant_colnames <- c("locus_id", "name")
        relevant_cols <- relevant_colnames %in% tolower(shared_cols)
        if (any(relevant_cols)) {
          color_scheme_column <- relevant_colnames[which(relevant_cols)[1]]
          color_scheme_column <- shared_cols[
            which(tolower(shared_cols) == color_scheme_column)
          ]
        } else {
          # Probably cannot happen as long as name is in relevant_cols
          color_scheme_column <- names_1[names_1 %in% shared_cols][1]
        }
      } else if (!all(sapply(dna_segs,
                             function(x) color_scheme_column %in% names(x)
                             ))
                 ) {
        stop("'color_scheme_column' must be either \"auto\" or a column ",
             "name. A column name was given, but not all dna_seg ",
             "objects contained the given column name"
             )
      }
    }

    # color_scheme_colors: choice of palette
    if (global_color_scheme == "gradient") {
      if (color_scheme_colors == "auto") {
        color_scheme_colors <- "red_blue"
      } else if (!any(color_scheme_colors == c("grey", "gray", "grays", "greys",
                                               "red_blue", "blue_red", "auto"
                                               ))
                 ) {
        stop("If 'global_color_scheme' is \"gradient\", then ",
             "'color_scheme_colors' must be one of: red_blue, greys, auto")
      }
    } else if (global_color_scheme == "uniform") {
      if (color_scheme_colors == "auto") {
        uniform_palette <- NULL
      } else {
        uniform_palette <- color_scheme_colors
      }
    }
    
    if (global_color_scheme == "gradient") {
      # gradient_scheme_direction: increasing or decreasing
      col_args_5 <- c("increasing", "decreasing", "auto")
      if (length(grep(gradient_scheme_direction, col_args_5)) != 1) {
        stop("'gradient_scheme_direction' must be one of: ",
             paste(col_args_5, collapse = ", "))
      } else {
        gradient_scheme_direction <- grep(gradient_scheme_direction,
                                          col_args_5,
                                          value = TRUE
                                          )
      }
      # gradient_scheme_direction: incr/decr. Turn it to decreasing
      if (gradient_scheme_direction == "auto") {
        gradient_scheme_direction <- FALSE
        if (color_scheme_column %in% c("mism", "gaps", "e_value")) {
          gradient_scheme_direction <- TRUE
        }
      } else if (gradient_scheme_direction == "decreasing") {
        gradient_scheme_direction <- TRUE
      } else if (gradient_scheme_direction == "increasing") {
        gradient_scheme_direction <- FALSE
      }
      # gather range of values from all comparisons
      range_col_from <- c(Inf, -Inf)
      for (i in 1:n_comparisons) {
        if (nrow(comparisons[[i]]) > 0) {
          range_col_from[1] <- min(c(range_col_from[1],
                                     comparisons[[i]][[color_scheme_column]]
                                     ))
          range_col_from[2] <- max(c(range_col_from[2],
                                     comparisons[[i]][[color_scheme_column]]
                                     ))
        }
      }
      # perform gradient_color_scheme
      for (i in 1:n_comparisons) {
        comparisons[[i]]$col <-
          gradient_color_scheme(x = comparisons[[i]][[color_scheme_column]],
                                direction = comparisons[[i]]$direction,
                                color_scheme = color_scheme_colors,
                                decreasing = gradient_scheme_direction,
                                rng = range_col_from,
                                alpha = alpha_comparisons
                                )
      }
    } else if (global_color_scheme == "uniform") {
      # gather all unique values that need to be colored 
      comp_ids <- character()
      seg_ids <- character()
      if ((color_scheme_dataset == "auto" & column_in_comps) | 
          color_scheme_dataset == "comparisons"
          ) {
        # comparisons were already trimmed, so gather all values
        comp_ids <- unlist(sapply(
          comparisons,
          function(x) x[, c(color_scheme_column), with = FALSE]
        ))
        comp_ids <- comp_ids[comp_ids != "NA"]
        id_count <- summary(as.factor(comp_ids), maxsum = 1000000)
        comp_ids <- unique(comp_ids)
      }
      if ((color_scheme_dataset == "auto" & column_in_segs) | 
          color_scheme_dataset == "dna_segs"
          ) {
        # use xlims to determine which values need to be taken from dna_segs
        for (i in 1:length(xlims)) {
          n_subsegs <- nrow(xlims[[i]])
          for (j in 1:n_subsegs) {
            seg_ids <- c(
              seg_ids,
              unlist(dna_segs[[i]][
                start >= xlims[[i]][j, "x0"] & end <= xlims[[i]][j, "x1"],
              ][, color_scheme_column, with = FALSE], use.names = FALSE)
            )
          }
        }
        seg_ids <- seg_ids[seg_ids != "NA"]
        id_count <- summary(as.factor(seg_ids), maxsum = 1000000)
        seg_ids <- unique(seg_ids)
      }
      id_vector <- c(comp_ids, seg_ids)
      id_vector <- sort(unique(id_vector))
      if (length(id_vector[id_vector != "-"]) > 34) {
        # Get rid of IDs that appear only once, if there are too many IDs
        id_count <- id_count[id_count != 1 | names(id_count) == "-"]
        if (length(id_count[names(id_count) != "-"]) >= 10) {
          # If after filtering one-offs 10 or less are left, don't filter at all
          # Otherwise, continue with possibly more filtering
          id_vector <- names(id_count)
          if (length(id_vector[id_vector != "-"]) > 50) {
            # Still too many IDs left, remove ones that appear only twice
            id_count <- id_count[id_count != 2 | names(id_count) == "-"]
            if (length(id_count[names(id_count) != "-"]) >= 20) {
              # When more than 20, move on with this filter, otherwise keep last
              id_vector <- names(id_count)
            }
          }
        }
      }
      
      cluster_ids <- if (color_scheme_column == "gene") TRUE else FALSE
      if (color_scheme_dataset == "comparisons") {
        comparisons <- uniform_color_scheme(
          comparisons = comparisons,
          id_column = color_scheme_column,
          ids = id_vector,
          colors = uniform_palette,
          cluster_ids = cluster_ids,
          alpha_comparisons = alpha_comparisons
        )
      } else if (color_scheme_dataset == "dna_segs") {
        dna_segs <- uniform_color_scheme(dna_segs = dna_segs,
                                         id_column = color_scheme_column,
                                         ids = id_vector,
                                         colors = uniform_palette,
                                         cluster_ids = cluster_ids,
                                         alpha_dna_segs = alpha_dna_segs
                                         )
      } else {
        if (column_in_comps & column_in_segs) {
          # Column is in both datasets, allows for simple uniform coloring
          full_data <- uniform_color_scheme(
            dna_segs = dna_segs,
            comparisons = comparisons,
            id_column = color_scheme_column,
            ids = id_vector,
            colors = uniform_palette,
            cluster_ids = cluster_ids,
            alpha_dna_segs = alpha_dna_segs,
            alpha_comparisons = alpha_comparisons
          )
          dna_segs <- full_data$dna_segs
          comparisons <- full_data$comparisons
        } else if (column_in_segs) {
          # Column is only in dna_segs, comparisons are made to match dna_segs
          dna_segs <- uniform_color_scheme(dna_segs = dna_segs,
                                           id_column = color_scheme_column,
                                           ids = id_vector,
                                           colors = uniform_palette,
                                           cluster_ids = cluster_ids,
                                           alpha_dna_segs = alpha_dna_segs
                                           )
          comparisons <- update_comparisons(dna_seg_input = dna_segs,
                                            comparison_input = comparisons,
                                            update_positions = FALSE,
                                            update_region_plot = FALSE,
                                            color_var = "fill"
                                            )
        } else if (column_in_comps) {
          comparisons <- uniform_color_scheme(
            comparisons = comparisons,
            id_column = color_scheme_column,
            ids = id_vector,
            colors = uniform_palette,
            cluster_ids = cluster_ids,
            alpha_comparisons = alpha_comparisons
          )
          dna_segs <- update_dna_segs(dna_seg_input = dna_segs,
                                      comparison_input = comparisons,
                                      update_region_plot = FALSE,
                                      color_var = "fill"
                                      )
        }
      }
    } else if (global_color_scheme == "sequential") {
      if (color_scheme_dataset == "auto") {
        full_data <- sequential_color_scheme(dna_segs = dna_segs, 
                                             comparisons = comparisons,
                                             seg_id = color_scheme_column,
                                             color_var = "fill",
                                             both_directions = TRUE
                                             )
        dna_segs <- full_data$dna_segs
        comparisons <- full_data$comparisons
      } else if (color_scheme_dataset == "comparisons") {
        comparisons <- update_comparisons(dna_seg_input = dna_segs,
                                          comparison_input = comparisons,
                                          seg_id = color_scheme_column,
                                          update_positions = FALSE,
                                          update_region_plot = FALSE,
                                          color_var = "fill"
                                          )
      } else if (color_scheme_dataset == "dna_segs") {
        dna_segs <- update_dna_segs(dna_seg_input = dna_segs,
                                    comparison_input = comparisons,
                                    seg_id = color_scheme_column,
                                    update_region_plot = FALSE,
                                    color_var = "fill"
                                    )
      }
    }
  }
  
  # check legend
  if (!is.null(legend_colors) & !is.null(legend_labels)) {
    # make legend based on these colors and labels
    legend_type <- "custom"
    if (length(legend_colors != length(legend_labels))) {
      stop("If 'legend_colors' and 'legend_labels' are both provided, then ",
           "they must be of equal length")
    }
    if (!is.character(legend_labels)) {
      stop("'legend_labels' must be a character vector")
    }
    if (all(sapply(legend_colors, 
                   function(x) tryCatch(is.matrix(col2rgb(x)),
                                        error = function(e) FALSE)
                   ))
        ) {
      stop("Not all elements of 'legend_colors' were valid colors")
    }
  } else if (!is.null(legend_column)) {
    legend_type <- "auto"
    if (legend_column == "auto") {
      if (global_color_scheme == "uniform") {
        legend_column <- color_scheme_column
      } else if (all(sapply(dna_segs, function(x) "gene" %in% names(x)))) {
        legend_column <- "gene"
      } else {
        legend_column <- "name"
      }
    }
    if (!all(sapply(dna_segs, function(x) legend_column %in% names(x)))) {
      stop("'legend_column' must be \"auto\" refer to a column that is ",
           "present in each dna_seg object")
    }
    legend_ids <- character()
    legend_colors <- character()
  } else {
    legend_type <- "none"
  }
  
  # deal with resetting symbols
  if (!is.null(gene_type)) {
    if (gene_type == "auto") {
      n_genes <- sapply(dna_segs, nrow)
      gene_type <- auto_gene_type(n_genes)
    }
    for (i in 1:n_dna_segs) {
      dna_segs[[i]]$gene_type <- gene_type
    }
  }
  
  # trim dna_segs & seg_plot # 
  # initiate new object: create subsegments by trimming original dna_seg
  seg_subplots <- list()
  dna_subsegs <- list()
  n_genes <- rep(0, n_dna_segs)
  for (i in 1:n_dna_segs) {
    n_subsegs <- nrow(xlims[[i]])
    dna_subsegs[[i]] <- list()
    seg_subplots[[i]] <- list()
    for (j in 1:n_subsegs) {
      dna_subsegs[[i]][[j]] <- trim.dna_seg(
        dna_segs[[i]],
        c(xlims[[i]]$x0[j], xlims[[i]]$x1[j])
      )
      if (seg_plot_h[[i]] > 0) {
        seg_subplots[[i]][[j]] <- trim.seg_plot(
          seg_plots[[i]],
          c(xlims[[i]]$x0[j], xlims[[i]]$x1[j])
        )
      }
      if (!is.null(alpha_dna_segs)) {
        # Convert fill into hexadecimal
        rgb_matrix <- col2rgb(dna_subsegs[[i]][[j]]$fill)
        hex_color <- rgb(red = rgb_matrix[1,], green = rgb_matrix[2,], 
                         blue = rgb_matrix[3,], maxColorValue = 255)
        tpc <- floor(alpha_dna_segs*256)
        tpc <- sprintf("%X", tpc)
        if (nchar(tpc) == 1) tpc <- paste("0", tpc, sep = "")
        dna_subsegs[[i]][[j]]$fill <- paste(hex_color, tpc, sep = "")
      }
      n_genes[i] <- n_genes[i] + length(
        dna_subsegs[[i]][[j]][gene_type != "boundaries", gene_type]
      )
    }
    if (legend_type == "auto") {
      legend_ids <- c(legend_ids, sapply(dna_subsegs[[i]], function(x) 
        unlist(x[, legend_column, with = FALSE], use.names = FALSE)))
      legend_colors <- c(legend_colors, sapply(dna_subsegs[[i]], function(x) 
        unlist(x[, "fill", with = FALSE], use.names = FALSE)))
    }
  }
  
  # calculate offsets #
  if (is.null(offsets)) {
    if (verbose) cat_time("Calculating offsets.")
    prel_offsets <- lapply(xlims, function(x)
                           c(0, rep(def_gap_length, nrow(x) - 1))
                           )
    offsets <- minimize_comps(comparisons, xlims, unpadded_lengths,
                              prel_offsets, fixed_gap_length)
    if (verbose) cat_time("Finished calculating offsets.")
  } else {
    offsets <- as.list(offsets)
    for (i in 1:n_dna_segs) {
      if (length(offsets[[i]]) == 1) {
        offsets[[i]] <- c(offsets[[i]],
                          rep(def_gap_length, nrow(xlims[[i]]) - 1)
                          )
      } else if (length(offsets[[i]]) != nrow(xlims[[i]])) {
        stop("The length of each element of 'offsets' must be either 1, or ",
             "equal to the number of subsegments in the corresponding dna_seg")
      }
    }
  }
  
  if (print_offsets) {
    formatted_offsets <- character(n_dna_segs)
    for (i in 1:n_dna_segs) {
      formatted_offsets[[i]] <- paste(offsets[[i]], collapse = " ")
    }
    cat_time("offsets:\n",
             paste0(formatted_offsets, collapse = ","),
             "\n", sep = ""
             )
  }
  
  # check if total length doesn't exceed max length
  if (limit_to_longest_dna_seg) {
    for (i in 1:n_dna_segs) {
      tot_length <- sum(c(xlims[[i]]$length, offsets[[i]]))
      if (tot_length > max_length) {
        excess <- tot_length - max_length
        # reduce offsets from the end
        for (j in length(offsets[[i]]):1) {
          # decrease offsets or set them to min_gap
          if ((offsets[[i]][j] - excess) < def_gap_length) {
            excess <- excess - offsets[[i]][j] + def_gap_length
            offsets[[i]][j] <- def_gap_length
          } else {
            offsets[[i]][j] <- offsets[[i]][j] - excess
            break
          }
        }
      }
    }
  } else {
    # recalculate lengths
    padded_lengths <- sapply(xlims, function(x) sum(x$length)) +
      sapply(offsets, sum)
    max_length <- max(padded_lengths)
  }
  
  # recalculate comps coordinates
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      comparisons[[i]] <- calc_comp_coor(offsets[[i]],
                                         xlims[[i]],
                                         comparisons[[i]],
                                         side = 1
                                         )
      comparisons[[i]] <- calc_comp_coor(offsets[[i+1]],
                                         xlims[[i+1]],
                                         comparisons[[i]],
                                         side = 2
                                         )
    }
  }
  # recalculate lengths
  padded_lengths <-
    sapply(xlims, function(x) sum(x$length)) + sapply(offsets, sum)
  max_length <- max(padded_lengths)
  longest_segment <- which.max(padded_lengths)
  
  if (verbose) cat_time("Collecting grobs.")
  
  #----------------------------------------------------------------------------#
  # collect grobs
  #----------------------------------------------------------------------------#
  # collect dna_seg, dna_seg_scale & seg_plot grobs #
  dna_seg_grobs <- list()
  dna_seg_scale_grobs <- list()
  for (i in 1:n_dna_segs) {
    dna_seg_grobs[[i]] <- list()
    dna_seg_scale_grobs[[i]] <- list()
    for (j in 1:length(dna_subsegs[[i]])) {
      # debug
      if (debug > 0 && debug < nrow(dna_subsegs[[i]][[j]])) {
        dna_subsegs[[i]][[j]] <- dna_subsegs[[i]][[j]][1:debug,]
      }
      # end debug
      dna_seg_grobs[[i]][[j]] <- dna_seg_grob(dna_subsegs[[i]][[j]],
                                              arrow_head_len,
                                              i,
                                              ...
                                              )
      if (dna_seg_scale[[i]]) {
        dna_seg_scale_grobs[[i]][[j]] <- dna_seg_scale_grob(
          range = xlims[[i]][j,c("x0","x1")],
          cex = scale_cex,
          unit = scale_unit,
          i = i,
          j = j
          )
      } else {
        dna_seg_scale_grobs[[i]][[j]] <- NULL
      }
    }
  }
  # collect seg_plot_grobs & ylims #
  seg_plot_grobs <- list()
  seg_plot_ylims <- list()
  for (i in 1:n_dna_segs) {
    seg_plot_grobs[[i]] <- list()
    xl_sg <- c(Inf, -Inf)
    for (j in 1:length(dna_subsegs[[i]])) {
      if (length(seg_plots[[i]]) > 0) {
        grb <- do.call(seg_subplots[[i]][[j]]$func, seg_subplots[[i]][[j]]$args)
        rng <- nice_ylim.seg_plot(seg_subplots[[i]][[j]])
        xl_sg[1] <- min(xl_sg[1], rng[1])
        xl_sg[2] <- max(xl_sg[2], rng[2])
        seg_plot_grobs[[i]][[j]] <- grb
      } else {
        seg_plot_grobs[[i]] <- list(NULL)
      }
    }
    if (is.null(seg_plots[[i]]$ylim)) {
      seg_plot_ylims[[i]] <- xl_sg
    } else {
      seg_plot_ylims[[i]] <- seg_plots[[i]]$ylim
    }
  }
  # collect seg_plot_yaxis_grobs #
  seg_plot_yaxis_grobs <- list()
  for (i in 1:n_dna_segs) {
    if (length(seg_plots[[i]]) > 0 && seg_plot_yaxis > 0) {
      seg_plot_yaxis_grobs[[i]] <- yaxis_grob(seg_plot_ylims[[i]],
                                              cex = seg_plot_yaxis_cex,
                                              n = seg_plot_yaxis,
                                              i
                                              )
    } else {
      seg_plot_yaxis_grobs[[i]] <- NULL
    }
  }
  # collect comparison grobs #
  comparison_grobs <- list()
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      # debug
      if (debug > 0 && debug < nrow(comparisons[[i]]))
        comparisons[[i]] <- comparisons[[i]][1:debug,]
      # end debug
      comparison_grobs[[i]] <- comparison_grob(comparisons[[i]],
                                               i,
                                               alpha = alpha_comparisons
                                               )
    }
  }
  
  # annotations #
  if (!is.null(annotations)) {
    annotation_grobs <- list()
    for (i in 1:n_dna_segs) {
      if (is.null(annotations[[i]])) {
        annotation_grobs[[i]] <- list(NULL)
      } else {
        annotation_grobs[[i]] <- list()
        for (j in 1:length(dna_subsegs[[i]])) {
          annot <- trim.annotation(annotations[[i]],
                                   xlims[[i]][j, c("x0", "x1")]
                                   )
          annotation_grobs[[i]][[j]] <- annotation_grob(annot, annotation_cex)
        }
      }
    }
  }
  
  # scale #
  if (scale) {
    scale_grob <- scale_grob(max_length)
    scale_h <- 1
  } else {
    scale_h <- 0
  }

  # main title #
  if (!is.null(main)) {
    main_grob <- textGrob(x = main_x,
                          y = 0.9,
                          label = main,
                          gp = gpar(cex = 1.2),
                          just = c(main_just, "top")
                          )
    main_h <- 1.8
  } else {
    main_h <- 0
  }
  # helps prevent annotations from being plotted out of bounds
  if (!is.null(annotations[[1]])) {
    main_h <- main_h + ifelse(annot_h[1] <= 3, 3 - annot_h[1], 0)
  } 

  # tree #
  if (!is.null(tree)) {
    # tree
    # check that a nice permutation is OK, return ys
    y <- permute_tree(tree, seg_labels)
    # feed tree grob with permutation transformed as y coords
    tree_grob <- phylog_grob(tree,
                             1 - ((y - 1) / (n_dna_segs - 1)),
                             clabel.leaves = dna_seg_label_cex,
                             labels.col = dna_seg_label_col,
                             clabel.nodes = tree_branch_labels_cex,
                             tree.scale = tree_scale
                             )
    tree_w <- unit(0.1, "npc") + tree_grob$width
  } else if (!is.null(seg_labels)) {
    # just labels
    tree_grob <- seg_label_grob(seg_labels,
                                cex = dna_seg_label_cex,
                                col = dna_seg_label_col
                                )
    tree_w <- tree_grob$width
  } else {
    # nothing
    tree_grob <- NULL
    tree_w <- unit(0, "npc")
  }
  if (!is.null(tree_width)) tree_w <- unit(tree_width, "inches")
  # reset scale_h if tree_scale is true
  if (tree_scale) scale_h <- 1
  
  # legend #
  if (legend_type == "auto") {
    legend_ids <- unlist(legend_ids, use.names = FALSE)
    legend_colors <- unlist(legend_colors, use.names = FALSE)
    
    # First filter out NA values and for IDs the user specified to put in
    if (is.null(legend_labels)) {
      ids_to_keep <- which(legend_ids != "NA")
    } else {
      ids_to_keep <- which(legend_ids != "NA" & legend_ids %in% legend_labels)
    }
    legend_labels <- legend_ids[ids_to_keep]
    legend_colors <- legend_colors[ids_to_keep]
    
    # Then try and find the first value for each color that is not just a "-"
    realname <- which(legend_labels != "-")
    unique_cols <- unique(legend_colors)
    ids_to_keep <- numeric(length(unique_cols))
    for (i in 1:length(unique_cols)) {
      col_inds <- which(legend_colors == unique_cols[i])
      named_inds <- which(col_inds %in% realname)
      if (length(named_inds) > 0) {
        ids_to_keep[i] <- col_inds[named_inds[1]]
      } else {
        ids_to_keep[i] <- col_inds[1]
      }
    }
    legend_labels <- legend_labels[ids_to_keep]
    legend_colors <- legend_colors[ids_to_keep]
    ### Put a limit on the legend in some fashion, i.e. based on approx_height var
  }
  if (legend_type != "none") {
    legend_grob <- legendGrob(
      legend_labels,
      pch = 15,
      default.units = "native",
      vgap = unit(0.3, "lines"),
      gp = gpar(col = legend_colors,cex = 0.9, name = "legendgrob")
    )
    longest_legend <- legend_labels[which.max(nchar(legend_labels))]
    legend_w <- unit(1, "strwidth", paste("OOO", longest_legend))
    if (any(dna_seg_scale)) {
      legend_w <- legend_w + unit(scale_cex, "lines")
    }
  } else {
    if (any(dna_seg_scale)) {
      legend_w <- unit(scale_cex, "lines")
    } else {
      legend_w <- unit(0, "null")
    }
  }
  
  # Calculate plot height based on some simple facts
  ### rework this, could be way more robust
  approx_height <- length(dna_segs)/2
  if (n_comparisons > 0) approx_height <- approx_height * 1.8
  if (!is.null(annotations)) approx_height = approx_height * 1.4
  
  # outfiles #
  if (!is.null(outfile)) {
    if (outfile_height == "auto") {
      outfile_height <- approx_height
    } else {
      if (suppressWarnings(is.na(as.double(outfile_height)))) {
        stop('The value "', outfile_height, '" for the \'outfile_height\' ',
             'argument could not be recognized as a numeric value')
      }
      outfile_height <- as.double(outfile_height)
    }
    
    if (outfile_width == "auto") {
      # Calculate plot width by taking maximum number of bases
      # maximum length of seg_labels, and existence of a tree
      outfile_width <- max_length / 4300
      if (!is.null(tree_grob)) {
        outfile_width <- outfile_width + max(nchar(seg_labels)) / 10
      }
      if (!is.null(tree)) {
        outfile_width <- outfile_width + tree_grob$xbase / 2.5
      }
      if (outfile_width < 5) {
        outfile_width <- 5
      }
      if (legend_type != "none") {
        outfile_width <- outfile_width + ((nchar(longest_legend) + 3) / 10)
      }
    } else {
      if (suppressWarnings(is.na(as.double(outfile_width)))) {
        stop('The value "', outfile_width, '" for the \'outfile_width\' ',
             'argument could not be recognized as a numeric value')
      }
      outfile_width <- as.double(outfile_width)
    }

    if (outfile_format == "pdf") {
      if (capabilities("cairo")) {
        cairo_pdf(outfile, height = outfile_height, width = outfile_width)
      } else {
        pdf(outfile, height = outfile_height, width = outfile_width)
      }
    } else if (outfile_format == "png") {
      png(outfile, 
          height = outfile_height,
          width = outfile_width,
          units = "in",
          res = 72
          )
    } else if (outfile_format == "bmp") {
      bmp(outfile,
          height = outfile_height,
          width = outfile_width,
          units = "in",
          res = 72
          )
    } else {
      stop("'outfile_format' must be one of: pdf, png, bmp")
    }
  }
  
  if (verbose) cat_time("Finished collecting grobs, now plotting.")
  
  #----------------------------------------------------------------------------#
  # plotting
  #----------------------------------------------------------------------------#
  # overall frame
  if (plot_new) grid.newpage()
  pushViewport(
    viewport(width = unit(1, "npc") - unit(1, "lines"),
             height = unit(1, "npc") - unit(1, "lines"),
             name = "oma"
             ),
    viewport(
      layout = grid.layout(2,
                           1,
                           heights = unit(c(main_h, 1), c("lines", "null"))
                           ),
      name = "oma_layout"
    )
  )
  # main title
  if (!is.null(main)) {
    pushViewport(viewport(layout.pos.row = 1, name = "main"))
    grid.draw(main_grob)
    upViewport()
  }

  # frame: columns=tree,maps,legend rows=maps+tree+legend,scale
  if (seg_plot_yaxis > 0 && !is.null(seg_plots[[longest_segment]])) {
    seg_plot_yaxis_w <- unit(3,
                             "grobwidth",
                             data = seg_plot_yaxis_grobs[[longest_segment]]
                             )
  } else {
    seg_plot_yaxis_w <- unit(0, "null")
  }
  pushViewport(viewport(
    layout.pos.row = 2,
    layout = grid.layout(
      2,
      4,
      heights = unit(c(1, scale_h), c("null", "lines")),
      widths = unit.c(tree_w, unit(1, "null"), seg_plot_yaxis_w, legend_w)
    ),
    name = "frame"
  ))
  # scale
  if (scale) {
    pushViewport(viewport(layout.pos.row = 2,
                          layout.pos.col = 2,
                          xscale = c(0, max_length),
                          name = "scale"
                          )
                 )
    grid.draw(scale_grob)
    upViewport()
  }

  # tree or labels. Height is 1-3 lines because margin is 2 in plotarea,
  # and 1 to center to the middle of each dna_seg (1/2 line top and bottom)
  if (!is.null(tree_grob)) {
    # extra margin if there is an annotation in the first dna_seg
    if (is.null(annotations[[1]])) {
      annot_margin <- unit(0, "lines")
    } else {
      annot_margin <- unit(annot_h[1], "lines")
    }
    # extra margin is there is a seg_plot in the first dna_seg
    if (is.null(seg_plot_grobs[[1]][[1]])) {
      seg_plot_margin <- unit(0, seg_plot_height_unit)
    } else {
      seg_plot_margin <- unit(seg_plot_height, seg_plot_height_unit)
    }
    # make a supplementary margin if there is a scale in the last
    # dna_seg to get labels facing the text. Hack.
    hli <- unit(0.5, "lines")
    dna_scale_margin <- unit(scale_cex * dna_seg_scale[[n_dna_segs]], "lines")
    margin_heights <- unit.c(seg_plot_margin,
                             annot_margin,
                             hli,
                             unit(n_dna_segs * (1 + seg_plot_height), "null"),
                             hli,
                             dna_scale_margin
                             )
    pushViewport(viewport(layout.pos.row = 1,
                          layout.pos.col = 1,
                          layout = grid.layout(6, 1, heights = margin_heights),
                          name = "tree_outer"
                          )
                 )
    pushViewport(viewport(layout.pos.row = 4,
                          width = unit(1, "npc")-unit(1, "lines"),
                          just = c("centre", "bottom"),
                          name = "tree"
                          )
                 )
    grid.draw(tree_grob$grob)
    upViewport(2) # up tree & tree_outer vp
  } 
  
  # plotting area
  pushViewport(viewport(layout.pos.row = 1,
                        layout.pos.col = 2,
                        name = "plotarea_outer",
                        clip = "on"
                        ),
               viewport(width = unit(1, "npc") - unit(1, "lines"),
                        height = unit(1, "npc") - unit(0, "lines"),
                        name = "plotarea", clip = "off"
                        )
               )

  # map grid
  pushViewport(viewport(
    layout = grid.layout(n_rows, 1, heights = dna_seg_heights),
    name = "map"
  ))
  # comparisons #
  if (n_comparisons > 0) {
    for (i in 1:n_comparisons) {
      pushViewport(viewport(layout.pos.row = 3 * i,
                            yscale = c(0,1),
                            xscale = c(0, max_length),
                            #clip = "on",
                            name = paste("comparison", i, sep = ".")
                            )
                   )
      # draw comparison grobs
      grid.draw(comparison_grobs[[i]])
      upViewport() # pop comparisons[[i]] vp
    }
  }
  # seg_plots, dna_segs, annotations, scales #
  for (i in 1:n_dna_segs) {
    n_dna_subsegs <- length(dna_subsegs[[i]])
    n_cols <- n_dna_subsegs * 2 + 1
    widths <- numeric(n_cols)
    widths[1:n_dna_subsegs * 2] <- xlims[[i]]$length
    widths[1:n_dna_subsegs * 2 - 1] <- offsets[[i]]
    widths_units <- unit(widths, rep("native", n_cols))
    heights <- unit(c(annot_h[i], 1, scale_cex*dna_seg_scale[i]),
                    c("lines", "null", "lines")
                    )
    # push seg_plot grid
    pushViewport(viewport(
      layout.pos.row = 3 * i - 2,
      layout = grid.layout(1,
                           n_cols,
                           widths = widths_units,
                           just = c("left", "centre")
                           ),
      xscale = c(0,max_length),
      name = paste("seg_plot", i, sep = ".")
    ))
    for (j in 1:n_dna_subsegs) {
      idx <- if (xlims[[i]]$strand[j] == 1) c("x0", "x1") else c("x1", "x0")
      xscale <- as.numeric(xlims[[i]][j,idx])
      if (!is.null(seg_plots[[i]])) {
        pushViewport(viewport(layout.pos.col = j * 2,
                              yscale = seg_plot_ylims[[i]],
                              xscale = xscale,
                              just = c("left", "centre"),
                              name = paste("seg_subplot", i, j, sep = ".")
                              )
                     )
        grid.draw(seg_plot_grobs[[i]][[j]])
        upViewport() # up seg_subplot vp
        
      }
    }
    ## Draw y axis
    if (!is.null(seg_plots[[i]]) && seg_plot_yaxis > 0) {
      pushViewport(viewport(
        layout.pos.col = n_cols,
        yscale = seg_plot_ylims[[i]],
        width = unit(1, "grobwidth", data = seg_plot_yaxis_grobs[[i]]),
        just = c("left", "centre"),
        name = paste("seg_plot_yaxis", i, sep = ".")
      ))
      grid.draw(seg_plot_yaxis_grobs[[i]])
      upViewport() # up seg plot yaxis
    }
    upViewport() # up seg_plot vp
    # push dna_seg grid (subsegments in cols, annotations, genes and
    # scales in rows)
    pushViewport(viewport(
      layout.pos.row = 3 * i - 1,
      layout = grid.layout(3,
                           n_cols,
                           heights = heights,
                           widths = widths_units,
                           just = c("left", "centre")
                           ),
      xscale = c(0, max_length),
      name = paste("scale_and_dna_seg", i, sep = ".")
    ))
    for (j in 1:n_dna_subsegs) {
      # calculate xscale
      idx <- if (xlims[[i]]$strand[j] == 1) c("x0", "x1") else c("x1", "x0")
      xscale <- as.numeric(xlims[[i]][j,idx])
      # annotation
      if (!is.null(annotations[[i]])) {
        pushViewport(viewport(layout.pos.row = 1,
                              layout.pos.col = j * 2,
                              yscale = c(0, 1),
                              xscale = xscale,
                              just = c("left", "centre"),
                              name = paste("annotation", i, j, sep = ".")
                              )
                     )
        grid.draw(annotation_grobs[[i]][[j]])
        upViewport() # up annotation vp
      }
      # dna_seg_scale
      if (dna_seg_scale[i]) {
        pushViewport(viewport(layout.pos.row = 3,
                              layout.pos.col = j * 2,
                              yscale = c(0, 1),
                              xscale = xscale,
                              just = c("left", "centre"),
                              name = paste("dna_seg_scale", i, j, sep = ".")
                              )
                     )
        grid.draw(dna_seg_scale_grobs[[i]][[j]])
        upViewport() # up dna_seg_scale vp
      }
      # dna_seg itself
      pushViewport(viewport(layout.pos.row = 2,
                            layout.pos.col = j * 2,
                            yscale = c(0, 1),
                            xscale = xscale,
                            just = c("left", "centre"),
                            name = paste("dna_seg", i, j, sep = ".")
                            )
                   )
      # draw segment line
      if (!dna_seg_line[i] == "FALSE") {
        grid.segments(x0 = unit(xlims[[i]]$x0[j], "native"),
                      y0 = unit(0.5, "native"),
                      x1 = unit(xlims[[i]]$x1[j], "native"),
                      y1 = unit(0.5, "native"),
                      name = paste("dna_seg_line", i, j, sep = "."),
                      gp = gpar(col = dna_seg_line[i])
                      )
      }
      # draw dna_seg grobs
      grid.draw(dna_seg_grobs[[i]][[j]])
      upViewport() # up dna_seg
      # draw gap, but not at pos 1
      if (j > 1) {
        pushViewport(viewport(layout.pos.row = 2,
                              layout.pos.col = j * 2 - 1,
                              yscale = c(0, 1),
                              xscale = c(0, widths[2 * j - 1]),
                              just = c("centre", "centre"),
                              name = paste("gap", i, j, sep = ".")
                              )
                     )
        grid.draw(gap_grob(w = def_gap_length, m = widths[2 * j - 1] / 2, i, j))
        upViewport() # up gap vp
      }
    }
    upViewport() # up scale_and_dna_seg vp
  }
  upViewport(2) # pop map viewports
  upViewport() # pop plotarea viewport
  
  # legend
  if (legend_type != "none") {
    pushViewport(viewport(layout.pos.row = 1,
                          layout.pos.col = 4,
                          just = c("left", "top"),
                          name = "legend"
                          )
                 )
    grid.draw(legend_grob)
    upViewport()
  }
  
  upViewport() # pop plotarea viewport
  upViewport(2) # pop frame+oma viewport
 
  if (!is.null(outfile)) dev.off()
  if (verbose) cat_time("Finished plotting.")
  invisible()
}

################################################################################
# Data set Documentation
################################################################################


#' Comparison of 4 Bartonella genomes
#' 
#' A comparison of 4 Bartonella genomes by BLAST. 
#' 
#' @format `barto`, a list of 3 data frame lists, representing the
#' four genomes and their pairwise comparisons:
#' \describe{
#'   \item{`dna_segs`}{A list of 4 `dna_seg` objects, containing all the 
#'   protein-coding genes for each genome. Obtained by reading ptt files
#'   downloaded from NCBI with [read_dna_seg_from_ptt].}
#'   \item{`comparisons`}{A list of 3 `comparison` objects, obtained by doing
#'   genome-to-genome (fasta files) BLASTS, and then reading the resulting tab
#'   files with [read_comparison_from_blast].}
#'   \item{`rnt_segs`}{A list of 4 `dna_seg` objects, containing all the RNA
#'   genes of the four genomes. Obtained by reading rnt files downloaded from
#'   NCBI with [read_dna_seg_from_ptt].}
#' }
#' 
#' @name barto
#' @usage data(barto)
#' 
#' @examples
#' data(barto)
#' plot_gene_map(barto$rnt_segs, barto$comparisons, gene_type = "blocks")
#' 
"barto"

#' Comparisons of subsegments of the Y chromosome in human and chimp
#' 
#' A subsegment of the Y chromosome in Homo sapiens and Pan troglodytes, to
#' illustrate support for exons and introns.
#' 
#' @format `chrY_subseg`, a list of two data frame lists, representing the Y
#' segment in the two species and a comparison between them:
#' \describe{
#'   \item{`dna_segs`}{A list of 2 `dna_seg` objects, containing 3 genes each.
#'   Each exon and intron is a seperate feature (row) in the `dna_seg`.}
#'   \item{`comparison`}{A list containing 1 `comparison` object.}
#' }
#' 
#' @name chrY_subseg
#' @usage data(chrY_subseg)
#' 
#' @examples
#' data(chrY_subseg)
#' plot_gene_map(chrY_subseg$dna_segs, chrY_subseg$comparison,
#'               dna_seg_scale = TRUE, scale = FALSE)
#' 
"chrY_subseg"

#' Mauve backbone of 4 Bartonella genomes
#' 
#' The result of a multiple genome alignment of 4 Bartonella genomes with Mauve 
#' 
#' @format `mauve_bbone`, a list of two data frame lists, representing the
#' regions which are conserved in at least 2 genomes:
#' \describe{
#'   \item{`dna_segs`}{A list of 4 `dna_seg` objects, containing the mauve 
#'   blocks for each genome.}
#'   \item{`comparison`}{A list of 3 `comparison` objects.}
#' }
#' A bash script to obtain the same file as in the data is available in the
#' extdata folder of the package. Find its location by running
#' `system.file('extdata/mauve.sh',package = 'genoPlotR')`.
#' 
#' The resulting backbone file can then be read with [read_mauve_backbone].
#' 
#' @name mauve_bbone
#' @usage data(mauve_bbone)
#' 
#' @references Mauve: https://darlinglab.org/mauve/mauve.html
#' 
#' @examples
#' data(mauve_bbone)
#' plot_gene_map(mauve_bbone$dna_segs, mauve_bbone$comparisons)
#' 
"mauve_bbone"

#' Three genes data set
#' 
#' A set of three made-up genes, compared in three chromosomes.
#' 
#' @format `three_genes`, a list of two data frame lists, representing three
#' genes in three DNA segments:
#' \describe{
#'   \item{`dna_segs`}{A list of 3 `dna_seg` objects, containing three genes
#' (rows) each.}
#'   \item{`comparisons`}{A list of 2 `comparison` objects.}
#' }
#' 
#' @name three_genes
#' @usage data(three_genes)
#' 
#' @examples
#' data(three_genes)
#' plot_gene_map(three_genes$dna_segs, three_genes$comparisons)
#' 
"three_genes"


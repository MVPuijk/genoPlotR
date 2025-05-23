#!/usr/bin/env Rscript

################################################################################
# Wrapper script for genoPlotR
################################################################################

################################################################################
# Parsing and checking options
################################################################################

# Define the options and their usage text
option_list <- list(
  optparse::make_option(
    c("-s", "--seg_files"), action = "store", default = NA, type = "character",
    metavar = '["file path1,file path2"]',
    help = paste(
      '(Mandatory): Path to the dna_seg files to read in, with',
      'file names separated by ",". GenBank files are recommended.',
      'Each file will become a single dna_seg, so genomes',
      'consisting of multiple contigs or chromosomes should be' ,
      'concatenated into a single file.',
      'Alternatively, a single file with the suffix ".RDS" can be',
      'provided, which will be read as an R object containing a list',
      'of dna_segs.',
      sep = "\n\t\t")),
  optparse::make_option(
    c("-i", "--ids"), action = "store", default = NA, type = "character",
    metavar = "[file path]",
    help = paste(
      'A tab-delimited file with IDs corresponding to features',
      'that should be plotted. The header should at least contain',
      'an "id" column, but can include other column names to add',
      'to the dna_seg objects. The IDs can refer to values from',
      'any dna_seg column, specified by the --id_tags option.',
      'IMPORTANT: If --comp_format is a blastp algorithm or',
      'DIAMOND, then the "region_plot" column will be carried over',
      'to the next dna_seg object, for any features linked by the',
      'comparisons.',
      'An example header:\n',
      'id\tseg_label\tfill\tregion_plot\n',
      '- id: An identifying feature ID',
      '- seg_label: The name of the dna_seg object to query',
      '- fill: Anything identifiable as a color by R',
      '  (e.g. "red", "blue4", "salmon", "#88CCEE")',
      '- region_plot: If the value is "plot", "true", "yes", or',
      '  "y", then the genome neighbourhood of that feature will',
      '  be plotted.', sep = "\n\t\t")),
  optparse::make_option(
    c("-I", "--id_tags"), action = "store", default = "name, locus_id", 
    type = "character", metavar = '["column1,column2"]',
    help = paste(
      'A list of dna_seg column names to use alongside --ids.',
      'The IDs provided by the --ids option will be used as',
      'queries in the columns provided by this option. Column',
      'names must be separated by commas.',
      '(default: "%default")', sep = "\n\t\t")),
  optparse::make_option(
    c("--dna_seg_mode"), action = "store", default = "all", type = "character",
    metavar = "[dna_seg parsing mode]",
    help = paste(
      'How should the dna_seg files supplied by --seg_files be',
      'read? Must be one of: "all", "ids", "fast"',
      '- all: Read all dna_seg files provided.',
      '- ids: Read only the dna_seg files mentioned in the --ids',
      '  file, using the "seg_label" column.',
      '- fast: Only allowed when making a regional plot using',
      '  blastp or DIAMOND. After reading in a dna_seg, an',
      '  alignment is made between it and the next dna_seg,',
      '  transferring over the "region_plot" column. If there',
      '  are no features to be plotted in this next dna_seg,',
      '  then it is skipped.', sep = "\n\t\t")),
  optparse::make_option(
    c("-c", "--comp_format"), action = "store", default = NA, 
    type = "character", metavar = "[comparison format]",
    help = paste(
      'The format to use for reading in comparison files.',
      'If the necessary files are not found in --comp_path, then',
      'BLAST or DIAMOND can be chosen here to create them.', 
      'Must be one of: "diamond", "blast", "tab", "orthofinder",',
      '"orthomcl, "mmseqs2", or one of the supported blast',
      'programs ("blastp", "blastp-fast", "blastp-short",',
      '"blastn", "blastn-short", "megablast", "dc-megablast").',
      '',
      '- diamond: tab-separated DIAMOND result file (--outfmt 6)',
      '- blast: tab-separated BLAST result file (-outfmt 6)',
      '- tab: a custom tab-delimited file with the following',
      '  structure: "start1 end1 start2 end2 col gene1 gene2"',
      '- orthofinder: An OrthoFinder format file containing',
      '  orthogroups for the provided dna_segs',
      '- orthomcl: A file containing orthogroups in the orthomcl',
      '  format',
      '- mmseqs2: A file containing clusters of genes in the ',
      '  "Cluster TSV format" from mmseqs2, or the cluster ',
      '  output from DIAMOND', sep = "\n\t\t")),
  optparse::make_option(
    c("-C", "--comp_path"), action = "store", default = NA, type = "character",
    metavar = "[directory/file path]",
    help = paste(
      'Path to the directory for the comparison files, or the path',
      'to a file containing orthogroups, if --comp_format is',
      '"orthofinder" or "orthomcl". This option is ignored if',
      '--comp_format is not provided. With a --comp_format of',
      'DIAMOND or any BLAST program, this script will attempt to',
      'create the necessary comparison files if it can not locate',
      'them in the directory provided. The file names should',
      'contain the names of the dna_segs joined by an',
      'underscore, e.g. "Genome1_Genome2.suffix"',
      'default: current directory',
      sep = "\n\t\t")),
  optparse::make_option(
    c("--comp_mode"), action = "store", default = "full", 
    type = "character", metavar = "[comparison mode]",
    help = paste(
      'The filter mode used when sequentially parsing the sequence',
      'alignments for each dna_seg object. This option is ignored',
      'unless a BLAST program or DIAMOND was chosen using',
      '--comp_format. It is recommended to use "besthit" and',
      '"bidirectional" only when doing comparisons on a list of',
      'features, like genes or proteins. Must be one of:',
      '"besthit", "bidirectional", "full". (default: %default)',
      '- besthit: Filters the sequence alignment results to',
      '  include only the best hit for each query dna_seg',
      '- bidirectional: Filters the sequence alignment results to',
      '  include only bidirectional best hits for each pair of',
      '  dna_segs',
      '- full: The sequence alignment results are not filtered for',
      '  best hits', sep = "\n\t\t")),
  optparse::make_option(
    c("--diamond_sensitivity"), action = "store", default = "default", 
    type = "character", metavar = '["sensitivity"]',
    help = paste(
      'The sensitivity option used when running DIAMOND. This',
      'option is ignored unless DIAMOND was chosen using',
      '--comp_format. Must be one of the following arguments,',
      'in order of least to most sensitive: "fast", "default",',
      '"mid-sensitive", "sensitive", "more-sensitive",',
      '"very-sensitive", "ultra-sensitive". (default: %default)',
      sep = "\n\t\t")),
  optparse::make_option(
    c("--update_positions"), action = "store", default = "auto",
    type = "character", metavar = "[auto|true|false]",
    help = paste(
      'Decides whether or not positions of comparison objects',
      'should be updated to match the positions in the dna_segs.',
      'Must be one of: "auto", "true", "false". If true, always',
      'update, if false, never update.',
      'The default is auto, which updates the positions when',
      '--comp_format is DIAMOND or any of the blastp algorithms.',
      sep = "\n\t\t")),
  optparse::make_option(
    c("--update_region_plot"), action = "store", default = "auto",
    type = "character", metavar = "[auto|true|false]",
    help = paste(
      'Decides whether or not the region_plot column of dna_segs',
      'should be updated. If "true", the region_plot attribute',
      'will be given to each feature that can be linked to a',
      'feature that already has the region_plot attribute.',
      'Must be one of: "auto", "true", "false". If true, always',
      'update, if false, never update.',
      'The default is auto, which updates region_plot when',
      '--comp_format is DIAMOND or any of the blastp algorithms.',
      sep = "\n\t\t")),
  optparse::make_option(
    c("-o", "--out"), action = "store", default = "out.pdf",
    type = "character", metavar = "[file path]",
    help = paste(
      'Name of output PDF file. (default: %default)',
      sep = "\n\t\t")),
  optparse::make_option(
    c("--height"), action = "store", default = "auto",
    metavar = "[num]", help = "PDF height, in inches. (default: %default)"),
  optparse::make_option(
    c("--width"), action = "store", default = "auto", 
    metavar = "[num]", help = "PDF width, in inches. (default: %default)"),
  optparse::make_option(
    c("-t", "--tree"), action = "store", default = NA, type = "character",
    metavar = "[file path]",
    help = paste(
      'Path to newick format tree file with the dna_seg labels',
      '(file names) as leaves.', sep = "\n\t\t")),
  optparse::make_option(
    c("--annotations"), action = "store_true", default = FALSE,
    help = paste(
      'Gene annotations from "gene name" and "product"',
      'will be plotted.', sep = "\n\t\t")),
  optparse::make_option(
    c("-l", "--legend_column"), action = "store", default = NA, 
    type = "character", metavar = "[column]",
    help = paste(
      'A dna_seg column to use as the basis for plotting a legend.',
      'Must be the name of a column that appears in all dna_segs or',
      '"auto", in which case it will pick the "gene" column if it is',
      'present, and the "name" column if it is not. A legend will',
      'then be plotted based on the colors of the dna_seg features,',
      'using all of the unique values found in the chosen column.',
      sep = "\n\t\t")),
  optparse::make_option(
    c("-g", "--global_color_scheme"), action = "store", default = NA, 
    type = "character", metavar = "[color scheme]",
    help = paste(
      'A color scheme to apply to the dna_segs and/or comparisons.',
      'Must be one of: "gradient", "uniform", "sequential".',
      '- gradient: Calculates color gradients based on numerical',
      '  columns. It uses red hues for direct comparisons, and',
      '  blue hues for reverse comparisons. It is applied only to',
      '  comparisons, and bases the coloring on comparison',
      '  metrics, namely "per_id", "bit_score", and "e_value"',
      '  columns, in that order.',
      '- uniform: Generates a color palette based on the possible',
      '  values a column can have for all the features that are',
      '  being displayed. Looks for a column related to homology,',
      '  followed by "gene" and "gene_type" if it cannot find any.',
      '  Applies itself to both dna_segs and comparisons.',
      '- sequential: Takes alterations made on the colors of',
      '  features in dna_segs, and transfers them over to the',
      '  other dna_segs and comparisons. Only works if changes',
      '  have been made to at least 1 dna_seg already. This can',
      '  be done by using --id_tags, for instance.',
      sep = "\n\t\t")),
  optparse::make_option(
    c("-G", "--color_scheme_dataset"), action = "store", default = "auto", 
    type = "character", metavar = "[dataset]",
    help = paste(
      'Which dataset to apply the color scheme to. Only applicable',
      'when --global_color_scheme is set to uniform or sequential.',
      'Must be one of: "auto", "dna_segs", "comparisons".',
      'The default is auto, which will attempt to automatically',
      'determine how best to apply the color scheme based on the',
      'available data.',
      sep = "\n\t\t")),
  optparse::make_option(
    c("-a", "--alpha_dna_segs"), action = "store", default = NA,
    type = "character", metavar = "[num]",
    help = paste(
      'A number between 0 (transparent) and 1 (opaque) for the',
      'alpha-transparency level that will be applied to dna_segs.',
      sep = "\n\t\t")),
  optparse::make_option(
    c("-A", "--alpha_comparisons"), action = "store", default = NA,
    type = "character", metavar = "[num]",
    help = paste(
      'A number between 0 (transparent) and 1 (opaque) for the',
      'alpha-transparency level that will be applied to',
      'comparisons.', sep = "\n\t\t")),
  optparse::make_option(
    c("--dna_seg_scale"), action = "store_true", default = FALSE,
    help = paste(
      'A scale will be plotted for each dna_seg.', sep = "\n\t\t")),
  optparse::make_option(
    c("-r", "--region_size"), action = "store", default = 10000, 
    type = "character", metavar = "[num]",
    help = paste(
      'Neighbourhood size to plot, unless a sequence boundary is',
      'found first, in bp. Setting this argument to 0 will plot',
      'the entirety of the dna_segs. This argument is ignored when',
      'xlims are provided via --xlims_in or --xlims_from_file',
      '(default: %default)', sep = "\n\t\t")),
  optparse::make_option(
    c("-x", "--xlims_in"), action = "store", default = NA, type = "character",
    metavar = '["xlims_seg1,xlims_seg2"]',
    help = paste(
      'Predefined positions to plot for each dna_seg (and',
      'subsegment). If --xlims_from_file is given, this option is',
      'ignored. dna_segs must be separated by commas, in plotting',
      'order, with subsegments separated by spaces (same format as',
      '--print_xlims). e.g.:\n',
      '"1 10000 155000 154000,500 10500"', sep = "\n\t\t")),
  optparse::make_option(
    c("-X", "--xlims_from_file"), action = "store", default = NA, 
    type = "character", metavar = "[file path]",
    help = paste(
      'Path to a tab-separated file with predefined positions to',
      'plot for each dna_seg (and subsegment) in the same format',
      'as that of --xlims_out. One line per subsegment, in',
      'plotting order, with both positions and their dna_seg',
      'label. e.g.:\n',
      'x0\tx1\tseg_label',
      '0\t13000\tdna_seg1',
      '25000\t35000\tdna_seg1',
      '500\t10500\tdna_seg2',
      '25901\t35901\tdna_seg2',
      'Note: If provided, only dna_segs present in this file',
      'will be plotted.', sep = "\n\t\t")),
  optparse::make_option(
    c("-p", "--print_xlims"), action = "store_true", default = FALSE,
    help = 'Prints out neighbourhood limit positions to the terminal.'),
  optparse::make_option(
    c("--xlims_out"), action = "store", default = NA, type = "character",
    metavar = "[file path]",
    help = paste(
      'Name of output xlim file, which will contain the',
      'neighbourhood limits that were used to plot with.', sep = "\n\t\t")),
  optparse::make_option(
    c("-f", "--offsets_in"), action = "store", default = NA, type = "character",
    metavar = '["offsets1,offsets2"]',
    help = paste(
      'A list of predefined offset values for each dna_seg and',
      'subsegment. If --offsets_from_file is given, this option is',
      'ignored. For each dna_seg 1 value must be provided, or as',
      'many values as there are subsegments for that dna_seg.',
      'dna_segs must be separated by commas, in plotting order,',
      'with subsegments separated by spaces. e.g.:\n',
      '"150 1300,300"', sep = "\n\t\t")),
  optparse::make_option(
    c("-F", "--offsets_from_file"), action = "store", default = NA,
    type = "character", metavar = "[file path]",
    help = paste(
      'Path to a file with predefined offset values for each',
      'dna_seg and subsegment. One line per dna_seg, in plotting',
      'order, with subsegments separated by spaces, e.g.:\n',
      '150 1300','300', sep = "\n\t\t")),
  optparse::make_option(
    c("--print_offsets"), action = "store_true", default = FALSE,
    help = 'Prints out offset positions to the terminal.'),
  optparse::make_option(
    c("-S", "--dna_segs_out"), action = "store", default = NA, 
    type = "character", metavar = "[file path]",
    help = paste(
      'Name of output dna_segs file. If provided, the dna_segs',
      'will be saved as an R object just after reading them in,',
      'before any alterations are made to it. This way, this',
      'script can be called again using the same arguments, except',
      'for providing this file to --seg_files, which should in',
      'most cases produce the same result, but without having to',
      'wait for the dna_seg files to be parsed.', sep = "\n\t\t")),
  optparse::make_option(
    c("--workspace_out"), action = "store", default = NA, 
    type = "character", metavar = "[file path]",
    help = paste(
      'Name of output file. If provided, all of the arguments that',
      'are passed to main plot function (plot_gene_map) are saved',
      'to this file as well, as an R object that can be loaded in',
      'R using the readRDS() function.', sep = "\n\t\t")),
  optparse::make_option(
    c("-n", "--threads"), action = "store", default = 1, type = "integer",
    metavar = "[NUM]",
    help = paste(
      'Number of threads for parallel processing, used by BLAST',
      'and DIAMOND. (default: %default)', 
      sep = "\n\t\t")),
  optparse::make_option(
    c("-v", "--verbose"), action = "store_true", default = FALSE,
    help = paste(
      'Print out extra information, mostly on writing and ',
      'reading in files', sep = "\n\t\t"))
)

# Parse options
opt_p <- optparse::OptionParser(
  option_list = option_list, add_help_option=TRUE,
  usage = '%prog -g ["dna_seg files"] [options]')
opt = optparse::parse_args(opt_p)

### DEBUG arguments

# input_dir <- ifelse(Sys.info()["sysname"] == "Windows",
#                     "D:/Documents", "/mnt/d/Documents")
# dataset <- "report_test"
# output_dir <- file.path(input_dir, "genoPlotR/out", dataset)

# opt$seg_files <- "out/report_test_raw.RDS"
# opt$dna_segs_out <- "out/report_test_raw2.RDS"
# opt$ids <- file.path(input_dir, "genoPlotR/data/ids.tsv")
# opt$id_tags <- "name, locus_id"
# opt$dna_seg_mode <- "fast"
# opt$comp_format <- "diamond"
# opt$comp_path <- output_dir
# opt$comp_mode <- "full"
# opt$update_positions <- "auto"
# opt$update_region_plot <- "auto"
# opt$out <- file.path(output_dir, "report_test.pdf")
# opt$height <- "auto"
# opt$width <- "auto"
# opt$tree <- NA
# opt$annotations <- TRUE
# opt$dna_seg_scale <- TRUE
# opt$global_color_scheme <- "sequential"
# opt$alpha_comparisons <- 0.5
# opt$alpha_dna_segs <- NA
# opt$region_size <- 10000
# opt$print_xlims <- TRUE
# opt$xlims_out <- NA
# opt$xlims_in <- NA
# opt$xlims_from_file <- NA
# opt$offsets_in <- NA
# opt$offsets_from_file <- NA
# opt$threads <- 1
# opt$verbose <- TRUE
# opt$help <- FALSE


################################################################################
# Parsing and checking of arguments
################################################################################
# Function for printing time stamped verbose messages
cat_time <- function(...) {
  cat(paste0(format(Sys.time()), ": ", ..., "\n"))
}

# Function for printing out help function and then stopping the program.
stop_help <- function(opt_p, ...) {
  optparse::print_help(opt_p)
  stop(..., call. = FALSE)
}

# Internal helper function to get dna_seg labels
# Is here to avoid making its copy in the package itself public
get_seg_labels <- function(dna_segs, seg_labels, required = FALSE) {
  if (is.null(seg_labels)) {
    if (length(dna_segs) == sum(names(dna_segs) != "", na.rm = TRUE)) {
      # If each element of dna_segs is named, use those
      seg_labels <- names(dna_segs)
    } else {
      # Does each element of the dna_segs list have a "seg_label" atrribute?
      seg_labels <- unlist(lapply(dna_segs, attr, "seg_label"))
      if (length(dna_segs) == length(seg_labels) & 
          !any(duplicated(seg_labels))
      ) {
        seg_labels <- unlist(lapply(dna_segs, attr, "seg_label"))
      }
    }
  }
  if (is.null(seg_labels) & required) {
    seg_labels <- paste(rep("dna_seg", length(dna_segs)),
                        seq(1:length(dna_segs)),
                        sep = "_"
    )
  }
  # check length of labels
  if (!is.null(seg_labels) && !(length(seg_labels) == length(dna_segs)))
    stop("'seg_labels' must be NULL or the same length as 'dna_segs'")
  
  seg_labels
}

# Check flags
# Flags end up as "NA" when the user enters an argument after the long flag
# It already errors out automatically when you attempt this with the short flag
# So the long flag functionality is made to error out to match
use_annots <- opt$annotations
if (is.na(use_annots)) {
  stop('The --annotations option does not support an argument')
}
print_xlims <- opt$print_xlims
if (is.na(print_xlims)) {
  stop('The --print_xlims option does not support an argument')
}
print_offsets <- opt$print_offsets
if (is.na(print_offsets)) {
  stop('The --print_offsets option does not support an argument')
}
dna_seg_scale <- opt$dna_seg_scale
if (is.na(dna_seg_scale)) {
  stop('The --dna_seg_scale option does not support an argument')
}
verbose <- opt$verbose
if (is.na(verbose)) {
  stop('The --verbose option does not support an argument')
}

if (verbose) cat_time("Checking arguments.")

# Check seg_files
if (is.na(opt$seg_files)) {
  stop_help(opt_p,
            'dna_seg files must be provided using the --seg_files option')
}
dna_seg_files <- unlist(strsplit(opt$seg_files, "\\, ?"))
if (length(dna_seg_files) == 1 & all(endsWith(dna_seg_files, ".RDS"))) {
  read_in_RDS <- TRUE
  if (!file.exists(dna_seg_files)) {
    stop_help(opt_p,
              "Could not find the specified dna_seg files (--seg_files)")
  }
} else {
  read_in_RDS <- FALSE
  seg_files <- c()
  for (i in dna_seg_files) {
    seg_files <- c(seg_files, Sys.glob(i))
  }
  dna_seg_files <- seg_files
}
if (! length(dna_seg_files)) {
  stop_help(opt_p, "Could not find the specified dna_seg files (--seg_files)")
}

# Check dna_segs_out
outfile_dna_segs <- opt$dna_segs_out
if (!is.na(outfile_dna_segs) & verbose) {
  if (file.exists(outfile_dna_segs)) {
    cat_time("dna_segs output file: ", outfile_dna_segs, " already ",
             "exists. File will be overwritten.")
  }
}

# Check workspace_out
workspace_out <- opt$workspace_out
if (!is.na(workspace_out) & verbose) {
  if (file.exists(workspace_out)) {
    cat_time("Workspace output file: ", workspace_out, " already ",
             "exists. File will be overwritten.")
  }
}

# Check ids
if (!is.na(opt$ids)) {
  if (!file.exists(opt$ids)) {
    stop_help(opt_p, "Could not find the specified ids file (--ids)")
  }
  use_ids <- TRUE
} else {
  use_ids <- FALSE
}

# Check id_tags
id_tags <- unlist(strsplit(opt$id_tags, ", ?"))

# Check comp format
if (!is.na(opt$comp_format)) {
  comp_format <- tolower(opt$comp_format)
  if (!any(comp_format == c("diamond", "blast", "tab", "orthofinder", 
                            "orthomcl", "mmseqs2", "blastp", "blastp-fast", 
                            "blastp-short", "blastn", "blastn-short", 
                            "megablast", "dc-megablast", "tblastx"))) {
    stop_help(opt_p, comp_format, " is not a valid argument for --comp_format")
  }
  if (comp_format == "blast") comp_format <- "blastp"
  if (comp_format == "diamond") {
    sensitivity <- opt$diamond_sensitivity
    if (!any(sensitivity == c("fast", "default", "mid-sensitive", 
                              "sensitive", "more-sensitive", 
                              "very-sensitive", "ultra-sensitive"))) {
      stop_help(opt_p, sensitivity, 
                " is not a valid argument for --diamond_sensitivity")
    }
  }
} else {
  comp_format <- "none"
}

# Check comp mode
if (any(comp_format == c("diamond", "blastp", "blastp-fast", "blastp-short",
                         "blastn", "blastn-short", "megablast", "dc-megablast", 
                         "tblastx"))) {
  comp_mode <- tolower(opt$comp_mode)
  if (!any(comp_mode == c("besthit", "bidirectional", "full",
                          "besthits", "bbh"))) {
    stop_help(opt_p, comp_mode, " is not a valid argument for --comp_mode")
  }
  if (comp_mode == "besthits") comp_mode <- "besthit"
  if (comp_mode == "bbh") comp_mode <- "bidirectional"
} else {
  comp_mode <- "none"
}

# Check comp files
if (comp_format != "none") {
  if (is.na(opt$comp_path)) {
    if (any(comp_format == c("orthofinder", "orthomcl", "mmseqs2"))) {
      stop_help(opt_p, 'If "orthofinder", "orthomcl", or "mmseqs2" is chosen ',
                'using --comp_format then a single file containing ',
                'orthogroups should be provided using --comp_path')
    }
    comp_path <- getwd()
  } else {
    comp_path <- opt$comp_path
    if (any(comp_format == c("orthofinder", "orthomcl", "mmseqs2"))) {
      if (!file.exists(comp_path)) {
        stop_help(opt_p, 'If "orthofinder", "orthomcl", or "mmseqs2" is ',
                  'chosen using --comp_format then a single file containing ',
                  'orthogroups should be provided using --comp_path'
                  )
      }
    } else {
      if (!dir.exists(comp_path)) {
        stop_help(opt_p, "Could not find the specified comparison directory ",
                  "(--comp_path)")
      }
      # Gets rid of trailing slashes if there are any
      comp_path <- file.path(comp_path)
    }
  }
}

# Checks the arguments specifying xlims
if (!is.na(opt$xlims_from_file)) {
  if (!file.exists(opt$xlims_from_file)) {
    stop_help(opt_p, "Could not find xlims file (--xlims_from_file)")
  }
  if (verbose) {
    cat_time("xlim file detected. Will plot positions as defined by: ",
             opt$xlims_from_file)
  }
  xlim_option <- "infile"
  region_size <- NULL
} else if (!is.na(opt$xlims_in)) {
  xlim_option <- "input"
  region_size <- NULL
} else {
  region_size <- suppressWarnings(as.numeric(opt$region_size))
  if (is.na(region_size)) {
    stop_help(opt_p, 'The argument "', opt$region_size,'" for the ',
              '--region_size option could not be recognized as a numeric value')
  }
  region_size <- round(region_size)
  if (region_size > 0) {
    xlim_option <- "length"
  } else {
    xlim_option <- "none"
  }
}

# Check dna_seg mode
seg_mode <- opt$dna_seg_mode
if (!any(seg_mode == c("all", "ids", "fast"))) {
  stop_help(opt_p, seg_mode, " is not a valid argument for --dna_seg_mode")
} else if (seg_mode == "fast") {
  if (xlim_option != "length" | !any(comp_format == c("blastp-fast",
                                                      "blastp-short",
                                                      "diamond",
                                                      "blastp"))
  ) {
    stop_help(opt_p, paste0('--dna_seg_mode can only be set to "fast" when ',
                            'making a regional plot using blastp or DIAMOND'))
  }
}

# Check update_positions and update_region_plot
if (comp_format != "none") {
  # Check update positions
  update_pos <- tolower(opt$update_positions)
  if (any(update_pos == c("yes", "y", "true"))) {
    update_pos <- TRUE
  } else if (any(update_pos == c("no", "n", "false"))) {
    update_pos <- FALSE
  } else if (update_pos == "auto") {
    if (any(comp_format == c("diamond", "blastp", 
                             "blastp-fast", "blastp-short"))) {
      update_pos <- TRUE
    } else {
      update_pos <- FALSE
    }
  } else {
    stop_help(opt_p, 
              '--update_positions must be one of: "auto", "true", "false"'
              )
  }
  
  # Check update_region_plot
  if (xlim_option == "length") {
    update_reg <- tolower(opt$update_region_plot)
    if (any(update_reg == c("yes", "y", "true"))) {
      update_reg <- TRUE
    } else if (any(update_reg == c("no", "n", "false"))) {
      update_reg <- FALSE
    } else if (update_reg == "auto") {
      if (any(comp_format == c("diamond", "blastp", 
                               "blastp-fast", "blastp-short"))) {
        update_reg <- TRUE
      } else {
        update_reg <- FALSE
      }
    } else {
      stop_help(opt_p,
                '--update_region_plot must be one of: "auto", "true", "false"')
    }
  } else {
    update_reg <- FALSE
  }
} else {
  update_pos <- FALSE
  update_reg <- FALSE
}

# Check out
outfile <- opt$out
if (verbose) {
  if (file.exists(outfile)) {
    cat_time("PDF output file: ",outfile," already exists. ",
                    "File will be overwritten.")
  }
}


# Check tree
if (!is.na(opt$tree)) {
  tree_input <- opt$tree
  if (!file.exists(tree_input)) {
    stop_help(opt_p, "Could not find the specified tree file (--tree)")
  }
  use_tree <- TRUE
} else {
  use_tree <- FALSE
}

# Check plotting dimension parameters
if (opt$height == "auto") {
  if (verbose) {
    cat_time("No plot height argument was provided. It will be ",
             "calculated based on the plotted data.")
  }
  height <- opt$height
} else {
  if (suppressWarnings(is.na(as.double(opt$height)))) {
    stop_help(opt_p, 'The argument "', opt$height, '" for the --height option ',
              'could not be recognized as a numeric value'
              )
  }
  height <- as.double(opt$height)
}
if (opt$width == "auto") {
  if (verbose) {
    cat_time("No plot width argument was provided. It will be ",
             "calculated based on the plotted data.")
  }
  width <- opt$width
} else {
  if (suppressWarnings(is.na(as.double(opt$width)))) {
    stop_help(opt_p, 'The argument "', opt$width, '" for the --width option ',
              'could not be recognized as a numeric value')
  }
  width = as.double(opt$width)
}

# Check global_color_scheme
if (!is.na(opt$global_color_scheme)) {
  global_color_scheme <- opt$global_color_scheme
  if (!any(global_color_scheme == c("gradient", "uniform", "sequential"))) {
    stop_help(opt_p,
              '--global_color_scheme must be one of: ',
              '"gradient", "uniform", "sequential"'
              )
  }
} else {
  global_color_scheme <- NULL
}

# Check color_scheme_dataset
if (any(global_color_scheme == c("uniform", "sequential"))) {
  color_scheme_dataset <- opt$color_scheme_dataset
  if (!any(color_scheme_dataset == c("auto", "dna_segs", "comparisons"))) {
    stop_help(opt_p,
              '--color_scheme_dataset must be one of: ',
              '"auto", "dna_segs", "comparisons"'
              )
  }
} else {
  color_scheme_dataset <- "auto"
}

# Check alpha_dna_segs
if (!is.na(opt$alpha_dna_segs)) {
  alpha_dna_segs <- opt$alpha_dna_segs
  if (suppressWarnings(is.na(as.double(alpha_dna_segs)))) {
    stop_help(opt_p,
              'The argument "', alpha_dna_segs, '" for the --alpha_dna_segs ',
              'option could not be recognized as a numeric value'
              )
  }
} else {
  alpha_dna_segs <- NULL
}

# Check alpha_comparisons
if (!is.na(opt$alpha_comparisons)) {
  alpha_comparisons <- opt$alpha_comparisons
  if (suppressWarnings(is.na(as.double(alpha_comparisons)))) {
    stop_help(opt_p,
              'The argument "', alpha_comparisons, 
              '" for the --alpha_comparisons option ',
              'could not be recognized as a numeric value'
              )
  }
} else {
  alpha_comparisons <- NULL
}

# Check the arguments specifying offsets
if (!is.na(opt$offsets_from_file)) {
  if (!file.exists(opt$offsets_from_file)) {
    stop_help(opt_p, "Could not find offsets file (--offsets_from_file)")
  }
  offset_input <- "infile"
} else if (!is.na(opt$offsets_in)) {
  offset_input <- "input"
} else {
  offset_input <- "none"
}

# Check xlims out
outfile_xlims <- opt$xlims_out
if (!is.na(outfile_xlims) & verbose) {
  if (file.exists(outfile_xlims)) {
    cat_time("xlim output file: ",outfile_xlims," already exists. ",
             "File will be overwritten.")
  } else {
    cat_time("Plotted positions will be written to file: ", outfile_xlims)
  }
} else {
  outfile_xlims <- NULL
}

# Check legend_column
if (!is.na(opt$legend_column)) {
  legend_column <- opt$legend_column
} else {
  legend_column <- NULL
}

# Check threads
threads <- suppressWarnings(as.numeric(opt$threads))
if (is.na(threads)) {
  stop_help(opt_p,
            'The argument "', opt$threads,'" for the --threads option ',
            'could not be recognized as a numeric value')
}
threads <- round(threads)

# Load necessary packages
library(genoPlotR)
library(data.table)

### DEBUG
# library(devtools)
# load_all(".")

################################################################################
# Parsing/ preparing dna_segs
################################################################################

# Read in dna_segs if RDS is provided
if (read_in_RDS) {
  dna_segs <- readRDS(dna_seg_files)
  seg_labels <- get_seg_labels(dna_segs, seg_labels = NULL, required = TRUE)
} else {
  seg_labels <- gsub("\\.[^\\.]*$", "", basename(dna_seg_files))
}
old_labels <- seg_labels

# Depending on seg_mode, filter dna_segs based on seg_labels from ids
region_matched_labels <- character()
if (use_ids) {
  ids <- fread( opt$ids, header = TRUE, sep = "auto", fill = TRUE)
  if (any("seg_label" == names(ids))) {
    ids_labels <- unique(ids[, seg_label])
    ids_labels <- ids_labels[ids_labels != ""]
    if (!is.null(ids$region_plot)) {
      ids$region_plot <- as.character(ids$region_plot)
      ids$region_plot[ids$region_plot %ilike% "^true$|^plot$|^y$|^yes$"] <- "TRUE"
      ids$region_plot[!ids$region_plot %ilike% "^true$|^plot$|^y$|^yes$|^start$|^end$"] <- "NA"
      region_labels <- unique(ids[region_plot == "TRUE", seg_label])
      region_labels <- region_labels[region_labels != ""]
    }
    matched_labels <- character()
    for (i in ids_labels) {
      query <- grep(i, seg_labels)
      if (length(query) == 0) {
        if (read_in_RDS) {
          query <- grep(gsub("\\.[^\\.]*$", "", i), seg_labels)
        } else {
          query <- grep(i, basename(dna_seg_files))
        }
      }
      if (length(query) == 0) {
        stop_help(opt_p, 'seg_label "', i, '" could not be matched to any ',
                  'file from --seg_files')
      } else if (length(query) > 1) {
        # If there are multiple matches, look for an exact match instead
        query <- grep(paste0("\\<", i, "\\>"), seg_labels)
        if (length(query) == 0 & !read_in_RDS) {
          query <- grep(paste0("\\<", i, "\\>"), basename(dna_seg_files))
        }
        if (length(query) != 1) {
          stop_help(opt_p,
                    'Multiple matches found for seg_label "', i, '" in ',
                    'the dna_seg file names, make sure the seg_labels match ',
                    'the filenames exactly'
                    )
        }
      }
      query <- as.numeric(query)
      matched_labels <- c(matched_labels, seg_labels[query])
      if (any(i == region_labels)) {
        region_matched_labels <- c(region_matched_labels, seg_labels[query])
      }
    }
    match_order <- sapply(matched_labels,
                          function(x) which(x == seg_labels),
                          USE.NAMES = FALSE
                          )
    if (seg_mode != "ids") {
      match_order <- c(match_order, seq(1:length(seg_labels))[-match_order])
    }
    seg_labels <- seg_labels[match_order]
    old_labels <- old_labels[match_order]
    if (read_in_RDS) {
      dna_segs <- dna_segs[seg_labels]
    } else {
      dna_seg_files <- dna_seg_files[match_order]
    }
  }
}

# Read in tree, reorder dna_segs based on tree
if (use_tree) {
  if (verbose) cat_time("Reading tree from: ", tree_input)
  tree_read <- ape::read.tree(tree_input)
  
  # Use the tree to reorder dna_segs, starting with labels
  permuted <- permute_dna_segs(seg_labels, tree_read, return_old_labels = TRUE)
  to_keep <- sapply(permuted$old_labels,
                    function(x) which(x == seg_labels),
                    USE.NAMES = FALSE
                    )
  seg_labels <- permuted$dna_segs
  old_labels <- permuted$old_labels

  # To reorder dna_segs we use indices given by the return_old_labels argument
  if (read_in_RDS) {
    dna_segs <- dna_segs[to_keep]
    names(dna_segs) <- seg_labels
  } else {
    dna_seg_files <- dna_seg_files[to_keep]
  }
  
} else {
  if (verbose) {
    cat_time("No tree given; will plot without a tree.")
  }
  tree_read <- NULL
}

# Determine xlims, potentially reordering/ filtering dna_segs
if (xlim_option == "infile") {
  xlim_data <- read_xlims_file(seg_labels, opt$xlims_from_file,
                               reorder_dna_segs = TRUE)
  xlims <- xlim_data$xlims
  to_keep <- sapply(xlim_data$seg_labels,
                    function(x) which(x == seg_labels),
                    USE.NAMES = FALSE
                    )
  seg_labels <- xlim_data$seg_labels
  old_labels <- old_labels[to_keep]
  if (read_in_RDS) {
    dna_segs <- dna_segs[to_keep]
  } else {
    dna_seg_files <- dna_seg_files[to_keep]
  }
} else if (xlim_option == "input") {
  xlims <- read_xlims(opt$xlims_in)
} else {
  xlims <- NULL
}
if (verbose) cat_time("Reading in dna_segs and potentially comparisons.")

################################################################################
# Making and/or parsing comparisons alongside dna_segs
################################################################################

if (seg_mode == "fast") {
  if (!read_in_RDS) {
    dna_segs <- list()
    unordered_dna_segs <- list()
    if (!is.na(outfile_dna_segs)) {
      raw_dna_segs <- list()
    }
    comparisons <- list()
    seg_order <- 1:length(seg_labels)
    if (use_tree) {
      seg <- read_dna_seg_from_file(
        dna_seg_files[1],
        read_sequence = TRUE,
        verbose = verbose,
        extra_fields = "note",
        tagsToParse = c("CDS", "tRNA", "rRNA", "repeat_region"),
        boundariesToParse = c("contig", "source", "chromosome")
      )
      if (!is.na(outfile_dna_segs)) raw_dna_segs[[seg_labels[1]]] <- seg
      if (use_ids) {
        seg <- edit_dna_segs(seg, ids, id_tags = id_tags,
                             seg_labels = seg_labels[1])
      }
      if (!any("TRUE" == unique(seg$region_plot))) {
        stop_help(opt_p, 'Region plotting was selected but no features were ',
                  'marked to be plotted in "', basename(dna_seg_files[1]), 
                  '", the first dna_seg in the provided tree. Use the --ids ',
                  'option to mark features to be plotted using the ',
                  '"region_plot" column, or turn off region plotting using ',
                  '--region_size, --xlims_in, or --xlims_from_file.'
                  )
      }
      dna_segs[[seg_labels[1]]] <- seg
    } else if (!use_tree) {
      start_found <- FALSE
      for (i in seg_order) {
        seg <- read_dna_seg_from_file(
          dna_seg_files[i],
          read_sequence = TRUE,
          verbose = verbose,
          extra_fields = "note",
          tagsToParse = c("CDS", "tRNA", "rRNA", "repeat_region"),
          boundariesToParse = c("contig", "source", "chromosome")
        )
        if (!is.na(outfile_dna_segs)) raw_dna_segs[[seg_labels[i]]] <- seg
        if (use_ids) {
          seg <- edit_dna_segs(seg, ids, id_tags = id_tags,
                               seg_labels = seg_labels[i])
        }
        
        unordered_dna_segs[[seg_labels[i]]] <- seg
        if (any("TRUE" == unique(seg$region_plot))) {
          dna_segs[[seg_labels[i]]] <- seg
          seg_order <- c(i, seg_order[-which(i == seg_order)])
          dna_seg_files <- dna_seg_files[seg_order]
          seg_labels <- seg_labels[seg_order]
          old_labels <- old_labels[seg_order]
          start_found <- TRUE
          break
        }
      }
      if (!start_found) {
        stop_help(opt_p, 'Region plotting was selected but no features in any ',
                  'dna_seg were marked to be plotted. Use the --ids option to ',
                  'mark features to be plotted using the "region_plot" ',
                  'column, or turn off region plotting using --region_size, ',
                  '--xlims_in, or --xlims_from_file.')
      }
    }
    for (i in 1:(length(seg_labels)-1)) {
      # Now load in remaining dna_segs and comparisons sequentially
      comp <- comparisons_from_dna_segs(files = c(dna_seg_files[i:(i+1)]),
                                        mode = comp_mode,
                                        tool = comp_format, 
                                        algorithm = comp_format,
                                        output_path = comp_path,
                                        seg_labels = seg_labels[i:(i+1)],
                                        sensitivity = sensitivity,
                                        num_threads = threads,
                                        verbose = TRUE
                                        )
      # Comparison is updated using previous dna_seg to check for region_plot
      comp <- update_comparisons(dna_segs[[seg_labels[i]]],
                                 comp,
                                 color_var = NULL,
                                 update_positions = FALSE,
                                 update_region_plot = TRUE
                                 )[[1]]
      if (any("TRUE" == unique(comp$region_plot)) | 
          any(old_labels[i+1] == region_matched_labels )) {
        # Load in dna_seg when there region_plot is present in this comparison
        # or if we know from the ids file that something should be there
        if (!any(seg_labels[i+1] == names(unordered_dna_segs))) {
          # Only when dna_seg has not been loaded already
          seg <- read_dna_seg_from_file(
            dna_seg_files[i+1],
            read_sequence = TRUE,
            verbose = verbose,
            extra_fields = "note",
            tagsToParse = c("CDS", "tRNA", "rRNA", "repeat_region"),
            boundariesToParse = c("contig", "source", "chromosome")
          )
          if (!is.na(outfile_dna_segs)) raw_dna_segs[[seg_labels[i+1]]] <- seg
          if (use_ids) {
            seg <- edit_dna_segs(seg, ids, id_tags = id_tags,
                                 seg_labels = seg_labels[i+1])
          }
          
          dna_segs[[seg_labels[i+1]]] <- seg
        } else {
          dna_segs[[seg_labels[i+1]]] <- unordered_dna_segs[[seg_labels[i+1]]]
        }
        if (update_reg) {
          dna_segs[[seg_labels[i+1]]] <- update_dna_segs(
            dna_seg_input = dna_segs[[seg_labels[i+1]]],
            comparison_input = comp, color_var = NULL,
            update_region_plot = update_reg
          )
        }
        
        comparisons[[length(comparisons)+1]] <- comp
      } else {
        # region_plot not here, remove dna_seg in variables
        if (!is.na(outfile_dna_segs)) raw_dna_segs[[seg_labels[i+1]]] <- NULL
        seg_labels[i+1] <- seg_labels[i]
        dna_seg_files[i+1] <- dna_seg_files[i]
        
      }
    }
    seg_labels <- names(dna_segs)
    dna_seg_files <- unique(dna_seg_files) # Might not be necessary
  } else {
    if (!is.na(outfile_dna_segs)) raw_dna_segs <- dna_segs
    comparisons <- list()
    new_order <- list()
    if (use_ids) {
      dna_segs <- edit_dna_segs(dna_segs, ids, id_tags = id_tags,
                                seg_labels = seg_labels)
    }
    seg_order <- 1:length(seg_labels)
    if (use_tree) {
      if (!any("TRUE" == unique(dna_segs[[1]]$region_plot))) {
        stop_help(opt_p, 'Region plotting was selected but no features were ',
                  'marked to be plotted in "', seg_labels[1], '", the first ',
                  'dna_seg in the provided tree. Use the --ids option to mark ',
                  'features to be plotted using the "region_plot" column, ',
                  'or turn off region plotting using --region_size, ',
                  '--xlims_in, or --xlims_from_file.')
      }
      new_order[[seg_labels[1]]] <- dna_segs[[1]]
    } else if (!use_tree) {
      start_found <- FALSE
      for (i in seg_order) {
        if (any("TRUE" == unique(dna_segs[[i]]$region_plot))) {
          new_order[[seg_labels[i]]] <- dna_segs[[i]]
          seg_order <- c(i, seg_order[-which(i == seg_order)])
          seg_labels <- seg_labels[seg_order]
          old_labels <- old_labels[seg_order]
          start_found <- TRUE
          break
        }
      }
      if (!start_found) {
        stop_help(opt_p, 'Region plotting was selected but no features in any ',
                  'dna_seg were marked to be plotted. Use the --ids option to ',
                  'mark features to be plotted using the "region_plot" ',
                  'column, or turn off region plotting using --region_size, ',
                  '--xlims_in, or --xlims_from_file.')
      }
    }
    for (i in 1:(length(seg_labels)-1)) {
      # Now load in remaining dna_segs and comparisons sequentially
      comp <- comparisons_from_dna_segs(
        dna_segs = list(new_order[[seg_labels[i]]],
                        dna_segs[[seg_labels[i+1]]]
                        ),
        mode = comp_mode,
        tool = comp_format,
        algorithm = comp_format,
        output_path = comp_path,
        seg_labels = old_labels[i:(i+1)],
        sensitivity = sensitivity,
        num_threads = threads,
        verbose = verbose
      )
      # Comparison is updated using previous dna_seg to check for region_plot
      comp <- update_comparisons(new_order[[seg_labels[i]]],
                                     comp,
                                     color_var = NULL,
                                     update_positions = FALSE,
                                     update_region_plot = TRUE)[[1]]
      if (any("TRUE" == unique(comp$region_plot)) | 
          any(old_labels[i+1] == region_matched_labels )
          ) {
        if (update_reg) {
          new_order[[seg_labels[i+1]]] <- update_dna_segs(
            dna_seg_input = dna_segs[[seg_labels[i+1]]],
            comparison_input = comp, color_var = NULL,
            update_region_plot = update_reg
          )
        } else {
          new_order[[seg_labels[i+1]]] <- dna_segs[[seg_labels[i+1]]]
        }
        
        comparisons[[length(comparisons)+1]] <- comp
      } else {
        # region_plot not here, remove dna_seg in variables with previous
        if (!is.na(outfile_dna_segs)) raw_dna_segs[[seg_labels[i+1]]] <- NULL
        seg_labels[i+1] <- seg_labels[i]
        old_labels[i+1] <- old_labels[i]
      }
    }
    dna_segs <- new_order
    seg_labels <- names(dna_segs)
    
  }
  if (update_reg) {
    full_data <- sequential_updates(dna_segs, comparisons, 
                                    color_var = NULL,
                                    update_region_plot = TRUE,
                                    update_positions = FALSE)
    dna_segs <- full_data$dna_segs
    comparisons <- full_data$comparisons
  }
} else {
  # Fast mode isn't enabled, read in dna_segs normally 
  if (!read_in_RDS) {
    dna_segs <- read_dna_seg_from_files(
      dna_seg_files, read_sequence = TRUE,
      verbose = verbose,
      extra_fields = "note",
      tagsToParse = c("CDS", "tRNA", "rRNA", "repeat_region"),
      boundariesToParse = c("contig", "source", "chromosome")
    )
    names(dna_segs) <- seg_labels
  }
  if (!is.na(outfile_dna_segs)) raw_dna_segs <- dna_segs
  
  ### Write out dna_segs for fast reuse
  ### DEBUG
  # if (!is.na(outfile_dna_segs)) {
  #   if (verbose) cat_time("Saving raw dna_segs to: ", outfile_dna_segs)
  #   saveRDS(raw_dna_segs, file = paste0(outfile_dna_segs, "_rawer.RDS"))
  # }
  
  # Edit dna_segs with the ids file 
  if (use_ids & !any(comp_format == c("orthofinder", "orthomcl", "mmseqs2"))) {
    dna_segs <- edit_dna_segs(dna_segs, ids, id_tags = id_tags,
                              seg_labels = seg_labels)
  }
  
  # Start reading in comparisons
  if (comp_format == "none") {
    comparisons <- NULL
    if (xlim_option == "length") {
      # Remove dna_segs that don't have region_plot features
      plot_seg <- sapply(dna_segs,
                         function(x) any("TRUE" == unique(x$region_plot)),
                         USE.NAMES = FALSE
                         )
      to_keep <- which(plot_seg)
      dna_segs <- dna_segs[to_keep]
      seg_labels <- seg_labels[to_keep]
    }
  } else {
    reading_comps <- TRUE
    while (reading_comps) {
      # Reading comps loop
      if (any(comp_format == c("diamond", "blast", "blastp", "blastp-fast", 
                               "blastp-short", "blastn", "blastn-short", 
                               "megablast", "dc-megablast", "tblastx"
                               ))) {
        if (read_in_RDS) {
          comparisons <- comparisons_from_dna_segs(dna_segs = dna_segs,
                                                   mode = comp_mode,
                                                   tool = comp_format,
                                                   algorithm = comp_format,
                                                   output_path = comp_path,
                                                   seg_labels = old_labels,
                                                   sensitivity = sensitivity,
                                                   num_threads = threads,
                                                   verbose = verbose
                                                   )
        } else {
          comparisons <- comparisons_from_dna_segs(files = dna_seg_files,
                                                   mode = comp_mode,
                                                   tool = comp_format,
                                                   algorithm = comp_format,
                                                   output_path = comp_path,
                                                   seg_labels = seg_labels,
                                                   sensitivity = sensitivity,
                                                   num_threads = threads,
                                                   verbose = verbose
                                                   )
        }
      } else if (comp_format == c("tab")) {
        comparisons <- read_comparison_from_files(
          files = file.path(comp_path, "*"),
          fileType = "tab",
          seg_labels = seg_labels
        )
      } else if (any(comp_format == c("orthofinder", "orthomcl", "mmseqs2"))) {
        full_data <- read_orthogroup_from_file(file = comp_path,
                                               dna_segs = dna_segs,
                                               fileType = comp_format,
                                               alter_dna_segs = TRUE
                                               )
        dna_segs <- full_data$dna_segs
        comparisons <- full_data$comparisons
        
        # Edit dna_segs with the ids file now because dna_segs were altered
        if (use_ids) {
          dna_segs <- edit_dna_segs(dna_segs, ids, id_tags = id_tags,
                                    seg_labels = seg_labels)
        }
      }
      
      # Done reading comps in this iteration, checking for region_plot
      if (xlim_option != "length") {
        # If we are not going for a region_plot, stop the comparison loop
        reading_comps <- FALSE
      } else {
        if (update_reg) {
          full_data <- sequential_updates(dna_segs, comparisons, 
                                          color_var = NULL,
                                          update_region_plot = TRUE,
                                          update_positions = FALSE)
          plot_seg <- sapply(full_data$dna_segs,
                             function(x) any("TRUE" == unique(x$region_plot)),
                             USE.NAMES = FALSE
                             )
          if (!all(plot_seg)) {
            # Remove first dna_seg without region_plot features, then try again
            to_rm <- which(!plot_seg)[1]
            dna_segs <- dna_segs[-to_rm]
            dna_seg_files <- dna_seg_files[-to_rm]
            seg_labels <- seg_labels[-to_rm]
            old_labels <- old_labels[-to_rm]
            if (!is.na(outfile_dna_segs)) raw_dna_segs <- raw_dna_segs[-to_rm]
            if (length(dna_segs) < 2) {
              stop_help(opt_p, 'Region plotting was selected but no features ',
                        'in any dna_seg were marked to be plotted. Use the ',
                        '--ids option to mark features to be plotted using ',
                        'the "region_plot" column, or turn off region ',
                        'plotting using --region_size, --xlims_in, or ',
                        '--xlims_from_file.')
            }
          } else {
            dna_segs <- full_data$dna_segs
            comparisons <- full_data$comparisons
            reading_comps <- FALSE
          }
        } else if (reading_comps) {
          plot_seg <- sapply(dna_segs,
                             function(x) any("TRUE" == unique(x$region_plot)),
                             USE.NAMES = FALSE
                             )
          if (!all(plot_seg)) {
            # In this case, remove all dna_segs without region_plot features
            to_rm <- which(!plot_seg)
            dna_segs <- dna_segs[-to_rm]
            dna_seg_files <- dna_seg_files[-to_rm]
            seg_labels <- seg_labels[-to_rm]
            old_labels <- old_labels[-to_rm]
            if (!is.na(outfile_dna_segs)) raw_dna_segs <- raw_dna_segs[-to_rm]
          } else {
            reading_comps <- FALSE
          }
        }
      }
    }
  }
}


################################################################################
# Final preparations before calling plot_gene_map
################################################################################

if (verbose) {
  cat_time("Running genoPlotR with ", length(dna_segs), " genomes.")
}


# Update comparison positions to match dna_segs
if (comp_format != "none" & update_pos) {
  comparisons <- update_comparisons(dna_segs,
                                    comparisons,
                                    color_var = NULL,
                                    update_region_plot = FALSE,
                                    update_positions = TRUE
                                    )
}

# Trim tree now that dna_segs is done
if (use_tree) {
  # Trim unused leaves and convert it to a format plot_gene_map understands
  tree_read <- trim_tree(seg_labels, tree_read)
}

### DEBUG
# saveRDS(dna_segs, file.path(dirname(outfile), paste0(substr(basename(outfile), 1, nchar(basename(outfile))-4), "_segs.RDS")))
# saveRDS(comparisons, file.path(dirname(outfile), paste0(substr(basename(outfile), 1, nchar(basename(outfile))-4), "_comps.RDS")))

# Determine offsets
if (offset_input == "infile") {
  offsets <- read_offsets_file(opt$offsets_from_file)
} else if (offset_input == "input") {
  offsets <- read_offsets(opt$offsets_in)
} else {
  offsets <- NULL
}

# Create annotations
if (use_annots) {
  annots <- auto_annotate(dna_segs, rot=60, keep_genes_only = FALSE)
} else {
  annots <- NULL
}

# Write out dna_segs for fast reuse
if (!is.na(outfile_dna_segs)) {
  if (verbose) cat_time("Saving raw dna_segs to: ", outfile_dna_segs)
  saveRDS(raw_dna_segs, file = outfile_dna_segs)
}

# Write out dna_segs for fast reuse
if (!is.na(workspace_out)) {
  if (verbose) {
    cat_time("Will save all plot arguments (workspace) as a single ",
             "R object to: ", workspace_out)
  }
  saveRDS(list(dna_segs = dna_segs,
               comparisons = comparisons,
               tree = tree_read,
               annotations = annots,
               legend_column = legend_column,
               region_size = region_size,
               xlims = xlims,
               print_xlims = print_xlims,
               outfile_xlims = outfile_xlims,
               offsets = offsets,
               print_offsets = print_offsets,
               seg_labels = seg_labels,
               outfile = outfile,
               outfile_height = height,
               outfile_width = width,
               global_color_scheme = global_color_scheme,
               color_scheme_dataset = color_scheme_dataset,
               alpha_dna_segs = alpha_dna_segs,
               alpha_comparisons = alpha_comparisons,
               dna_seg_scale = dna_seg_scale,
               verbose = verbose
               ),
          file = workspace_out
          )
}



# Plotting
if (verbose) cat_time("Plotting to: ", outfile)


'
                          reorder_xlims=FALSE, # unimplemented
                          gene_type=NULL,      # if not null, resets gene_type
                          arrow_head_len=200,  # force arrow head length
                          dna_seg_line=TRUE,   # draw a line on each dna_seg
                          scale=TRUE,          # scale in the bottom right
                          dna_seg_scale=FALSE, # scale on each dna_seg
                          n_scale_ticks=7,     # number of tick marks for these
                          scale_cex=0.6,       # size of text on scale
'

plot_gene_map(
  dna_segs = dna_segs,
  comparisons = comparisons,
  tree = tree_read,
  annotations = annots,
  legend_column = legend_column,
  # legend_labels = notes,
  region_size = region_size,
  xlims = xlims,
  print_xlims = print_xlims,
  outfile_xlims = outfile_xlims,
  offsets = offsets,
  print_offsets = print_offsets,
  seg_labels = seg_labels,
  outfile = outfile,
  outfile_height = height,
  outfile_width = width,
  global_color_scheme = global_color_scheme,
  color_scheme_dataset = color_scheme_dataset,
  # color_scheme_column = "gene",
  alpha_dna_segs = alpha_dna_segs,
  alpha_comparisons = alpha_comparisons,
  dna_seg_scale = dna_seg_scale,
  verbose = verbose
)

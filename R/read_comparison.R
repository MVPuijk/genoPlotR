################################################################################
# File reading functions: read comparison
################################################################################

#' Creating comparisons between genomic sequences from files
#' 
#' Functions to parse `comparison` objects from tabular data or BLAST output.
#' 
#' @details
#' To make `comparisons` from tabular files, the columns in these files must 
#' include `start1`, `end1`, `start2`, and `end2`. If no header is specified,
#' these columns will be expected in that order, and an optional fifth column
#' will be parsed as the color of the comparison under the `col` column.
#' Additional columns can be included regardless of whether a header was 
#' specified.
#' 
#' BLAST output must be provided in the default tabular format (`-outfmt 6`),
#' which contains the following columns in this order:
#' `qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`,
#' `qstart`, `qend`, `sstart`, `send`, `evalue`, and `bitscore`.
#' 
#' The resulting `comparison` object will contain these columns under these
#' names: `name1`, `name2`, `per_id`, `aln_len`, `mism`, `gaps`, `start1`,
#' `start2`,`end1`, `end2`, `e_value`, and `bit_score`.
#' 
#' @name read_comparison
#' @returns A list of `comparison` objects for `read_comparison_from_files`, 
#' and a single `comparison` object otherwise.
#' @export
#' 
#' @param file A character string containing a file path, or a file connection.
#' @param fileType A character string containing the file format to parse.
#' Must be one of: `"blast"`, or `"tab"`.
#' @param ... Further arguments to pass to `read_comparison_from_blast` or 
#' `read_comparison_from_tab`, see arguments below.
#' 
#' @author Lionel Guy, Jens Roat Kultima, Mike Puijk
#' 
#' @seealso [comparison], [run_sequence_alignment]
#' 
#' @examples
#' 
#' ## Read comparison from tabular data
#' tab_file <- system.file('extdata/comparison2.tab', package = 'genoPlotR')
#' tab_comp <- read_comparison_from_tab(tab_file)
#' data.frame(tab_comp)
#' 
#' ## Load dna_segs
#' data(barto)
#' bh <- barto$dna_segs[[3]]
#' bq <- barto$dna_segs[[4]]
#' 
#' ## Read comparison from BLAST output
#' blast <- system.file('extdata/BH_vs_BQ.blastn.tab', package = 'genoPlotR')
#' blast_comp <- read_comparison_from_blast(blast, color_scheme = "red_blue")
#' 
#' ## Plot
#' plot_gene_map(dna_segs = list(BH = bh, BQ = bq),
#'               comparisons = list(blast_comp),
#'               xlims = list(c(1,50000), c(1, 50000)))
#' 
read_comparison_from_file <- function(file, fileType, ...) {
  
  if (tolower(fileType) == "tab") {
    comparison <- read_comparison_from_tab(file, ...)
  } else if (tolower(fileType) == "blast") {
    comparison <- read_comparison_from_blast(file, ...)
  } else {
    stop("'fileType' has to be either: blast or tab")
  }
  return(comparison)
}

#' @name read_comparison
#' @export
#' 
#' @param files A list or character vector containing file paths. Supports
#' wildcard expansion (e.g. *.txt).
#' @param seg_labels A character vector containing `dna_seg` labels. When
#' provided, this function will search for file names that match these `dna_seg`
#' labels. For example, with `seg_labels = c("seg1", "seg2", "seg3")`, it will
#' look for `"seg1_seg2"` and `"seg2_seg3"` among the file names provided by the
#' `files` argument, parse those files, and ignore all other files provided.
#' 
# Added for support of reading many files at once, including globs
read_comparison_from_files <- function(
  files,
  fileType,
  seg_labels = NULL,
  ...
) {
  
  filenames <- c()
  for (i in files) {
    filenames <- c(filenames, Sys.glob(i))
  }
  if (! length(filenames)) {
    stop("Could not find any files at the location(s) specified.")
  }
  comparisons <- list()
  
  if (is.null(seg_labels)) {
    for (i in 1:length(filenames)) {
      file <- filenames[i]
      comparisons[[basename(file)]] <- read_comparison_from_file(file,
                                                                 fileType,
                                                                 ...)
    }
  } else {
    for (i in 1:(length(seg_labels) -1)) {
      file <- grep(paste0(seg_labels[[i]], "_", seg_labels[[i+1]]),
                   filenames, value = TRUE)[1]
      if (file.exists(file)) {
        comparisons[[basename(file)]] <- read_comparison_from_file(file,
                                                                   fileType,
                                                                   ...)
      } else {
        stop('Could not find a comparison file for ', seg_labels[[i]],
             ' and ', seg_labels[[i+1]], ', using the search string "',
             seg_labels[[i]], "_", seg_labels[[i+1]], '"')
      }
    }
  }
  comparisons
}

#' @name read_comparison
#' @export
#' 
#' @param sort_by A character string containing the name of a column to sort by.
#' Must be a numeric column.
#' @param filt_high_evalue A numerical, filters out all comparisons with an 
#' e-value higher than this value (unfiltered when left as `NULL`).
#' @param filt_low_per_id A numerical, filters out all comparisons with a
#' percentage identity lower than this value (unfiltered when left as `NULL`).
#' @param filt_length A numerical, filters out all comparisons that have
#' alignments shorter than this value (unfiltered when left as `NULL`).
#' @param color_scheme A color scheme to apply. Possible values include `grey`,
#' `red_blue`, and `NULL` (which applies no color scheme). See 
#' [gradient_color_scheme] for more details.
#' @param reverse When provided, the `comparison` will be reversed.
#' If `reverse = 1`, the first side will be reversed. 
#' If `reverse = 2`, the second side will be reversed.
#' If `reverse < 1`, no side is reversed.
#' If `reverse > 2`, both sides are reversed.
#' 
read_comparison_from_blast <- function(
  file,
  sort_by = NULL,
  filt_high_evalue = NULL,
  filt_low_per_id = NULL,
  filt_length = NULL,
  color_scheme = NULL,
  reverse = 0
) {
  # read table
  table <- read.table(file, as.is = TRUE, header = FALSE,
                      sep = "\t", quote = "")
  # from blast output
  col_names <- c("name1", "name2", "per_id", "aln_len", "mism", "gaps",
                 "start1", "end1", "start2", "end2", "e_value", "bit_score")
  names(table) <- col_names
  # sort
  if (!is.null(sort_by)) {
    if (!sort_by %in% col_names) {
      stop("'sort_by' must refer to a column from ",
           "this set:\n", paste(col_names, collapse = ", "))
    }
    if (!is.numeric(table[[sort_by]])) {
      stop("'sort_by' must designate a numeric column")
    }
    decr <- if (sort_by %in% c("mism", "gaps", "e_value")) FALSE else TRUE
    # reorder from weakest to strongest
    table <- table[do.call(order, list(table[[sort_by]], decreasing = decr)), ]
  }
  # filter evalue
  if (!is.null(filt_high_evalue)) {
    if (!is.numeric(filt_high_evalue)) {
      stop("'filt_high_evalue' must be numeric")
    }
    table <- table[table$e_value <= filt_high_evalue, ]
  }
  # filter per id
  if (!is.null(filt_low_per_id)) {
    if (!is.numeric(filt_low_per_id)) stop("'filt_low_per_id' must be numeric")
    table <- table[table$per_id >= filt_low_per_id,]
  }
  # reorder columns
  table <- table[,c(col_names[c(7:10,1:6,11:12)])]
  
  # send further
  comparison <- .read_comparison(table, filt_length = filt_length)
  # apply color scheme if requested
  if (!is.null(color_scheme)) {
    comparison$col <- gradient_color_scheme(comparison$per_id,
                                            direction = comparison$direction,
                                            color_scheme = color_scheme
                                            )
  }
  comparison
}

#' @name read_comparison
#' @export
#' 
#' @param header Logical. If `TRUE`, the first line will be parsed as a header
#' containing column names.
#' 
read_comparison_from_tab <- function(
  file,
  header = TRUE,
  filt_length = NULL,
  reverse = 0
) {
  
  col_names <-  c("start1", "end1", "start2", "end2", "col")
  # from tabular data
  table <- read.table(file, as.is = TRUE, header = header, 
                      sep = "\t", quote = "")
  if (ncol(table) < 4) stop("Insufficient number of columns in table")
  if (!header) {
    if (ncol(table) == 4) {
      names(table) <- col_names[1:4]
    } else {
      names(table)[1:length(col_names)] <- col_names[1:length(col_names)]
    }
  }
  
  .read_comparison(table)
}

#' Create comparisons between dna_segs by reading in a file of groupings
#' 
#' Functions to create `comparisons` between `dna_seg` objects by parsing 
#' (ortho)groups from a file. 
#' 
#' @details
#' `read_orthogroup_from_orthomcl`, `read_orthogroup_from_orthofinder`,
#' `read_orthogroup_from_mmseqs2`, and `read_orthogroup_from_diamond` are
#' all just convenience functions for `read_orthogroup_from_file`.
#' 
#' This function was created to create a list of `comparisons` from a list of
#' `dna_segs` and a file that contains orthologous groups of genes
#' (orthogroups). However, it could theoretically be used for any grouping of
#' (genetic) elements on a genomic track. For instance, a group could represent
#' an operon, pathway, or general function. This function creates the
#' `comparisons` by linking together columns from `dna_segs` (specified by 
#' `id_tag`). 
#' 
#' Because `"orthofinder"` format files contain columns representing the
#' different genomes as input, this function will attempt to 
#' query only the `dna_seg` whose label matches the column. `dna_seg` labels 
#' will be determined automatically, but they can also be provided using the
#' `seg_labels` argument. If the labels cannot be matched, it will continue
#' without matching `dna_seg` names, querying each `dna_seg` for each column.
#' 
#' @name read_orthogroup
#' @returns With `alter_dna_segs = TRUE`, a list with 2
#' named elements: `dna_segs` and `comparisons`, which are both lists containing
#' the `dna_seg` and `comparison` objects, respectively.
#' @returns With `alter_dna_segs = FALSE`, a list of `comparison` objects.
#' @export
#' 
#' @param file A character string containing a file path.
#' @param dna_segs A list of `dna_seg` objects.
#' @param fileType A character string containing the file format to parse.
#' Must be one of: `"orthomcl"`, `"orthofinder"`, `"mmseqs2"`, or `"diamond"`.
#' @param seg_labels Only for use with `"orthofinder"` format files, a character
#' string of `dna_seg` labels. See details.
#' @param id_tag A character string with a `dna_seg` column. The (gene) names
#' taken from the (ortho)groups file will be matched to the names found in this
#' `dna_seg` column.
#' @param group_name A character string containing a column name. This column
#' will contain the group names found in the (ortho)groups file and will be
#' added to the `comparisons`, as well as the `dna_segs` when `alter_dna_segs`
#' is set to `TRUE`.
#' @param alter_dna_segs Logical. If `TRUE`, a group column will be added to
#' each `dna_seg` containing the groupings found in the (ortho)groups file.
#' @param verbose Logical. If `TRUE`, will report a warning when the column
#' specified by `group_name` is already present in the `dna_segs` and will 
#' therefore be overwritten. Has no effect unless `alter_dna_segs` is set to 
#' `TRUE`.
#' @param ... Arguments to pass to [fread] and `read_orthogroup_from_file`.
#' 
#' @author Mike Puijk
#' 
#' @seealso [comparison], [read_comparison_from_file]
#' 
#' @examples
#' ## Generate data
#' names1 <- c("1_FeatA1", "1_FeatA2", "1_FeatB")
#' names2 <- c("2_FeatA", "2_FeatB")
#' names3 <- c("3_FeatA", "3_FeatB")
#' df1 <- data.frame(name = names1, start = c(1, 501, 1501),
#'                   end = c(400, 900, 2200), strand = c(1, 1, 1))
#' df2 <- data.frame(name = names2, start = c(1, 501),
#'                   end = c(400, 1200), strand = c(1, 1))
#' df3 <- data.frame(name = names3, start = c(1, 501),
#'                   end = c(400, 1200), strand = c(1, 1))
#' 
#' ## Create list of dna_segs
#' dna_segs <- list(dna_seg(df1), dna_seg(df2), dna_seg(df3))
#' 
#' ## Read feature groups from OrthoFinder (Orthogroups.tsv) format
#' file <- system.file('extdata/OrthoFinder_format.tsv', package = 'genoPlotR')
#' full_data <- read_orthogroup_from_file(file = file, dna_segs = dna_segs,
#'                                        fileType = "orthofinder", 
#'                                        id_tag = "name")
#' 
#' ## Plot data
#' plot_gene_map(dna_segs = full_data$dna_segs,
#'               comparisons = full_data$comparisons,
#'               global_color_scheme = "uniform",
#'               alpha_comparisons = 0.5)
#' 
#' ## Examples of these groups in the different supported formats:
#' OrthoFinder_file <- system.file('extdata/OrthoFinder_format.tsv',
#'                                 package = 'genoPlotR')
#' OrthoMCL_file <- system.file('extdata/OrthoMCL_format.txt',
#'                              package = 'genoPlotR')
#' MMSeqs2_or_DIAMOND_file <- system.file('extdata/MMseqs2_DIAMOND_format.tsv',
#'                                        package = 'genoPlotR')
#' 
#' cat(readLines(OrthoFinder_file), sep = "\n")
#' cat(readLines(OrthoMCL_file), sep = "\n")
#' cat(readLines(MMSeqs2_or_DIAMOND_file), sep = "\n")
#' 
read_orthogroup_from_file <- function(
  file,
  dna_segs,
  fileType,
  seg_labels = NULL,
  id_tag = "locus_id",
  group_name = "orthogroup",
  alter_dna_segs = TRUE,
  verbose = FALSE,
  ...
) {
  
  # Check mandatory arguments
  if (missing(dna_segs)) stop("'dna_segs' must be provided")
  if (!is.list(dna_segs) || !all(sapply(dna_segs, is.dna_seg))) {
    stop("'dna_segs' must be a list of dna_seg objects")
  }

  fileType <- tolower(fileType)
  if (!any(fileType == c("orthomcl", "orthofinder", "mmseqs2", "diamond"))) {
    stop("'fileType' must be one of: orthomcl, orthofinder, mmseqs2, diamond")
  }
  
  if (fileType == "diamond") fileType <- "mmseqs2"
  
  if (file.exists(file)) {
    if (fileType == "orthomcl") {
      ortho_file <- fread(file, header = FALSE, sep = ":", ...)
    } else if (fileType == "orthofinder") {
      ortho_file <- fread(file, header = TRUE, sep = "\t", ...)
    } else if (fileType == "mmseqs2") {
      melted <- fread(file, header = FALSE, sep = "\t", ...)
    }
  } else {
    stop('Could not find file at "', file ,'".')
  }
  
  if (alter_dna_segs) dna_segs <- copy(dna_segs)
  
  for (i in 1:length(dna_segs)) {
    if (! any(id_tag == names(dna_segs[[i]]))) {
      stop("Not all dna_seg objects contained the column specified by 'id_tag'")
    }
  }
  
  
  # Determining how to query for genomes, pre-preparing ortho_file accordingly
  if (fileType == "orthomcl") {
    use_names <- FALSE
    maxSplits <- max(lengths(strsplit(ortho_file[["V2"]], " ")))
    ortho_file[, 
               paste0("featureID", 1:maxSplits) := tstrsplit(V2, " ", fixed = T)
               ]
    setnames(ortho_file, old = "V1", new = "temp000")
    melted <- melt(ortho_file,
                   id.vars = "temp000",
                   measure.vars = patterns("^featureID"),
                   value.factor = TRUE
                   )
  } else if (fileType == "orthofinder") {
    seg_labels <- get_seg_labels(dna_segs, seg_labels)
    use_names <- if (is.null(seg_labels)) FALSE else TRUE
    setnames(ortho_file, old = "Orthogroup", new = "temp000")
    if (use_names) {
      cols <- names(ortho_file)
      order <- unlist(lapply(seg_labels, function(x) which(x == cols)))
      if (length(order) != length(seg_labels)) {
        warning(as.character(length(seg_labels)-length(order)),
                ' dna_seg objects could not be matched to the genomes in "',
                file, '" by name. Proceeding without name and quering every ',
                'dna_seg for each ortholog.'
                )
        use_names <- FALSE
      }
    }
    # If the dna_seg objects have no names, query all dna_segs
    if (!use_names) {
      # Get column names to be combined and construct a call to combine them
      cols <- names(ortho_file)[2:ncol(ortho_file)]
      comb <- as.call(c(
        list(as.name("paste")),
        lapply(cols, as.name),
        list(sep = quote(", "))
      ))
      # Combine the columns in such a way that they can be split in 1 command
      maxSplits <- max(lengths(strsplit(
        ortho_file[, gsub("(, $|^, )", "", gsub("(, )+", ", ", eval(comb)))],
        ", "
      )))
      ortho_file[, paste0("featureID", 1:maxSplits) := tstrsplit(
        gsub("(, $|^, )", "", gsub("(, )+", ", ", eval(comb))), ", "
      )]
      melted <- melt(ortho_file,
                     id.vars = "temp000",
                     measure.vars = patterns("^featureID"),
                     value.factor = TRUE
                     )
    }
  } else if (fileType == "mmseqs2") {
    # mmseqs2 files have a format that matches that of melted almost exactly already
    use_names <- FALSE
    setnames(melted, new = c("temp000", "value"))
  }
  
  comparisons <- list()
  for (i in 1:(length(dna_segs)-1)) {
    if (use_names) {
      # Obtain matching genome name from orthogroup file
      # Weird var names are chosen to avoid using a real column name,
      # as that can cause issues
      label_000 <- cols[order[i]]
      # Determine the maximum number of genes a grouping can have
      maxSplits <- max(lengths(strsplit(ortho_file[[label_000]], ", ")))
      # Take current genome and split so that there is only 1 gene per column
      splitted <- ortho_file[,
                             c(label_000, "temp000"),
                             with = FALSE
                             ][,
                               paste0("featureID",1:maxSplits) := tstrsplit(
                                 eval(as.name(label_000)), ", ", fixed = T)
                               ]
      # Melt so that there is only 1 gene per row
      melted <- melt(splitted,
                     id.vars = "temp000",
                     measure.vars = patterns("^featureID"),
                     value.factor = TRUE
                     )
    }
    
    # Merge the result with the dna_seg object to get the start and end values
    merged1 <- merge(dna_segs[[i]],
                     melted,
                     by.y = "value",
                     by.x = id_tag,
                     sort = FALSE,
                     allow.cartesian = TRUE
                     )
    merged1[strand == -1, c("start", "end") := .(end, start)]
    # Filter the results to include only the necessary columns
    merged1 <- merged1[, c(id_tag, "start", "end", "temp000"), with = FALSE]
    
    if (use_names) {
      # Now repeat some of the steps with the next genome in the list
      label_000 <- cols[order[i+1]]
      maxSplits <- max(lengths(strsplit(ortho_file[[label_000]], ", ")))
      splitted <- ortho_file[,
                             c(label_000, "temp000"),
                             with = FALSE
                             ][,
                               paste0("featureID", 1:maxSplits) := tstrsplit(
                                 eval(as.name(label_000)), ", ", fixed = T)
                               ]
      melted <- melt(splitted,
                     id.vars = "temp000",
                     measure.vars = patterns("^featureID"),
                     value.factor = TRUE
                     )
    }
    
    merged2 <- merge(dna_segs[[i+1]],
                     melted,
                     by.y = "value",
                     by.x = id_tag,
                     sort = FALSE,
                     allow.cartesian = TRUE
                     )
    merged2[strand == -1, c("start", "end") := .(end, start)]
    merged2 <- merged2[, c(id_tag, "start", "end", "temp000"), with = FALSE]
    
    if (alter_dna_segs) {
      # Add group_name column to dna_segs
      seg_merge <- merge(dna_segs[[i]],
                         merged1[, c(id_tag, "temp000"), with = FALSE],
                         all.x = TRUE,
                         sort = FALSE,
                         allow.cartesian = TRUE
                         )
      setcolorder(seg_merge, neworder = names(dna_segs[[i]]))
      set(seg_merge,
          i = which(is.na(seg_merge[["temp000"]])),
          j = "temp000",
          value = "NA"
          )
      if (any(names(seg_merge) == group_name)) {
        if (verbose) {
          warning('One of the dna_segs already contains a "', 
                  group_name, '" column, so it will be overwritten.')
        }
        seg_merge[[group_name]] <- NULL 
      }
      setnames(seg_merge, old = c("temp000"), new = c(group_name))
      dna_segs[[i]] <- seg_merge
      
      if (i == (length(dna_segs) - 1)) {
        # On the last iteration of the loop, also edit the final dna_seg
        seg_merge <- merge(dna_segs[[i+1]],
                           merged2[, c(id_tag, "temp000"), with = FALSE],
                           all.x = TRUE,
                           sort = FALSE,
                           allow.cartesian = TRUE
                           )
        setcolorder(seg_merge, neworder = names(dna_segs[[i+1]]))
        set(seg_merge,
            i = which(is.na(seg_merge[["temp000"]])),
            j = "temp000",
            value = "NA"
            )
        if (any(names(seg_merge) == group_name)) {
          if (verbose) {
            warning('One of the dna_segs already contains a "', 
                    group_name, '" column, so it will be overwritten.')
          }
          seg_merge[[group_name]] <- NULL 
        }
        setnames(seg_merge, old = c("temp000"), new = c(group_name))
        dna_segs[[i+1]] <- seg_merge
      }
    }
    
    # Merge the 2 results to create a comparison object and then clean it up
    final_merge <- merge(merged1, merged2, by = "temp000",
                         sort = FALSE, allow.cartesian = TRUE)
    setnames(final_merge,
             old = c("start.x", "end.x",
                     "start.y", "end.y",
                     paste0(id_tag, ".x"), paste0(id_tag, ".y"),
                     "temp000"
                     ),
             new = c("start1", "end1",
                     "start2", "end2",
                     paste0(id_tag, "1"), paste0(id_tag, "2"),
                     group_name
                     )
             )
    setcolorder(final_merge, c(
      "start1", "end1",
      "start2", "end2",
      paste0(id_tag, "1"), paste0(id_tag, "2"),
      group_name
    ))
    class(final_merge) <- c("data.table", "data.frame")
    comparisons[[i]] <- as.comparison(final_merge)
  }
  if (alter_dna_segs) {
    return(list(dna_segs = dna_segs, comparisons = comparisons))
  }
  
  comparisons
}

#' @name read_orthogroup
#' @export
#' 
# Added for support, to read an orthoMCL format file directly
read_orthogroup_from_orthomcl <- function(file, dna_segs, ...) {
  read_orthogroup_from_file(file, dna_segs, fileType = "orthomcl", ...)
}

#' @name read_orthogroup
#' @export
#' 
# Added for support, to read an OrthoFinder format file directly
read_orthogroup_from_orthofinder <- function(file, dna_segs, ...) {
  read_orthogroup_from_file(file, dna_segs, fileType = "Orthofinder", ...)
}

#' @name read_orthogroup
#' @export
#' 
# Added for support, to read an mmseqs2 format file directly
read_orthogroup_from_mmseqs2 <- function(file, dna_segs, ...) {
  read_orthogroup_from_file(file, dna_segs, fileType = "mmseqs2", ...)
}

#' @name read_orthogroup
#' @export
#' 
# Added for support, to read a DIAMOND orthogroup format file directly
read_orthogroup_from_diamond <- function(file, dna_segs, ...) {
  read_orthogroup_from_file(file, dna_segs, fileType = "diamond", ...)
}

# Internal function, final checks before making comparison
.read_comparison <- function(table, filt_length = NULL, reverse = 0) {
  # check arguments
  if (ncol(table) < 4) stop("Insufficient number of columns in table")
  if (!all(c("start1", "end1", "start2", "end2") %in% names(table))) {
    stop("Table must contain columns 'start1', 'end1', 'start2', 'end2'")
  }
  
  # filter length
  if (!is.null(filt_length)) {
    if (!is.numeric(filt_length)) stop("'filt_length' must be numeric")
    av_len <- apply(
      cbind(abs(table$end1 - table$start1), abs(table$end2 - table$start2)),
      1,
      mean
    )
    table <- table[av_len >= filt_length,]
  }
  # reverse
  if (is.numeric(reverse) && reverse > 0) {
    table <- reverse.comparison(table, side = reverse)
  }
  # make comparison object
  as.comparison(table)
}

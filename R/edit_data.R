################################################################################
# Edit and update functions
################################################################################


#' Edit dna_seg features
#' 
#' Edit the features from `dna_segs` by supplying a set of IDs, and a 
#' set of new values for features that match these IDs.
#' 
#' @details
#' If `ids` is a character string, it is assumed to be a file path, and the file
#' will be read using [fread]. If not, it has to be a `data.frame`
#' or `data.table` object, with a mandatory `id` column. It will then update
#' the `dna_segs` by querying each value (the IDs) from the `id` 
#' column, updating each matching row. It will look for the IDs in the columns
#' provided by the `id_tags` argument. This can be constrained so that it only
#' looks in a specific `dna_seg` for each ID by including a `seg_label` column
#' in `ids`.
#' 
#' This function can be used to alter `dna_seg` attributes on mass, by providing
#' IDs that match to very general attributes, like their color or the presence
#' of a certain word in their functions. But, it can also be used to modify very
#' specific features by making use of attributes with locus tags or the like.
#' 
#' @returns Either a single `dna_seg` object or a list of `dna_seg` objects,
#' matching the input given using `dna_seg_input`.
#' @export
#' 
#' @param dna_seg_input Either a single `dna_seg` or a list of `dna_seg`
#' objects.
#' @param ids Either a character string (specifying a file path), or a
#' data.frame object. Contains the information on how to edit the `dna_segs`,
#' and must contain at least an `id` column. See details.
#' @param seg_labels Either `NULL` or a character vector of the same length
#' as `dna_seg_input`. If `ids` contains a `seg_label` column, then changes
#' will be made only in the `dna_segs` with the corresponding labels. These 
#' labels will be determined both from the `seg_labels` argument, but also from
#' the `dna_segs` themselves. As such, `seg_labels` can be used to provide an 
#' alternate set of names.
#' @param id_tags A character vector of `dna_seg` column names to match to the
#' `id` column from `ids`.
#' @param fixed Logical. If `TRUE`, values from the `id` column have to match
#' exactly to the values found in the `dna_segs`. If `FALSE`, [grep] is used 
#' to search instead, allowing for regular expressions to be used.
#' @param verbose Logical. If `TRUE`, generates warnings when no `dna_seg` could
#' be found for the labels found in the `seg_label` column from `ids`.
#' Additionally, generates a warning when 
#' the columns provided by `id_tags` could not be found in the `dna_segs`.
#' @param ... Arguments to pass to [fread], which is used when the `ids` 
#' argument refers to a file.
#' 
#' @author Mike Puijk
#' 
#' @seealso [dna_seg]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("1A", "1B", "1C")
#' names2 <- c("2A", "2C", "2B")
#' names3 <- c("3B", "3A", "3C")
#' 
#' ## Make dna_segs
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
#' dna_segs <- list("Genome 1" = dna_seg1,
#'                  "Genome 2" = dna_seg2,
#'                  "Genome 3" = dna_seg3)
#' 
#' ## Colors before using edit_dna_segs
#' lapply(dna_segs, function(x) x[, .(name, fill)])
#' 
#' ## Add colors based on exact feature names
#' id_fixed <- c("1A", "1B", "2A", "2B")
#' fill_fixed <- c("red", "blue", "red", "blue")
#' dna_segs1 <- edit_dna_segs(dna_seg_input = dna_segs,
#'                            ids = data.frame(id = id_fixed,
#'                                             fill = fill_fixed),
#'                            fixed = TRUE)
#' lapply(dna_segs1, function(x) x[, .(name, fill)])
#' 
#' ## Add colors based on the presence of a string in the feature names
#' id_grep <- c("A", "B")
#' fill_grep <- c("red", "blue")
#' dna_segs2 <- edit_dna_segs(dna_seg_input = dna_segs,
#'                            ids = data.frame(id = id_grep, fill = fill_grep))
#' lapply(dna_segs2, function(x) x[, .(name, fill)])
#' 
#' ## Use seg_labels to add colors only to specific dna_segs
#' id_labels <- c("Genome 1", "Genome 2")
#' dna_segs3 <- edit_dna_segs(dna_seg_input = dna_segs,
#'                            ids = data.frame(id = id_grep, fill = fill_grep))
#' lapply(dna_segs3, function(x) x[, .(name, fill)])
#' 
edit_dna_segs <- function(
  dna_seg_input,
  ids,
  seg_labels = NULL,
  id_tags = c("name", "locus_id"),
  fixed = FALSE,
  verbose = FALSE,
  ...
) {
  
  # Check dna_seg_input argument
  if (missing(dna_seg_input)) stop("'dna_seg_input' must be provided")
  # If a single dna_seg object was supplied, turn it into a list
  if (is.dna_seg(dna_seg_input)) {
    dna_segs <- list(dna_seg_input)
    return_list <- FALSE
  } else if (!is.list(dna_seg_input) || 
             !all(sapply(dna_seg_input, is.dna_seg))
             ) {
    stop("'dna_seg_input' must be a dna_seg object or a list of dna_seg objects")
  } else {
    dna_segs <- dna_seg_input
    return_list <- TRUE
  }
  
  # This copy is made since this function changes dna_segs by reference
  dna_segs <- copy(dna_segs)
  
  if (missing(ids)) stop("'ids' must be provided")
  if (is.character(ids) && file.exists(ids)) {
    # If ids argument is a string, it is assumed to be a file path
    ids_frame <- fread( ids, header=TRUE, fill = TRUE, ...)
  } else if (is.data.frame(ids)) {
    # If not a file path, then assume it is a data.frame
    ids_frame <- as.data.table(ids)
  } else if (is.character(ids)) {
    stop('The file "', id, " as specified by 'ids' could not be found")
  } else {
    stop("'ids' must be either a character string (specifying a file path), ",
         "or a data.frame")
  }
  
  # Check if the required id column is there
  if (! any("id" == names(ids_frame))) {
    stop("'ids' must contain an 'id' column")
  }
  
  # Remove rows with no value in the id column
  ids_frame <- ids_frame[id != "",]

  # Fill empty strings with NA, but not for the seg_label and id columns
  for (col in names(ids_frame)) {
    if (!any(col == c("seg_label", "id"))) {
      set(ids_frame,
          i = which(ids_frame[[col]] %in% c("", "NA")),
          j = col,
          value = NA
          )
    }
  }
  
  # Remove completely empty columns
  empty <- ids_frame[, names(which(sapply(.SD, function(x) all(is.na(x)))))]
  ids_frame[, (empty) := NULL ]

  # Check region_plot if present
  if (!is.null(ids_frame$region_plot)) {
    ids_frame$region_plot <- as.character(ids_frame$region_plot)
    ids_frame$region_plot[ids_frame$region_plot %ilike%
                            "^true$|^plot$|^y$|^t$|^yes$"] <- "TRUE"
    ids_frame$region_plot[!ids_frame$region_plot %ilike% 
                            "^true$|^plot$|^y$|^t$|^yes$|^start$|^end$" & 
                            !is.na(ids_frame$region_plot)] <- "NA"
    if (!is.null(ids_frame$gene_type)) {
      ids_frame$gene_type[ids_frame$region_plot %ilike% 
                            "^start$|^end$"] <- "boundaries"
    }
  }
  
  # Check seg_labels
  if (!is.null(seg_labels)) {
    seg_labels <- get_seg_labels(dna_segs, seg_labels)
    alt_labels <- get_seg_labels(dna_segs, seg_labels = NULL)
  } else {
    seg_labels <- get_seg_labels(dna_segs, seg_labels)
    alt_labels <- seg_labels
  }
  
  # This code block is here to setup the upcoming loop by attempting to match
  # the seg_labels from the ids file with the dna_seg labels
  if (is.null(ids_frame$seg_label)) {
    # No seg labels, so turn upcoming loop into a simple loop over numbers
    ids_seg_labels <- seq(1,(length(dna_segs)))
    seg_labels <- ids_seg_labels
  } else if (is.null(seg_labels) & is.null(alt_labels)) {
    if (verbose) {
      warning("'ids' has a 'seg_label' column, but dna_seg labels could ",
              "not be determined, so the 'seg_label' column will be ignored")
    }
    # seg_labels could not be determined, so behave as if they were not provided
    ids_frame$seg_label <- NULL
    ids_seg_labels <- seq(1,(length(dna_segs)))
    seg_labels <- ids_seg_labels
  } else {
    ids_seg_labels <- unique(ids_frame$seg_label)
    ids_ori_labels <- ids_seg_labels
    # Try all the different possible variations of matching labels
    label_test1 <- sum(ids_seg_labels %in% seg_labels)
    label_test2 <- sum(ids_seg_labels %in% alt_labels)
    label_test3 <- sum(gsub("\\.[^\\.]*$", "", ids_seg_labels) %in% seg_labels)
    label_test4 <- sum(gsub("\\.[^\\.]*$", "", ids_seg_labels) %in% alt_labels)
    label_tests <- c(label_test1, label_test2, label_test3, label_test4)
    label_test_result <- which.max(label_tests)
    if ((label_test_result %% 2) == 0) {
      # alt_labels matched instead of seg_labels
      seg_labels <- alt_labels
    }
    if (label_test_result >= 3) {
      # ids_seg_labels matched when their "file extensions" were removed
      ids_frame[, seg_label := gsub("\\.[^\\.]*$", "", seg_label)][]
      ids_seg_labels <- gsub("\\.[^\\.]*$", "", ids_seg_labels)
    }
    names(dna_segs) <- seg_labels
    if (any(ids_seg_labels == "")) {
      # Some rows had no seg_label, extend the loop to account for this
      extra_loop <- 0
      ids_seg_labels <- ids_seg_labels[-which(ids_seg_labels == "")]
      ids_seg_labels <- c(ids_seg_labels, rep("", length(dna_segs)))
      seg_labels <- c(seg_labels, "")
    }
  }
  if (verbose) {
    tag_warnings <- rep(0, length(id_tags))
    label_warnings <- character()
  }
  
  # Begin the loop through the seg_labels, or dna_segs, depending on prior setup
  for (i in ids_seg_labels) {
    if (any(i == seg_labels)) {
      if (is.null(ids_frame$seg_label)) {
        # No seg_label columns, so take entire IDs table
        ids_i <- copy(ids_frame)
      } else if (i == "") {
        # seg_label column is there, but some values for it were empty, take 
        # those rows. Loop has to be repeated entirely for these rows.
        extra_loop <- extra_loop + 1
        i <- extra_loop
        ids_i <- ids_frame[seg_label == ""]
        ids_i$seg_label <- NULL
      } else {
        # seg_label column was there
        ids_i <- ids_frame[seg_label == i]
        ids_i$seg_label <- NULL
      }
      
      idnames <- names(ids_i)
      ifelse("id" %in% c(names(dna_segs[[i]]), id_tags),
             remove_id <- FALSE,
             remove_id <- TRUE)
      
      # Search dna_seg object for each user-defined tag in ids_tags
      for (id_tag000 in id_tags) {
        # id_tag000 is named to avoid using an actual column name, as that could create problems
        if (!id_tag000 %in% names(dna_segs[[i]])) {
          indice <- which(id_tag000 == id_tags)
          if (verbose) tag_warnings[indice] <- tag_warnings[indice] + 1
          next
        }
        if (fixed) {
          # NA values in ids are replaced if the column matches a dna_seg column
          idx <- dna_segs[[i]][ids_i,
                               which = TRUE,
                               on = c(paste0(id_tag000, "==id")), mult = "all"
                               ]
          for (col in idnames) {
            nas <- which(is.na(ids_i[[col]]))
            this_idx <- idx[nas]
            if (length(nas) != 0 && col %in% names(dna_segs[[i]])) {
              set(ids_i,
                  i = nas,
                  j = col,
                  value = dna_segs[[i]][[col]][this_idx]
                  )
            }
          }
          
          # dna_seg object is modified in place by merging with ids_i
          dna_segs[[i]][ids_i,
                        on = c(paste0(id_tag000, "==id")),
                        (idnames) := mget(paste0("i.", idnames)),
                        mult = "all"
                        ]
          
          # If the merge resulted in an unnecessary "id" column, remove it
          if (remove_id) dna_segs[[i]][, id := NULL]
        } else {
          # Take only columns that are in ids_i
          matched <- names(
            dna_segs[[i]])[which(names(dna_segs[[i]]) %in% idnames)]
          seg <- copy(dna_segs[[i]][, c(id_tag000, matched), with = FALSE])
          setnames(seg,
                   old = c(id_tag000, matched),
                   new = c(id_tag000, paste0(matched, ".i"))
                   )
          setkeyv(seg, id_tag000)
          setkey(ids_i, id)

          # First merge, only has matched ids
          merged <- ids_i[, seg[grep(id, get(id_tag000))], by = mget(idnames)]
          merged <- unique(merged)
          
          # Replace NA values from ids_i with their dna_seg values
          for (col in matched) {
            merged[is.na(get(col)), (col) := get(paste0(col, ".i"))]
          }
          # Get rid of duplicate columns
          merged[, c(paste0(matched, ".i")) := NULL]
          if (remove_id) merged[, id := NULL]
          
          # Final merge
          m_names <- names(merged)
          dna_segs[[i]][merged,
                        on = c(id_tag000),
                        (m_names) := mget(paste0("i.", m_names)),
                        mult = "all"
                        ][]
        }
      }
    } else if (verbose) {
      label_warnings <- c(label_warnings,
                          ids_ori_labels[which(i == ids_seg_labels)])
    }
  }
  
  if (verbose) {
    for (i in 1:length(tag_warnings)) {
      if (tag_warnings[i] > 0) {
        warning("The \"", id_tags[i], "\" column as supplied by the 'id_tags' ",
                "argument could not be found in ", tag_warnings[i],
                " of the dna_seg object(s)")
      }
    }
    if (length(label_warnings) > 0) {
      warning('No dna_seg object could be found matching: ', 
              paste(label_warnings, collapse = ", "))
    }
  }
  
  
  for (seg in dna_segs) {
    for (col in names(seg)) {
      # Replacing missing values with "NA" in dna_segs
      set(seg, i=which(is.na(seg[[col]])), j=col, value="NA")
    }
    # Check region_plot again to be sure
    if (!is.null(seg$region_plot)) {
      seg$region_plot <- as.character(seg$region_plot)
      seg$region_plot[seg$region_plot %ilike% 
                        "^true$|^plot$|^y$|^t$|^yes$"] <- "TRUE"
      seg$region_plot[!seg$region_plot %ilike% 
                        "^true$|^plot$|^y$|^t$|^yes$|^start$|^end$"] <- "NA"
      seg$gene_type[seg$region_plot %ilike% "^start$|^end$"] <- "boundaries"
    }
  }
  if (!return_list) dna_segs <- dna_segs[[1]]
  dna_segs
}

#' Filter dna_seg features by looking at a maximum within groups
#' 
#' Takes a `dna_seg` or list of `dna_seg` objects. It groups them based on 
#' `group_by`, and per group takes the feature with the maximum value in the
#' column given by `longest`.
#' 
#' @details
#' This was intended to take the longest transcript per gene, although it can
#' be used for other purposes. If `group_by`
#' points to a column with gene IDs, it intentionally mimics the output
#' of the `primary_transcript.py` script from OrthoFinder, so that `dna_segs`
#' can be loaded in from FASTA files before `primary_transcript.py` is used.
#' This preserves the metadata from the FASTA files, since
#' `primary_transcript.py` will remove this metadata.
#' 
#' @returns Either a single `dna_seg` object or a list of `dna_seg` objects,
#' matching the input given using `dna_seg_input`.
#' @export
#' 
#' @param dna_seg_input Either a single `dna_seg` or a list of `dna_seg`
#' objects.
#' @param group_by A character string, representing a `dna_seg` attribute that
#' the features will be grouped by.
#' @param longest A character string, representing a `dna_seg` attribute. After
#' grouping, features will be taken with the maximum value in the column given
#' by this argument.
#' @param ignore_boundaries Logical. If `TRUE`, any features with 
#' `"boundaries"` as their `gene_type` will be kept regardless.
#' 
#' @author Mike Puijk
#' 
#' @seealso [dna_seg]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("1A", "1B", "2A", "2B", "2C")
#' genes1 <- c("1", "1", "2", "2", "2")
#' starts1 <- c(1, 1, 101, 101, 101)
#' ends1 <- c(30, 60, 160, 130, 160)
#' lengths1 <- abs(starts1 - ends1)+1
#' 
#' ## Make dna_seg
#' dna_seg_raw <- dna_seg(data.frame(name=names1, start=starts1, end=ends1,
#'                                   strand=rep(1, 5), length=lengths1,
#'                                   gene=genes1))
#' dna_seg_raw
#' 
#' ## Take longest feature per gene name
#' dna_seg_edit <- max_by_group(dna_seg_input = dna_seg_raw, group_by = "gene")
#' dna_seg_edit
#' 
max_by_group <- function(
  dna_seg_input,
  group_by,
  longest = "length",
  ignore_boundaries = TRUE
) {
  
  # Check dna_seg_input argument
  if (missing(dna_seg_input)) {
    stop("'dna_seg_input' must be provided")
  }
  # If a single dna_seg object was supplied, turn it into a list
  if (is.dna_seg(dna_seg_input)) {
    dna_segs <- list(dna_seg_input)
    return_list <- FALSE
  } else if (!is.list(dna_seg_input) ||
             !all(sapply(dna_seg_input, is.dna_seg))
             ) {
    stop("'dna_seg_input' must be a dna_seg object or a list of dna_seg objects")
  } else {
    dna_segs <- dna_seg_input
    return_list <- TRUE
  }
  
  for (i in 1:length(dna_segs)) {
    if (!is.character(group_by) || !all(group_by %in% names(dna_segs[[i]]))) {
      stop("Each element of 'group_by' must be a character string referring ",
           "to a column in the dna_seg object(s)")
    }
    if (!is.character(longest) || length(longest) != 1 || 
        !all(longest %in% names(dna_segs[[i]]))
        ) {
      stop("'longest' must be a single character string referring to a column ",
           " that is present in all dna_seg objects")
    }
    # Create indices for rows to keep
    # longest_000 is named to avoid using an actual column name, as that would create problems
    longest_000 <- longest
    to_keep <- dna_segs[[i]][,
                             .I[which.max(eval(as.name(longest_000)))],
                             by = c(group_by)
                             ]$V1
    # Extend rows to keep with boundary features
    if (ignore_boundaries) {
      to_keep <- unique(sort(c(to_keep,
                               dna_segs[[i]][, .I[gene_type == "boundaries"]]
                               )))
    } 
    dna_segs[[i]] <- dna_segs[[i]][to_keep]
  }
  if (!return_list) dna_segs <- dna_segs[[1]]
  
  dna_segs
}

#' Make unique IDs for dna_segs
#' 
#' Generates unique identifiers (IDs) for each `dna_seg` features. They can be
#' based on the values from existing columns, or generated from scratch.
#' 
#' @details
#' This function generates unique identifiers for `dna_segs`. Having unique
#' identifiers is necessary for certain other functions, like converting a
#' `dna_seg` into a FASTA file, as most tools that make use of FASTA files
#' require unique headers for each sequence in the FASTA file.
#' 
#' If `old_id` is left as `NULL`, the generated IDs are simply row numbers for 
#' each feature. If `old_id` refers to one or multiple `dna_seg` columns,
#' then the values of those columns are concatenated, separated by `"_"`. Then,
#' a number is added to these values, which starts at 1 for each combination
#' of values, and goes up each time the same combination is found. See the
#' examples below.
#' 
#' @returns Either a single `dna_seg` object or a list of `dna_seg` objects,
#' matching the input given using `dna_seg_input`.
#' @export
#' 
#' @param dna_seg_input Either a single `dna_seg` or a list of `dna_seg`
#' objects.
#' @param old_id Either a character vector representing `dna_seg`
#' columns, or `NULL`. The IDs will be generated based on the vector of 
#' `dna_seg` columns provided, or generated from scratch if this argument is 
#' `NULL`.
#' @param new_id A character string, the generated IDs will be stored in the
#' `dna_seg` column given by this argument. Will create a new column if it does
#' not exist in the `dna_segs`.
#' 
#' @author Mike Puijk
#' 
#' @seealso [dna_seg]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("A", "A", "B", "B", "B", "C")
#' types1 <- c("gene", "gene", "gene", "protein", "gene", "gene")
#' 
#' ## Make dna_seg
#' dna_seg_raw <- dna_seg(data.frame(name = names1,
#'                                   start = (1:6) * 3,
#'                                   end = (1:6) * 3 + 2,
#'                                   strand = rep(1, 6),
#'                                   type = types1))
#' 
#' ## Generate IDs based on 1 column 
#' dna_seg_edit <- make_unique_ids(dna_seg_input = dna_seg_raw,
#'                                 old_id = "name")
#' dna_seg_edit[, .(name, type, id)]
#' 
#' ## Generate IDs based on multiple columns
#' dna_seg_edit <- make_unique_ids(dna_seg_input = dna_seg_raw, 
#'                                 old_id = c("name", "type"))
#' dna_seg_edit[, .(name, type, id)]
#' 
make_unique_ids <- function(dna_seg_input, old_id = NULL, new_id = "id") {
  
  # To avoid assigning by reference
  dna_seg_input <- copy(dna_seg_input)
  
  # Check dna_seg_input argument
  if (missing(dna_seg_input)) stop("'dna_seg_input' must be provided")
  
  # If a single dna_seg object was supplied, turn it into a list
  if (is.dna_seg(dna_seg_input)) {
    dna_segs <- list(dna_seg_input)
    return_list <- FALSE
  } else if (!is.list(dna_seg_input) ||
             !all(sapply(dna_seg_input, is.dna_seg))
             ) {
    stop("'dna_seg_input' must be a dna_seg object or a list of dna_seg objects")
  } else {
    dna_segs <- dna_seg_input
    return_list <- TRUE
  }
  if (!is.null(old_id) && 
      !all(sapply(dna_segs, function(x) old_id %in% names(x)))
      ) {
    stop("Each element of 'old_id' must refer to a column that is present in ",
         "all dna_seg objects or it must be left as NULL")
  }
  
  for (i in 1:length(dna_segs)) {
    # If old_id is left empty, use row numbers as a unique ID
    if (is.null(old_id)) {
      dna_segs[[i]][, as.character(new_id)] <- 1:nrow(dna_segs[[i]])
    } else {
      # Combine the columns supplied by the old_id argument
      id_1 <- dna_segs[[i]][, do.call(paste, c(.SD, sep = "_")), .SDcols = old_id]
      # Add a number representing how often the values of old_id were present
      id_2 <- rowidv(dna_segs[[i]], cols = old_id)
      dna_segs[[i]][, as.character(new_id)] <- paste0(id_1, "_", id_2)
    }
  }
  if (!return_list) dna_segs <- dna_segs[[1]]
  dna_segs
}


#' Update dna_segs using comparisons
#' 
#' A `dna_seg` or list of `dna_segs` is updated using a matching list of
#' comparisons. This can be used to update the `region_plot` or color attributes
#' of `dna_seg` features, by linking them with comparisons.
#' 
#' @details
#' If `dna_seg_input` is a single `dna_seg`, it will be updated using the first
#' (or only) `comparison` from `comparison_input`. 
#' 
#' When updating colors, the colors will be taken from the `col` column from
#' the `comparison` objects, regardless of whether the `col` or `fill` 
#' column is updated in the `dna_segs`, because the `fill` column is absent
#' in `comparison` objects.
#' 
#' The objects are linked together through shared
#' values. The columns for these shared values are determined by the `seg_id`
#' and `comparison_id` arguments, for the `dna_segs` and comparisons, 
#' respectively. `comparison_id` refers to 2 columns, and defaults to `"auto"`, 
#' which will attempt to determine which columns to use automatically.
#' If for example, `comparison_id` is set as `"name"`, it will look for the
#' `"name1"` and `"name2"` columns to match to the `seg_id` in the `dna_segs`
#' above, and under it, respectively.
#' 
#' When `update_from_top` is `TRUE`, it assumes the list of `dna_segs` and
#' comparisons are in plotting order. The first `dna_seg` will not be updated
#' as there is no `comparison` above it to update it with. If instead 
#' `update_from_top` is `FALSE`, the `dna_segs` and comparisons would have to 
#' be supplied in reverse plot order, which is why this is not recommended.
#' In this case the last `dna_seg` in plot order will not be updated, as there 
#' was no `comparison` below it to update it with.
#' 
#' @returns Either a single `dna_seg` object or a list of `dna_seg` objects,
#' matching the input given using `dna_seg_input`.
#' @export
#' 
#' @param dna_seg_input Either a single `dna_seg` or a list of `dna_seg`
#' objects.
#' @param comparison_input Either a single `comparison` or a list of 
#' `comparison` objects.
#' @param seg_id The name of a `dna_seg` column, whose values will be used to
#' make the links to the comparisons.
#' @param comparison_id The shared name of the `comparison` columns, whose 
#' values will be used to make the links to the `dna_segs`. See details.
#' @param update_region_plot Logical. If `TRUE`, updates the `region_plot`
#' attribute of `dna_segs`, which determines whether the neighborhood of these
#' features is plotted in a regional plot.
#' @param color_var The color column to update in the `dna_segs`. Must be either
#' `col`, `fill`, or left as the default `NULL`, which will result in no color
#' column being updated.
#' @param default_color A character string providing the default color of the
#' comparisons, must be either `NULL` or a valid color. The color given by
#' this argument will be ignored when updating, never overwriting any color 
#' in the `dna_segs`.
#' @param update_from_top Logical. If `TRUE`, updates the `dna_segs` with
#' the comparisons above them in the plotting order. Setting this to `FALSE`
#' will make it update the `dna_segs` using the comparisons below them in the
#' plotting order, but this also means that the `dna_segs` and `comparisons`
#' need to be provided in reverse plotting order. Therefore, this is recommended
#' to be used only when providing a single `dna_seg` and `comparison`.
#' 
#' @author Mike Puijk
#' 
#' @seealso [update_comparisons], [sequential_updates], [dna_seg], [comparison]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("1A", "1B", "1C")
#' names2 <- c("2A", "2C", "2B")
#' names3 <- c("3B", "3A", "3C")
#' 
#' ## Make dna_segs
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
#'                                strand = rep(1, 3),
#'                                fill = c("grey80", "grey80", "green")))
#' dna_segs <- list("Genome 1" = dna_seg1,
#'                  "Genome 2" = dna_seg2,
#'                  "Genome 3" = dna_seg3)
#' 
#' ## Make comparisons
#' comp1 <- comparison(data.frame(start1 = c(3, 6, 9), end1 = c(5, 8, 11),
#'                                start2 = c(3, 9, 6), end2 = c(5, 11, 8),
#'                                name1 = c("1A", "1B", "1C"),
#'                                name2 = c("2A", "2B", "2C"),
#'                                col = c("red", "blue", "grey80"),
#'                                direction = c(1, 1, 1)))
#' comp2 <- comparison(data.frame(start1 = c(3, 9, 6), end1 = c(5, 11, 8),
#'                                start2 = c(6, 3, 9), end2 = c(8, 5, 11),
#'                                name1 = c("2A", "2B", "2C"),
#'                                name2 = c("3A", "3B", "3C"),
#'                                col = c("red", "blue", "grey80"),
#'                                direction = c(1, 1, 1)))
#' 
#' ## Before using update_dna_segs
#' plot_gene_map(dna_segs = dna_segs, 
#'               comparisons = list(comp1, comp2),
#'               alpha_comparisons = 0.6)
#' 
#' ## Apply update_dna_segs to update the colors
#' dna_segs <- update_dna_segs(dna_seg_input = dna_segs,
#'                             comparison_input = list(comp1, comp2),
#'                             seg_id = "name",
#'                             update_region_plot = FALSE,
#'                             color_var = "fill",
#'                             default_color = "grey80")
#' 
#' ## After using update_dna_segs
#' ## Because default_color was "grey80", that color will not overwrite dna_segs
#' plot_gene_map(dna_segs = dna_segs, 
#'               comparisons = list(comp1, comp2),
#'               alpha_comparisons = 0.6)
#' 
update_dna_segs <- function(
  dna_seg_input,
  comparison_input,
  seg_id = "locus_id",
  comparison_id = "auto",
  update_region_plot = TRUE,
  color_var = NULL,
  default_color = "grey80",
  update_from_top = TRUE
) {
  
  # Check color_var argument
  if (any(color_var == c("col", "fill"))) {
    update_color <- TRUE
  } else if (is.null(color_var)) {
    update_color <- FALSE
  } else {
    stop("'color_var' must be NULL, or a character string consisting ",
         "of either: col, fill")
  }
  
  # Check update_region_plot argument
  update_region_plot <- as.logical(update_region_plot)
  if (is.na(update_region_plot)) {
    stop("'update_region_plot' could not be coerced to type 'logical'")
  }
  
  # Check dna_seg_input argument
  if (missing(dna_seg_input)) stop("'dna_seg_input' must be provided")
  # If a single dna_seg object was supplied, turn it into a list
  if (is.dna_seg(dna_seg_input)) {
    dna_segs <- list(dna_seg_input)
    return_list <- FALSE
  } else if (!is.list(dna_seg_input) ||
             !all(sapply(dna_seg_input, is.dna_seg))
             ) {
    stop("'dna_seg_input' must be a dna_seg object or ",
         "a list of dna_seg objects")
  } else {
    dna_segs <- dna_seg_input
    return_list <- TRUE
  }
  if (!all(sapply(dna_segs, function(x) seg_id %in% names(x)))) {
    stop("'seg_id' must refer to a column that is present ",
         "in each dna_seg object")
  }
  
  # Check comparison_input argument
  if (missing(comparison_input)) stop("'comparison_input' must be provided")
  # If a single comparison object was supplied, turn it into a list
  if (is.comparison(comparison_input)) {
    comparisons <- list(comparison_input)
  } else if (!is.list(comparison_input) || 
             !all(sapply(comparison_input, is.comparison))
             ) {
    stop("'comparison_input' must be a comparison object or ",
         "a list of comparison objects")
  } else {
    comparisons <- comparison_input
  }
  
  # Determine comparison_id
  if (comparison_id == "auto") {
    cols <- lapply(comparisons, function(x) names(x))
    shared_cols <- names(which(table(unlist(cols)) == length(cols)))
    shared_cols <- shared_cols[!shared_cols %in% c("start1", "end1", 
                                                   "start2", "end2")]
    shared_cols <- shared_cols[which(endsWith(shared_cols, "1") | 
                                       endsWith(shared_cols, "2"))]
    if (all(c("name1", "name2") %in% shared_cols)) {
      comparison_id1 <- "name1"
      comparison_id2 <- "name2"
    } else if (all(c("locus_id1", "locus_id2") %in% shared_cols)) {
      comparison_id1 <- "locus_id1"
      comparison_id2 <- "locus_id2"
    } else {
      # Remove last character of each shared column, then find first duplicate
      numberless <- sapply(shared_cols,
                           function(x) substr(x, 1, nchar(x) - 1),
                           USE.NAMES = FALSE
                           )
      numberless <- duplicated(numberless)
      if (any(numberless)) {
        comparison_id1 <- shared_cols[which(numberless)[1] - 1]
        comparison_id2 <- shared_cols[which(numberless)[1]]
      } else {
        stop('The comparisons must have a set of columns whose names end ',
             'with "1" and "2" that can link them to the dna_seg objects, ',
             'but none could be found')
      }
    }
  } else {
    comparison_id1 <- paste0(comparison_id, "1")
    comparison_id2 <- paste0(comparison_id, "2")
    if (!all(sapply(comparisons,
                    function(x) c(comparison_id1, comparison_id2) %in% names(x)
                    ))
        ) {
      stop("'", comparison_id1, "' and '", comparison_id2, "' (as provided by ",
           "'comparison_id') must be present in all comparison objects")
    }
  }
  
  if (update_color & 
      !all(sapply(comparisons, function(x) "col" %in% names(x)))
      ) {
    stop("Each comparison object must have a 'col' column in order to ", 
         "change the colors of the dna_seg objects")
  }
  
  if (update_region_plot & 
      !all(sapply(comparisons, function(x) "region_plot" %in% names(x)))
      ) {
    stop("Each comparison object must have a 'region_plot' column in ",
         "order to change the region_plot of the dna_seg objects")
  }
  
  # check that there are enough comparisons compared to dna segments
  if (length(dna_segs) - length(comparisons) > 1) {
    stop("There must be at least 1 comparison object for each dna_seg object")
  }
  
  # Check default_color
  if (is.null(default_color)) default_color <- NA
  
  # To prevent accidental overwriting of the original
  dna_segs <- copy(dna_segs)
  
  for (i in 1:length(dna_segs)) {
    if (update_color | update_region_plot) {
      # First look at how and if we are updating this dna_seg
      if (length(dna_segs) > 1 & i == 1) {
        # For the first dna_seg of the list, there is nothing to update it with
        next
      } else {
        # If only 1 dna_seg is provided, use the only provided one
        if (length(dna_segs) == 1) comp_i <- 1
        else comp_i <- i-1
        comp_id <- ifelse(update_from_top, comparison_id2, comparison_id1)
        
        seg <- copy(dna_segs[[i]])
        # Deal with duplicate values of seg_id
        ori_ids <- unlist(seg[[seg_id]], use.names = FALSE)
        dups <- rle(sort(ori_ids))
        dups_i <- which(dups$lengths > 1)
        if (length(dups_i) > 0) {
          temp_ids <- ori_ids
          dups$lengths <- dups$lengths[dups_i]
          dups$values <- dups$values[dups_i]
          for (j in 1:length(dups_i)) {
            dup_id <- which(temp_ids == dups$values[j])
            temp_ids[dup_id] <- paste0(rep(dups$values[j], dups$lengths[j]),
                                       "_",
                                       seq(1:dups$lengths[j])
                                       )
          }
          seg[[seg_id]] <- temp_ids
        }
        
        
      }
      if (update_color) {
        merged <- merge(
          seg[, c(seg_id, color_var, "region_plot"), with = FALSE],
          comparisons[[comp_i]][, c(comp_id, "col"), with = FALSE],
          by.x = seg_id, by.y = comp_id,
          all.x = TRUE,
          sort = FALSE,
          allow.cartesian = TRUE
        )
        setnames(merged, new = c(seg_id, "seg_col", "region_plot", "comp_col"))
        # Prevent default color from overwriting anything
        merged[comp_col == default_color | comp_col == "NA", comp_col := NA]
        # Prevent boundaries from getting overwritten
        merged[region_plot == "start" | region_plot == "end", comp_col := NA]
        merged[, thisone := fifelse(is.na(comp_col), 0, 1)]
        # Takes per seg_id, the first row which is not NA for the comparison color, or first row in general
        to_keep <- merged[, .I[which.max(thisone)], by=c(seg_id)]$V1
        merged <- merged[to_keep, .(ifelse(is.na(comp_col), seg_col, comp_col))]
        set(dna_segs[[i]],
            j = color_var,
            value = unlist(merged, use.names = FALSE)
            )
      }
      if (update_region_plot) {
        merged <- merge(
          seg[, c(seg_id, "region_plot"), with = FALSE],
          comparisons[[comp_i]][, c(comp_id, "region_plot"), with = FALSE],
          by.x = seg_id, by.y = comp_id, 
          all.x = TRUE,
          sort = FALSE,
          allow.cartesian = TRUE
        )
        merged[, region_plot := fifelse(
          region_plot.y == "TRUE" | region_plot.x == "TRUE",
          1,
          0
        )]
        merged[is.na(region_plot), region_plot := 0]
        # Takes per seg_id, the first row which has TRUE for region_plot, or the first row in general
        to_keep <- merged[, .I[which.max(region_plot)], by = c(seg_id)]$V1
        # Also add boundaries to the list of things to keep
        boundaries_to_keep <- which(merged$region_plot.x %in% c("start", "end"))
        to_keep <- sort(unique(c(to_keep, boundaries_to_keep)))
        merged <- merged[to_keep,
                         c(seg_id, "region_plot", "region_plot.x"),
                         with = FALSE
                         ]
        merged[, region_plot := fifelse(region_plot == 1, "TRUE", "NA")]
        merged[region_plot.x %in% c("start", "end"),
               region_plot := region_plot.x
               ]
        set(dna_segs[[i]],
            j = "region_plot",
            value = unlist(merged[, region_plot], use.names = FALSE)
            )
      }
    } else {
      # Nothing is meant to be updated?
      warning("Nothing was updated. Check if this is correct.")
      if (!return_list) dna_segs <- dna_segs[[1]]
      return(dna_segs)
    }
  }
  if (!return_list) dna_segs <- dna_segs[[1]]
  
  dna_segs
}

#' Update comparisons using dna_segs
#' 
#' A `comparison` or list of `comparisons` is updated using a matching
#' list of `dna_segs`. This can be used to update the `region_plot`, `col`, 
#' or position columns of the comparisons, by linking them with `dna_seg`
#' features.
#' 
#' @details
#' If `dna_seg_input` is a single `dna_seg`, it will be updated using the first
#' (or only) `comparison` from `comparison_input`.
#' 
#' The objects are linked together through shared
#' values. The columns for these shared values are determined by the `seg_id`
#' and `comparison_id` arguments, for the `dna_segs` and comparisons, 
#' respectively. `comparison_id` refers to 2 columns, and defaults to `"auto"`, 
#' which will attempt to determine which columns to use automatically.
#' If for example, `comparison_id` is set as `"name"`, it will look for the
#' `"name1"` and `"name2"` columns to match to the `seg_id` in the `dna_segs`
#' above, and under it, respectively.
#' 
#' When `update_from_top` is `TRUE`, it assumes the list of `dna_segs` and
#' comparisons are in plotting order. The first `dna_seg` will not be updated
#' as there is no `comparison` above it to update it with. If instead 
#' `update_from_top` is `FALSE`, the `dna_segs` and comparisons would have to 
#' be supplied in reverse plot order, which is why this is not recommended.
#' In this case the last `dna_seg` in plot order will not be updated, as there 
#' was no `comparison` below it to update it with.
#' 
#' @returns Either a single `comparison` object or a list of 
#' `comparison` objects, matching the input given using `comparison_input`.
#' @export
#' 
#' @param dna_seg_input Either a single `dna_seg` or a list of `dna_seg`
#' objects.
#' @param comparison_input Either a single `comparison` or a list of 
#' `comparison` objects.
#' @param seg_id The name of a `dna_seg` column, whose values will be used to
#' make the links to the comparisons.
#' @param comparison_id The shared name of the `comparison` columns, whose 
#' values will be used to make the links to the `dna_segs`. See details.
#' @param update_positions Logical. If `TRUE`, updates the plotted positions
#' of the comparisons to match the `dna_segs`. `start1` and `end1` will be
#' updated using the `dna_seg` above the `comparison` in plotting order,
#' while `start2` and `end2` will be updated using the `dna_seg` under the 
#' `comparison` in plotting order.
#' @param update_region_plot Logical. If `TRUE`, updates the `region_plot`
#' column of the comparisons. This only serves the purpose of being able to
#' transfer this information to other `dna_segs`, the `region_plot` column has
#' no other function in `comparison` objects.
#' @param color_var The color column to update the comparisons with. If it is
#' `col` or `fill`, then that column will be taken from the `dna_segs` to update
#' the `col` column in the comparisons. Must be either `col`, `fill`, or left as
#' the default `NULL`, which will result in no color column being updated. 
#' @param default_color A character string providing the default color of the
#' `dna_segs`, must be either `NULL` or a valid color. The color given by
#' this argument will be ignored when updating, never overwriting any color 
#' in the comparisons.
#' @param update_from_top Logical. If `TRUE`, updates the comparisons with
#' the `dna_segs` above them in the plotting order. Setting this to `FALSE`
#' will make it update the comparisons using the `dna_segs` below them in the
#' plotting order, but this also means that the `dna_segs` and comparisons
#' need to be provided in reverse plotting order. Therefore, it is recommended
#' to be used only when providing a single `dna_seg` and `comparison`.
#' 
#' @author Mike Puijk
#' 
#' @seealso [update_dna_segs], [sequential_updates], [dna_seg], [comparison]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("1A", "1B", "1C")
#' names2 <- c("2A", "2C", "2B")
#' names3 <- c("3B", "3A", "3C")
#' 
#' ## Make dna_segs
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
#' dna_segs <- list("Genome 1" = dna_seg1,
#'                  "Genome 2" = dna_seg2,
#'                  "Genome 3" = dna_seg3)
#' 
#' ## Add colors based on the presence of a string in the feature names
#' id_grep <- c("A", "B")
#' fill_grep <- c("red", "blue")
#' dna_segs <- edit_dna_segs(dna_seg_input = dna_segs,
#'                           ids = data.frame(id = id_grep, fill = fill_grep))
#' 
#' ## Make comparisons
#' comp1 <- comparison(data.frame(start1 = c(3, 6, 8), end1 = c(5, 8, 10),
#'                                start2 = c(3, 9, 7), end2 = c(5, 11, 9),
#'                                name1 = c("1A", "1B", "1C"),
#'                                name2 = c("2A", "2B", "2C"),
#'                                direction = c(1, 1, 1)))
#' comp2 <- comparison(data.frame(start1 = c(3, 9, 7), end1 = c(5, 11, 9),
#'                                start2 = c(6, 3, 8), end2 = c(8, 5, 10),
#'                                name1 = c("2A", "2B", "2C"),
#'                                name2 = c("3A", "3B", "3C"),
#'                                direction = c(1, 1, 1)))
#' 
#' ## Before using update_comparisons
#' plot_gene_map(dna_segs = dna_segs, 
#'               comparisons = list(comp1, comp2),
#'               alpha_comparisons = 0.6)
#' 
#' ## Apply update_comparisons to fix the positions and update the colors
#' comps <- update_comparisons(dna_seg_input = dna_segs,
#'                             comparison_input = list(comp1, comp2),
#'                             seg_id = "name",
#'                             color_var = "fill",
#'                             update_positions = TRUE)
#' 
#' ## After using update_comparisons
#' plot_gene_map(dna_segs = dna_segs, 
#'               comparisons = comps,
#'               alpha_comparisons = 0.6)
#' 
update_comparisons <- function(
  dna_seg_input,
  comparison_input,
  seg_id = "locus_id",
  comparison_id = "auto",
  update_positions = FALSE,
  update_region_plot = TRUE,
  color_var = NULL,
  default_color = "grey80",
  update_from_top = TRUE
) {
  
  # Check mandatory arguments
  if (missing(dna_seg_input)) stop("'dna_segs' must be provided")
  if (is.dna_seg(dna_seg_input)) {
    dna_segs <- list(dna_seg_input)
  } else if (!is.list(dna_seg_input) || 
             !all(sapply(dna_seg_input, is.dna_seg))
             ) {
    stop("'dna_seg_input' must be a dna_seg object or ",
         "a list of dna_seg objects")
  } else {
    dna_segs <- dna_seg_input
  }
  if (!all(sapply(dna_segs, function(x) seg_id %in% names(x)))) {
    stop("'seg_id' must refer to a column that is present in each dna_seg")
  }
  
  # Check comparison_input argument
  if (missing(comparison_input)) stop("'comparison_input' must be provided")
  # If a single comparison object was supplied, turn it into a list
  if (is.comparison(comparison_input)) {
    comparisons <- list(comparison_input)
    return_list <- FALSE
  } else if (!is.list(comparison_input) || 
             !all(sapply(comparison_input, is.comparison))
             ) {
    stop("'comparison_input' must be a comparison object or ",
         "a list of comparison objects")
  } else {
    comparisons <- comparison_input
    return_list <- TRUE
  }
  # check that there are enough comparisons compared to dna segments
  if (update_positions) {
    if (length(dna_segs) - length(comparisons) != 1) {
      stop("To update positions, the amount of comparison objects ",
           "provided has to be 1 less than the amount of dna_seg objects")
    }
  } else {
    if (!length(dna_segs) - length(comparisons) >= 0) {
      stop("There must be at least 1 dna_seg for each provided comparison ",
           "object")
    }
  }
  
  # Determine comparison_id
  if (comparison_id == "auto") {
    cols <- lapply(comparisons, function(x) names(x))
    shared_cols <- names(which(table(unlist(cols)) == length(cols)))
    shared_cols <- shared_cols[!shared_cols %in% c("start1", "end1", 
                                                   "start2", "end2")]
    shared_cols <- shared_cols[which(endsWith(shared_cols, "1") | endsWith(shared_cols, "2"))]
    if (all(c("name1", "name2") %in% shared_cols)) {
      comparison_id1 <- "name1"
      comparison_id2 <- "name2"
    } else if (all(c("locus_id1", "locus_id2") %in% shared_cols)) {
      comparison_id1 <- "locus_id1"
      comparison_id2 <- "locus_id2"
    } else {
      # Remove last character of each shared column, then find first duplicate
      numberless <- sapply(shared_cols, function(x) substr(x, 1, nchar(x)-1), USE.NAMES = FALSE)
      numberless <- duplicated(numberless)
      if (any(numberless)) {
        comparison_id1 <- shared_cols[which(numberless)[1]-1]
        comparison_id2 <- shared_cols[which(numberless)[1]]
      } else {
        stop('The comparisons must have a set of columns whose names end ',
             'with "1" and "2" that can link them to the dna_seg objects, ',
             'but none could be found')
      }
    }
  } else {
    comparison_id1 <- paste0(comparison_id, "1")
    comparison_id2 <- paste0(comparison_id, "2")
    if (!all(sapply(comparisons,
                    function(x) c(
                      comparison_id1, comparison_id2) %in% names(x)))
        ) {
      stop("'", comparison_id1, "' and '", comparison_id2, "' (as provided by ",
           "'comparison_id') must be present in all comparison objects")
    }
  }
  
  # To prevent accidental overwriting of the original
  comparisons <- copy(comparisons)
  
  # Check update_positions argument
  update_positions <- as.logical(update_positions)
  if (is.na(update_positions)) {
    stop("'update_positions' could not be coerced to type 'logical'")
  }
  
  # Check update_region_plot argument
  update_region_plot <- as.logical(update_region_plot)
  if (is.na(update_region_plot)) {
    stop("'update_region_plot' could not be coerced to type 'logical'")
  }
  region_var <- if (update_region_plot) "region_plot" else NULL
  
  # Check color_var argument
  if (any(color_var == c("col", "fill"))) {
    update_color <- TRUE
  } else if (is.null(color_var)) {
    update_color <- FALSE
  } else {
    stop("'color_var' must be NULL, or a character string consisting ",
         "of either: col, fill")
  }
  
  # Check default_color
  if (is.null(default_color)) default_color <- NA
  
  start_length <- numeric(length(comparisons))
  end_length <- numeric(length(comparisons))
  for (i in 1:length(comparisons)) {
    start_length[i] <- nrow(comparisons[[i]])
    if (update_positions) {
      # Merge dna_segs and comparisons to update positions and possibly colors
      if (update_from_top) {
        merged <- merge(comparisons[[i]],
                        dna_segs[[i]][, c(seg_id, "start", "end", "strand", color_var, region_var), with = FALSE],
                        by.x = comparison_id1, by.y = seg_id, sort = FALSE, allow.cartesian = TRUE)
        merged <- merge(merged,
                        dna_segs[[i+1]][, c(seg_id, "start", "end", "strand"), with = FALSE],
                        by.x = comparison_id2, by.y = seg_id, sort = FALSE, allow.cartesian = TRUE)
      } else {
        merged <- merge(comparisons[[i]],
                        dna_segs[[i]][, c(seg_id, "start", "end", "strand"), with = FALSE],
                        by.x = comparison_id1, by.y = seg_id, sort = FALSE, allow.cartesian = TRUE)
        merged <- merge(merged,
                        dna_segs[[i+1]][, c(seg_id, "start", "end", "strand", color_var, region_var), with = FALSE],
                        by.x = comparison_id2, by.y = seg_id, sort = FALSE, allow.cartesian = TRUE)
      }
      merged[, c("start1", "end1", "start2", "end2") := 
               .(ifelse(strand.x == 1, start.x, end.x),
                 ifelse(strand.x == 1, end.x, start.x),
                 ifelse(strand.y == 1, start.y, end.y),
                 ifelse(strand.y == 1, end.y, start.y))]
      merged[, c("start.x", "end.x", "strand.x",
                 "start.y", "end.y", "strand.y", "direction") := 
               .(NULL, NULL, NULL, NULL, NULL, NULL,
                 ifelse(sign(start1-end1) * sign(start2-end2) > 0, 1, -1))]
    } else if (update_region_plot | update_color) {
      # Merge dna_segs and comparisons without updating positions
      if (update_from_top) {
        merged <- merge(comparisons[[i]],
                        dna_segs[[i]][, c(seg_id, color_var, region_var), with = FALSE],
                        by.x = comparison_id1, by.y = seg_id, sort = FALSE,
                        allow.cartesian = TRUE)
      } else {
        merged <- merge(comparisons[[i]],
                        dna_segs[[i+1]][, c(seg_id, color_var, region_var), with = FALSE],
                        by.x = comparison_id2, by.y = seg_id, sort = FALSE,
                        allow.cartesian = TRUE)
      }
    } else {
      # Nothing is meant to be updated?
      warning("Nothing was updated. Check if this is correct.")
      if (!return_list) comparisons <- comparisons[[1]]
      return(comparisons)
    }
    if (update_color & any(paste0(color_var, ".x") == names(merged))) {
      # Clean up color columns when the chosen column were already present
      if (color_var == "col") {
        merged[col.y != default_color, col := col.y]
        merged[col.y == default_color, col := col.x]
        merged[, c("col.x", "col.y") := .(NULL, NULL)]
      } else if (color_var == "fill") {
        merged[fill.y != default_color, col := fill.y]
        merged[, fill := fill.x]
        merged[, c("fill.x", "fill.y") := .(NULL, NULL)]
      }
    } else if (update_color) {
      if (color_var == "col") {
        merged[col == default_color, col := "NA"]
      } else if (color_var == "fill") {
        merged[fill != default_color, col := fill]
        merged[, fill := NULL]
      }
    }
    if (update_region_plot & any("region_plot.x" == names(merged))) {
      merged[, region_plot := fifelse(region_plot.y == "TRUE" | region_plot.x == "TRUE", "TRUE", "NA")]
      merged[, c("region_plot.x", "region_plot.y") := NULL]
    }
    
    # Clean up missing values 
    if (any("col" == names(merged))) {
      merged[is.na(col), col := "NA"]
    }
    
    # Reorder mandatory comparison columns back to their original positions
    setcolorder(merged, neworder = c("start1", "end1", "start2", "end2",
                                     comparison_id1, comparison_id2))
    comparisons[[i]] <- merged
    end_length[i] <- nrow(comparisons[[i]])
  }
  
  if (any(start_length > end_length)) {
    culled <- which(start_length>end_length)
    warning("Some identifiers could not be matched to anything in dna_segs, ",
            "leading to the removal of rows in these comparisons:\n", 
            paste(culled, collapse = ", "), "\n", "  Resulting in this many ",
            "rows being discarded for these comparisons:\n", 
            paste((start_length-end_length)[culled], collapse = ", "))
  }
  
  if (!return_list) comparisons <- comparisons[[1]]
  comparisons
}


#' Update dna_segs and comparisons sequentially in plotting order
#' 
#' Takes a list of `dna_seg` and a list of `comparison` objects, and then
#' updates both sequentially, in plotting order. This can be used to update the
#' positions of the comparisons, as well as the color and `region_plot` 
#' attributes of both the `dna_segs` and comparisons. It does this by taking the
#' values of features for these attributes and transferring those over to
#' comparisons directly connected to it, and then to the features connected to
#' those comparisons. For example, if a feature from a single `dna_seg` has red
#' as its `fill` attribute, the comparisons that can be linked to this feature
#' will become red as well. This is then followed up by updating any `dna_seg`
#' features linked to those comparisons, and so on.
#' 
#' @details
#' When updating colors, the existing colors from the input `dna_seg` and 
#' `comparison` objects are transferred over to the next object in the
#' plotting order, with the exception of their default colors, provided by
#' `default_color`. As `comparison` objects only have a single color
#' attribute `col`, those will be updated using the column provided by
#' `color_var` from the `dna_segs`, while the `dna_segs` themselves will be
#' updated using the `col` column from the comparisons regardless of `color_var`
#' (unless it is left as `NULL` to avoid updating colors entirely).
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
#' There are some circumstances with very interconnected comparisons where
#' you might want to set `both_directions` to `FALSE` to avoid transferring
#' over the `region_plot` attribute to too many features.
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
#' @param color_var A character string denoting which color attribute to update,
#' must be one of: `"fill"`, `"col"`, or left as `NULL`, which means no color
#' will be updated.
#' @param default_color A character string providing a default color, must be 
#' either `NULL` or a valid color. The color given by this argument will be
#' ignored when updating, never overwriting any other color.
#' @param update_region_plot Logical. If `TRUE`, updates the `region_plot`
#' attribute of `dna_segs` (and of `comparisons` but only to pass this
#' information to the `dna_segs` that follow), which determines whether the 
#' neighborhood of these features is plotted in a regional plot.
#' @param update_positions Logical. If `TRUE`, updates the plotted positions
#' of the comparisons to match the `dna_segs`. `start1` and `end1` will be
#' updated using the `dna_seg` above the `comparison` in plotting order,
#' while `start2` and `end2` will be updated using the `dna_seg` under the 
#' `comparison` in plotting order.
#' @param both_directions Logical. If `FALSE`, updates are applied
#' sequentially in plotting order, starting from the first `dna_seg`. When 
#' `both_directions` is `TRUE`, it will then additionally update each `dna_seg`
#' and `comparison` in reverse plotting order.
#' 
#' @author Mike Puijk
#' 
#' @seealso [update_comparisons], [update_dna_segs], [dna_seg], [comparison]
#' 
#' @examples
#' ## Prepare dna_seg
#' names1 <- c("1A", "1B", "1C")
#' names2 <- c("2A", "2C", "2B")
#' names3 <- c("3B", "3A", "3C")
#' 
#' ## Make dna_segs with some alternate colors
#' dna_seg1 <- dna_seg(data.frame(name=names1, start=(1:3)*3, end=(1:3)*3+2,
#'                                strand=rep(1, 3),
#'                                fill=c("darkred", "grey80", "darkblue")))
#' dna_seg2 <- dna_seg(data.frame(name=names2, start=(1:3)*3, end=(1:3)*3+2,
#'                                strand=rep(1, 3),
#'                                fill=c("grey80", "grey80", "darkgreen")))
#' dna_seg3 <- dna_seg(data.frame(name=names3, start=(1:3)*3, end=(1:3)*3+2,
#'                                strand=rep(1, 3)))
#' 
#' ## Make comparisons
#' comp1 <- comparison(data.frame(start1=c(3,6,9), end1=c(5,8,11),
#'                                start2=c(3,9,6), end2=c(5,11,8),
#'                                name1=c("1A", "1B", "1C"), 
#'                                name2=c("2A", "2B", "2C"),
#'                                direction=c(1,1,1)))
#' comp2 <- comparison(data.frame(start1=c(3,9,6), end1=c(5,11,8),
#'                                start2=c(6,3,9), end2=c(8,5,11),
#'                                name1=c("2A", "2B", "2C"), 
#'                                name2=c("3A", "3B", "3C"),
#'                                direction=c(1,1,1)))
#' 
#' ## Before applying sequential_updates
#' plot_gene_map(dna_segs=list(dna_seg1, dna_seg2, dna_seg3), 
#'               comparisons=list(comp1, comp2),
#'               alpha_comparisons=0.6)
#' 
#' ## Applying sequential_updates without updating in both directions
#' full_data <- sequential_updates(list(dna_seg1, dna_seg2, dna_seg3),
#'                                 comparisons=list(comp1, comp2),
#'                                 seg_id = "name",
#'                                 color_var = "fill",
#'                                 both_directions = FALSE)
#' plot_gene_map(dna_segs=full_data$dna_segs,
#'               comparisons=full_data$comparisons,
#'               alpha_comparisons=0.6)
#' 
#' ## Applying sequential_updates using both directions
#' full_data <- sequential_updates(list(dna_seg1, dna_seg2, dna_seg3),
#'                                 comparisons=list(comp1, comp2),
#'                                 seg_id = "name",
#'                                 color_var = "fill")
#' plot_gene_map(dna_segs=full_data$dna_segs, 
#'               comparisons=full_data$comparisons,
#'               alpha_comparisons=0.6)
#' 
sequential_updates <- function(dna_segs,
                               comparisons,
                               seg_id = "locus_id",
                               comparison_id = "auto",
                               color_var = NULL,
                               default_color = "grey80",
                               update_region_plot = TRUE,
                               update_positions = FALSE,
                               both_directions = TRUE) {
  
  # Check mandatory arguments
  if (missing(dna_segs)) stop("'dna_segs' must be provided")
  if (!is.list(dna_segs) || !all(sapply(dna_segs, is.dna_seg))) {
    stop("'dna_segs' must be a list of dna_seg objects")
  }
  
  if (!is.list(comparisons) || !all(sapply(comparisons, is.comparison))) {
    stop("'comparisons' must be a list of comparison objects")
  }
  
  # check that there are enough comparisons compared to dna segments
  if (length(dna_segs) - length(comparisons) != 1) {
    stop("The amount of comparison objects must be 1 less than the amount of dna_seg objects")
  }
  
  for (i in 1:length(dna_segs)) {
    # There is nothing to update the first dna_seg with
    if (i != 1) {
      # Update each dna_seg using the colors from the comparison above it
      dna_segs[[i]] <- update_dna_segs(dna_segs[[i]],
                                       comparisons[[i-1]],
                                       seg_id = seg_id,
                                       comparison_id = comparison_id,
                                       color_var = color_var,
                                       update_region_plot = update_region_plot)
    }
    # There is 1 less comparison, so final iteration will not update comparison
    if (i != length(dna_segs)) {
      comparisons[[i]] <- update_comparisons(list(dna_segs[[i]],dna_segs[[i+1]]),
                                                 comparisons[[i]],
                                                 seg_id = seg_id,
                                                 comparison_id = comparison_id,
                                                 color_var = color_var,
                                                 update_positions = update_positions,
                                                 update_region_plot = update_region_plot)
    }
  }
  
  if (both_directions) {
    for (i in length(dna_segs):1) {
      # Last dna_seg will already have been updated
      if (i != length(dna_segs)) {
        # Update each dna_seg using the colors from the comparison under it
        dna_segs[[i]] <- update_dna_segs(dna_segs[[i]],
                                         comparisons[[i]],
                                         seg_id = seg_id,
                                         comparison_id = comparison_id,
                                         color_var = color_var,
                                         update_region_plot = update_region_plot,
                                         update_from_top = FALSE)
      }
      # There is 1 less comparison, so final iteration will not update comparison
      if (i != 1) {
        comparisons[[i-1]] <- update_comparisons(list(dna_segs[[i-1]],dna_segs[[i]]),
                                                     comparisons[[i-1]],
                                                     seg_id = seg_id,
                                                     comparison_id = comparison_id,
                                                     color_var = color_var,
                                                     update_positions = update_positions,
                                                     update_region_plot = update_region_plot,
                                                     update_from_top = FALSE)
      }
    }
  }
  
  list(dna_segs = dna_segs, comparisons = comparisons)
}

#' Filter a comparison to include only the best hits for each query
#' 
#' Takes a `comparison` object and filters it to include only the best hit
#' for each query. The query names are provided by `group_by`, and the best hit
#' is determined by sorting each set of queries based on a set of columns
#' provided by `sort_order`.
#' 
#' @details
#' Designed to find the best hits per query from `comparison` objects
#' generated from tabular BLAST or DIAMOND results. For each column from
#' `sort_order`, the maximum (or minimum) value is taken for each unique query
#' name from `group_by`, in the order provided by `sort_order`, until each query
#' only 1 hit left. The maximum values are taken for each `sort_order` column,
#' unless the column provided is one of: `"mism"`, `"gaps"`, `"e_value"`,
#' `"name1"`, `"name2"`, `"start1"`, `"start2"`, `"end1"`, or `"end2"`.
#' 
#' @returns A `comparison` object.
#' @export
#' 
#' @author Mike Puijk
#' 
#' @seealso [bidirectional_best_hit], [comparison], [run_blast], [run_diamond]
#' 
#' @param comparison A `comparison` object to filter.
#' @param group_by A character string referring to a column in `comparison`
#' that holds the query names.
#' @param sort_order A character vector of column names that determine what the
#' best hit is, in order of importance. The first column is used to determine 
#' the best hits per query, and each subsequent column is only used in case of 
#' ties.
#' 
#' @examples
#' ## Read example blastp results
#' infile <- system.file('extdata/blastp_example1.tab', package = 'genoPlotR')
#' 
#' ## comparison before filtering for best hits
#' blast_comparison <- read_comparison_from_blast(infile)
#' print_comparison(blast_comparison)
#' 
#' ## Filter for best hits and print results
#' bh_comparison <- best_hit(blast_comparison)
#' print_comparison(bh_comparison)
#' 
best_hit <- function(
    comparison,
    group_by = "name1",
    sort_order = c("bit_score", "aln_len", "per_id", "direction")
    ) {
  if (!is.character(group_by) || !all(group_by %in% names(comparison))) {
    stop("'group_by' must be a character string referring ",
         "to a column in 'comparison")
  }
  if (!is.character(sort_order) || !all(sort_order %in% names(comparison))) {
    stop("Each element of 'sort_order' must be a character string referring ",
         "to a column in 'comparison")
  }
  
  comparison <- copy(comparison)
  
  n_ids <- length(unique(unlist(comparison[, group_by, with = FALSE])))
  if (nrow(comparison) == n_ids) {
    return(comparison)
  }
  
  for (i in 1:length(sort_order)) {
    if (any(sort_order[i] == c("mism", "gaps", "e_value", "name1", "name2",
                               "start1", "start2", "end1", "end2"))
        ) {
      to_keep <- comparison[, .I[which(
        min(eval(as.name(sort_order[i]))) == eval(as.name(sort_order[i])))],
        by=c(group_by)]$V1
    } else {
      to_keep <- comparison[, .I[which(
        max(eval(as.name(sort_order[i]))) == eval(as.name(sort_order[i])))],
        by=c(group_by)]$V1
    }
    comparison <- comparison[to_keep]
    if (length(to_keep) == n_ids) {
      return(comparison)
    }
  }
  # If a "best" hit was still not selected, take the first "best" hit
  to_keep <- comparison[, .I[which.max(eval(as.name(sort_order[i])))], by=c(group_by)]$V1
  comparison <- comparison[to_keep]
  
  comparison
}

#' Filter a comparison to include only bidirectional best hits
#' 
#' Takes a `comparison` object and filters it to include only
#' bidirectional best hits, with the use of a second `comparison`,
#' provided by `other_direction`. Both `comparison` objects must be
#' filtered for best hits already (see [best_hit]).
#' 
#' @details
#' The best hits from the first
#' `comparison` are only kept when their query-subject combinations can also
#' be found in the best hits from the second `comparison`. For example, take
#' a best hit in the first `comparison` with the query name `"geneA"` and
#' subject name `"geneB"`. In the `comparison` provided by
#' `other_direction`, the best hit for the query name `"geneB"` has to be the
#' subject name `"geneA"`. Only then is the hit kept as a bidirectional best
#' hit. The query names are provided by `group_by1` and `group_by2`, for the
#' `comparison` to filter and the other `comparison` respectively.
#' 
#' @returns A `comparison` object.
#' @export
#' 
#' @author Mike Puijk
#' 
#' @seealso [best_hit], [comparison], [run_blast], [run_diamond]
#' 
#' @param comparison A `comparison` object to filter.
#' @param other_direction A `comparison` object that the `comparison` to
#' filter is compared with.
#' @param group_by1 A character string referring to a column in `comparison`
#' that holds its query names.
#' @param group_by2 A character string referring to a column in
#' `other_direction` that holds its query names.
#' 
#' @examples
#' ## Read example blastp results
#' infile1 <- system.file('extdata/blastp_example1.tab', package = 'genoPlotR')
#' 
#' ## comparison before filtering for best hits
#' blast_comparison1 <- read_comparison_from_blast(infile1)
#' print_comparison(blast_comparison1)
#' 
#' ## Filter for best hits and print results
#' bh_comparison1 <- best_hit(blast_comparison1)
#' print_comparison(bh_comparison1)
#' 
#' ## Repeat steps BLAST results in the other direction
#' 
#' infile2 <- system.file('extdata/blastp_example2.tab', package = 'genoPlotR')
#' blast_comparison2 <- read_comparison_from_blast(infile2)
#' bh_comparison2 <- best_hit(blast_comparison2)
#' 
#' ## Filter for bidirectional best hits and print results
#' bbh_comparison1 <- bidirectional_best_hit(comparison = bh_comparison1,
#'                                          other_direction = bh_comparison2)
#' print_comparison(bbh_comparison1)
#' 
bidirectional_best_hit <- function(comparison,
                                   other_direction,
                                   group_by1 = "name1",
                                   group_by2 = "name2") {
  
  names1 <- paste(unlist(comparison[, group_by1, with = FALSE], use.names = FALSE),
                  unlist(comparison[, group_by2, with = FALSE], use.names = FALSE))
  
  names2 <- paste(unlist(other_direction[, group_by2, with = FALSE], use.names = FALSE),
                  unlist(other_direction[, group_by1, with = FALSE], use.names = FALSE)) 
  
  comparison <- comparison[which(names1 %in% names2)]
  comparison
}



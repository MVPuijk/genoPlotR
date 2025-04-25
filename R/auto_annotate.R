#' Auto-annotate dna_segs
#' 
#' Annotate `dna_segs` in a smart way. This is especially designed for 
#' `dna_segs` read from GenBank or EMBL files, but can be used for other
#' purposes as well. With the default arguments, it produces `annotation`
#' objects based on the `gene` attribute of `dna_segs`.
#' 
#' @details
#' `keep_genes_only` is intended to be used with the `gene` column. When `TRUE`,
#' it will only make annotations for `names` that are not 'empty' (`"-"`, 
#' `"NA"`, or `""`). When it is `FALSE` however, `locus_tag_pattern` becomes 
#' relevant. For any element of `names` that is 'empty', it will take the `name`
#' attribute of the `dna_seg` and remove the `locus_tag_pattern` from it 
#' (e.g. `Eco003456` becomes `003456`, with `Eco` as the `locus_tag_pattern`).
#' If `locus_tag_pattern` is left as `NULL` it will attempt to determine a
#' common prefix automatically. To turn this behavior off, use
#' `locus_tag_pattern = ""`.
#' 
#' If `names` refers to gene names, it will create spanning annotations for
#' operons or sequences of genes. For this to work, gene names have to be 
#' consecutive and end with a number or capital letter.
#' 
#' @returns If `dna_seg_input` is a single `dna_seg`, then a single
#' `annotation` object will be returned.
#' @returns If `dna_seg_input` is a list of `dna_seg` objects, a list of 
#' `annotation` objects will be returned of equal length.
#' @export
#' 
#' @param dna_seg_input Either a single `dna_seg` or a list of `dna_seg`
#' objects.
#' @param names A character string with the name of a `dna_seg` attribute, to
#' base the annotations on. If `dna_seg_input` is a single `dna_seg`, it
#' can also be a character vector with as many elements as there are rows in
#' that `dna_seg`.
#' @param basic_mode Logical. If `TRUE`, annotate the `dna_seg_input` using the
#' `names` argument, ignoring all other arguments and never generating spanning
#' annotations.
#' @param locus_tag_pattern A character string giving a pattern, used to
#' simplify names. To turn this off, use `locus_tag_pattern = ""`. See details. 
#' @param keep_genes_only Logical. If `TRUE`, then for any row where `names` is
#' `"-"`, `"NA"`, or empty (`""`), no annotations will be made. See details.
#' @param dna_seg Deprecated, included for backwards compatibility.
#' When provided, replaces `dna_seg_input`.
#' it.
#' @param ... Further arguments to be passed to the `annotation` function,
#' like `rot` or `col`.
#' 
#' @author Lionel Guy, Mike Puijk
#' 
#' @seealso [annotation], [dna_seg]
#' 
#' @examples
#' ## Prepare dna_seg
#' names <- paste("Eco", sprintf("%04d", 1:20), sep = "")
#' gene <- c("-", "atpC", "atpB", "atpA", "atp2", 
#'           "-", "-", "cda1", "cda2", "cda3",
#'           "vcx23", "vcx22", "vcx21", "cde20",
#'           "-", "gfrU", "gfrT", "gfrY", "gfrX", "gfrW")
#' ds <- dna_seg(data.frame(name = names,
#'                          start = (1:20) * 3,
#'                          end = (1:20) * 3 + 2,
#'                          strand = rep(1, 20),
#'                          gene = gene,
#'                          stringsAsFactors = FALSE))
#' 
#' ## Original annotation
#' annot1 <- annotation(x1 = middle(ds), text = ds$gene, rot = 45)
#' 
#' ## auto_annotate with various options
#' annot2 <- auto_annotate(ds)
#' annot3 <- auto_annotate(ds, keep_genes_only = FALSE,
#'                         locus_tag_pattern = "", rot = 75)
#' annot4 <- auto_annotate(ds, keep_genes_only = FALSE, col = "red")
#' 
#' ## Plot
#' plot_gene_map(dna_segs = list(ds, ds, ds, ds),
#'               annotations = list(annot1, annot2, annot3, annot4))
#' 
# For a dna_seg, auto-annotates:
# Takes gene name if possible, nothing or locus_tag for the rest.
# For operons or sequences of genes, creates spanning annotations
# Gene names have to be consecutive and end with a number or a capital letter
auto_annotate <- function(
  dna_seg_input,
  names = "gene",
  basic_mode = FALSE,
  locus_tag_pattern = NULL,
  keep_genes_only = TRUE,
  dna_seg = NULL, # backwards compatibility
  ...
) {
  
  # Check dna_seg_input argument
  if (missing(dna_seg_input)) {
    if (is.null(dna_seg)) {
      stop("'dna_seg_input' must be provided")
    }
    dna_seg_input <- dna_seg
  }
  # If a single dna_seg object was supplied, turn it into a list
  if (is.dna_seg(dna_seg_input)) {
    dna_segs <- list(dna_seg_input)
    return_list <- FALSE
  } else if (!is.list(dna_seg_input) || 
             !all(sapply(dna_seg_input, is.dna_seg))
             ) {
    stop("'dna_seg_input' must be a dna_seg object or", 
         "a list of dna_seg objects")
  } else {
    dna_segs <- dna_seg_input
    return_list <- TRUE
  }
  
  if (return_list & length(names) != 1 & 
      !all(sapply(dna_segs, function(x) any(names == names(x))))
      ) {
    stop("If 'dna_seg_input' is a list of dna_seg objects, then 'names' must ",
         "be a single character string referring to a column that is present ",
         "in all dna_segs")
  }
  
  annots <- list()
  for (seg in 1:length(dna_segs)) {
    dna_seg <- dna_segs[[seg]][gene_type != "boundaries" & name != "NA", ]
    if (length(names) == 1) {
      names_seg <- unlist(dna_seg[, names, with = FALSE], use.names = FALSE)
    } else {
      names_seg <- names
    }
    nr <- nrow(dna_seg)
    if (is.null(names_seg) || length(names_seg) != nr) {
      stop("'names' must be of the same length as 'dna_seg', or a character ",
           "string referring to a column that is present in all dna_segs")
    }
      
    idx <- rep(TRUE, nr)
    x1 <- middle(dna_seg)
    x2 <- rep(NA, nr)
    if (basic_mode) {
      # If basic mode is enabled, skip all the rest, only remove empty labels
      annots[[seg]] <- annotation(x1 = x1, x2 = x2, text = names_seg, ...)
      annots[[seg]] <- annots[[seg]][!text %in% c("", "-", "NA") & !is.na(text)]
      next
    }
    last_num <- numeric(nr)
    last_char_idx <- numeric(nr)
    prefix <- character(nr)
    if (is.null(locus_tag_pattern)) {
      if (!keep_genes_only) {
        # For all features with no gene info, get common prefix
        sorted_names <- sort(
          dna_seg$name[which(names_seg %in% c("-", "", "NA"))]
        )
        sorted_names <- sorted_names[which(sorted_names != "NA")]
        splitted <- strsplit(sorted_names[c(1, length(sorted_names))], "")
        ### Potentially creates a lot of warnings when lengths don't match
        prefix_indice <- match(FALSE, do.call("==", splitted)) - 1
        if (prefix_indice != 0) {
          # If a common prefix was found, use it as locus_tag_pattern
          seg_pattern <- gsub("0*$",
                              "",
                              substr(sorted_names[1], 1, prefix_indice)
                              )
        } else {
          seg_pattern <- NULL
        }
      }
    } else {
      seg_pattern <- locus_tag_pattern
    }
    for (i in 1:nr) {
      # If gene is not named
      if (is.na(names_seg[i]) || any(names_seg[i] == c("-", "", "NA"))) {
        if (keep_genes_only) {
          idx[i] <- FALSE
        } else if (!is.null(seg_pattern) && is.character(seg_pattern)) {
          names_seg[i] <- gsub(paste("^", seg_pattern, "(.*)", sep = ""),
                               "\\1",
                               dna_seg$name[i],
                               ignore.case = TRUE
                               )
        } else {
          names_seg[i] <- dna_seg$name[i]
        }
      } else {
        # Else, which "type" is the name
        type <- NA
        if (length(grep("[A-Z]$", names_seg[i])) > 0) {
          last_char_idx[i] <- which(LETTERS == gsub(".*(.$)",
                                                    "\\1",
                                                    names_seg[i]
                                                    ))
          prefix[i] <- gsub("^(.*).$", "\\1", names_seg[i])
        } else if (length(grep("\\d+$", names_seg[i])) > 0) {
          ln <- gsub("^\\D*(\\d+$)", "\\1", names_seg[i])
          last_num[i] <- if (ln != names_seg[i]) ln else 0
          prefix[i] <- gsub("^(\\D+)\\d+$", "\\1", names_seg[i])
        }
      }
    }
    # Make sure last_* are numeric
    last_num <- as.numeric(last_num)
    # Letter series
    series_idx <- numeric(0)
    series_dir <- 0
    series_type <- "" # num or char
    for (i in 1:(nr+1)) {
      break_me <- FALSE
      # Has pref + last_num/last_char?
      if (i <= nr && 
          nchar(prefix[i]) > 0 && 
          (last_char_idx[i] > 0 || last_num[i] > 0)
          ) {
        valid <- TRUE
        curr_type <- if (last_num[i] > 0) "num" else "char"
        # Is series started
        if (length(series_idx) > 0) {
          # Is type + prefix compatible?
          if (prefix[i] == prefix[series_idx[1]] && curr_type == series_type) {
            # Has series direction?
            if (series_dir != 0) {
              # Can I increment?
              if ((series_type == "num" && 
                   last_num[i-1] + series_dir == last_num[i]
                   ) ||
                  (series_type == "char" &&
                   last_char_idx[i-1] + series_dir == last_char_idx[i]
                   )
                  ) {
                # Increment (do nothing)
              } else {
                # Else break
                break_me <- TRUE
              }
            } else {
              # No direction. Can I give one?
              dif <- ifelse(series_type == "num",
                            last_num[i] - last_num[i-1],
                            last_char_idx[i] - last_char_idx[i-1]
                            )
              if (abs(dif) == 1) {
                series_dir <- dif
              } else {
                # Else break
                break_me <- TRUE
              }
            }
          } else {
            # Type + prefix not compatible
            break_me <- TRUE
          }
        } else {
          # Series not started, start it
          series_type <- if (last_num[i] > 0) "num" else "char"
        }
      } else {
        # No pref + last_num/last_char
        break_me <- TRUE
        valid <- FALSE
      }
      # Break or increment series
      if (break_me || i > nr) {
        if (length(series_idx) > 1) {
          # Record series
          idx[series_idx[2:length(series_idx)]] <- FALSE
          x1[series_idx[1]] <- dna_seg$start[series_idx[1]]
          x2[series_idx[1]] <- dna_seg$end[series_idx[length(series_idx)]]
          if (series_type == "num") {
            suffices <- last_num[series_idx]
          } else {
            suffices <- LETTERS[last_char_idx[series_idx]]
          }
          
          names_seg[series_idx[1]] <- paste(prefix[series_idx[1]],
                                            suffices[1],
                                            "-",
                                            suffices[length(suffices)],
                                            sep = ""
                                            )
        } else {
          # No series to record, thus do nothing
        } 
        # Reset series
        series_dir <- 0
        # Start a new one at once if valid
        if (valid) {
          series_idx <- i 
          series_type <- curr_type
        } else {
          series_idx <- numeric(0)
          series_type <- ""
        }
      } else {
        # Else increment
        series_idx <- c(series_idx, i)
      }
    }
    annots[[seg]] <- annotation(x1 = x1, x2 = x2, text = names_seg, ...)
    annots[[seg]] <- annots[[seg]][idx,]
  }
  
  if (return_list) return(annots)
  else return(annots[[1]])
}

################################################################################
# Performing sequence alignment within R to create comparison objects
################################################################################

#' Perform sequence alignments using DIAMOND or BLAST
#' 
#' Executes DIAMOND, or a blastn, blastp, or tblastx algorithm.
#' 
#' @details
#' This function requires that the command-line implementation of the tool of
#' choice is installed. In the case of BLAST, this is `blast+`. Uses every
#' sequence in the FASTA file provided by `query` as queries, aligning them
#' against the reference database provided by `subject`.
#' 
#' `evalue`, `max_target_seqs`, and
#' `num_threads` are all directly passed to DIAMOND/ BLAST as arguments. See
#' their respective documentations for further details.
#' 
#' If out is left as `NULL`, the current working directory will be used
#' instead. The file name will be generated based on the query and subject 
#' file names.
#' 
#' @name run_sequence_alignment
#' @returns Returns nothing (invisible `NULL`).
#' @export
#' 
#' @param query Path to a FASTA file containing sequences to query.
#' @param subject Path to a FASTA file or database. If a FASTA file is
#' given, a database will be created first, in the same folder as `out`. File
#' paths to databases should be without the file extension.
#' @param out Output file/ folder path. If the path is an existing folder, a
#' file name will be generated to save the hits to, based on the query
#' and subject file names.
#' @param tool_path Path to the folder containing the blast+ executables. 
#' Only necessary if they are not in the `$PATH` variable.
#' @param evalue Expectation value (E) threshold for hits
#' @param max_target_seqs The maximum amount of hits per query
#' @param num_threads Number of threads (CPUs) to use in the BLAST search
#' @param always_make_db Logical. If `TRUE`, always treats `subject` as
#' a FASTA file and make a database from it first. If `FALSE`, it 
#' will check for the presence of a database first, and only creates one 
#' when necessary.
#' @param verbose Logical. If `TRUE`, reports timings when
#' it starts making the subject database and when it executes the actual 
#' DIAMOND/ BLAST command.
#' @param algorithm Choice of BLAST algorithm, one of: `"blastp"`, 
#' `"blastp-fast"`, `"blastp-short"`, `"blastn"`, `"blastn-short"`, 
#' `"megablast"`, `"dc-megablast"`, or `"tblastx"`
#' 
#' @author Mike Puijk
#' 
#' @seealso [dna_seg_to_fasta], [genbank_to_fasta], [comparisons_from_dna_segs]
#' 
#' @examples
#' \dontrun{
#' ## Running BLAST with the blastp algorithm
#' run_blast(query = "query.fa",
#'           subject = "subject.fa",
#'           out = "query_subject",
#'           algorithm = "blastp",
#'           tool_path = "path/to/blast/executable")
#' 
#' ## Running DIAMOND with the "very-sensitive" option
#' run_diamond(query = "query.fa",
#'             subject = "subject.fa",
#'             out = "query_subject",
#'             sensitivity = "very-sensitive",
#'             tool_path = "path/to/diamond/executable")
#' }
#' 
run_blast <- function(
  query,
  subject,
  out = NULL,
  tool_path = NULL,
  evalue = "1E-3",
  max_target_seqs = 10000,
  num_threads = 1,
  always_make_db = FALSE,
  verbose = FALSE,
  algorithm = "blastp"
) {
  
  # Check BLAST program
  if (any(algorithm == c("blastp", "blastp-fast", "blastp-short"))) {
    db_type <- "prot"
  } else if (any(algorithm == c("blastn", "blastn-short", 
                                "megablast", "dc-megablast", "tblastx"))
             ) {
    db_type <- "nucl"
  } else {
    stop("'algorithm' must be one of the following BLAST programs or ",
         "tasks:\nblastp\n blastp-fast\n blastp-short\nblastn\n blastn-short",
         "\n megablast\n dc-megablast\ntblastx")
  }
  
  # Check always_make_db
  always_make_db <- as.logical(always_make_db)
  if (is.na(always_make_db)) {
    stop("'always_make_db' could not be coerced to type logical")
  }
  
  # Check query
  if (!is.character(query) || !file.exists(query)) {
    stop("The file \"", query,
         "\", as specified by 'query', could not be found")
  }

  # Check subject, and to see whether a database should be made
  if (is.character(subject)) {
    if (always_make_db) {
      if (!file.exists(subject)) {
        stop("The FASTA file \"", subject,
             "\", as specified by 'subject', could not be found")
      }
      make_db <- TRUE
    } else {
      path_str <- ifelse(is.null(tool_path),
                         "",
                         paste0("export PATH=", tool_path, "; ")
                         )
      db_exists <- system(
        paste0(path_str, "blastdbcmd -info", " -db ", subject),
        ignore.stdout = TRUE,
        ignore.stderr = TRUE
      )
      if (db_exists == 0) {
        make_db <- FALSE
      } else if (!file.exists(subject)) {
        stop('No FASTA file or database could be found at "', subject, 
             "\", as specified by 'subject'")
      } else {
        make_db <- TRUE
      }
    }
  } else {
    stop("'subject' must be a character string referring to a FASTA file ",
         "or BLAST database")
  }
  
  # Check out
  if (is.null(out)) {
    outfile_path <- file.path(getwd(), paste0(
      gsub("\\.[^\\.]*$", "", basename(query)), "_",
      gsub("\\.[^\\.]*$", "", basename(subject))
    ))
  } else if (is.character(out)) {
    # If a directory was provided, only construct a file name
    if (dir.exists(out)) {
      outfile_path <- file.path(file.path(out), paste0(
        gsub("\\.[^\\.]*$", "", basename(query)), "_",
        gsub("\\.[^\\.]*$", "", basename(subject))
      ))
    } else {
      outfile_path <- file.path(out)
    }
  } else {
    stop("'out' must be either NULL or a character string referring to a ",
         "folder or file location")
  }
  
  # Create subject database if necessary
  if (make_db) {
    subject_db <- file.path(dirname(outfile_path), basename(subject))
    if (verbose) {
      cat_time('Creating BLAST database for: "', subject, '"')
    }
    tryCatch({
      system(paste0(
        ifelse(is.null(tool_path), "", paste0("export PATH=", tool_path, "; ")),
        "makeblastdb -in ", subject,
        " -input_type fasta",
        " -dbtype ", db_type,
        " -out ", subject_db,
        " -hash_index"
      ))
    }, error = function(e) { 
      stop("makeblastdb was not executed succesfully, check your inputs")
    })
  } else {
    subject_db <- subject
  }
  
  # Converting algorithm into both the BLAST program and task
  algorithm_string <- switch(algorithm,
                             "blastp" = {"blastp -task blastp"},
                             "blastp-fast" = {"blastp -task blastp-fast"},
                             "blastp-short" = {"blastp -task blastp-short"},
                             "blastn" = {"blastn -task blastn"},
                             "blastn-short" = {"blastn -task blastn-short"},
                             "megablast" = {"blastn -task megablast"},
                             "dc-megablast" = {"blastn -task dc-megablast"},
                             "tblastx" = {"tblastx"}
                             )
  
  # Execute BLAST command
  command <- paste0(
    ifelse(is.null(tool_path), "", paste0("export PATH=", tool_path, "; ")), 
    algorithm_string,
    " -query ", query, 
    " -db ", subject_db,
    " -evalue ", evalue,
    " -max_target_seqs ", max_target_seqs,
    " -num_threads ", num_threads, 
    " -outfmt 6",
    " -out ", outfile_path
  )
  if (verbose) {
    cat_time('Running BLAST using the following command:\n', command, "\n")
  }
  tryCatch({
    system(command)
  }, error = function(e) { 
    stop("The BLAST program was not executed succesfully, check your inputs")
  })
  invisible()
}

#' @name run_sequence_alignment
#' @export
#' 
#' @param iterate Logical. If `TRUE`, DIAMOND will start looking for hits using
#' the fast sensitivity setting, increasing the sensitivity if no hits are found
#' until the target sensitivity is reached. Also causes it to stop searching 
#' after the first hit is found for each query.
#' @param sensitivity Choice of the sensitivity option used by DIAMOND, one of:
#' `"fast"`, `"default"`, `"mid-sensitive"`, `"sensitive"`, `"more-sensitive"`,
#' `"very-sensitive"`, or `"ultra-sensitive"`.
#' 
run_diamond <- function(
  query,
  subject,
  out = NULL,
  sensitivity = "default",
  iterate = FALSE,
  tool_path = NULL,
  evalue = "1E-3",
  max_target_seqs = 10000,
  num_threads = 1,
  always_make_db = FALSE,
  verbose = FALSE
) {
  
  # Check sensitivity
  if (sensitivity == "default") {
    sensitivity <- ""
  } else if (any(sensitivity == c("fast", "default", "mid-sensitive", 
                                  "sensitive", "more-sensitive", 
                                  "very-sensitive", "ultra-sensitive"))
             ) {
    sensitivity <- paste0(" --", sensitivity)
  } else {
    stop("'sensitivity' must be one of: fast, default, mid-sensitive, ",
         "sensitive, more-sensitive, very-sensitive, ultra-sensitive")
  }
  
  # Check always_make_db
  always_make_db <- as.logical(always_make_db)
  if (is.na(always_make_db)) {
    stop("'always_make_db' could not be coerced to type logical")
  }
  
  # Check iterate
  iterate <- as.logical(iterate)
  if (is.na(iterate)) stop("'iterate' could not be coerced to type logical")
  if (iterate) {
    iterate <- " --iterate"
  } else {
    iterate <- ""
  }
  
  # Check evalue
  evalue <- as.numeric(evalue)
  if (is.na(evalue)) stop("'evalue' could not be recognized as numeric")
  
  # Check query
  if (!is.character(query) || !file.exists(query)) {
    stop('The file "', query, "\" as specified by 'query', could not be found")
  }
  
  # Check subject, and to see whether a database should be made
  if (is.character(subject)) {
    if (always_make_db) {
      if (!file.exists(subject)) {
        stop('The FASTA file "', subject,
             "\" as specified by 'subject', could not be found")
      }
      make_db <- TRUE
    } else if (file.exists(paste0(subject, ".dmnd"))) {
      make_db <- FALSE
    } else if (!file.exists(subject)) {
      stop('No FASTA file or database could be found at "', subject, 
           "\", as specified by 'subject'")
    } else {
      make_db <- TRUE
    }
  } else {
    stop("'subject' must be a character string referring to a FASTA file ",
         "or DIAMOND database")
  }
  
  # Check out
  if (is.null(out)) {
    outfile_path <- file.path(getwd(), paste0(
      gsub("\\.[^\\.]*$", "", basename(query)), "_",
      gsub("\\.[^\\.]*$", "", basename(subject))
    ))
  } else if (is.character(out)) {
    # If a directory was provided, only construct a file name
    if (dir.exists(out)) {
      outfile_path <- file.path(file.path(out), paste0(
        gsub("\\.[^\\.]*$", "", basename(query)), "_",
        gsub("\\.[^\\.]*$", "", basename(subject))
      ))
    } else {
      outfile_path <- file.path(out)
    }
  } else {
    stop("'out' must be either NULL or a character string referring to a ",
         "folder or file location")
  }
  # Create subject database if necessary
  if (make_db) {
    subject_db <- file.path(dirname(outfile_path), basename(subject))
    if (verbose) {
      cat_time('Creating DIAMOND database for: "', basename(subject), '"\n')
    }
    tryCatch({
      system(paste0(
        ifelse(is.null(tool_path), "", paste0("export PATH=", tool_path, "; ")),
        "diamond makedb --in ", subject,
        " --db ", subject_db
      ))
    }, error = function(e) { 
      stop("diamond makedb was not executed succesfully, check your inputs")
    })
  } else {
    subject_db <- subject
  }
  # Execute DIAMOND command
  command <- paste0(
    ifelse(is.null(tool_path), "", paste0("export PATH=", tool_path, "; ")), 
    "diamond blastp",
    " --query ", query,
    " --db ", subject_db,
    sensitivity,
    iterate,
    " --evalue ", evalue,
    " --max-target-seqs ", max_target_seqs,
    " --threads ", num_threads, 
    " --outfmt 6",
    " --out ", outfile_path
  )
  if (verbose) {
    cat_time('Running DIAMOND using the following command:\n', command, "\n")
  }
  tryCatch({
    system(command)
  }, error = function(e) { 
    stop("The DIAMOND program was not executed succesfully, check your inputs")
  })
  
  invisible()
}

#' FASTA File check
#' 
#' Checks to see if a FASTA file exists for the given inputs. If it does not,
#' creates a FASTA file.
#' 
#' @details
#' This function takes two different FASTA paths to try and cover more ground,
#' mostly for use with the wrapper script. 
#' 
#' @returns The file path to the FASTA file.
#' @noRd
#' 
#' @param fasta_path A file path, it will look for the FASTA file here. If none
#' is found, one will be made and saved in this location.
#' @param fasta_path_alt A file path, it will look for the FASTA file here if
#' none can be found at `fasta_path`.
#' @param use_cache Logical. If `FALSE`, it will always create a FASTA file
#' instead of checking if one already exists.
#' @param file A file path containing the source for the FASTA file. Must be a
#' FASTA or GenBank file. 
#' @param dna_seg A `dna_seg` object. Will be converted to a FASTA file if none
#' is foumd. This argument is ignored if `file` is provided.
#' @param verbose Logical. If `TRUE`, reports timings when creating the FASTA 
#' file.
#' @param take_origin Logical. If `TRUE`, takes the ORIGIN field(s) from the
#' GenBank file provided by `file` instead of the CDS. Fails if a `dna_seg` is
#' provided instead of `file`, since those does not contain the source sequence.
#' 
#' @author Mike Puijk
#' 
#' @seealso [comparisons_from_dna_segs]
#' 
check_fasta <- function(
  fasta_path,
  fasta_path_alt,
  use_cache = TRUE,
  file = NULL,
  dna_seg = NULL,
  verbose = FALSE,
  take_origin = FALSE
) {
  
  # Is a file supplied?
  if (!is.null(file)) {
    peek <- readLines(file, n = 2)
    if (any(grepl("^>", peek))) {
      # Is this a FASTA?
      return(file)
    } else if (any(grepl("^LOCUS", peek))) {
      # If it's a GenBank file, make a FASTA file at the expected path
      if (! use_cache | 
          (! file.exists(fasta_path) & !file.exists(fasta_path_alt))
          ) {
        genbank_to_fasta(file,
                         out = fasta_path,
                         verbose = verbose,
                         take_origin = take_origin
                         )
      } else if (!file.exists(fasta_path) & file.exists(fasta_path_alt)) {
        return(fasta_path_alt)
      }
      return(fasta_path)
    } else {
      stop("Could not recognize file as either FASTA nor GenBank")
    }
  } else if (! use_cache |
             (! file.exists(fasta_path) & !file.exists(fasta_path_alt))
             ) {
    # dna_seg was supplied, but a FASTA of it didn't exist yet, so make one
    if (take_origin) {
      stop("dna_seg objects do not contain the whole source sequence, so a ",
           "FASTA cannot be made when nucleotides are required instead of ",
           "amino acids.")
    } 
    dna_seg_to_fasta(dna_seg,
                     output_path = dirname(fasta_path),
                     file_names = basename(fasta_path),
                     verbose = verbose
                     )
    return(fasta_path)
  } else if (!file.exists(fasta_path) & file.exists(fasta_path_alt)) {
    return(fasta_path_alt)
  } else {
    return(fasta_path)
  }
}

#' Create comparisons between DNA segments
#' 
#' Create a list of `comparison` objects from a list of `dna_seg` objects or
#' files by parsing (and executing) sequence alignments. If these files already
#' exist, then those will be parsed. If not, DIAMOND or a BLAST program can be
#' executed to generate the sequence alignment results between the `dna_seg`,
#' with respect to the order of `dna_segs`. Executing DIAMOND or BLAST requires
#' that the command-line implementations of these tools are installed.
#' 
#' @details
#' Unless `use_cache` is set to `FALSE`, this function will look for the files 
#' required using a combination of the `seg_labels` (if these are provided), and
#' the names of the `dna_segs` or `files` that were provided as input. If it
#' cannot find sequence alignment results in the form of `"query_subject"`
#' (or to put it differently, `"dna_seg1_dna_seg2"`), then it will run DIAMOND
#' or BLAST to generate these results. Using this system, it also looks for the
#' FASTA files required as input for the sequence alignment.
#' 
#' If `output_path` is left as `NULL`, the current working directory will be 
#' used instead.
#' 
#' @returns A list of `comparison` objects.
#' @export
#' 
#' @param dna_segs A list of `dna_seg` objects to create comparisons between.
#' Either `dna_segs` or `files` must be provided.
#' @param seg_labels A character vector containing DNA segment labels.
#' @param files A character vector, containing file paths to the FASTA or
#' GenBank files. The comparisons will be made between these files. Either
#' `dna_segs` or `files` must be provided.
#' @param mode Determines how the comparisons will be filtered.
#' `"besthit"`, `"bidirectional"`, or `"full"`. If mode is `"besthit"`, only the
#' best hit will be taken from each input query (see [best_hit]). If mode is
#' `"bidirectional"`, then hits are only kept if they are the best hits for
#' their query in both directions (see [bidirectional_best_hit]). `"full"` means
#' that all sequence alignment results are considered.
#' @param tool Choice of sequence alignment tool. Either `"blast"` or 
#' `"diamond"`.
#' @param algorithm Choice of BLAST algorithm to run. One of: `"blastp"`,
#' `"blastp-fast"`, `"blastp-short"`, `"tblastx"`, `"blastn"`, `"blastn-short"`,
#' `"megablast"`, or `"dc-megablast"`.
#' @param sensitivity Choice of sensitivity option when running DIAMOND. One of:
#' `"fast"`, `"default"`, `"mid-sensitive"`, `"sensitive"`, `"more-sensitive"`,
#' `"very-sensitive"`, or `"ultra-sensitive"`.
#' @param output_path Path to the folder that will contain the output files.
#' Both the sequence alignment result and the FASTA files used to make them will
#' be stored here.
#' @param all_vs_all Logical. If `TRUE`, sequence alignments will be performed
#' for every combination of the inputs, instead of just the ones necessary for
#' plotting. Note that this can take a long time, so use with caution.
#' @param filt_high_evalue
#' @param filt_low_per_id
#' @param filt_length A number indicating the minimum length required for hits,
#' or `"auto"`. If `"auto"`, it will be determined based on the choice of `tool`
#' and `algorithm` (150 for DIAMOND or any blastp algorithm, 450 for tblastx,
#' 900 for any blastn algorithm).
#' @param use_cache Logical. If `FALSE`, it will never check for existing files.
#' This includes the FASTA files used as input for sequcence alignment, the 
#' database files used by DIAMOND and BLAST, and the sequence alignment results
#' themselves.
#' @param verbose Logical. If `TRUE`, reports timings when creating new files.
#' @param ... Arguments to pass to other functions (the functions executing
#' the sequence alignments tools, [run_blast], and [run_diamond]).
#' 
#' @author Mike Puijk
#' 
#' @seealso [run_blast], [run_diamond], [dna_seg_to_fasta], [genbank_to_fasta]
#' 
#' @examples
#' \dontrun{
#' ## Comparisons from a vector of GenBank files using DIAMOND
#' comparisons <- comparisons_from_dna_segs(
#'   files = c("genome1.gb", "genome2.gb", "genome2.gb"),
#'   tool = "diamond",
#'   output_path = "output/diamond",
#'   sensitivity = "very-sensitive",
#'   verbose = TRUE
#' )
#' }
#' 
comparisons_from_dna_segs <- function(
  dna_segs = NULL,
  seg_labels = NULL,
  files = NULL,
  mode = "full",
  tool = "blast",
  algorithm = "blastp",
  sensitivity = "default",
  output_path = NULL,
  all_vs_all = FALSE,
  filt_high_evalue = NULL,
  filt_low_per_id = NULL,
  filt_length = "auto",
  use_cache = TRUE,
  verbose = FALSE,
  ...
) {
  
  # Check dna_seg
  if (is.null(files)) {
    if (is.null(dna_segs)) {
      stop("Either 'dna_segs' or 'files' must be provided")
    }
    if (!is.list(dna_segs) || !all(sapply(dna_segs, is.dna_seg))) {
      stop("'dna_segs' must be a list of dna_seg objects")
    }
    if (length(dna_segs) < 2) stop("At least 2 dna_segs must be provided")
    # Check seg_labels
    seg_labels <- get_seg_labels(dna_segs, seg_labels)
    if (is.null(seg_labels)) {
      seg_labels <- get_seg_labels(dna_segs, seg_labels, required = TRUE)
      alt_labels <- seg_labels
    } else {
      alt_labels <- get_seg_labels(dna_segs, seg_labels = NULL)
      if (is.null(alt_labels)) {
        alt_labels <- get_seg_labels(dna_segs, seg_labels = NULL,
                                     required = TRUE)
      }
    }
    filenames <- NULL
  } else {
    filenames <- c()
    for (i in files) {
      filenames <- c(filenames, Sys.glob(i))
    }
    if (! length(filenames)) {
      stop("Could not find any files at the location(s) specified.")
    }
    if (is.null(seg_labels)) {
      seg_labels <- gsub("\\.[^\\.]*$", "", basename(filenames))
    }
    alt_labels <- gsub("\\.[^\\.]*$", "", basename(filenames))
    dna_segs <- NULL
  }
  
  mode_options <- c("besthit", "bidirectional", "full")
  tool_options <- c("blast", "diamond",
                    "blastp-fast", "blastp-short",
                    "blastp", "tblastx",
                    "blastn", "blastn-short", 
                    "megablast", "dc-megablast"
                    )
  # Check mode
  if (!any(tolower(mode) == mode_options)) {
    stop("'mode' must be one of: ", paste(mode_options, collapse = ", "))
  } else {
    mode <- tolower(mode)
  }
  # Check tool
  if (!any(tolower(tool) == tool_options)) {
    stop("'tool' must be one of: blast, diamond")
  } else {
    tool <- tolower(tool)
    if (any(tool == c("blastp", "blastp-fast", "blastp-short",
                      "blastn", "blastn-short", 
                      "megablast", "dc-megablast"))
        ) {
      tool <- "blast"
    }
    if (tool == "diamond") {
      if (!any(sensitivity == c("fast", "default", "mid-sensitive", 
                               "sensitive", "more-sensitive", 
                               "very-sensitive", "ultra-sensitive"))
          ) {
        stop("'sensitivity' must be one of: fast, default, mid-sensitive, ",
             "sensitive, more-sensitive, very-sensitive, ultra-sensitive")
      }
    }
  }
  
  # Check output_path 
  if (!is.null(output_path) && !dir.exists(output_path)) {
    stop('The directory "', output_path, '" as specified by ', "'output_path' ",
         'could not be found')
  }
  if (is.null(output_path)) output_path <- getwd()
  output_path <- file.path(output_path)

  # Check use_cache 
  use_cache <- as.logical(use_cache)
  if (is.na(use_cache)) {
    stop("'use_cache' could not be coerced to type logical")
  }
  
  # Check use_cache 
  verbose <- as.logical(verbose)
  if (is.na(verbose)) {
    stop("'verbose' could not be coerced to type logical")
  }
  
  # Check filt_length
  blastp_options <- c("blastp", "blastp-fast", "blastp-short")
  blastn_options <- c("blastn", "blastn-short", "megablast", "dc-megablast")
  
  if (filt_length == "auto") {
    if (any(algorithm == blastp_options) | tool == "diamond") {
      filt_length <- 100
    } else if (any(algorithm == blastn_options) & tool != "diamond") {
      filt_length <- 500
    } else if (any(algorithm == "tblastx") & tool != "diamond") {
      filt_length <- 300
    }
  } else {
    filt_length <- as.numeric(filt_length)
  }
  
  if (any(algorithm == blastn_options)) {
    take_origin <- TRUE
  } else {
    take_origin <- FALSE
  }
  
  fasta_path <- character(length(seg_labels))
  fasta_path_alt <- character(length(alt_labels))
  for (i in 1:length(seg_labels)) {
    fasta_path[i] <- file.path(output_path, seg_labels[[i]])
    fasta_path_alt[i] <- file.path(output_path, alt_labels[[i]])
  }
  
  # Run alignment
  comparisons <- list()
  file_type <- switch(tool, "blast", "diamond" = {"blast"})
  for (i in 1:(length(seg_labels) -1)) {
    full_path <- file.path(output_path,
                           paste0(seg_labels[[i]], "_", seg_labels[[i+1]])
                           )
    full_path_alt <- file.path(output_path,
                               paste0(alt_labels[[i]], "_", alt_labels[[i+1]])
                               )
    if (!use_cache | (!file.exists(full_path) & !file.exists(full_path_alt))) {
      fasta_path[i] <- check_fasta(fasta_path = fasta_path[i],
                                   fasta_path_alt = fasta_path_alt[i],
                                   use_cache = use_cache,
                                   file = filenames[i],
                                   dna_seg = dna_segs[[i]],
                                   verbose = verbose,
                                   take_origin = take_origin
                                   )
      fasta_path[i+1] <- check_fasta(fasta_path = fasta_path[i+1],
                                     fasta_path_alt = fasta_path_alt[i+1],
                                     use_cache = use_cache,
                                     file = filenames[i+1],
                                     dna_seg = dna_segs[[i+1]],
                                     verbose = verbose,
                                     take_origin = take_origin
                                     )
      if (tool == "blast") {
        run_blast(query = fasta_path[i],
                  subject = fasta_path[i+1],
                  algorithm = algorithm,
                  out = full_path,
                  verbose = verbose,
                  ...
                  )
      } else if (tool == "diamond") {
        run_diamond(query = fasta_path[i],
                    subject = fasta_path[i+1],
                    out = full_path,
                    verbose = verbose,
                    sensitivity = sensitivity,
                    ...
                    )
      }
    } else if (!file.exists(full_path) & file.exists(full_path_alt)) {
      # File was found under an alternate name, adjust full_path to match
      full_path <- full_path_alt
    }
    
    comparisons[[i]] <- read_comparison_from_blast(
      full_path,
      fileType = file_type,
      filt_length = filt_length,
      filt_high_evalue = filt_high_evalue,
      filt_low_per_id = filt_low_per_id
    )
    if (mode == "besthit" | mode == "bidirectional") {
      comparisons[[i]] <- best_hit(comparisons[[i]])
    }
    
    if (mode == "bidirectional") {
      # Repeat steps in the opposite direction 
      full_path <- file.path(output_path,
                             paste0(seg_labels[[i+1]], "_", seg_labels[[i]])
                             )
      full_path_alt <- file.path(output_path,
                                 paste0(alt_labels[[i+1]], "_", alt_labels[[i]])
                                 )
      if (!use_cache |
          (!file.exists(full_path) & !file.exists(full_path_alt))
          ) {
        fasta_path[i] <- check_fasta(fasta_path = fasta_path[i],
                                     fasta_path_alt = fasta_path_alt[i],
                                     use_cache = use_cache,
                                     file = filenames[i],
                                     dna_seg = dna_segs[[i]],
                                     verbose = verbose,
                                     take_origin = take_origin
                                     )
        fasta_path[i+1] <- check_fasta(fasta_path = fasta_path[i+1],
                                       fasta_path_alt = fasta_path_alt[i+1],
                                       use_cache = use_cache,
                                       file = filenames[i+1],
                                       dna_seg = dna_segs[[i+1]],
                                       verbose = verbose,
                                       take_origin = take_origin
                                       )
        if (tool == "blast") {
          run_blast(query = fasta_path[i+1],
                    subject = fasta_path[i],
                    algorithm = algorithm,
                    out = full_path,
                    verbose = verbose,
                    ...
                    )
        } else if (tool == "diamond") {
          run_diamond(query = fasta_path[i+1],
                      subject = fasta_path[i],
                      out = full_path,
                      verbose = verbose,
                      sensitivity = sensitivity,
                      ...
                      )
        }
      } else if (!file.exists(full_path) & file.exists(full_path_alt)) {
        # File was found under an alternate name, adjust full_path to match
        full_path <- full_path_alt
      }
      other_direction <- read_comparison_from_blast(
        full_path,
        fileType = file_type,
        filt_length = filt_length,
        filt_high_evalue = filt_high_evalue,
        filt_low_per_id = filt_low_per_id
      )
      other_direction <- best_hit(other_direction)
      comparisons[[i]] <- bidirectional_best_hit(comparisons[[i]],
                                                 other_direction)
    }
    
  }
  if (all_vs_all) {
    # Seperate loop since this needs to include the last dna_seg as well
    for (i in 1:(length(seg_labels) )) {
      # Most of the steps above are repeated, but results are not read.
      other_dna_segs <- 1:(length(seg_labels))
      other_dna_segs <- other_dna_segs[ other_dna_segs != i ]
      for (j in other_dna_segs) {
        full_path <- file.path(output_path,
                               paste0(seg_labels[[i]], "_", seg_labels[[j]])
                               )
        full_path_alt <- file.path(output_path,
                                   paste0(alt_labels[[i]], "_", alt_labels[[j]])
                                   )
        if (!use_cache |
            (!file.exists(full_path) & !file.exists(full_path_alt))
            ) {
          fasta_path[i] <- check_fasta(fasta_path = fasta_path[i],
                                       fasta_path_alt = fasta_path_alt[i],
                                       use_cache = use_cache,
                                       file = filenames[i],
                                       dna_seg = dna_segs[[i]],
                                       verbose = verbose,
                                       take_origin = take_origin
                                       )
          fasta_path[j] <- check_fasta(fasta_path = fasta_path[j],
                                       fasta_path_alt = fasta_path_alt[j],
                                       use_cache = use_cache,
                                       file = filenames[j],
                                       dna_seg = dna_segs[[j]],
                                       verbose = verbose,
                                       take_origin = take_origin
                                       )
          if (tool == "blast") {
            run_blast(query = fasta_path[i],
                      subject = fasta_path[j],
                      algorithm = algorithm,
                      out = full_path,
                      verbose = verbose,
                      ...
                      )
          } else if (tool == "diamond") {
            run_diamond(query = fasta_path[i],
                        subject = fasta_path[j], 
                        out = full_path,
                        verbose = verbose,
                        sensitivity = sensitivity,
                        ...
                        )
          }
        } else if (!file.exists(full_path) & file.exists(full_path_alt)) {
          # File was found under an alternate name, adjust full_path to match
          # Currently unnecessary, is here in case all_vs_all is expanded upon
          full_path <- full_path_alt
        }
      }
    }
  }
  comparisons
}

#' Convert GenBank files to FASTA files
#' 
#' Takes the `CDS` tags with `translation` fields and creates a FASTA file from
#' them, or takes the ORIGIN field(s) to create a FASTA file containing
#' nucleotide sequences instead.
#' 
#' @details
#' If `out` is left as `NULL`, the current working directory will be used
#' instead.
#' 
#' @returns Returns nothing (invisible `NULL`).
#' @export
#' 
#' @param file A GenBank file to convert.
#' @param id_tag A character string, denotes which CDS field to use as the
#' header of the FASTA file sequences. This field must exist in all CDS tags
#' present in the GenBank file.
#' @param take_origin Logical. If `TRUE`, will take the ORIGIN field(s) to 
#' create a FASTA file containing nucleotide sequences instead.
#' @param out Output file path.
#' @param verbose Logical. If `TRUE`, report the time when it starts writing
#' to the output file.
#' 
#' @author Mike Puijk
#' 
#' @seealso [dna_seg_to_fasta], [run_blast], [run_diamond]
#' 
#' @examples
#' \dontrun{
#' genbank_to_fasta(file = "genome1.gb", out = "genome1")
#' }
#' 
genbank_to_fasta <- function(file,
                             id_tag = "locus_tag",
                             take_origin = FALSE,
                             out = NULL,
                             verbose = FALSE) {
  if (!file.exists(file)) {
    stop('Could not find GenBank file "', file, '"')
  }
  
  if (is.null(out)) {
    out <- paste0(gsub("\\.[^\\.]*$", "", basename(file)), ".fa")
  }
  
  imported_data <- readLines(file)
  
  
  main_inds <- grep("^[[:alnum:]\\/]", imported_data)
  main_names <- gsub("*| .*", "",
                     grep("^[[:alnum:]\\/]", imported_data, value = TRUE)
                     )
  if (take_origin) {
    deflines <- imported_data[main_inds[which(main_names == "VERSION")]]
    if (length(deflines) == 0) {
      # Probably not necessary
      deflines <- main_inds[which(main_names == "LOCUS")]
    }
    deflines <- strsplit(gsub(" +", " ", deflines), " ")
    deflines <- sapply(deflines, function(x) paste0(">", unlist(x)[2]))
    
    origin_start <- main_inds[which(main_names == "ORIGIN")]+1
    origin_end <- main_inds[which(main_names == "ORIGIN")+1]-1
    if (any(is.na(origin_end))) {
      stop('Could not detect the end of an ORIGIN field. Does the GenBank ',
           'file properly end with a "//" line?')
    }
    
    if (length(origin_start) == 0 ) {
      stop("Could not detect any ORIGIN field in this file.")
    }
    origins <- character(length(origin_start))
    for (i in 1:length(origin_start)) {
      origins[i] <- paste(
        gsub("( |[0-9]+)", "", imported_data[origin_start[i]:origin_end[i]]),
        collapse = ""
      )
    }
    
    format <- paste0(deflines, "\n", gsub("(.{60})", "\\1\n", origins))
    format <- trimws(format, which = "right")
    
  } else {
    
    features_end <- main_inds[which(main_names == "FEATURES")+1]
    all_indices <- grep(
      "^( {5}|\\t|FT {3})[[:alnum:]'_-]+[[:space:]]+(complement|join|order|[[:digit:]<,])",
      imported_data
    )
    all_indices <- sort(c(all_indices, features_end))
    all_indices_end <- numeric(length(all_indices)-1)
    for (i in 1:(length(all_indices)-1)) {
      all_indices_end[i] <- all_indices[i+1]-1
    }
    
    # Find the start and end of CDS tags
    tag_indices <- grep(
      "^( {5}|\\t|FT {3})CDS[[:space:]]+(complement|join|order|[[:digit:]<,])",
      imported_data
    )
    tag_indices_end <- all_indices_end[which(all_indices %in% tag_indices)]
    # Find the translation fields
    seq_indices <- grep("/translation=", imported_data, fixed = TRUE)
    
    # Find CDS tags without translation fields and remove them from the list
    to_keep <- which(sapply(
      1:length(tag_indices),
      function(x) any(tag_indices[x]:tag_indices_end[x] %in% seq_indices)
    ))
    tag_indices <- tag_indices[to_keep]
    tag_indices_end <- tag_indices_end[to_keep]
    
    # Now remove any translation fields that are not in a CDS tag
    # (should not be able to occur, but just in case)
    if (length(seq_indices) > length(tag_indices)) {
      i_var <- 1
      for (i in 1:length(seq_indices)) {
        # Look for translation fields not in CDS tags
        in_cds <- between(seq_indices[i_var],
                          tag_indices[i_var],
                          tag_indices_end[i_var]
                          )
        if (!in_cds) {
          # When one is found, remove it and don't increment i_var
          seq_indices <- seq_indices[-i_var]
        } else if (i_var == length(seq_indices)) {
          # If end of seq_indices is reached, end loop
          # break
        } else {
          i_var <- i_var + 1
        }
      }
    } else if (length(seq_indices) < length(tag_indices)) {
      # This should not be able to happen if the code so far is correct
      stop('CDS features (and only CDS features) should always include a ',
           '"/translation" qualifier.')
    }
    
    
    deflines <- character(length(tag_indices))
    translations <- character(length(tag_indices))
    for (i in 1:length(tag_indices)) {
      deflines[i] <- grep(
        paste0("/", id_tag, "="),
        imported_data[tag_indices[i]:tag_indices_end[i]],
        fixed = TRUE,
        value = TRUE
      )
      translations[i] <- strsplit(paste(
        imported_data[seq_indices[i]:tag_indices_end[i]],
        collapse = ""
        ), "   /"
      )[[1]][2]
    }
    deflines <- paste0(">", gsub(paste0("(/", id_tag, "=|\"| )"), "", deflines))
    translations <- gsub("( |translation=|\")", "", translations)
    
    format <- paste0(deflines, "\n", gsub("(.{60})", "\\1\n", translations))
    format <- trimws(format, which = "right")
  }
  
  if (verbose) {
    cat_time('Writing "', file, '" as a FASTA file to: ', out)
  }
  conn <- file(out, open = "wb")
  writeLines(format, con = conn, sep = "\n")
  close(conn)
  invisible()
}

#' Convert dna_seg objects to FASTA files
#' 
#' Takes the (non-boundary) features from `dna_seg_input` with sequences
#' and makes a FASTA file from them.
#' 
#' @details
#' If `output_path` is left as `NULL`, the current working directory will be
#' used instead.
#' 
#' @returns Returns nothing (invisible `NULL`).
#' @export
#' 
#' @param dna_seg_input Either a single `dna_seg` object or a list of `dna_seg`
#' objects,
#' @param output_path Path to the folder that will contain the output files.
#' @param id The `dna_seg` column to use as the header for the FASTA sequences.
#' @param file_names A character vector of file names to use. If provided,
#' must be the same length as the amount of provided `dna_seg` objects.
#' @param unique_ids Logical. If `TRUE`, all values in the column provided by
#' `id` have to be unique. If `FALSE` it will use only the first sequence for
#' each unique value of `id`.
#' @param verbose Logical. If `TRUE`, report the time when it starts writing
#' to the output file, as well as a potential warning when the `id` column has
#' duplicate values.
#' @param ... Arguments to pass to [fwrite].
#' 
#' @author Mike Puijk
#' 
#' @seealso [genbank_to_fasta], [run_blast], [run_diamond]
#' 
#' @examples
#' \dontrun{
#' dna_seg_to_fasta(dna_seg_input = list(dna_seg1, dna_seg2), 
#'                  output_path = "path/to/output/folder")
#' }
#' 
dna_seg_to_fasta <- function(dna_seg_input,
                             output_path = NULL,
                             id = "locus_id",
                             file_names = NULL, 
                             unique_ids = FALSE,
                             verbose = FALSE, 
                             ...) {
  # Check dna_seg_input argument
  if (missing(dna_seg_input)) stop("'dna_seg_input' must be provided")
  # If a single dna_seg object was supplied, turn it into a list
  if (is.dna_seg(dna_seg_input)) {
    dna_segs <- list(dna_seg_input)
  } else if (!is.list(dna_seg_input) ||
             !all(sapply(dna_seg_input, is.dna_seg))
             ) {
    stop("'dna_seg_input' must be a dna_seg object or a list of dna_seg objects")
  } else {
    dna_segs <- dna_seg_input
  }
  dna_segs <- copy(dna_segs)
  
  # Check unique_ids argument
  unique_ids <- as.logical(unique_ids)
  if (is.na(unique_ids)) {
    stop("'unique_ids' could not be coerced to type 'logical'")
  }
  
  # Check id argument
  if (!all(sapply(dna_segs, function(x) id %in% names(x)))) {
    stop("'id' must refer to a column that is present in all dna_seg objects")
  }
  
  # Check sequence argument
  if (!all(sapply(dna_segs, function(x) "sequence" %in% names(x)))) {
    stop("All dna_seg objects must contain a 'sequence' column.")
  }
  
  # Check output_path 
  if (!is.null(output_path) && !dir.exists(output_path)) {
    stop('The directory "', output_path, '" as specified by ',
         "'output_path' ", 'could not be found')
  }
  if (is.null(output_path)) output_path <- getwd()
  output_path <- file.path(output_path)
  
  # Check region_plot argument
  # region_plot <- as.logical(region_plot)
  # if (is.na(region_plot)) {
  #   stop("'region_plot' could not be coerced to type 'logical'")
  # }
  
  # Get seg labels 
  if (is.null(file_names)) {
    seg_labels <- get_seg_labels(dna_segs, seg_labels = NULL, required = TRUE)
  } else {
    if (!is.character(file_names) | length(file_names) != length(dna_segs)) {
      stop("'file_names' must be a vector with as many character strings as ", 
           " the amount of supplied dna_seg objects")
    }
    seg_labels <- file_names
  }
  
  print_warning <- FALSE
  # Format and then write out dna_seg objects
  for (i in 1:length(dna_segs)) {
    
    # Filter out features without sequences or an id value of NA
    dna_segs[[i]] <- dna_segs[[i]][gene_type != "boundaries" & sequence != "-"]
    idlist <- unlist(dna_segs[[i]][, id, with = FALSE], use.names = FALSE)
    dna_segs[[i]] <- dna_segs[[i]][which(idlist != "NA")]
    # Check to see if the list of ids is unique
    if (length(idlist) != length(unique(idlist))) {
      if (unique_ids) {
        stop("All values in the '", id, "' column (as provided by 'id'), must ",
             "be unique in all dna_seg objects. You can use the ",
             "make_unique_ids  function to generate unique ids, or you can ",
             "turn this behaviour off by setting 'unique_ids' to FALSE, in ",
             "which case only the first feature will be used for each ",
             "distinct id value."
             )
      } else {
        print_warning <- TRUE
      }
    }
    
    dups <- rowidv(dna_segs[[i]], cols = id)
    format <- dna_segs[[i]][dups == 1, c(id, "sequence"), with = FALSE
                            ][, c("V1", "V2") := .SD
                              ][, .(paste0(">", V1),
                                    trimws(gsub("(.{60})", "\\1\n", V2),
                                           which = "right")
                                    )
                                ]
    
    if (verbose) {
      if (print_warning) {
        warning("Not all values in the '", id, "' column (as provided by ",
                "'id'), are unique in all dna_seg objects, so only the first ",
                "feature will be used for each distinct id value. If the ",
                "'dna_seg_input' contained genome data with introns, you ",
                "safely ignore this warning.")
      }
      cat_time('Writing "', seg_labels[[i]], '" as a FASTA file to: ',
               file.path(output_path, seg_labels[[i]]))
    }
    fwrite(format,
           file = file.path(output_path, seg_labels[[i]]),
           sep = "\n",
           quote = FALSE,
           col.names = FALSE,
           eol = "\n",
           ...
           )
  }
  
  invisible()
}

################################################################################
# File reading functions: read dna_seg
################################################################################
# SUPPORT FUNCTIONS FOR READ_DNA_SEG_FROM_GENBANK
# Get start and end position for a currentFeature Object
# Please note, the function makes sure that only start and end of the
# following two types are ok:
# 'XXX [numeric]..[numeric]' or 'XXX complement([numeric]..[numeric])',
# where XXX may contain a-Z, 0-9 and the [numeric] parts can include < and >
get_start <- function(line) {
  if (length(grep("complement", line)) > 0) {
    start <- as.numeric(gsub(
      "_|[[:blank:]]|[[:alpha:]]|<|\\(|\\)|\\.\\..*", "",
      grep("^[[:graph:]]+ complement\\(<?[[:digit:]]+\\.\\.>?[[:digit:]]+\\)$",
      line, value = TRUE))
    )
  } else {
    start <- as.numeric(gsub(
      "_|[[:blank:]]|[[:alpha:]]|<|\\(|\\)|\\.\\..*", "",
      grep("^[[:graph:]]+ <?[[:digit:]]+\\.\\.>?[[:digit:]]+$",
      line, value = TRUE))
    )
  }
  
  start
}

get_end <- function(line) {
  if (length(grep("complement", line)) > 0) {
    end <- as.numeric(gsub(
      "_|[[:blank:]]|[[:alpha:]]|>|\\(|\\)|.*\\.\\.", "",
      grep("^[[:graph:]]+ complement\\(<?[[:digit:]]+\\.\\.>?[[:digit:]]+\\)$",
      line, value = TRUE))
    )
  } else {
    end <- as.numeric(gsub(
      "_|[[:blank:]]|[[:alpha:]]|>|\\(|\\)|.*\\.\\.", "",
      grep("^[[:graph:]]+ <?[[:digit:]]+\\.\\.>?[[:digit:]]+$",
      line, value = TRUE))
    )
  }
  
  end
}

# Extracts data from feature lines.
extract_data <- function(extract, cF) {
  extract <- gsub(extract,
                  "",
                  grep(paste("^", extract, sep = ""), cF, value = TRUE)
                  )
  if (length(extract) == 0) extract <- "NA"
  extract[1]
}

#' Creating dna_segs from files
#' 
#' Functions to parse `dna_seg` objects from different file formats. Support
#' is included for the following file formats: GenBank, EMBL, ptt, FASTA, and 
#' tabular data.
#' 
#' @details
#' GenBank and EMBL files are two commonly used file types that often contain a
#' great variety of information. To properly extract data from these files, the
#' user has to choose which features to extract. Commonly 'CDS' features are of
#' interest, but other feature tags such as `"gene"`, `"misc_feature"`, or
#' `"tRNA"` may be of interest. Should a feature contain an inner "pseudo" tag
#' indicating this CDS or gene is a pseudo gene, this will be presented as a
#' `"CDS_pseudo"` or a `"gene_pseudo"` feature type respectively in the 
#' resulting table. In these two file types, the following fields are parsed
#' (in addition to the mandatory name, start, end, and strand): gi
#' (from db_xref=GI), uniprot_id, gene, locus_id (from locus_tag=),
#' proteinid, product, color, and region_plot. The sequence itself from CDS tags
#' can also be read using `read_sequence = TRUE`. In addition, extra tags can be 
#' parsed with the argument `extra_fields`. If there is more than one field for
#' a given name, only the first one is parsed.
#' 
#' Tab files representing DNA segments should have at least the following
#' columns: name, start, end, and strand. If these column names are not all
#' present, or if there is no header at all, then these 4 columns are assumed
#' to be the first 4 columns of the file, in that order. If the tab file does
#' have headers, then any additional columns can be supplied as needed, like
#' line width and type, pch and/ or cex. See [dna_seg] for more information.
#' An example:
#' 
#' \tabular{lllll}{
#' name   \tab start \tab end   \tab strand \tab fill  \cr
#' feat1A \tab 2     \tab 1345  \tab 1      \tab blue  \cr
#' feat1B \tab 1399  \tab 2034  \tab 1      \tab red   \cr
#' feat1C \tab 2101  \tab 2932  \tab -1     \tab grey  \cr
#' feat1D \tab 2800  \tab 3120  \tab 1      \tab green \cr
#' }
#' 
#' FASTA files are parsed differently, depending on the first sequence header
#' (defline) found in the file. When parsing FASTA files from UniprotKB,
#' metadata is parsed fromthe deflines, creating extra columns in the resulting
#' `dna_seg` object, including locus_id, gene, and product. Alternatively,
#' positional information of features will be parsed when using FASTA files from
#' ensembl when they contain this information in their headers. In all other
#' cases, each entry in a FASTA file will result in a single feature, all
#' concatenated behind each other in the resulting `dna_seg`. Some support is
#' included for parsing metadata from FASTA files from NCBI.
#' 
#' Ptt (or protein table) files are a tabular format providing information on
#' each protein of a genome (or plasmid, or virus, etc).
#' 
#' @name read_dna_seg
#' @returns A list of `dna_seg` objects for `read_dna_seg_from_files`, and a
#' single `dna_seg` object otherwise.
#' @export
#' 
#' @param file A character string containing a file path, or a file connection.
#' @param files A list or character vector containing file paths. Supports
#' wildcard expansion (e.g. *.txt).
#' @param tagsToParse A character vector of tags to parse for GenBank or EMBL
#' files. Common examples include `"CDS"`, `"gene"`, `"tRNA"`,
#' `"repeat_region"`, and `"misc_feature"`.
#' @param fileType A character string containing the file format to parse.
#' Must be one of: `"genbank"`, `"embl"`, `"ptt"`, `"fasta"`, or `"detect"`. If
#' `"detect"` is chosen, then it will attempt to determine the file format
#' automatically.
#' @param meta_lines The number of lines in the ptt file that represent 
#' metadata, not counting the header lines. Standard for NCBI files is 2 (name 
#' and length, number of proteins).
#' @param gene_type A character string, determines how genes are visualized.
#' Must be a valid gene type (see [gene_types]). For GenBank and EMBL files, if
#' this argument is `"auto"`, the genes will appear as arrows if there are no
#' introns, and as exons (blocks) when there are introns present.
#' @param header Logical. If `TRUE`, parses the first line of the tabular file
#' as a header containing column names.
#' @param extra_fields A character vector of extra fields to parse for GenBank
#' or EMBL files. These fields will be added as columns in the resulting
#' `dna_seg` object, unless the field was always empty.
#' @param boundariesToParse A character vector of tags to parse as sequence
#' boundaries for GenBank or EMBL files. Common examples include `"source"`,
#' `"contig"`, `"chromosome"`, and `"scaffold"`.
#' @param read_sequence Logical. If `TRUE`, will add a sequence column to the
#' `dna_seg` containing the DNA or amino acid sequence of the features.
#' @param verbose Logical. If `TRUE`, reports timings whenever it starts parsing
#' a file.
#' @param ... Further arguments to pass to [as.dna_seg].
#' 
#' @author Lionel Guy, Jens Roat Kultima, Mike Puijk
#' 
#' @seealso [dna_seg], [gene_types]
#' 
#' @examples
#' ## Read DNA segment from tab
#' dna_seg3_file <- system.file('extdata/dna_seg3.tab', package = 'genoPlotR')
#' dna_seg3 <- read_dna_seg_from_tab(dna_seg3_file)
#' 
#' ## From GenBank file
#' bq_genbank <- system.file('extdata/BG_plasmid.gbk', package = 'genoPlotR')
#' bq <- read_dna_seg_from_file(bq_genbank, fileType = "detect")
#' 
#' ## Parsing extra fields in the GenBank file
#' bq <- read_dna_seg_from_file(bq_genbank,
#'                              extra_fields = c("db_xref", "transl_table"))
#' names(bq)
#' 
#' ## From embl file
#' bq_embl <- system.file('extdata/BG_plasmid.embl', package = 'genoPlotR')
#' bq <- read_dna_seg_from_embl(bq_embl)
#' 
#' ## From ptt files
#' bq_ptt <- system.file('extdata/BQ.ptt', package = 'genoPlotR')
#' bq <- read_dna_seg_from_ptt(bq_ptt)
#' 
# MAIN FUNCTION
read_dna_seg_from_file <- function(
  file,
  tagsToParse = c("CDS"),
  fileType = "detect",
  meta_lines = 2,
  gene_type = "auto",
  header = TRUE,
  extra_fields = NULL,
  boundariesToParse = NULL,
  read_sequence = FALSE,
  verbose = FALSE,
  ...
) {
  
  if (verbose) {
    if (is.character(file)) {
      cat_time("Reading in: ", file)
    } else if (any(class(file) == "connection")) {
      cat_time("Reading in: ", summary(file)$description)
    }
  }
  
  # Import data from file into variable
  importedData <- readLines(file)
  
  # Find file type
  if (tolower(fileType) == "detect") {
    if (length(grep(">", importedData[1]))) TYPE <- "Fasta"
    if (length(grep("^ID   ", importedData))) TYPE <- "EMBL"
    if (length(grep("^LOCUS", importedData))) TYPE <- "Genbank"
  } else if (tolower(fileType) == "embl") {
    TYPE <- "EMBL"
  } else if (tolower(fileType) == "genbank") {
    TYPE <- "Genbank"
  } else if (tolower(fileType) == "ptt") {
    TYPE <- "PTT"
  } else if (tolower(fileType) == "fasta") {
    TYPE <- "Fasta"
  } else {
    stop('\'fileType\' has to be either "detect", "embl", "genbank", ',
         '"fasta", or "ptt"\nNote: if file type is .ptt, please specify this ',
         'rather than using "detect"'
         )
  }
  
  # If type corresponds to separate functions
  if (TYPE == "PTT") {
    dna_seg <- read_dna_seg_from_ptt(file, meta_lines, header, ...)
    return(dna_seg)
  } else if (TYPE == "Fasta") {
    dna_seg <- read_dna_seg_from_fasta(file, read_sequence, ...)
    return(dna_seg)
  } else {
  # If type does not correspond to a separate function
    
    # Extract and name main segments
    if (TYPE == "Genbank") {
      mainSegments <- grep("^[[:alnum:]]", importedData)
      names(mainSegments) <- gsub(
        "*| .*",
        "",
        grep("^[[:alnum:]]", importedData, value = TRUE)
      )
    }
    
    # SIMPLE ERROR HANDLING
    if (TYPE == "Genbank") {
      if (length(grep("FEATURES|DEFINITION", names(mainSegments))) < 2) {
        stop("FEATURES or DEFINITION segment missing in GBK File.")
      }
      if (length(grep("LOCUS", names(mainSegments))) == 0) {
        stop("The number of LOCUS lines has to be at least 1.")
      }
    }
    if (TYPE == "EMBL") {
      if (length(grep("^ID|^FH", importedData)) < 2) {
        stop("ID or FH segment missing in EMBL File.")
      }
      if (length(grep("^ID", importedData)) == 0) {
        stop("The number of ID lines has to be at least 1.")
      }
    }
    
    # Extract META data
    if (TYPE == "EMBL") {
      seg_label <- gsub(
        "^DE[[:blank:]]+",
        "",
        grep("^DE", importedData, value = TRUE)
      )
    }
    if (TYPE == "Genbank") {
      seg_label <- gsub(
        "DEFINITION {1,}",
        "",
        importedData[mainSegments["DEFINITION"]]
      )
    }
    
    # Extract features only, handles whether FEATURES is the last (or not)
    # entry in GBK file
    if (TYPE == "Genbank") {
      featureLines <- which(names(mainSegments) == "FEATURES")
      dataFeatures <- c()
      for (i in featureLines) {
        ifelse(
          i == length(mainSegments),
          dataFeatures <- c(
            dataFeatures,
            importedData[mainSegments[i]:(length(importedData) - 1)]
          ),
          dataFeatures <- c(
            dataFeatures,
            importedData[mainSegments[i]:(mainSegments[i+1] - 1)]
          )
        )
      }
    }
    if (TYPE == "EMBL") {
      dataFeatures <- grep("(^FT|^//)", importedData, value = TRUE)
    }
    
    # SIMPLE ERROR HANDLING
    if (TYPE == "Genbank") {
      if (length(dataFeatures) < 2) stop("No FEATURES in GBK file.")
    }
    if (TYPE == "EMBL") {
      if (length(dataFeatures) < 1) stop("No FEATURES in GBK file.")
    }
    
    
    # Extract each start line for each feature
    if (TYPE == "Genbank") {
      startLineOfFeature <- 
        c(1:length(dataFeatures))[- grep("^ {6,}",dataFeatures)]
    }  
    if (TYPE == "EMBL") {
      startLineOfFeature <- grep("(FT   [[:alnum:]]|^//)", dataFeatures)
    }
    startLineOfFeature <- c(startLineOfFeature, length(dataFeatures)+1)
    
    
    # This might require more testing, but this seems to accurately calculate
    # the eventual length of the dna_seg, which allows us to avoid making new
    # vectors constantly, speeding up this function.
    tagsToParse <- c(tagsToParse, boundariesToParse)
    tagQuery <- paste0(tagsToParse, collapse = "|")
    tagQuery <- paste0("^( {5}|\\t|FT {3})(",
                       tagQuery,
                       ")+[[:space:]]+(complement|join|order|[[:digit:]<,])"
                       )
    tagCount <- sum(grepl(tagQuery, importedData)
    )
    
    boundaryQuery <- paste0(boundariesToParse, collapse = "|")
    if (boundaryQuery != "") {
      boundaryQuery <- paste0(
        "^( {5}|\\t|FT {3})(",
        boundaryQuery,
        ")+[[:space:]]+(complement|join|order|[[:digit:]<,])"
      )
      boundaryCount <- sum(grepl(boundaryQuery, importedData))
    } else {
      boundaryCount <- 0
    }
    
    tagResult <- grep(tagQuery, importedData, value = TRUE)
    intronCount <- sum(lengths(regmatches(tagResult, gregexpr(",", tagResult))))
    
    segLength <- tagCount + boundaryCount + (intronCount * 2)
    
    
    # Define variables for storage
    nF <- length(startLineOfFeature) - 1
    name <- character(segLength)
    start <- numeric(segLength)
    end <- numeric(segLength)
    strand <- numeric(segLength)
    featLength <- numeric(segLength)
    gi <- character(segLength)
    uniprot_id <- character(segLength)
    gene <- character(segLength)
    locus_id <- character(segLength)
    product <- character(segLength)
    color <- character(segLength)
    proteinid <- character(segLength)
    feature <- character(segLength)
    geneType <- character(segLength)
    region_plot <- character(segLength)
    seq_origin <- character(segLength)
    sequence <- character(segLength)
    extra <- list()
    for (tag in extra_fields) {
      extra[[tag]] <- character(segLength)
    }
    
    curSequenceStart <- 0
    featCounter <- 0
    seqCounter <- 0
    
    # Loop over all features
    for (counter in 1:nF) {
      
      # Get feature, normally 20ish lines.
      currentFeature <- (dataFeatures[
        startLineOfFeature[counter] : (startLineOfFeature[counter + 1] - 1)
      ])
      # Clean up feature, decreases number of lines etc.
      if (TYPE == "Genbank") {
        currentFeature <- gsub(
          "^ |:|\"| $",
          "",
          gsub(
            "[[:blank:]]+|[[:space:]]+",
            " ",
            strsplit(paste(currentFeature, collapse = ""), "   /")[[1]]
          )
        )
      }
      if (TYPE == "EMBL") {
        currentFeature <- gsub(
          "^ |:|\"| $",
          "",
          gsub(
            "[[:blank:]]+|[[:space:]]+",
            " ",
            strsplit(paste(gsub("FT", "", currentFeature), collapse = ""),
                     "   /"
                     )[[1]]
          )
        )
      }
      # If the current feature is the start of a new sequence
      if (TYPE == "Genbank" &
          startsWith(currentFeature[1], "FEATURES") &
          counter != 1
          ) {
        curSequenceStart <- max(start, end)
      } else if (TYPE == "EMBL" & startsWith(currentFeature[1], "//")) {
        curSequenceStart <- max(start, end)
      }
      
      # If feature is of a type to parse. Default is only CDS tags.
      if (length(grep(gsub(" [[:print:]]+", "", currentFeature[1]),
                      tagsToParse)) > 0) {
        # Create list with exons to parse
        tag <- gsub(" [[:graph:]]+", "", currentFeature[1])
        exonVector <- strsplit(
          gsub("[[:alpha:]]|_| |\\(|\\)|", "", currentFeature[1]),
          ","
        )
        if (length(grep("complement", currentFeature[1])) > 0) {
          exonVector <- paste(tag,
                              " complement(", exonVector[[1]], ")",
                              sep = ""
                              )
        }
        if (length(grep("complement", currentFeature[1])) == 0) {
          exonVector <- paste(tag, " ", exonVector[[1]], sep = "")
        }

        # If there is more than one exon, insert introns
        if (length(exonVector) > 1) {
          exonVector2 <- 0
          for (i in 1:(length(exonVector) - 1)) {
            if (length(grep("complement", currentFeature[1])) > 0) {
              exonVector2 <- c(exonVector2,
                               exonVector[i],
                               paste(tag,
                                     "_intron complement(",
                                     get_end(exonVector[i]) + 1,
                                     "..",
                                     get_start(exonVector[i+1]) - 1,
                                     ")",
                                     sep = ""
                                     )
                               )
            }
            if (length(grep("complement", currentFeature[1])) == 0) {
              exonVector2 <- c(exonVector2,
                               exonVector[i],
                               paste(tag,
                                     "_intron ",
                                     get_end(exonVector[i]) + 1,
                                     "..",
                                     get_start(exonVector[i+1]) - 1,
                                     sep = ""
                                     )
                               )
            }
          }
          exonVector2[length(exonVector) * 2] <- exonVector[length(exonVector)]
          exonVector2 <- exonVector2[2:length(exonVector2)]
          exonVector <- exonVector2
        }
        
        # For each exon in currentFeature
        for (currentExon in exonVector) {
          
          # Set currentExon to currentFeature line 1
          currentFeature[1] <- currentExon
          # SIMPLE ERROR HANDLING AND Only continue parsing if start and stop
          # is valid...
          # Extract gene name or ID AND start and end, THEN, check if it's ok.
          ifelse(length(grep("gene=", currentFeature)) > 0,
                 nameTEMP <- extract_data("gene=", currentFeature),
                 nameTEMP <- extract_data("locus_tag=", currentFeature)
                 )
          startTEMP <-get_start(currentFeature[1])
          endTEMP <- get_end(currentFeature[1])
          if (length(startTEMP) == 0 || length(endTEMP) == 0) {
            warning('Start and stop position invalid for "', nameTEMP,
                    '", this entry has been excluded')
          }
          
          # Continue if start and end is ok... Otherwise, skip this feature
          if (length(startTEMP) > 0 && length(endTEMP) > 0) {
            startTEMP <- startTEMP + curSequenceStart
            endTEMP <- endTEMP + curSequenceStart
            
            
            # Special case for features parsed as sequence boundaries
            if (length(grep(gsub(" [[:print:]]+", "", currentFeature[1]),
                            boundariesToParse
                            )
                       ) > 0
                ) {
              featCounter <- featCounter + 1
              seqCounter <- seqCounter + 1
              name[featCounter:(featCounter+1)] <- c(
                paste0("start_seq_", as.character(seqCounter)),
                paste0("end_seq_", as.character(seqCounter))
              )
              start[featCounter:(featCounter+1)] <- c(startTEMP, endTEMP)
              end[featCounter:(featCounter+1)] <- c(startTEMP, endTEMP)
              featLength[featCounter:(featCounter+1)] <- c(1, 1)
              strand[featCounter:(featCounter+1)] <- c(1, 1)
              gi[featCounter:(featCounter+1)] <- c("NA", "NA")
              uniprot_id[featCounter:(featCounter+1)] <- c("NA", "NA")
              gene[featCounter:(featCounter+1)] <- c("NA", "NA")
              seq_origin[featCounter:(featCounter+1)] <- c(
                as.character(seqCounter),
                as.character(seqCounter)
              )
              locus_id[featCounter:(featCounter+1)] <- c(
                paste0("start_seq_", as.character(seqCounter)),
                paste0("end_seq_", as.character(seqCounter))
              )
              proteinid[featCounter:(featCounter+1)] <- c("NA", "NA")
              product[featCounter:(featCounter+1)] <- c("NA", "NA")
              sequence[featCounter:(featCounter+1)] <- c("NA", "NA")
              region_plot[featCounter:(featCounter+1)] <- c("start", "end")
              color[featCounter:(featCounter+1)] <- c("black", "black")
              for (tag in names(extra)) {
                extra[[tag]][featCounter:(featCounter+1)] <- c(
                  extract_data(paste(tag, "=", sep = ""), currentFeature),
                  extract_data(paste(tag, "=", sep = ""), currentFeature)
                )
              }
              geneType[featCounter:(featCounter+1)] <- c("boundaries", 
                                                         "boundaries")
              feature[featCounter:(featCounter+1)] <- c(
                gsub(" [[:print:]]+", "", currentFeature[1]),
                gsub(" [[:print:]]+", "", currentFeature[1])
              )
              
              # Two features are created, hence it it increments here again
              featCounter <- featCounter + 1
              
            } else {
              # Save name, start, end, length
              featCounter <- featCounter + 1
              name[featCounter] <- nameTEMP
              start[featCounter] <- startTEMP
              end[featCounter] <- endTEMP
              featLength[featCounter] <- (
                get_end(currentFeature[1]) - get_start(currentFeature[1]) + 1
              ) / 3 - 1
              
              # Set strand to 1 or -1
              ifelse(length(grep("complement", currentFeature[1])) > 0,
                     strand[featCounter] <- -1,
                     strand[featCounter] <- 1
                     )
              
              # Extract GI
              gi[featCounter] <- extract_data("db_xref=GI", currentFeature)

              # Extract UniprotKB protein accession
              uniprotTEMP <- extract_data("db_xref=UniProtKB/Swiss-Prot",
                                          currentFeature)
              if (uniprotTEMP == "NA") {
                uniprotTEMP <- extract_data("db_xref=UniProtKB/TrEMBL",
                                            currentFeature)
              }
              if (uniprotTEMP == "NA") {
                uniprotTEMP <- extract_data("db_xref=UniProt/Swiss-Prot",
                                            currentFeature)
              }
              if (uniprotTEMP == "NA") {
                uniprotTEMP <- extract_data("db_xref=UniProt/TrEMBL",
                                            currentFeature)
              }
              uniprot_id[featCounter] <- uniprotTEMP
              
              # Extract gene
              ifelse(length(grep("gene=", currentFeature)) > 0,
                     gene[featCounter] <- extract_data("gene=", currentFeature),
                     gene[featCounter] <- "-"
                     )
              
              # Extract locus_id
              locus_id[featCounter] <- extract_data("locus_tag=",
                                                    currentFeature)
              
              # Extract protein ID
              proteinid[featCounter] <- extract_data("protein_id=",
                                                     currentFeature)
              
              # Extract product
              product[featCounter] <- extract_data("product=", currentFeature)
              
              # Extract sequence
              if (read_sequence) {
                extract <- gsub(
                  "translation=",
                  "",
                  grep(paste0("^", "translation="), currentFeature, value = T)
                )
                extract <- gsub(" ", "", extract)
                if (length(extract) == 0) sequence[featCounter] <- "-"
                else sequence[featCounter] <- extract[1]
              } else {
                sequence[featCounter] <- "NA"
              }
              
              # Extract color
              color[featCounter] <- extract_data("(color|colour)=",
                                                 currentFeature)
              
              # Extract whether to plot this features region
              region_plot[featCounter] <- extract_data("region_plot=",
                                                       currentFeature)
              
              # Extract extra
              for (tag in names(extra)) {
                extra[[tag]][featCounter] <- extract_data(
                  paste(tag, "=", sep = ""),
                  currentFeature
                )
              }
              
              # Set geneType
              if (length(grep("intron", currentFeature[1])) > 0) {
                geneType[featCounter] <- "introns"
              }

              if (length(grep("intron", currentFeature[1])) == 0) {
                geneType[featCounter] <- gene_type
              }
              
              # Return tag feature info, with or without added _pseudo tag
              if (length(grep("^pseudo", currentFeature)) > 0) {
                feature[featCounter] <- paste(
                  gsub(" [[:print:]]+","", currentFeature[1]),
                  "_pseudo",
                  sep = ""
                )
              }
              if (length(grep("^pseudo", currentFeature)) == 0) {
                feature[featCounter] <- gsub(" [[:print:]]+",
                                             "",
                                             currentFeature[1]
                                             )
              }
              # Set sequence name
              ifelse(seqCounter == 0,
                     seq_origin[featCounter] <- "1",
                     seq_origin[featCounter] <- as.character(seqCounter)
                     )
            }
            
          # end of parse if start and end are ok
          }
          
        # End of exon loop
        }
        
      # End of CDS loop
      }

    # End of loop over all features
    }
    
    # SIMPLE ERROR HANDLING
    if (! is.numeric(start)) stop("'start' is not numeric.")
    if (! is.numeric(end)) stop("'end' is not numeric.")
    if (! is.numeric(featLength)) stop("'length' is not numeric.")
    if (! is.numeric(strand)) stop("'strand' is not numeric.")
    if (! is.character(gi)) stop("'gi' is not character.")
    if (! is.character(uniprot_id)) stop("'uniprot_id' is not character.")
    if (! is.character(name)) stop("'name' is not character.")
    if (! is.character(gene)) stop("'gene' is not character.")
    if (! is.character(locus_id)) stop("'locus_id' is not character.")
    if (! is.character(product)) stop("'product' is not character.")
    if (! is.character(proteinid)) stop("'proteinid' is not character.")
    if (! is.character(feature)) stop("'feature' is not character.")
    if (! is.character(sequence)) stop("'sequence' is not character.")
    if (! is.character(region_plot)) stop("'region_plot' is not character.")
    
    # Check color
    # Eventually, change Artemis colors to their RGB equivalent
    artCol <- artemisColors()
    if (length(color) > 0 && all(color %in% c("NA", artCol$n))) {
      for (i in 1:length(color)) {
        if (color[i] != "NA") color[i] <- artCol$colors[artCol$n == color[i]]
      }
    }

    # If gene_type is auto, arrows, unless exons and introns are present
    if (length(grep("intron", geneType)) >= 1 && gene_type == "auto") {
      geneType[geneType == "auto"] <- "exons"
    }
    if (length(grep("intron", geneType)) == 0 && gene_type == "auto") {
      geneType[geneType == "auto"] <- "arrows"
    }
    
    # Check for NA locus_ids
    if (any(locus_id == "NA")) {
      locus_id[which(locus_id == "NA")] <- sapply(
        which(locus_id == "NA"),
        function(x) paste(feature[x], start[x], sep = "_")
      )
    }
    
    # Check for NA names
    # if (any(name == "NA")) {
    #   name[which(name == "NA")] <- sapply(
    #     which(name == "NA"),
    #     function(x) paste(feature[x], start[x], sep = "_")
    #   )
    # }

    # Create table of features
    table <- data.table(name = name, start = start, 
                        end = end, strand = strand,
                        length = featLength, gi = gi,
                        uniprot_id = uniprot_id, gene = gene,
                        locus_id = locus_id, product = product,
                        proteinid = proteinid, feature = feature,
                        gene_type = geneType, region_plot = region_plot,
                        seq_origin = seq_origin, sequence = sequence
                        )
    
    # Add extra fields
    for (tag in extra_fields) {
      if (any(tag == c("lty", "lwd", "pch", "cex"))) {
        # Convert to numerical, replace NAs with 1
        extra[[tag]] <- suppressWarnings(as.numeric(extra[[tag]]))
        extra[[tag]][is.na(extra[[tag]])] <- 1
      }
      table[[tag]] <- extra[[tag]]
    }
    
    # Remove empty columns 
    empty <- table[, names(which(sapply(.SD, function(x) all(x == "NA"))))]
    table[, (empty) := NULL]
    
    ### Arguably fill is way more important than color, so maybe change this
    if (!all(color == "NA")) {
      color[color == "NA"] <- "grey80"
      table$fill <- color
    }
    
    # SIMPLE ERROR HANDLING
    if (dim(table)[1] == 0) {
      warning("No features were found, returning NULL")
      return (NULL)
    }
    
    # Go to next function
    .read_dna_seg(table, seg_label, ...)
    
  # End of not PTT
  }

}

#' @name read_dna_seg
#' @export
#' 
# Added for support of reading many files at once, including globs
read_dna_seg_from_files <- function(
  files,
  tagsToParse = c("CDS"),
  fileType = "detect",
  meta_lines = 2,
  gene_type = "auto",
  header = TRUE,
  extra_fields = NULL,
  boundariesToParse = NULL,
  read_sequence = FALSE,
  verbose = FALSE,
  ...
) {
  
  filenames <- c()
  for (i in files) {
    filenames <- c(filenames, Sys.glob(i))
  }
  if (! length(filenames)) {
    stop("Could not find any files at the location(s) specified")
  }
  dna_segs <- list()
  name_list <- list()
  for (i in 1:length(filenames)) {
    file <- filenames[i]
    file_base <- gsub("\\.[^\\.]*$", "", basename(file))
    if (any(names(name_list) == file_base)) { # To deal with duplicate names
      if (verbose) {
        cat_time('Reading in: "', file, '" using the seg_label: "',
                 file_base, "_", name_list[[file_base]], '"'
        )
      }
      name_list[[file_base]] <- name_list[[file_base]] + 1
      dna_segs[[paste0(file_base, "_", name_list[[file_base]])]] <-
        read_dna_seg_from_file(file = file,
                               tagsToParse = tagsToParse,
                               fileType = fileType,
                               meta_lines = meta_lines,
                               gene_type = gene_type,
                               header = header,
                               extra_fields = extra_fields,
                               boundariesToParse = boundariesToParse,
                               read_sequence = read_sequence,
                               verbose = FALSE,
                               ...)
    } else {
      if (verbose) {
        cat_time('Reading in: "', file, '" using the seg_label: "', 
                 file_base, '"'
        )
      }
      name_list[[file_base]] <- 1
      dna_segs[[file_base]] <- read_dna_seg_from_file(
        file = file,
        tagsToParse = tagsToParse,
        fileType = fileType,
        meta_lines = meta_lines,
        gene_type = gene_type,
        header = header,
        extra_fields = extra_fields,
        boundariesToParse = boundariesToParse,
        read_sequence = read_sequence,
        verbose = FALSE,
        ...
      )
    }
  }
  
  dna_segs
}

read_dna_seg_from_gff <- function(
  file,
  tagsToParse = c("CDS"),
  fileType = "detect",
  meta_lines = 2,
  gene_type = "auto",
  header = TRUE,
  extra_fields = NULL,
  boundariesToParse = NULL,
  read_sequence = FALSE,
  verbose = FALSE,
  ...
) {
  
  if (verbose) {
    if (is.character(file)) {
      cat_time("Reading in: ", file)
    } else if (any(class(file) == "connection")) {
      cat_time("Reading in: ", summary(file)$description)
    }
  }
  
  # Import data from file, get rid of empty lines
  raw_lines <- readLines(file)
  raw_lines <- raw_lines[raw_lines != ""]
  
  # Search for sequence header lines
  seqstart_inds <- grep("^>", raw_lines)
  if (length(seqstart_inds) > 0) {
    lines <- raw_lines[1:(seqstart_inds[1] - 1)]
    
    if (read_sequence) {
      # Take the sequence portion of the file, then get rid of any comment lines
      fasta_lines <- grep("^[^#]",
                          raw_lines[seqstart_inds[1]:length(raw_lines)],
                          value = TRUE
                          )
      seqstart_inds <- grep("^>", fasta_lines)
      seqend_inds <- numeric(length(seqstart_inds))
      for (i in 1:(length(seqstart_inds) - 1)) {
        seqend_inds[i] <- seqstart_inds[i+1] - 1
      }
      seqend_inds[length(seqend_inds)] <- length(fasta_lines)
    }
  } else {
    if (read_sequence) {
      stop("'read_sequence' was set to TRUE, but no sequences could be found")
    }
    lines <- raw_lines
    fasta_lines <- character(0)
  }
  
  comment_lines <- grep("^##", lines, value = TRUE)
  lines <- grep("^[^#]", lines, value = TRUE)
  
  # SIMPLE ERROR HANDLING
  if (length(lines) == 0) {
    stop("No data is left after removing any comments or FASTA sequences.")
  }
  
  # Create a table from the file, then start processing the attributes column
  table <- fread(text = lines, sep = "\t", fill = TRUE, header = FALSE,
                 col.names = c("seq_origin", "source", "feature", "start",
                               "end", "score", "strand", "phase", "attributes"
                               )
                 )
  attributes_raw <- strsplit(table[, attributes], ";", fixed = TRUE)
  all_fields <- unique(sub("=.*", "", unlist(attributes_raw)))
  query_fields <- c("ID", "Name", "Alias", "Parent", "Target", "Gap",
                    "Derives_from", "Note", "Dbxref", "Ontology_term",
                    "Is_circular", "locus_tag", "gene", extra_fields)
  # query_fields <- c("ID", "Name", "Alias", "Parent", "Target", "Gap",
  #                   "Derives_from", "Note", "Dbxref", "Ontology_term",
  #                   "Is_circular", "locus_tag", "gene", "function",
  #                   "product", extra_fields)
  parsed_fields <- list()
  for (field in query_fields) {
    if (any(field == all_fields)) {
      # These are the fields to parse
      # Make a list of vectors that belong to this field
      field_list <- lapply(
        attributes_raw,
        function(x) sub(
          ".*=",
          "",
          grep(paste0("^", field, "\\="), x, fixed = FALSE, value = TRUE)
        )
      )
      # Fill up empty character vectors with "NA", then URLdecode 
      field_list <- unlist(
        lapply(field_list,
               function(x) ifelse(length(x) == 0, "NA", URLdecode(x))
               )
      )
      parsed_fields[[field]] <- field_list
    }
  }
  
  ### Currently no support for exons and introns yet
  
  ### Boundaries
  
  ### Strand to 1 or -1 
  
 
  # Go to next function
  .read_dna_seg(table, seg_label, ...)
  
}

#' @name read_dna_seg
#' @export
#' 
# Added for support, to read an embl file directly
read_dna_seg_from_embl <- function(
  file,
  tagsToParse = c("CDS"),
  boundariesToParse = NULL,
  extra_fields = NULL,
  read_sequence = FALSE,
  gene_type = "auto",
  verbose = FALSE,
  ...
) {
  read_dna_seg_from_file(file = file,
                         tagsToParse = tagsToParse,
                         boundariesToParse = boundariesToParse,
                         extra_fields = extra_fields,
                         read_sequence = read_sequence,
                         gene_type = gene_type,
                         fileType = "embl",
                         verbose = verbose,
                         ...
                         )
}

#' @name read_dna_seg
#' @export
#' 
# Added for support, to read a genbank file directly
read_dna_seg_from_genbank <- function(
  file,
  tagsToParse = c("CDS"),
  boundariesToParse = NULL,
  extra_fields = NULL,
  read_sequence = FALSE,
  gene_type = "auto",
  verbose = FALSE,
  ...
) {
  read_dna_seg_from_file(file = file,
                         tagsToParse = tagsToParse,
                         boundariesToParse = boundariesToParse,
                         extra_fields = extra_fields,
                         read_sequence = read_sequence,
                         gene_type = gene_type,
                         fileType = "genbank",
                         verbose = verbose,
                         ...
                         )
}

#' @name read_dna_seg
#' @export
#' 
read_dna_seg_from_fasta <- function(file, read_sequence = FALSE, ...) {
  
  # read data
  data <- readLines(file)
  if (length(grep("^>", data)) < 1) {
    stop("'file' does not seem to be a valid fasta file")
  }
  entries <- grep("^>", data)
  # Define variables for storage
  name <- character(length(entries))
  start <- numeric(length(entries))
  end <- numeric(length(entries))
  strand <- numeric(length(entries))
  gi <- character(length(entries))
  uniprot_id <- character(length(entries))
  gene <- character(length(entries))
  locus_id <- character(length(entries))
  product <- character(length(entries))
  proteinid <- character(length(entries))
  feature <- character(length(entries))
  featLength <- numeric(length(entries))
  sequence <- character(length(entries))
  seq_type <- character(length(entries))
  seq_name <- character(length(entries))
  metadata <- character(length(entries))

  curEnd <- 0
  
  
  # Under the assumption that all sequences in the FASTA file follow the same
  # format, first line is tested to check for supported formats, which are
  # UniProtKB, ensembl, and some basic support for NCBI.
  
  # Matches the start of a UniProtKB entry and its unique accession code
  uniprot_regex <- paste0("^>(sp|tr)\\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z]",
                          "[0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})(\\-[0-9]+)?\\|")
  if (grepl(uniprot_regex, data[entries[1]])) {
    TYPE <- "uniprot"
  # Matches the coordinate system from ensembl
  } else if (grepl(
    "[^:[:^graph:]]+\\:[^:[:^graph:]]+\\:[[:alnum:]]+\\:[0-9]+\\:[0-9]+\\:(1|\\-1)",
    data[entries[1]],
    perl = TRUE
  )) {
    TYPE <- "ensembl"
  # Other fasta format files. Mostly tailored to NCBI.
  } else {
    TYPE <- "unknown"
  }

  for (i in 1:length(entries)) {
    # Extract name of the feature by removing ">" and whitespace directly after
    header <- data[entries[i]]
    header <- trimws(substr(header, 2, nchar(header)))
    if (read_sequence) {
      ifelse (
        i == length(entries),
        sequence[i] <- paste0(data[(entries[i]+1):length(data)], collapse = ""),
        sequence[i] <- paste0(data[(entries[i]+1):(entries[i+1]-1)],
                              collapse = ""
                              )
      )
      featLength[i] <- nchar(sequence[i])
    } else {
      ifelse (
        i == length(entries),
        featLength[i] <- sum(nchar(data[(entries[i]+1):length(data)])),
        featLength[i] <- sum(nchar(data[(entries[i]+1):(entries[i+1]-1)]))
      )
    }
    
    
    if (TYPE == "uniprot") {
      header_split <- unlist(strsplit(header, "|", fixed = TRUE))
      uniprot_id[i] <- header_split[2]
      locus_id[i] <- header_split[2]
      name[i] <- unlist(strsplit(header_split[3], " ", fixed = TRUE))[1]
      
      # Take remains of header, keep removing the final part after parsing it
      remaining <- sub(paste0(name[i], " "), "", header_split[3], fixed = TRUE)
      for (j in 1:length(gregexpr("=" , header_split[3], fixed = TRUE)[[1]])) {
        extract <- regmatches(remaining,
                              gregexpr("[A-Z][A-Z]\\=[^\\=]+$", remaining)
                              )[[1]]
        remaining <- sub(paste0(" ", extract), "", remaining, fixed = TRUE)
        if (i == 1 && grepl("OS=", extract, fixed = TRUE)) {
          seg_label <- unlist(strsplit(extract, "=", fixed = TRUE))[2]
        } else if (grepl("GN=", extract, fixed = TRUE)) {
          gene[i] <- unlist(strsplit(extract, "=", fixed = TRUE))[2]
        }
      }
      # What is left now is only the protein product
      product[i] <- remaining
      
      start[i] <- curEnd + 1
      end[i] <- curEnd + featLength[i]
      strand[i] <- 1
      
      curEnd <- curEnd + featLength[i]
      

    } else if (TYPE == "ensembl") {
      
      loc_regex <- paste0("[^:[:^graph:]]+\\:[^:[:^graph:]]+\\:[[:graph:]]+",
                          "\\:[0-9]+\\:[0-9]+\\:(1|\\-1)")
      loc <- regmatches(header, gregexpr(loc_regex, header, perl = TRUE))[[1]]
      split_loc <- unlist(strsplit(loc, ":", fixed = TRUE))
      if (i == 1) seg_label <- split_loc[2]
      seq_type[i] <- split_loc[1]
      seq_name[i] <- split_loc[3]
      start[i] <- as.numeric(split_loc[4])
      end[i] <- as.numeric(split_loc[5])
      strand[i] <- as.numeric(split_loc[6])
      
      split_header <- unlist(strsplit(header, " ", fixed = TRUE))

      # True for more elaborate FASTA headers, who start with a name and type
      if (length(split_header) >= 3 && split_header[1] != split_loc[3] ) {
        name[i] <- split_header[1]
        feature[i] <- split_header[2]
       
        # Loop over remaining metadata
        metaQuery <- grepl(
          "^(gene|transcript|gene_biotype|transcript_biotype|gene_symbol|description)\\:",
          split_header
        )
        while (any(metaQuery)) {
          # Finds last piece of metadata
          indice <- max(which(metaQuery)):length(metaQuery)
          meta <- split_header[indice]
          if (grepl("^description\\:", meta[1])[1]) {
            meta[1] <- gsub("^description\\:", "", meta[1])
            product[i] <- paste(meta, collapse = " ")
          } else if (grepl("^gene_symbol\\:", meta[1])[1]) {
            meta[1] <- gsub("^gene_symbol\\:", "", meta[1])
            gene[i] <- paste(meta, collapse = " ")
          } else if (grepl("^gene\\:", meta[1])[1]) {
            meta[1] <- gsub("^gene\\:", "", meta[1])
            locus_id[i] <- paste(meta, collapse = " ")
          }
          split_header <- split_header[-indice]
          metaQuery <- metaQuery[-indice]
        }
        
      } else {
        name[i] <- paste(split_loc[3], split_loc[4], split_loc[5],
                         split_loc[6], sep = ":")
        locus_id[i] <- name[i]
      }
      

    } else {
      brackets <- regmatches(header, gregexpr("\\[[^\\[]+\\]", header))[[1]]
      name[i] <- unlist(strsplit(header, " ", fixed = TRUE))[1]

      # Attempts to read organism name
      if (i == 1) {
        # Matches "[organism=Organism name]"
        extract <- regmatches(header,
                              gregexpr("\\[organism\\=[^\\[]+\\]", header)
                              )[[1]]
        if (length(extract > 0)) {
          seg_label <- gsub("(\\[organism\\=|\\])", "", extract)
        } else {
          # Matches "[anything]", could go wrong, but for NCBI proteins, when
          # there is only one set of square brackets, they contain organism name
          if (length(brackets) == 1) {
            seg_label <- gsub("(\\[|\\])", "", brackets)
          }
        }
      }
      # Attempt to parse information in square brackets
      if (length(brackets) > 0) {
        for (j in 1:length(brackets)) {
          if (grepl("(name\\=|protein\\=)", brackets[j]) &&
              length(product[i]) != 0
              ) {
            product[i] <- gsub("(\\[name\\=|\\[protein\\=|\\])",
                               "",
                               brackets[j]
                               )
          } else if (grepl("(protein_accession\\=|protein_id\\=)",
                           brackets[j]
                           ) && length(proteinid[i]) != 0
                     ) {
            proteinid[i] <- gsub(
              "(\\[protein_accession\\=|\\[protein_id\\=|\\])",
              "",
              brackets[j]
            )
          } else if (grepl("locus_tag\\=", brackets[j]) &&
                     length(locus_id[i]) != 0
                     ) {
            locus_id[i] <- gsub("(\\[locus_tag\\=|\\])", "", brackets[j])
          } else if (grepl("GeneID\\=", brackets[j]) && 
                     length(locus_id[i]) != 0
                     ) {
            gi[i] <- gsub("(\\[GeneID\\=|\\])", "", brackets[j])
          }
        }
      } else {
        # Arriving here means the FASTA format was unrecognized
        # So all the remaining information is just stored in 1 column
        metadata[i] <- sub(paste0(name[i], " "), "", header, fixed = TRUE)
      }
      
      if (locus_id[i] == "") locus_id[i] <- name[i]
      
      start[i] <- curEnd + 1
      end[i] <- curEnd + featLength[i]
      strand[i] <- 1
      
      curEnd <- curEnd + featLength[i]
    }

  }
  # Fallback option in case seg label could not be identified
  if (!exists("seg_label")) {
    seg_label <- gsub("\\.[^\\.]*$", "", basename(file))
  }

  if (TYPE == "ensembl") {
    table <- data.table(name = name, start = start,
                        end = end, strand = strand,
                        locus_id = locus_id, gene = gene,
                        product = product, proteinid = proteinid,
                        uniprot_id = uniprot_id, gi = gi,
                        feature = feature, length = featLength,
                        sequence = sequence, seq_type = seq_type,
                        seq_name = seq_name
                        )
    # Create sequence boundaries, and attach them to the table
    
    # Currently this has sequences starting and ending by the minimum and 
    # maximum values for start and end, but they should arguably start at
    # 1 unless minimum is negative
    boundaries <- table[,
                        range(start, end),
                        by = c("seq_type", "seq_name")
                        ][,
                          .(name = ifelse(.I %% 2 == 1,
                                          paste0("start_seq_", seq_name),
                                          paste0("end_seq_", seq_name)
                                          ),
                            start = V1, end = V1, strand = 1,
                            locus_id = ifelse(.I %% 2 == 1,
                                              paste0("start_seq_", seq_name),
                                              paste0("end_seq_", seq_name)
                                              ),
                            gene = "NA", product = "NA",
                            proteinid = "NA", uniprot_id = "NA",
                            gi = "NA", feature = seq_type,
                            length = 1, sequence = "NA",
                            seq_type, seq_name,
                            region_plot = ifelse(.I %% 2 == 1, "start", "end")
                            )
                          ]
    table <- rbindlist(list(table, boundaries), use.names = TRUE, fill = TRUE)
    table[, seq_origin := paste0(seq_type, ":", seq_name)]
    
    # Order on sequence name, numerics first
    # to get chromosomes before scaffolds/ contigs/ etc.
    order <- gtools::mixedorder(table$seq_origin)
    table <- table[order][, let(seq_name = NULL, seq_type = NULL)]
    # Finally, update start and end values for each sequence
    curEnd <- 0
    for (i in unique(table$seq_origin)) {
      table[seq_origin == i, let(start = start + curEnd, end = end + curEnd)]
      curEnd <- table[seq_origin == i, max(end)]
    }
    
  } else {
    table <- data.table(name = name, start = start,
                        end = end, strand = strand,
                        metadata = metadata, locus_id = locus_id,
                        gene = gene, product = product,
                        proteinid = proteinid, uniprot_id = uniprot_id,
                        gi = gi, sequence = sequence
                        )
  }
  
  # Fill empty sequences with "-" since NA could theoretically be a sequence
  if (read_sequence & !is.null(table$sequence)) {
    set(
      table,
      i = which(table$sequence == "" & 
                  !table$region_plot %ilike% "(^start$|^end$)"),
      j = "sequence",
      value = "-"
    )
  }
  
  # Fill empty string with "NA"
  for (col in names(table)) {
    set(table,
        i = which(table[[col]] %in% c("", "NA") | is.na(table[[col]])),
        j = col,
        value = "NA"
        )
  }
  
  # Cut table to include only added features
  empty <- table[, names(which(sapply(.SD, function(x) all(x == "NA"))))]
  table[, (empty) := NULL ]
  
  .read_dna_seg(table, seg_label, ...)
}

#' @name read_dna_seg
#' @export
#' 
# reading genes from a file. Use source=tab or ptt to specify type
read_dna_seg_from_ptt <- function(file, meta_lines = 2, header = TRUE, ...) {
  # reads meta info
  seg_label <- readLines(file, n = 1)
  seg_label <- strsplit(seg_label, "/,|-/", fixed = TRUE)[[1]][1]
  # reads ptt table
  ptt <- read.table(file,
                    skip = meta_lines,
                    as.is = TRUE,
                    header = header,
                    sep = "\t",
                    quote = ""
                    )
  if (header) {
    names(ptt) <- tolower(names(ptt))
  } else {
    names(ptt) <- c("location", "strand", "length", "pid", "gene",
                    "synonym", "code", "cog", "product")
  }
  # parse location
  location <- strsplit(ptt$location, "..", fixed = TRUE)
  start <- as.numeric(sapply(location, function(x) x[[1]]))
  end <- as.numeric(sapply(location, function(x) x[[2]]))
  # parse strand
  strand <- ptt$strand
  strand[strand == "-"] <- -1
  strand[strand == "+"] <- 1
  strand <- as.numeric(strand)
  # parse gene name from name or synonym if not present
  name <- ifelse(ptt$gene == "-", ptt$synonym, ptt$gene)
  table <- data.frame(name = name, start = start,
                      end = end, strand = strand,
                      length = ptt$length, pid = ptt$pid,
                      gene = ptt$gene, synonym = ptt$synonym,
                      code = ptt$code, cog = ptt$cog,
                      product = ptt$product,
                      stringsAsFactors = FALSE
                      )
  .read_dna_seg(table, seg_label, ...)
}

#' @name read_dna_seg
#' @export
#' 
read_dna_seg_from_tab <- function(file, header = TRUE, ...) {
  table <- read.table(file,
                      as.is = TRUE,
                      header = header,
                      sep = "\t",
                      quote = ""
                      )
  if (ncol(table) < 4) stop("Insufficient number of columns in table")
  col_names <-  c("name", "start", "end", "strand")
  if (all(col_names %in% names(table))) {
    table <- as.data.table(table)
    setcolorder(table, neworder = col_names)
  } else {
    names(table)[1:length(col_names)] <- col_names
  }
  # parse name from file name by default
  seg_label <- gsub("\\.[^\\.]*$", "", basename(file))
  .read_dna_seg(table, seg_label, ...)
}

# Internal function, final checks before making dna_seg
.read_dna_seg <- function(table, seg_label, reverse = FALSE, xlim = NULL, ...) {
  # check args
  if (ncol(table) < 4) stop("Insufficient number of columns in table")
  if (nrow(table) < 1) stop("No lines in table")
  col_names <-  c("name", "start", "end", "strand")
  if (!all(col_names %in% names(table))) {
    stop("'table' must contain at least columns ",
         "'name', 'start', 'end', 'strand'")
  }
    
  # make dna_seg object, set seg_label attribute
  dna_seg <- as.dna_seg(table, ...)
  dna_seg <- trim.dna_seg(dna_seg, xlim)
  if (reverse) dna_seg <- reverse.dna_seg(dna_seg)
  setattr(dna_seg, "seg_label", seg_label)
  dna_seg
}

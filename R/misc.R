################################################################################
# Misc. internal functions
################################################################################

# Internal helper function to get dna_seg labels
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

# Prints time, used for verbose arguments
cat_time <- function(...) {
  cat(paste0(format(Sys.time()), ": ", ..., "\n"))
}

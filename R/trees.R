################################################################################
# Phylo tree
################################################################################
# a function to plot phylogenetic trees ported to grid from ade4.
# to be rewritten...
phylog_grob <- function (
  x,
  y = NULL,
  f.phylog = 1,
  cleaves = 0,
  cnodes = 0,
  labels.leaves = names(x$leaves),
  clabel.leaves = 1,
  labels.nodes = names(x$nodes),
  clabel.nodes = 0.8,
  labels.col,
  leave.seg.col = grey(0.7),
  tree.scale = FALSE,
  sub = "",
  csub = 1.25,
  possub = "bottomleft",
  draw.box = FALSE,
  ...
) {
  
  if (!inherits(x, "phylog")) stop("'x' must be a phylog object")
  leaves.number <- length(x$leaves)
  leaves.names <- names(x$leaves)
  nodes.number <- length(x$nodes)
  nodes.names <- names(x$nodes)
  if (length(labels.leaves) != leaves.number) labels.leaves <- names(x$leaves)
  if (length(labels.nodes) != nodes.number) labels.nodes <- names(x$nodes)
  leaves.car <- gsub("[_]", " ", labels.leaves)
  nodes.car <- gsub("[_]", " ", labels.nodes)
  if (f.phylog < 0.05) f.phylog <- 0.05
  if (f.phylog > 0.95) f.phylog <- 0.95
  maxx <- max(x$droot)
  x.leaves <- x$droot[leaves.names]
  x.nodes <- x$droot[nodes.names]
  annot_nodes <- grep("^[^I]", names(x.nodes), value = TRUE)
  if (!is.null(y)) {
    # check that it is constrained between 0 and 1
    if (!all(range(y) %in% c(0,1))) stop("'y' must be a number between 0 and 1")
  } else {
    y <- seq(1, 0, len = leaves.number)
  }
  names(y) <- leaves.names
  xcar <- maxx * 1.05
  xx <- c(x.leaves, x.nodes)
  # prepare grobs
  labelGrobs <- gList()
  labelSegGrobs <- gList()
  branchLabelGrobs <- gList()
  branchesGrobs <- gList()
  # leaves labels and segments leading to it
  if (clabel.leaves > 0) {
    for (i in 1:leaves.number) {
      labelGrobs[[i]] <- textGrob(
        x = 0, y = y[i],
        label = leaves.car[i], just = "left",
        gp = gpar(cex = clabel.leaves, col = labels.col[i]),
        default.units = "native",
        name = paste("tree.leave.label", i, sep = ".")
      )
      labelSegGrobs[[i]] <- segmentsGrob(
        xcar, y[i], xx[i], y[i],
        gp = gpar(col = leave.seg.col), default.units = "native",
        name = paste("tree.leave.segment", i, sep = ".")
      )
    }
  }
  yleaves <- y[1:leaves.number]
  xleaves <- xx[1:leaves.number]
  ## if (cleaves > 0) {
  ##     for (i in 1:leaves.number) {
  ##         points(xx[i], y[i], pch = 21, bg = 1, cex = par("cex") *
  ##             cleaves)
  ##     }
  ## }
  yn <- rep(0, nodes.number)
  names(yn) <- nodes.names
  y <- c(y, yn)
  # plot tree
  for (i in 1:length(x$parts)) {
    w <- x$parts[[i]]
    but <- names(x$parts)[i]
    y[but] <- mean(y[w])
    b <- range(y[w])
    # vertical branches
    branchesGrobs[[i]] <- segmentsGrob(
      xx[but], b[1], xx[but], b[2],
      default.units = "native",
      name = paste("tree.branch.vert.segment", i, sep = ".")
    )
    x1 <- xx[w]
    y1 <- y[w]
    x2 <- rep(xx[but], length(w))
    # horizontal branches
    branchesGrobs[[i+length(x$parts)]] <- segmentsGrob(
      x1, y1, x2, y1,
      default.units = "native",
      name = paste("tree.branch.horiz.segment", i, sep = ".")
    )
    idx_node_label <- names(x1) %in% annot_nodes
    if (clabel.nodes > 0 && sum(idx_node_label) > 0) {
      label <- names(x1)[idx_node_label]
      label <- gsub("^X?([^_]+)(_[0-9]+)?$", "\\1", label, perl = TRUE)
      branchLabelGrobs[[i]] <- textGrob(
        label = label,
        x = (x1[idx_node_label] + x2[idx_node_label]) / 2,
        y = unit(y1[idx_node_label], "native") + unit(0.5, "line"),
        gp = gpar(cex = clabel.nodes), default.units = "native",
        name = paste("tree.branch.label", i, sep = ".")
      )
    }
  }
  # scale
  if (tree.scale) {
    rng <- diff(pretty(c(0, maxx), n = 2))[1]
    scaleGrob <- segmentsGrob(maxx, unit(-1.5, "lines"),
                              maxx - rng,  unit(-1.5, "lines"),
                              default.units = "native",
                              name = ("tree.scale.seg")
                              )
    scaleLabelGrob <- textGrob(label = rng,
                               x = maxx-(rng/2), y = unit(-1, "lines"),
                               default.units = "native",
                               name = ("tree.scale.label")
                               )
  } else {
    scaleGrob <- nullGrob()
    scaleLabelGrob <- nullGrob()
  }
  ## if (cnodes > 0) {
  ##     for (i in nodes.names) {
  ##         points(xx[i], y[i], pch = 21, bg = "white", cex = cnodes)
  ##     }
  ## }
  ## if (clabel.nodes > 0) {
  ##     scatterutil.eti(xx[names(x.nodes)], y[names(x.nodes)],
  ##         nodes.car, clabel.nodes)
  ## }
  x <- x.leaves
  y <- y[leaves.names]
  xbase <- xcar
  ## if (csub > 0)
  ##     scatterutil.sub(sub, csub = csub, possub = possub)
  ## if (draw.box)
  ##     box()
  ## if (cleaves > 0)
  ##     points(xleaves, yleaves, pch = 21, bg = 1, cex = par("cex") *
  ##         cleaves)
  # creating gTree for branches
  branchesTree <- gTree(
    children = gList(branchesGrobs, labelSegGrobs,
                     branchLabelGrobs, scaleGrob, scaleLabelGrob),
    vp = viewport(xscale = c(0, xcar), yscale = c(0, 1),
                  name = "tree.branches"),
    name = "tree.branchesTree"
  )
  labelTree <- gTree(
    children = labelGrobs,
    vp = viewport(xscale = c(0, 1), yscale = c(0, 1), name = "tree.labels"),
    name = "tree.labelsTree"
  )
  # finally plotting
  label_width <- unit(1, "grobwidth",
                      labelGrobs[[which.max(nchar(leaves.car))]])
  layout <- grid.layout(1, 2, widths = unit.c(unit(1, "null"), label_width))
  fg <- frameGrob(layout = layout, name = "treeFrameGrob",
                  vp = viewport(name = "treeFrame"))
  fg <- placeGrob(fg, branchesTree, col = 1)
  fg <- placeGrob(fg, labelTree, col = 2)
  return(invisible(list(xy = data.frame(x = x, y = y), xbase = xbase,
                        cleaves = cleaves, grob = fg, width = label_width)))
}


# permute tree leaves to match labels
permute_tree <- function(tree, labels, no.over = 100000) {
  # Convert "." to "_", as ade4 tree objects already do this automatically
  labels <- gsub(".", "_", labels, fixed = TRUE)
  ref <- names(tree$leaves)
  n <- length(ref)
  wanted_permut <- rep(NA, n)
  for (i in 1:n) {
    idx <- which(ref[i] == labels)
    if (length(idx) != 1) stop("Non-unique or non-matching labels")
    wanted_permut[i] <- idx
  }
  if (identical(as.numeric(1:n), as.numeric(wanted_permut))) {
    res <- 1:n
    names(res) <- labels
    return(res)
  }
  permuts <- ade4::enum.phylog(tree, no.over)
  if (is.null(permuts)) {
    stop("Number of permutations too large, use a tree that has labels in the ",
         "same order as dna_segs")
  }
  equals <- apply(permuts, 1, function(x)
                  identical(as.numeric(x), as.numeric(wanted_permut)))
  if (!any(equals)) {
    stop("No tree permutation compatible with label order. Change input order")
  }
  if (!sum(equals)) {
    stop("Several permutations matching. Something went wrong, ",
         "contact the author")
  }
    
  return(permuts[equals,])
}

#' Reorder dna_segs or labels to match a tree
#' 
#' Takes a list of `dna_seg` objects or `dna_seg` labels and reorganizes them
#' based on a given (phylogenetic) tree. 
#' 
#' @details
#' This function takes a character vector of `dna_seg` labels, either directly,
#' or by extracting them from a list of `dna_segs`, through the `dna_segs`
#' argument. Each of the labels is queried to find matching tree tip labels,
#' sorting them to match the order in which they are found in the tree. If
#' exactly 1 match is found, the `dna_seg` label is updated to match the tree
#' tip label, unless `exact_match = TRUE`. If multiple matches are found, an
#' error is returned that shows the offending `dna_seg` label.
#' 
#' @returns A list of `dna_seg` objects or a character vector of `dna_seg`
#' labels, matching the input given in the `dna_segs` argument.
#' @returns If `return_old_labels = TRUE`, a list with 2 named elements will 
#' be returned instead (`dna_segs`, the same return value as above, 
#' and `old_labels`, a character vector of the original labels that is now
#' sorted).
#' @export
#' 
#' @param dna_segs Either a character vector containing `dna_seg` labels, or a
#' list of `dna_seg` objects.
#' @param tree A (phylogenetic) tree, in the form of a `phylo` or `phylog`
#' object, or a character string containing a file path to a Newick tree format
#' file.
#' @param exact_match Logical. If `TRUE`, `dna_seg` labels will need to match
#' the labels of the tree exactly. If `exact_match = FALSE`, tree tip labels
#' only need to contain the `dna_seg` labels for a match to be found (e.g. the 
#' `dna_seg` label `"seq_1"` will match tree tip label `"E_coli_seq_1.fa"`).
#' @param return_old_labels Logical. If `TRUE`, then the `dna_seg` labels will
#' be returned using the original names provided by the `dna_segs` argument.
#' Only relevant when `exact_match = FALSE`, as this option can cause `dna_seg`
#' labels to be changed to match the tree tip labels.
#' 
#' @author Mike Puijk
#' 
#' @seealso [trim_tree], [plot_gene_map]
#' 
#' @examples
#' ## Generate data
#' seg_labels <- c("seq_2", "seq_3", "seq_1", "seq_4")
#' tree_str <- paste0("(seq_1_B_bacilliformis:0.5,",
#'                    "(seq_2_B_grahamii:0.1,",
#'                    "(seq_3_B_henselae:0.1,",
#'                    "seq_4_B_quintana:0.2):0.1):0.1);")
#' tree <- ade4::newick2phylog(tree_str)
#' 
#' ## Reorder and rename dna_seg labels to match tree
#' seg_labels
#' seg_labels <- permute_dna_segs(dna_segs = seg_labels, tree = tree)
#' seg_labels
#' 
permute_dna_segs <- function(
  dna_segs,
  tree,
  exact_match = FALSE,
  return_old_labels = FALSE
) {
  
  # Check mandatory arguments
  if (missing(dna_segs)) stop("'dna_segs' must be provided")
  if (is.character(dna_segs)) {
    seg_labels <- dna_segs
    return_list <- FALSE
  } else if (is.list(dna_segs) && all(sapply(dna_segs, is.dna_seg))) {
    seg_labels <- get_seg_labels(dna_segs, seg_labels = NULL)
    return_list <- TRUE
    if (is.null(seg_labels)) {
      stop("dna_seg labels could not be determined from this dna_seg list")
    }
  } else {
    stop("'dna_segs' must be a list of dna_seg objects or a character ",
         "vector containing the seg labels (i.e. names(dna_segs))")
  }
  
  # Read tree to determine the order of dna_segs and comparisons
  # If tree argument is a character string, it is assumed to be a file path
  if (missing(tree)) stop("'tree' must be provided")
  if (is.character(tree) && file.exists(tree)) {
    treeNames <- names(ade4::newick2phylog(readLines(tree))$leaves)
    # If not a file path, then check for known tree objects
  } else if (inherits(tree, "phylog")) {
    treeNames <- names(tree$leaves)
  } else if (inherits(tree, "phylo")) {
    treeNames <- tree$tip.label
  } else {
    stop("'tree' must be either a character string (specifying a file path) ",
         "or an object of class phylo or phylog")
  }
  
  new_labels <- numeric()
  if (return_list | return_old_labels) {
    old_labels <- numeric()
  }
  
  for (i in seg_labels) {
    if (exact_match) {
      # query will only go through if there is only 1 exact match
      query <- grep(paste0("\\<", i, "\\>"), treeNames)
      if (length(query) == 0) {
        stop('"', i, '" could not be found in the tree tip labels')
      } else if (length(query) > 1) {
        stop('Multiple matches found for "', i, '" in the tree tip labels, ',
             'make sure the names of dna_segs match the names of the provided ',
             'tree exactly'
             )
      }
    } else {
      query <- grep(i, treeNames)
      if (length(query) == 0) {
        stop(i, " could not be found in the tree tip labels")
        # If there are multiple matches, look for an exact match instead
      } else if (length(query) > 1) {
        query <- grep(paste0("\\<", i, "\\>"), treeNames)
        if (length(query) != 1) {
          stop('Multiple matches found for "', i, '" in the tree tip labels, ',
               'make sure the names of dna_segs match the names of the ',
               'provided tree exactly')
        }
      }
    }
    # Fill vector with indices, with their tree leaves as names
    new_labels[treeNames[as.numeric(query)]] <- as.numeric(query)
    
    if (return_list | return_old_labels) {
      old_labels[i] <- as.numeric(query)
    }
  }
  
  # Sort the labels based on the indices, which is the order of the tree leaves
  sorted_labels <- sort(new_labels)
  sorted_labels <- names(sorted_labels)
  if (return_list | return_old_labels) {
    sorted_old <- sort(old_labels)
    sorted_old <- names(sorted_old)
  }
  
  if (return_list) {
    # Fill a new list of dna_segs using the tree order
    sorted_return <- list()
    for (i in 1:length(sorted_labels)) {
      sorted_return[[sorted_labels[i]]] <- dna_segs[[sorted_old[i]]]
    }
  } else {
    # Sort the tree leaves by their indices and return ordered dna_seg labels
    sorted_return <- sorted_labels
  }
  if (return_old_labels) {
    return(list(dna_segs = sorted_return, old_labels = sorted_old))
  } else {
    return(sorted_return)
  }
}

#' Trim a tree to remove unused sequences
#' 
#' Takes a (phylogenetic) tree and removes all tree tips that are not found
#' in a given set of names. 
#' 
#' @details
#' This function takes a character vector of `dna_seg` labels, either directly,
#' or by extracting them from a list of `dna_segs`, through the `dna_segs`
#' argument. Each of the labels is queried to find matching tree tip labels,
#' and any tree tip label without a match will removed. If multiple matches are
#' found for a single `dna_seg` label, an
#' error is returned that shows the offending `dna_seg` label.
#' 
#' @returns A phylog object containing the filtered tree.
#' @export
#' 
#' @param dna_segs Either a character vector containing `dna_seg` labels, or a
#' list of `dna_seg` objects.
#' @param tree A (phylogenetic) tree, in the form of a `phylo` or `phylog`
#' object, or a character string containing a file path to a Newick tree format
#' file.
#' @param exact_match Logical. If `TRUE`, `dna_seg` labels will need to match
#' the labels of the tree exactly. If `exact_match = FALSE`, tree tip labels
#' only need to contain the `dna_seg` labels for a match to be found (e.g. the 
#' `dna_seg` label `"seq_1"` will match tree tip label `"E_coli_seq_1.fa"`).
#' 
#' @author Mike Puijk
#' 
#' @seealso [permute_dna_segs], [plot_gene_map]
#' 
#' @examples
#' ## Generate data
#' names <- c("seq_1", "seq_2", "seq_3")
#' tree_str <- paste0("(seq_1_B_bacilliformis:0.5,",
#'                    "(seq_2_B_grahamii:0.1,",
#'                    "(seq_3_B_henselae:0.1,",
#'                    "seq_4_B_quintana:0.2):0.1):0.1);"
#'                    )
#' tree <- ade4::newick2phylog(tree_str)
#' 
#' ## Filter tree
#' tree$tre
#' tree <- trim_tree(dna_segs = names, tree = tree)
#' tree$tre
#' 
trim_tree <- function(dna_segs, tree, exact_match = FALSE) {
  # Check mandatory arguments
  if (missing(dna_segs)) stop("'dna_segs' must be provided")
  if (is.character(dna_segs)) {
    seg_labels <- dna_segs
  } else if (is.list(dna_segs) && all(sapply(dna_segs, is.dna_seg))) {
    seg_labels <- get_seg_labels(dna_segs, seg_labels = NULL)
    if (is.null(seg_labels)) {
      stop("dna_seg labels could not be determined from this dna_seg list")
    }
  } else {
    stop("'dna_segs' must be a list of dna_seg objects or a character ",
         "vector containing the names (i.e. names(dna_segs))")
  }
  
  # Check tree argument, convert to phylo object if necessary
  if (missing(tree)) stop("'tree' must be provided")
  if (is.character(tree) && file.exists(tree)) {
    tree <- ape::read.tree(tree)
    # If not a file path, then check for known tree objects
  } else if (inherits(tree, "phylog")) {
    tree <- ape::read.tree(text = tree$tre)
  } else if (inherits(tree, "phylo")) {
    tree <- tree
  } else {
    stop("'tree' must be either a character string (specifying a file path) ",
         "or an object of class phylo or phylog")
  }
  
  toKeep <- character()
  for (i in seg_labels) {
    if (exact_match) {
      # query will only go through if there is exactly 1 exact match
      query <- grep(paste0("\\<", i, "\\>"), tree$tip.label, value = TRUE)
      if (length(query) == 0) {
        stop('"', i, '" could not be found in the tree tip labels')
      } else if (length(query) > 1) {
        stop('Multiple matches found for "', i, '" in the tree tip labels, ',
             'make sure the names of dna_segs match the names of the provided ',
             'tree exactly'
             )
      }
    } else {
      query <- grep(i, tree$tip.label, value = TRUE)
      if (length(query) == 0) {
        stop(i, " could not be found in the tree tip labels")
        # If there are multiple matches, look for an exact match instead
      } else if (length(query) > 1) {
        query <- grep(paste0("\\<", i, "\\>"), tree$tip.label, value = TRUE)
        if (length(query) != 1) {
          stop('Multiple matches found for "', i, '" in the tree tip labels, ',
               'make sure the names of dna_segs match the names of the ',
               'provided tree exactly'
               )
        }
      }
    }
    # Fill character vector with tree tip names to keep
    toKeep[query] <- query
  }
  if (length(tree$tip.label) > length(toKeep)) {
    remove <- tree$tip.label[grep(
      paste(toKeep, collapse = "|"),
      tree$tip.label,
      invert = TRUE
    )]
    tree <- ape::drop.tip(tree, remove)
  }
  
  # convert to phylog object
  tree <- paste(ape::write.tree(tree), collapse = "\n")
  tree <- ade4::newick2phylog(tree)
  tree
}

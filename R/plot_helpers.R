################################################################################
# Plot helpers
################################################################################
# calculate arrow coordinates from gene coordinates
arrow_coord <- function(x1,
                        x2,
                        y = 0.5,
                        strand = NULL,
                        width = 1,
                        head_len = 100) {
  # take care of strand, to get x1 as bottom and x2 as tip of arrow
  if (!is.null(strand) && strand == -1) {
    x_temp <- x2
    x2 <- x1
    x1 <- x_temp
  }
  w2 <- width / 4
  # if the head of the arrow is larger than half of the gene, reduce to half
  if (head_len > abs(x1 - x2) / 2) {
    head_len <- abs(x1 - x2) / 2
  }
  # calculate xi, x "internal"
  if (x2 > x1) {
    xi <- x2 - head_len
  } else {
    xi <- x2 + head_len
  }
  list(x = c(x1,     xi,     xi,         x2, xi,         xi,     x1),
       y = c(y - w2, y - w2, y - w2 * 2, y,  y + w2 * 2, y + w2, y + w2)
       )
}
# calculate arrow coordinates from gene coordinates
headless_arrow_coord <- function(x1,
                                 x2,
                                 y = 0.5,
                                 strand = NULL,
                                 width = 0.7,
                                 head_len = 100) {
  # take care of strand, to get x1 as bottom and x2 as tip of arrow
  if (!is.null(strand) && strand == -1) {
    x_temp <- x2
    x2 <- x1
    x1 <- x_temp
  }
  w2 <- width/2
  # if the head of the arrow is larger than half of the gene, reduce to half
  if (head_len > abs(x1 - x2) / 2) {
    head_len <- abs(x1 - x2) / 2
  }
  # calculate xi, x "internal"
  if (x2 > x1) {
    xi <- x2-head_len
  } else {
    xi <- x2+head_len
  }
  list(x = c(x1,     xi,     x2, xi,     x1),
       y = c(y - w2, y - w2, y,  y + w2, y + w2)
  )
}
# coords for a block
block_coord <- function(start, end, strand, y = 0.5) {
  x <- c(rep(start, 2), rep(end, 2))
  y <- c(y, y + strand / 2, y + strand / 2, y)
  list(x = x, y = y)
}
# exon coord
exon_coord <- function(start, end, strand) {
  x <- c(rep(start, 2), rep(end, 2))
  if (strand == 0 ) {
    y <- c(0.2, 0.8, 0.8, 0.2)
  }
  if (strand == 1 ) {
    y <- c(0.5, 0.8, 0.8, 0.5)
  }
  if (strand == -1 ) {
    y <- c(0.2, 0.5, 0.5, 0.2)
  }
  list(x = x, y = y)
}
# coords for a zone annotation
bracket_coord <- function(start, end, y = 0, w = 0.1) {
  x <- c(rep(start, 2), rep(end, 2))
  y <- c(y, rep(y + w, 2), y)
  list(x = x, y = y)
}
# axis coords
yaxis_coords <- function(at, x0 = 0, x1 = 0.5) {
  n <- length(at)
  list(x0 = c(rep(x0, n), x1),
       x1 = rep(x1, n+1),
       y0 = c(at, at[1]),
       y1 = c(at, at[n]))
}
# calculate comparison coordinates
calc_comp_coor <- function(gap, xlim, comp, side) {
  if (length(gap) != nrow(xlim))
    stop("'gap' must have the same length as 'xlim'")
  if (side < 1 || side > 2) stop("'side' must be 1 or 2")
  # x is the moving cursor
  x <- 0
  old_start <- if (side == 1) comp$start1 else comp$start2
  old_end <- if (side == 1) comp$end1 else comp$end2
  start <- old_start
  end <- old_end
  for (i in 1:nrow(xlim)) {
    # increment by the gap length
    x <- x + gap[i]
    # select comps
    idx <- old_start >= xlim$x0[i] & old_end <= xlim$x1[i]
    # re-number by substracting the xlim and adding x
    if (xlim$strand[i] == 1) {
      start[idx] <- old_start[idx] - xlim$x0[i] + x
      end[idx] <- old_end[idx] - xlim$x0[i] + x
    } else {
      start[idx] <- xlim$x1[i] - old_start[idx] + x
      end[idx] <- xlim$x1[i] - old_end[idx] + x
    }
    # increment x by the length of the segment
    x <- x + xlim$length[i]
  }
  # reattribute start and stop
  if (side == 1) comp$start1 <- start else comp$start2 <- start
  if (side == 1) comp$end1 <- end else comp$end2 <- end
  # return the modified comp
  comp
}

#' Find the middle point of dna_seg features
#' 
#' Returns a vector containing the middle point of each feature in a `dna_seg`
#' object. Useful to prepare annotations, for example.
#' 
#' @returns A numeric vector.
#' @export
#' 
#' @param dna_seg A `dna_seg` object.
#' 
#' @author Lionel Guy
#' 
#' @seealso [annotation], [dna_seg]
#' 
#' @examples
#' ## Load data
#' data(barto)
#' 
#' ## Get middles of the first dna_seg
#' mid <- middle(barto$dna_segs[[1]])
#' 
middle <- function(dna_seg) {
  if (!is.dna_seg(dna_seg)) stop("'dna_seg' must be a dna_seg object")
  apply(dna_seg[, c("start", "end")], 1, mean)
}

#' Human-readable nucleotide scale
#' 
#' Returns a human readable list from a nucleotide position or length.
#' 
#' @details
#' Return a nucleotide value in nt, kb, Mb or Gb, according to the value given.
#' This is particularly useful to display nice scales without too many 
#' trailing zeros.
#' 
#' @returns Returns a list with 4 elements.
#' \item{n}{A numeric value corresponding to `nt` divided by `mult` 
#' (see below).}
#' \item{tag}{A character, giving the multiplier used in text.}
#' \item{mult}{The muliplier used, in numeric value.}
#' \item{text}{A character, giving the value in a human readable format.}
#' 
#' @export
#' 
#' @param nt A nucleotide position.
#' @param signif Either a numeric or logical. If `FALSE`, `nt` is not rounded.
#' If this argument is numeric, it returns that amount of significant digits.
#' 
#' @author Lionel Guy
#' 
#' @examples
#' human_nt(123456)
#' human_nt(123456, signif = 2)
#' human_nt(123456890, signif = 2)
#' 
human_nt <- function(nt, signif = FALSE) {
  tag <- "nt"
  mult <- 1
  med <- median(nt)
  if (med >= 1e9) {
    nt <- nt/1e9
    tag <- "Gb"
    mult <- 1e9
  } else if (med >= 1e6) {
    nt <- nt/1e6
    tag <- "Mb"
    mult <- 1e6
  } else if (med >= 1e3) {
    nt <- nt/1e3
    tag <- "kb"
    mult <- 1e3
  }
  if (signif) nt <- signif(nt, signif)
  list(n = nt, tag = tag, mult = mult, text = paste(nt, tag))
}

#' Artemis Colors
#' 
#' Returns a data frame with the standard artemis colors.
#' 
#' @returns
#' A `data.frame` with the following columns: `n`, `names`, `colors`, `r`, `g`,
#' and `b`. The 3 first columns give the Artemis color number, its name, and
#' its equivalent in R. The final 3 columns give the r, g, and b values.
#' @export
#' 
#' @author Lionel Guy
#' 
#' @references Artemis website:
#' https://www.sanger.ac.uk/tool/artemis/
#' 
#' @examples
#' artCol <- artemisColors()
#' plot(rep(1, nrow(artCol)), artCol$n, xlim = c(1, 2), type = "n")
#' text(rep(1, nrow(artCol)), artCol$n, labels = artCol$n, col = artCol$colors)
#' text(rep(1, nrow(artCol)), artCol$n, labels = artCol$names,
#'      col = artCol$colors, pos = 4, offset = 1)
#' 
## Emulate artemis colors
## 0  white          (RGB values: 255 255 255)
## 1  dark grey      (RGB values: 100 100 100)
## 2  red            (RGB values: 255   0   0)
## 3  green          (RGB values:   0 255   0)
## 4  blue           (RGB values:   0   0 255)
## 5  cyan           (RGB values:   0 255 255)
## 6  magenta        (RGB values: 255   0 255)
## 7  yellow         (RGB values: 255 255   0)
## 8  pale green     (RGB values: 152 251 152)
## 9  light sky blue (RGB values: 135 206 250)
## 10 orange         (RGB values: 255 165   0)
## 11 brown          (RGB values: 200 150 100)
## 12 pale pink      (RGB values: 255 200 200)
## 13 light grey     (RGB values: 170 170 170)
## 14 black          (RGB values:   0   0   0)
## 15 mid red:       (RGB values: 255  63  63)
## 16 light red      (RGB values: 255 127 127)
## 17 pink           (RGB values: 255 191 191)
artemisColors <- function() {
  names <- c("white", "dark grey", "red", "green",
             "blue", "cyan", "magenta", "yellow", "pale green",
             "light sky blue", "orange", "brown", "pale pink",
             "light grey", "black", "mid red", "light red", "pink")
  numbers <- 0:(length(names) - 1)
  r <- c(255, 100, 255, 0, 0, 0, 255, 255, 152, 135, 255, 200, 255,
         170, 0, 255, 255, 255)
  g <- c(255, 100, 0, 255, 0, 255, 0, 255, 251, 206, 165, 150, 200,
         170, 0, 63, 127, 191)
  b <- c(255, 100, 0, 0, 255, 255, 255, 0, 152, 250, 0, 100, 200,
         170, 0, 63, 127, 191)
  colors <- rgb(r, g, b, maxColorValue = 255)
  data.frame(n = numbers, names = names, colors = colors, r = r, g = g, b = b,
             stringsAsFactors = FALSE)
}

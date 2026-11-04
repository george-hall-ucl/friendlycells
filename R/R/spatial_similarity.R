#' Spatial similarity between spatial transcriptomics samples
#'
#' @param samples Named list of data frames, one per sample, with columns
#'   `array_row`, `array_col` (hexagonal grid array coordinates) and `label`.
#' @param normalize Divide each sample's connection counts by their total.
#' @param pools Optional named list mapping pool names to sample names;
#'   connection counts are summed within each pool.
#' @param only_transitions If `TRUE`, only count connections between spots with
#'   distinct labels (i.e. borders); does not count connections between spots
#'   with the same label.
#' @return Symmetric similarity matrix with values in \[0, 1\], or `NA` for
#'   samples with no neighboring spots when `normalize = TRUE`.
#' @examples
#' samples <- list(
#'   a = data.frame(array_row = c(0, 0, 1), array_col = c(0, 2, 1),
#'                  label = c("x", "y", "x")),
#'   b = data.frame(array_row = c(0, 0, 0), array_col = c(0, 2, 4),
#'                  label = c("x", "x", "y")),
#'   c = data.frame(array_row = c(0, 0), array_col = c(0, 2),
#'                  label = c("y", "y"))
#' )
#' spatial_similarity(samples)
#' spatial_similarity(samples, pools = list(ab = c("a", "b"), c = "c"))
#' @importFrom stats dist
#' @export
spatial_similarity <- function(samples, normalize = TRUE, pools = NULL,
                               only_transitions = FALSE) {
  if (anyNA(unlist(lapply(samples, function(s) as.character(s$label))))) {
    # e.g. read.csv() reads a label literally called "NA" as missing
    warning("Some labels are 'NA'. These will be treated as a label called ",
            "'NA' and used in the comparison, rather than being treated as ",
            "the absence of a label and dropped.")
  }
  samples <- lapply(samples, function(s) {
    s$label <- ifelse(is.na(s$label), "NA", as.character(s$label))
    s
  })
  # Sorted so floating-point results don't depend on sample order
  labels <- sort(unique(unlist(lapply(samples, function(s) s$label))))
  counts <- do.call(rbind, lapply(samples, connection_counts, labels = labels,
                                  only_transitions = only_transitions))
  if (!is.null(pools)) {
    counts <- do.call(rbind, lapply(pools, function(p) {
      colSums(counts[p, , drop = FALSE])
    }))
  }
  if (normalize) {
    empty <- rowSums(counts) == 0
    if (any(empty)) {
      warning("Similarities are NA for samples with no neighboring spots",
              if (only_transitions) " with different labels (no transitions)",
              ": ", paste(rownames(counts)[empty], collapse = ", "))
    }
    counts <- counts / rowSums(counts)
  }
  d <- as.matrix(dist(counts))
  # dist() puts 0 on the diagonal even for NA rows
  diag(d)[is.na(rowSums(counts))] <- NA
  # Identical samples give a max distance of 0, so don't divide by it
  d_max <- max(0, d, na.rm = TRUE)
  if (d_max > 0) d <- d / d_max
  1 - d
}

# Hexagonal Visium neighbors in array coordinates.
hex_offsets <- list(c(0, -2), c(0, 2), c(-1, -1), c(-1, 1), c(1, -1), c(1, 1))

connection_counts <- function(sample, labels, only_transitions = FALSE) {
  label <- factor(sample$label, levels = labels)
  key <- paste(sample$array_row, sample$array_col)
  counts <- matrix(0, length(labels), length(labels))
  for (o in hex_offsets) {
    nbr <- match(paste(sample$array_row + o[1], sample$array_col + o[2]), key)
    hit <- !is.na(nbr)
    counts <- counts + unclass(table(label[hit], label[nbr[hit]]))
  }
  if (only_transitions) diag(counts) <- 0
  as.vector(t(counts))
}

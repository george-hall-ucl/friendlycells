#' Spatial similarity between spatial transcriptomics samples
#'
#' Computes spatial similarity between spatial transcriptomics samples by
#' counting connections between neighboring spots/cells.
#'
#' @param samples Named list of data frames, one per sample, with a row per
#'   spot/cell.
#' @param x,y Names of the columns holding the positions of spots/cells, in any
#'   unit (e.g. pixels or micrometres) consistent between samples.
#' @param labels Name of the column holding the labels (e.g. cell types).
#' @param sample_normalize Divide each sample's connection counts by their
#'   total to control for sample size.
#' @param label_normalize Control for the size of each label: add one
#'   connection of every label to itself, divide each count by the square root
#'   of the two labels' connections to themselves (so each label's connections
#'   to itself become 1), then divide each label's counts by their sum (its
#'   degree), so each label's counts become the proportions of its neighbors'
#'   labels, as in one step of a random walk. Applied after pooling. Use
#'   instead of `sample_normalize`: set `sample_normalize = FALSE`.
#' @param pools Optional named list mapping pool names to sample names;
#'   connection counts are summed within each pool.
#' @param only_transitions If `TRUE`, only count connections between spots with
#'   distinct labels (i.e. borders); does not count connections between spots
#'   with the same label.
#' @param neighborhood Definition of the neighbourhood of each spot/cell:
#'   `"radius"` selects all within `radius`; `"knn"` selects the `k` nearest
#'   (symmetric: both are included in each other's neighborhood if either is
#'   among the other's nearest); or `"delaunay"` selects those surrounding it by
#'   Delaunay triangulation, dropping links longer than `radius`. The default,
#'   `"auto"`, uses `"radius"` if the positions lie on a grid (at least 75% of
#'   nearest-neighbor distances within 10% of their median) and `"delaunay"`
#'   otherwise, with a message saying which was used.
#' @param radius Neighborhood radius for `"radius"`, or maximum link length for
#'   `"delaunay"`, in the units of the positions. By default, set from the
#'   spacing of the data: set to the median distance from each spot/cell to its
#'   nearest neighbor over all samples, multiplied by 1.2 for `"radius"` or 3
#'   for `"delaunay"`. Used to make neighborhoods wider or narrower.
#' @param k Number of nearest neighbors for `neighborhood = "knn"`. Ties are
#'   broken arbitrarily.
#' @details Under default settings, the neighbourhoods of each sample's
#'   spots/cells depend on technology:
#'   * Visium spots: six immediate neighbors.
#'   * Visium HD bins: four side neighbors. Set `radius` to 1.5 times the
#'     bin size to add the diagonals.
#'   * Segmented cells (Xenium, or Visium HD with cell segmentation): the cells
#'     surrounding them, ignoring links across gaps.
#'
#'   Positions can be the pixel or micrometre coordinates from Space Ranger or
#'   Xenium outputs, as long as all samples use the same units. As the spacing
#'   is taken over all samples, don't mix samples in different units or from
#'   different technologies.
#' @return Symmetric similarity matrix with values in \[0, 1\], or `NA` for
#'   samples with no neighboring spots when `sample_normalize = TRUE`.
#' @examples
#' samples <- list(
#'   a = data.frame(
#'     x = c(0, 100, 50), y = c(0, 0, 86.6),
#'     label = c("x", "y", "x")
#'   ),
#'   b = data.frame(
#'     x = c(0, 100, 200), y = c(0, 0, 0),
#'     label = c("x", "x", "y")
#'   ),
#'   c = data.frame(
#'     x = c(0, 100), y = c(0, 0),
#'     label = c("y", "y")
#'   )
#' )
#' spatial_similarity(samples, "x", "y", "label")
#' spatial_similarity(samples, "x", "y", "label",
#'   pools = list(ab = c("a", "b"), c = "c")
#' )
#' spatial_similarity(samples, "x", "y", "label", neighborhood = "delaunay")
#' @importFrom stats dist median
#' @export
spatial_similarity <- function(samples, x, y, labels,
                               sample_normalize = TRUE,
                               label_normalize = FALSE, pools = NULL,
                               only_transitions = FALSE,
                               neighborhood = c(
                                 "auto", "radius", "knn",
                                 "delaunay"
                               ),
                               radius = NULL, k = 6) {
  if (sample_normalize && label_normalize) {
    stop(
      "Use either sample_normalize or label_normalize, not both: ",
      "set sample_normalize = FALSE to use label_normalize."
    )
  }
  settings <- neighborhood_settings(
    samples, list(x = x, y = y, labels = labels), match.arg(neighborhood),
    radius
  )
  neighborhood <- settings$neighborhood
  radius <- settings$radius
  if (anyNA(unlist(lapply(samples, function(s) as.character(s[[labels]]))))) {
    # e.g. read.csv() reads a label literally called "NA" as missing
    warning(
      "Some labels are 'NA'. These will be treated as a label called ",
      "'NA' and used in the comparison, rather than being treated as ",
      "the absence of a label and dropped."
    )
  }
  # Sorted so floating-point results don't depend on sample order
  all_labels <- sort(unique(unlist(lapply(samples, function(s) {
    na_as_label(s[[labels]])
  }))))
  counts <- do.call(rbind, lapply(samples, connection_counts,
    x = x, y = y, labels = labels, all_labels = all_labels,
    only_transitions = only_transitions,
    neighborhood = neighborhood, radius = radius,
    k = k
  ))
  if (!is.null(pools)) {
    counts <- do.call(rbind, lapply(pools, function(p) {
      colSums(counts[p, , drop = FALSE])
    }))
  }
  if (label_normalize) {
    n <- length(all_labels)
    counts <- t(apply(counts, 1, function(v) {
      m <- matrix(v, n, n, byrow = TRUE) + diag(n)
      m <- m / sqrt(outer(diag(m), diag(m)))
      as.vector(t(m / rowSums(m)))
    }))
  }
  if (sample_normalize) {
    empty <- rowSums(counts) == 0
    if (any(empty)) {
      warning(
        "Similarities are NA for samples with no neighboring spots",
        if (only_transitions) " with different labels (no transitions)",
        ": ", paste(rownames(counts)[empty], collapse = ", ")
      )
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

# Checks each sample has the named columns (a list of x, y and maybe labels),
# then resolves "auto" and the default radius from all samples
neighborhood_settings <- function(samples, columns, neighborhood, radius) {
  for (arg in names(columns)) {
    column <- columns[[arg]]
    if (!is.character(column) || length(column) != 1) {
      stop("`", arg, "` must be the name of a column.")
    }
    missing <- !vapply(samples, function(s) column %in% names(s), logical(1))
    if (any(missing)) {
      where <- if (is.null(names(samples))) {
        which(missing)
      } else {
        names(samples)[missing]
      }
      stop(
        "Column \"", column, "\" (`", arg, "`) is missing from samples: ",
        paste(where, collapse = ", ")
      )
    }
  }
  if (neighborhood == "auto" || (is.null(radius) && neighborhood != "knn")) {
    nn <- nn_distances(samples, columns$x, columns$y)
  }
  if (neighborhood == "auto") {
    neighborhood <- if (on_grid(nn)) "radius" else "delaunay"
    message(
      "Using neighborhood = \"", neighborhood, "\", as the positions ",
      if (neighborhood == "radius") "lie" else "don't lie", " on a grid."
    )
  }
  if (is.null(radius) && neighborhood != "knn") {
    radius <- default_radius(nn, neighborhood)
  }
  list(neighborhood = neighborhood, radius = radius)
}

# Distance from each spot or cell to its nearest neighbor, over all samples
nn_distances <- function(samples, x, y) {
  unlist(lapply(samples, function(s) {
    if (nrow(s) < 2) {
      return(numeric(0))
    }
    dbscan::kNN(cbind(s[[x]], s[[y]]), k = 1)$dist
  }))
}

# On a grid, nearly all spots are the same distance from their nearest
# neighbor, while segmented cells vary widely
on_grid <- function(nn) {
  isTRUE(mean(abs(nn / median(nn) - 1) < 0.1) >= 0.75)
}

# The median nearest-neighbor distance, times 1.2 for "radius" (Visium's six
# neighbors, or HD bins' four, allowing for jitter) or 3 for "delaunay"
# (dropping only links across gaps)
default_radius <- function(nn, neighborhood) {
  median(nn) * if (neighborhood == "radius") 1.2 else 3
}

connection_counts <- function(sample, x, y, labels, all_labels,
                              only_transitions = FALSE,
                              neighborhood = "radius", radius = NULL, k = 6) {
  codes <- factor(na_as_label(sample[[labels]]), levels = all_labels)
  pairs <- neighbor_pairs(
    cbind(sample[[x]], sample[[y]]), neighborhood, radius, k
  )
  counts <- unclass(table(codes[pairs[, 1]], codes[pairs[, 2]]))
  if (only_transitions) diag(counts) <- 0
  as.vector(t(counts))
}

# Missing labels are counted as a label called "NA"
na_as_label <- function(labels) {
  ifelse(is.na(labels), "NA", as.character(labels))
}

# Pairs (i, j) of neighboring rows of xy, each pair listed in both directions
neighbor_pairs <- function(xy, neighborhood, radius = NULL, k = 6) {
  n <- nrow(xy)
  if (n < 2) {
    return(matrix(integer(0), 0, 2))
  }
  if (neighborhood == "radius") {
    nbr <- dbscan::frNN(xy, eps = radius)$id
    pairs <- cbind(rep(seq_len(n), lengths(nbr)), unlist(nbr))
  } else if (neighborhood == "knn") {
    nbr <- dbscan::kNN(xy, k = min(k, n - 1))$id
    pairs <- cbind(rep(seq_len(n), ncol(nbr)), as.vector(nbr))
  } else {
    # An explicit window lets deldir link collinear points too
    window <- c(range(xy[, 1]) + c(-1, 1), range(xy[, 2]) + c(-1, 1))
    tri <- deldir::deldir(xy[, 1], xy[, 2], rw = window, suppressMsge = TRUE)
    pairs <- as.matrix(tri$delsgs[, c("ind1", "ind2")])
    if (!is.null(radius)) {
      len <- sqrt(rowSums((xy[pairs[, 1], , drop = FALSE] -
        xy[pairs[, 2], , drop = FALSE])^2))
      pairs <- pairs[len <= radius, , drop = FALSE]
    }
  }
  # kNN links can be one-way, so add the reverse of every pair, then dedupe
  pairs <- rbind(pairs, pairs[, 2:1])
  pairs <- unname(pairs[!duplicated((pairs[, 1] - 1) * n + pairs[, 2]), ,
    drop = FALSE
  ])
  pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
}

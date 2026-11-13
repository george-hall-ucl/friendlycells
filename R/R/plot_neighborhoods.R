#' Plot example neighborhoods
#'
#' Plot example neighborhoods on your samples. Draws a patch of each sample,
#' highlighting a few spots/cells and their neighbors, to check that the
#' neighborhood settings give the neighborhoods intended.
#'
#' @inheritParams spatial_similarity
#' @param samples_to_plot Names of the samples to plot. By default, plot all
#'   samples.
#' @param n Number of neighborhoods to highlight in each sample. Colors repeat
#'   after 7.
#' @param seed Random seed, so the same patches and spots/cells are plotted
#'   each time. `NULL` for different ones each time.
#' @return Invisibly, a named list with an element per plotted sample: a list
#'   of `patch` (row numbers of the spots/cells drawn) and `neighborhoods` (a
#'   list of row numbers, the highlighted spot/cell first, then its
#'   neighbors).
#' @examples
#' grid <- expand.grid(x = 1:40, y = 1:40)
#' samples <- list(a = data.frame(grid, label = "x"))
#' plot_neighborhoods(samples, "x", "y")
#' @export
plot_neighborhoods <- function(samples, x, y, samples_to_plot = NULL,
                               neighborhood = c(
                                 "auto", "radius", "knn",
                                 "delaunay"
                               ),
                               radius = NULL, k = 6, n = 7, seed = 123) {
  if (is.null(samples_to_plot)) samples_to_plot <- names(samples)
  unknown <- setdiff(samples_to_plot, names(samples))
  if (length(unknown) > 0) {
    stop("Unknown samples_to_plot: ", paste(unknown, collapse = ", "))
  }
  if (!is.null(seed)) {
    # Restore the user's random number stream afterwards
    old_seed <- get0(".Random.seed", globalenv(), inherits = FALSE)
    on.exit(
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old_seed, globalenv())
      },
      add = TRUE
    )
    set.seed(seed)
  }
  settings <- neighborhood_settings(
    samples, list(x = x, y = y), match.arg(neighborhood), radius
  )
  cols <- ceiling(sqrt(length(samples_to_plot)))
  old_par <- graphics::par(
    mfrow = c(ceiling(length(samples_to_plot) / cols), cols),
    mar = c(0.5, 0.5, 2, 0.5)
  )
  on.exit(graphics::par(old_par), add = TRUE)
  # Okabe-Ito, which suits color blindness, without black and grey
  colors <- rep_len(grDevices::palette.colors(palette = "Okabe-Ito")[2:8], n)
  # Halfway to white
  light <- grDevices::rgb((t(grDevices::col2rgb(colors)) + 255) / 2,
    maxColorValue = 255
  )
  result <- lapply(stats::setNames(nm = samples_to_plot), function(name) {
    xy <- cbind(samples[[name]][[x]], samples[[name]][[y]])
    hoods <- example_neighborhoods(
      xy, settings$neighborhood, settings$radius, k, n
    )
    fill <- rep("grey85", nrow(xy))
    for (i in seq_along(hoods$neighborhoods)) {
      hood <- hoods$neighborhoods[[i]]
      fill[hood] <- c(colors[i], rep(light[i], length(hood) - 1))
    }
    patch <- xy[hoods$patch, , drop = FALSE]
    graphics::plot(patch,
      type = "n", asp = 1, axes = FALSE, xlab = "",
      ylab = "", main = name
    )
    # Circles in data units, so spots/cells keep their spacing at any size
    spacing <- stats::median(nn_distances(list(data.frame(
      x = patch[, 1], y = patch[, 2]
    )), "x", "y"))
    if (is.na(spacing)) spacing <- 1
    graphics::symbols(patch,
      circles = rep(0.45 * spacing, nrow(patch)),
      inches = FALSE, add = TRUE, bg = fill[hoods$patch], fg = NA
    )
    hoods
  })
  invisible(result)
}

# A patch of the patch_size spots/cells nearest a random one, and up to n
# separate neighborhoods of spots/cells in its inner half, where they
# match the neighborhoods in the full sample
example_neighborhoods <- function(xy, neighborhood, radius, k, n,
                                  patch_size = 1000) {
  start <- xy[sample.int(nrow(xy), 1), ]
  patch <- order((xy[, 1] - start[1])^2 + (xy[, 2] - start[2])^2)
  patch <- patch[seq_len(min(patch_size, nrow(xy)))]
  pairs <- neighbor_pairs(xy[patch, , drop = FALSE], neighborhood, radius, k)
  neighbors <- split(pairs[, 2], factor(pairs[, 1], levels = seq_along(patch)))
  inner <- seq_len(ceiling(length(patch) / 2))
  used <- logical(length(patch))
  hoods <- list()
  for (i in inner[sample.int(length(inner))]) {
    hood <- c(i, neighbors[[i]])
    if (length(hoods) == n) break
    if (any(used[hood])) next
    # Their neighbors too, so neighborhoods don't touch
    used[c(hood, unlist(neighbors[hood]))] <- TRUE
    hoods <- c(hoods, list(patch[hood]))
  }
  list(patch = patch, neighborhoods = hoods)
}

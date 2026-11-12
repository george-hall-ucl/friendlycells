# Visium spots, from array coordinates to micrometres
visium <- function(array_row, array_col, label) {
  data.frame(x = 50 * array_col, y = 50 * sqrt(3) * array_row, label = label)
}
samples <- list(
  A = visium(c(0, 0, 1), c(0, 2, 1), c("a", "b", "a")),
  B = visium(c(0, 0, 0, 2), c(0, 2, 4, 0), c("a", "a", "b", "b")),
  C = visium(c(0, 0), c(0, 2), c("b", "b"))
)

expect_sim <- function(sim, ab, ac, bc) {
  n <- rownames(sim)
  expected <- matrix(c(1, ab, ac, ab, 1, bc, ac, bc, 1), 3,
    dimnames = list(n, n)
  )
  expect_equal(sim, expected)
}

test_that("hexagonal neighbor counts give similarity 1 - d / max(d)", {
  sim <- spatial_similarity(samples, "x", "y", "label", normalize = FALSE)
  expect_sim(sim, 1 - sqrt(2) / 4, 0, 1 - sqrt(10) / 4)
})

test_that("output row/column names respect input sample order", {
  sim <- spatial_similarity(rev(samples), "x", "y", "label", normalize = FALSE)
  expect_equal(rownames(sim), c("C", "B", "A"))
})

test_that("identical samples are seen as fully similar", {
  sim <- spatial_similarity(
    list(A = samples$A, B = samples$A), "x", "y", "label"
  )
  expect_equal(sim, matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))))
})

test_that("factor labels yield the same output as string labels", {
  sim <- spatial_similarity(samples, "x", "y", "label", normalize = FALSE)
  factors <- lapply(samples, function(s) {
    s$label <- factor(s$label)
    s
  })
  expect_equal(
    spatial_similarity(factors, "x", "y", "label", normalize = FALSE), sim
  )
})

test_that("labels read as missing (NA) become a label 'NA', with a warning", {
  sim <- spatial_similarity(samples, "x", "y", "label", normalize = FALSE)
  # read.csv() turns a label written as NA into a missing value
  csv <- read.csv(text = "sample,x,y,label
A,0,0,a
A,100,0,NA
A,50,86.60254037844386,a
B,0,0,a
B,100,0,a
B,200,0,NA
B,0,173.20508075688772,NA
C,0,0,NA
C,100,0,NA")
  from_file <- split(csv[-1], csv$sample)
  expect_warning(
    expect_equal(
      spatial_similarity(from_file, "x", "y", "label", normalize = FALSE),
      sim
    ),
    "treated as a label called 'NA'"
  )
})

test_that("normalize divides counts by their total", {
  dmax <- sqrt(11 / 8)
  sim <- spatial_similarity(samples, "x", "y", "label")
  expect_sim(sim, 1 - (sqrt(6) / 12) / dmax, 1 - sqrt(4 / 3) / dmax, 0)
})

test_that("samples where no spot has any neighbor yield NA and a warning", {
  sim <- spatial_similarity(samples, "x", "y", "label")
  # D's two spots don't touch, so it has no connections to normalize
  d <- visium(0, c(0, 4), "a")
  expect_warning(
    with_d <- spatial_similarity(c(samples, list(D = d)), "x", "y", "label",
      neighborhood = "radius"
    ),
    "no neighboring spots: D"
  )
  expect_equal(with_d[1:3, 1:3], sim)
  expect_true(all(is.na(with_d["D", ])) && all(is.na(with_d[, "D"])))
})

test_that("pools sum counts of their samples", {
  pools <- list(p1 = c("A", "B"), p2 = "C", p3 = "A")
  sim <- spatial_similarity(samples, "x", "y", "label",
    normalize = FALSE, pools = pools
  )
  expect_equal(rownames(sim), c("p1", "p2", "p3"))
  expect_sim(sim, 0, 1 - sqrt(6) / sqrt(38), 1 - 4 / sqrt(38))
})

test_that("reproduces the pig esophagus paper's similarities", {
  # Durkin et al. (2026)
  spots <- read.csv(test_path("pig_esophagus_spots.csv.xz"))
  spots <- cbind(visium(spots$array_row, spots$array_col, spots$label),
    sample_name = spots$sample_name
  )
  pools <- list(
    "Pre-implant" = "Pre-implant",
    "1 month" = "1 month (D)",
    "3 months" = c("3 months (E)", "3 months (F)"),
    "6 months" = c("6 months (B)", "6 months (C)"),
    "Native" = c("Native DP1", "Native DP2", "Native DP3", "Native (C)")
  )
  sim <- spatial_similarity(split(spots, spots$sample_name), "x", "y", "label",
    pools = pools
  )
  paper <- matrix(
    c(
      1.00, 0.01, 0.03, 0.09, 0.00,
      0.01, 1.00, 0.90, 0.68, 0.27,
      0.03, 0.90, 1.00, 0.70, 0.28,
      0.09, 0.68, 0.70, 1.00, 0.52,
      0.00, 0.27, 0.28, 0.52, 1.00
    ),
    5,
    dimnames = list(names(pools), names(pools))
  )
  # The paper's 7.5 px neighbor cutoff missed a few neighbors in Native DP2
  paper[1, 2] <- paper[2, 1] <- 0.02
  expect_equal(round(sim, 2), paper)
})

test_that("only_transitions drops same-label contacts, so only borders count", {
  sim <- spatial_similarity(samples, "x", "y", "label",
    normalize = FALSE, only_transitions = TRUE
  )
  expect_sim(sim, 0.5, 0, 0.5)
})

test_that("one-label samples have no transitions, giving NA and a warning", {
  expect_warning(
    sim <- spatial_similarity(samples, "x", "y", "label",
      only_transitions = TRUE
    ),
    "with different labels \\(no transitions\\): C"
  )
  expect_equal(sim["A", "B"], 1)
  expect_true(all(is.na(sim["C", ])) && all(is.na(sim[, "C"])))
})

test_that("only_transitions = FALSE gives the default result", {
  expect_equal(
    spatial_similarity(samples, "x", "y", "label", only_transitions = FALSE),
    spatial_similarity(samples, "x", "y", "label")
  )
})

test_that("pools sum the border counts of their samples", {
  pools <- list(p1 = c("A", "B"), p2 = "C", p3 = "A")
  sim <- spatial_similarity(samples, "x", "y", "label",
    normalize = FALSE, pools = pools,
    only_transitions = TRUE
  )
  expect_sim(sim, 0, 2 / 3, 1 / 3)
})

test_that("radius links cells up to and including the radius", {
  xy <- cbind(c(0, 10, 20.5), 0)
  expect_equal(
    neighbor_pairs(xy, "radius", radius = 10),
    rbind(c(1, 2), c(2, 1))
  )
})

test_that("knn links cells both ways if either is among the other's nearest", {
  # 3's nearest is 2, but 2's nearest is 1
  xy <- cbind(c(0, 1, 3), 0)
  expect_equal(
    neighbor_pairs(xy, "knn", k = 1),
    rbind(c(1, 2), c(2, 1), c(2, 3), c(3, 2))
  )
})

test_that("delaunay links surrounding cells, dropping links over radius", {
  # A triangle, and a far cell that Delaunay links to it
  xy <- cbind(c(0, 10, 5, 100), c(0, 0, 8, 0))
  expect_equal(
    neighbor_pairs(xy, "delaunay", radius = 20),
    rbind(c(1, 2), c(1, 3), c(2, 1), c(2, 3), c(3, 1), c(3, 2))
  )
})

test_that("delaunay links collinear cells to the next along the line", {
  xy <- cbind(c(0, 20, 10), 0)
  expect_equal(
    neighbor_pairs(xy, "delaunay"),
    rbind(c(1, 3), c(2, 3), c(3, 1), c(3, 2))
  )
})

test_that("columns are chosen by name, with an error naming missing ones", {
  renamed <- lapply(samples, function(s) {
    data.frame(col = s$x, row = s$y, cell_type = s$label, label = "z")
  })
  expect_equal(
    spatial_similarity(renamed, "col", "row", "cell_type"),
    spatial_similarity(samples, "x", "y", "label")
  )
  expect_error(
    spatial_similarity(samples, "x", "y", "cell_type"),
    '"cell_type" \\(`labels`\\) is missing from samples: A, B, C'
  )
})

test_that("default radius scales the median spacing over all samples", {
  # Nearest-neighbor distances 10, 10 and 22, 22, 28, so the spacing is 22
  s <- list(
    a = data.frame(x = c(0, 10), y = 0),
    b = data.frame(x = c(0, 22, 50), y = 0)
  )
  nn <- nn_distances(s, "x", "y")
  expect_equal(default_radius(nn, "radius"), 1.2 * 22)
  expect_equal(default_radius(nn, "delaunay"), 3 * 22)
})

test_that("auto uses radius on a grid and delaunay otherwise, saying which", {
  expect_message(
    spatial_similarity(samples, "x", "y", "label"), 'neighborhood = "radius"'
  )
  cells <- list(A = data.frame(
    x = c(0, 4, 11, 13, 30), y = c(0, 9, 2, 20, 6),
    label = "a"
  ))
  expect_message(
    spatial_similarity(cells, "x", "y", "label"), 'neighborhood = "delaunay"'
  )
})

test_that("with the default radius, results don't depend on units", {
  pixels <- lapply(samples, function(s) transform(s, x = x / 0.3, y = y / 0.3))
  expect_equal(
    spatial_similarity(pixels, "x", "y", "label"),
    spatial_similarity(samples, "x", "y", "label")
  )
})

test_that("example neighborhoods don't touch and match the full sample's", {
  set.seed(1)
  xy <- cbind(runif(3000, 0, 100), runif(3000, 0, 100))
  hoods <- example_neighborhoods(xy, "delaunay", radius = 6, k = 6, n = 10)
  expect_length(hoods$neighborhoods, 10)
  expect_false(anyDuplicated(unlist(hoods$neighborhoods)) > 0)
  pairs <- neighbor_pairs(xy, "delaunay", radius = 6)
  hood_of <- integer(nrow(xy))
  hood_of[unlist(hoods$neighborhoods)] <- rep(
    seq_along(hoods$neighborhoods), lengths(hoods$neighborhoods)
  )
  linked <- hood_of[pairs[, 1]] > 0 & hood_of[pairs[, 2]] > 0
  expect_equal(hood_of[pairs[linked, 1]], hood_of[pairs[linked, 2]])
  for (hood in hoods$neighborhoods) {
    expect_setequal(hood[-1], pairs[pairs[, 1] == hood[1], 2])
  }
})

test_that("unknown samples_to_plot give an error", {
  expect_error(
    plot_neighborhoods(samples, "x", "y", samples_to_plot = "D"), "Unknown"
  )
})

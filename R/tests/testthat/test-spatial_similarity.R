samples <- list(
  A = data.frame(array_row = c(0, 0, 1), array_col = c(0, 2, 1), label = c("a", "b", "a")),
  B = data.frame(array_row = c(0, 0, 0, 2), array_col = c(0, 2, 4, 0), label = c("a", "a", "b", "b")),
  C = data.frame(array_row = c(0, 0), array_col = c(0, 2), label = c("b", "b"))
)

expect_sim <- function(sim, ab, ac, bc) {
  n <- rownames(sim)
  expected <- matrix(c(1, ab, ac, ab, 1, bc, ac, bc, 1), 3, dimnames = list(n, n))
  expect_equal(sim, expected)
}

test_that("hexagonal neighbor counts give similarity 1 - d / max(d)", {
  sim <- spatial_similarity(samples, normalize = FALSE)
  expect_sim(sim, 1 - sqrt(2) / 4, 0, 1 - sqrt(10) / 4)
})

test_that("output row/column names respect input sample order", {
  sim <- spatial_similarity(rev(samples), normalize = FALSE)
  expect_equal(rownames(sim), c("C", "B", "A"))
})

test_that("identical samples are seen as fully similar", {
  sim <- spatial_similarity(list(A = samples$A, B = samples$A))
  expect_equal(sim, matrix(1, 2, 2, dimnames = list(c("A", "B"), c("A", "B"))))
})

test_that("factor labels yield the same output as string labels", {
  sim <- spatial_similarity(samples, normalize = FALSE)
  factors <- lapply(samples, function(s) {
    s$label <- factor(s$label)
    s
  })
  expect_equal(spatial_similarity(factors, normalize = FALSE), sim)
})

test_that("labels read as missing (NA) become a label 'NA', with a warning", {
  sim <- spatial_similarity(samples, normalize = FALSE)
  # read.csv() turns a label written as NA into a missing value
  csv <- read.csv(text = "sample,array_row,array_col,label
A,0,0,a
A,0,2,NA
A,1,1,a
B,0,0,a
B,0,2,a
B,0,4,NA
B,2,0,NA
C,0,0,NA
C,0,2,NA")
  from_file <- split(csv[-1], csv$sample)
  expect_warning(
    expect_equal(spatial_similarity(from_file, normalize = FALSE), sim),
    "treated as a label called 'NA'"
  )
})

test_that("normalize divides counts by their total", {
  dmax <- sqrt(11 / 8)
  sim <- spatial_similarity(samples)
  expect_sim(sim, 1 - (sqrt(6) / 12) / dmax, 1 - sqrt(4 / 3) / dmax, 0)
})

test_that("samples where no spot has any neighbor yield NA and a warning", {
  sim <- spatial_similarity(samples)
  # D's two spots don't touch, so it has no connections to normalize
  d <- data.frame(array_row = 0, array_col = c(0, 4), label = "a")
  expect_warning(with_d <- spatial_similarity(c(samples, list(D = d))),
                 "no neighboring spots: D")
  expect_equal(with_d[1:3, 1:3], sim)
  expect_true(all(is.na(with_d["D", ])) && all(is.na(with_d[, "D"])))
})

test_that("pools sum counts of their samples", {
  pools <- list(p1 = c("A", "B"), p2 = "C", p3 = "A")
  sim <- spatial_similarity(samples, normalize = FALSE, pools = pools)
  expect_equal(rownames(sim), c("p1", "p2", "p3"))
  expect_sim(sim, 0, 1 - sqrt(6) / sqrt(38), 1 - 4 / sqrt(38))
})

test_that("reproduces the pig esophagus paper's similarities", {
  # Durkin et al. (2026)
  spots <- read.csv(test_path("pig_esophagus_spots.csv.xz"))
  pools <- list(
    "Pre-implant" = "Pre-implant",
    "1 month" = "1 month (D)",
    "3 months" = c("3 months (E)", "3 months (F)"),
    "6 months" = c("6 months (B)", "6 months (C)"),
    "Native" = c("Native DP1", "Native DP2", "Native DP3", "Native (C)")
  )
  sim <- spatial_similarity(split(spots, spots$sample_name), pools = pools)
  paper <- matrix(c(1.00, 0.01, 0.03, 0.09, 0.00,
                    0.01, 1.00, 0.90, 0.68, 0.27,
                    0.03, 0.90, 1.00, 0.70, 0.28,
                    0.09, 0.68, 0.70, 1.00, 0.52,
                    0.00, 0.27, 0.28, 0.52, 1.00),
                  5, dimnames = list(names(pools), names(pools)))
  # The paper's 7.5 px neighbor cutoff missed a few neighbors in Native DP2
  paper[1, 2] <- paper[2, 1] <- 0.02
  expect_equal(round(sim, 2), paper)
})

test_that("only_transitions drops same-label contacts, so only borders count", {
  sim <- spatial_similarity(samples, normalize = FALSE, only_transitions = TRUE)
  expect_sim(sim, 0.5, 0, 0.5)
})

test_that("one-label samples have no transitions, giving NA and a warning", {
  expect_warning(sim <- spatial_similarity(samples, only_transitions = TRUE),
                 "with different labels \\(no transitions\\): C")
  expect_equal(sim["A", "B"], 1)
  expect_true(all(is.na(sim["C", ])) && all(is.na(sim[, "C"])))
})

test_that("only_transitions = FALSE gives the default result", {
  expect_equal(spatial_similarity(samples, only_transitions = FALSE),
               spatial_similarity(samples))
})

test_that("pools sum the border counts of their samples", {
  pools <- list(p1 = c("A", "B"), p2 = "C", p3 = "A")
  sim <- spatial_similarity(samples, normalize = FALSE, pools = pools,
                            only_transitions = TRUE)
  expect_sim(sim, 0, 2 / 3, 1 / 3)
})

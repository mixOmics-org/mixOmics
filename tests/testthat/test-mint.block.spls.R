context("mint.block.spls")

.mint_block_spls_data <- function(nvar = NULL) {
  data(breast.TCGA)
  X <- lapply(c("mrna", "mirna"), function(block) {
    x <- rbind(breast.TCGA$data.train[[block]], breast.TCGA$data.test[[block]])
    if (is.null(nvar)) x else x[, seq_len(nvar)]
  })
  names(X) <- c("mrna", "mirna")
  Y <- as.matrix(X$mrna[, 1])
  colnames(Y) <- "Y"
  X$mrna <- X$mrna[, -1]
  list(X = X, Y = Y, study = rep(c("study1", "study2"), c(150, 70)))
}

test_that("(mint.block.spls:basic): breast.TCGA", {
  data <- .mint_block_spls_data()
  res <- mint.block.spls(data$X, data$Y, study = data$study, ncomp = 2,
                         keepX = list(mrna = c(10, 10), mirna = c(20, 20)),
                         keepY = c(1, 1))

  expect_true("mint.block.spls" %in% class(res))
  expect_equal(res$loadings$mrna[5], 0, tolerance = 1e-5)
  expect_equal(res$mode, "regression")
})

test_that("(mint.block.spls:data): Y supplied through 'indY'", {
  data <- .mint_block_spls_data(nvar = 20)
  res <- mint.block.spls(c(data$X, list(Y = data$Y)), indY = 3,
                         study = data$study, ncomp = 1,
                         keepX = list(mrna = 2, mirna = 3))

  expect_true("mint.block.spls" %in% class(res))
  expect_equal(res$keepX$comp1$Y, 1)
  expect_null(res$keepY)
})

test_that("(mint.block.spls:parameter): partially specified 'keepX' and 'keepY'", {
  data <- .mint_block_spls_data(nvar = 20)
  res <- mint.block.spls(data$X, data$Y, study = data$study, ncomp = 2,
                         keepX = list(mrna = 2), keepY = 1)

  expect_equal(res$keepX$comp1$mrna, 2)
  expect_equal(res$keepX$comp2$mrna, ncol(data$X$mrna))
  expect_equal(res$keepX$comp1$Y, ncol(data$Y))
})

test_that("(mint.block.spls:error): invalid model inputs", {
  data <- .mint_block_spls_data(nvar = 20)

  expect_error(mint.block.spls(data$X, study = data$study),
               "Either 'Y' or 'indY' is needed", fixed = TRUE)
  expect_error(mint.block.spls(data$X, data$Y[, 1], study = data$study),
               "'Y' must be a numeric matrix.", fixed = TRUE)
  expect_error(mint.block.spls(data$X, data$Y, study = data$study,
                               keepX = list(bad = c(1, 1)), ncomp = 2),
               "Each entry of 'keepX' must have a unique name", fixed = TRUE)
})

test_that("(mint.block.spls:edge.case): single component", {
  data <- .mint_block_spls_data(nvar = 20)
  res <- mint.block.spls(data$X, data$Y, study = data$study, ncomp = 1,
                         keepX = list(mrna = 2, mirna = 3), keepY = 1)

  expect_true("mint.block.spls" %in% class(res))
  expect_equal(ncol(res$variates$mrna), 1L)
})

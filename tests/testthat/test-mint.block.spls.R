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
  expect_equal(names(res$keepX), c("mrna", "mirna"))
  expect_equal(as.numeric(res$keepX$mrna), 2)
  expect_equal(as.numeric(res$keepX$mirna), 3)
  expect_equal(as.numeric(res$keepY), 1)
})

test_that("(mint.block.spls:parameter): partially specified 'keepX' and 'keepY'", {
  data <- .mint_block_spls_data(nvar = 20)
  res <- mint.block.spls(data$X, data$Y, study = data$study, ncomp = 2,
                         keepX = list(mrna = 2), keepY = 1)

  # unspecified components and blocks default to keeping all variables
  expect_equal(as.numeric(res$keepX$mrna), c(2, ncol(data$X$mrna)))
  expect_equal(as.numeric(res$keepX$mirna), rep(ncol(data$X$mirna), 2))
  expect_equal(as.numeric(res$keepY), c(1, 1))
})

test_that("(mint.block.spls:output): keepX/keepY, design and weights returned", {
  data <- .mint_block_spls_data()
  res <- mint.block.spls(data$X, data$Y, study = data$study, ncomp = 2,
                         keepX = list(mrna = c(10, 10), mirna = c(20, 20)))

  .expect_numerically_close(res$variates$mrna[1, 1], -1.8696)

  # keepX returned as supplied, keepY defaults to keeping all of Y
  expect_equal(res$keepX, list(mrna = c(10, 10), mirna = c(20, 20)))
  expect_equal(as.numeric(res$keepY), c(1, 1))

  # the selection matches keepX
  expect_equal(unname(colSums(res$loadings$mrna != 0)), c(10, 10))
  expect_equal(unname(colSums(res$loadings$mirna != 0)), c(20, 20))

  # X contains all blocks including the outcome, located by indY
  expect_equal(names(res$X), c("mrna", "mirna", "Y"))
  expect_equal(res$indY, 3)

  # design and weights are returned
  expect_equal(dim(res$design), c(3, 3))
  expect_equal(rownames(res$weights), c("mrna", "mirna"))
})

test_that("(mint.block.spls:predict): predict works on the fitted model", {
  data <- .mint_block_spls_data()
  res <- mint.block.spls(data$X, data$Y, study = data$study, ncomp = 2,
                         keepX = list(mrna = c(10, 10), mirna = c(20, 20)))
  pred <- predict(res, newdata = data$X, study.test = data$study)

  expect_true(all(c("predict", "variates", "B.hat") %in% names(pred)))
  expect_equal(unname(dim(pred$predict$mrna)), c(220, 1, 2))
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

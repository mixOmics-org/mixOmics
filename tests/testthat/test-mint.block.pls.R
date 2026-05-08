context("mint.block.pls")

.mint_block_pls_data <- function(nvar = NULL) {
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

test_that("(mint.block.pls:basic): breast.TCGA", {
  data <- .mint_block_pls_data()
  res <- mint.block.pls(data$X, data$Y, study = data$study, ncomp = 2)

  expect_true("mint.block.pls" %in% class(res))
  expect_equal(res$loadings$mrna[5], -0.08196021, tolerance = 1e-5)
  expect_equal(res$mode, "regression")
})

test_that("(mint.block.pls:data): Y supplied through 'indY'", {
  data <- .mint_block_pls_data(nvar = 20)
  res <- mint.block.pls(c(data$X, list(Y = data$Y)), indY = 3,
                        study = data$study, ncomp = 1)

  expect_true("mint.block.pls" %in% class(res))
  expect_equal(res$ncomp, c(mrna = 1, mirna = 1, Y = 1))
})

test_that("(mint.block.pls:parameter): block-only 'design'", {
  data <- .mint_block_pls_data(nvar = 20)
  design <- matrix(c(0, 0.5, 0.5, 0), nrow = 2)

  expect_message(
    res <- mint.block.pls(data$X, data$Y, study = data$study,
                          design = design, ncomp = 1),
    "Design matrix has changed to include Y"
  )
  expect_true("mint.block.pls" %in% class(res))
  expect_equal(res$ncomp["Y"], c(Y = 1))
})

test_that("(mint.block.pls:error): invalid model inputs", {
  data <- .mint_block_pls_data(nvar = 20)

  expect_error(mint.block.pls(data$X, study = data$study),
               "Either 'Y' or 'indY' is needed", fixed = TRUE)
  expect_error(mint.block.pls(data$X, data$Y[, 1], study = data$study),
               "'Y' must be a numeric matrix.", fixed = TRUE)
})

test_that("(mint.block.pls:edge.case): missing 'study'", {
  data <- .mint_block_pls_data(nvar = 20)
  res <- mint.block.pls(data$X, data$Y, ncomp = 1)

  expect_true("mint.block.pls" %in% class(res))
  expect_equal(table(res$study)[["1"]], nrow(data$X$mrna))
})

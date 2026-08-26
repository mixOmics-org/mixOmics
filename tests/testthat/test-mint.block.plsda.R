context("mint.block.plsda")

.mint_block_plsda_data <- function(nvar = NULL) {
  data(breast.TCGA)
  X <- lapply(c("mrna", "mirna"), function(block) {
    x <- rbind(breast.TCGA$data.train[[block]], breast.TCGA$data.test[[block]])
    if (is.null(nvar)) x else x[, seq_len(nvar)]
  })
  names(X) <- c("mrna", "mirna")
  Y <- factor(c(as.character(breast.TCGA$data.train$subtype),
                as.character(breast.TCGA$data.test$subtype)))
  list(X = X, Y = Y, study = rep(c("study1", "study2"), c(150, 70)))
}

test_that("(mint.block.plsda:basic): breast.TCGA", {
  data <- .mint_block_plsda_data()
  res <- mint.block.plsda(data$X, data$Y, study = data$study, design = "full")

  expect_true("mint.block.plsda" %in% class(res))
  .expect_numerically_close(res$variates$mrna[1, 1], -7.503)
  .expect_numerically_close(res$variates$Y[220, 2], -0.796)
  .expect_numerically_close(res$loadings$mrna[1, 1], 0.100)
  .expect_numerically_close(res$loadings$mirna[100, 2], -0.019)
})

test_that("(mint.block.plsda:data): Y supplied through 'indY'", {
  data <- .mint_block_plsda_data(nvar = 20)
  res <- mint.block.plsda(c(data$X, list(Y = data$Y)), indY = 3,
                          study = data$study, ncomp = 1)

  expect_true("mint.block.plsda" %in% class(res))
  expect_equal(res$Y, data$Y)
})

test_that("(mint.block.plsda:data): near-zero variance predictors", {
  data <- .mint_block_plsda_data(nvar = 10)
  data$X$mrna <- cbind(zero = 1, data$X$mrna)

  res <- mint.block.plsda(data$X, data$Y, study = data$study,
                          ncomp = 1, near.zero.var = TRUE)
  expect_false("zero" %in% colnames(res$X$mrna))
})

test_that("(mint.block.plsda:parameter): block-only 'design'", {
  data <- .mint_block_plsda_data(nvar = 20)
  design <- matrix(c(0, 0.5, 0.5, 0), nrow = 2)

  expect_message(
    res <- mint.block.plsda(data$X, data$Y, study = data$study,
                            design = design, ncomp = 1),
    "Design matrix has changed to include Y"
  )
  expect_equal(res$indY, 3L)
})

test_that("(mint.block.plsda:parameter): both 'Y' and 'indY' supplied", {
  data <- .mint_block_plsda_data(nvar = 20)

  expect_warning(
    res <- mint.block.plsda(data$X, data$Y, indY = 1,
                            study = data$study, ncomp = 1),
    "'Y' and 'indY' are provided, 'Y' is used.", fixed = TRUE
  )
  expect_equal(res$Y, data$Y)
})

test_that("(mint.block.plsda:output): DA structure and weights, predict works", {
  data <- .mint_block_plsda_data()
  res <- mint.block.plsda(data$X, data$Y, study = data$study, ncomp = 2)

  expect_true(all(c("mint.block.plsda", "mint.block.pls", "block.plsda", "block.pls")
                  %in% class(res)))

  # X excludes the outcome block; Y is the original factor
  expect_equal(names(res$X), c("mrna", "mirna"))
  expect_true(is.factor(res$Y))
  expect_equal(dim(res$ind.mat), c(220, 3))
  expect_equal(res$indY, 3)
  expect_equal(rownames(res$weights), c("mrna", "mirna"))

  pred <- predict(res, newdata = data$X, study.test = data$study)
  expect_true(all(c("predict", "class", "MajorityVote", "WeightedVote") %in% names(pred)))
})

test_that("(mint.block.plsda:error): invalid inputs", {
  data <- .mint_block_plsda_data(nvar = 20)
  Y.with.na <- data$Y
  Y.with.na[1] <- NA

  expect_error(mint.block.plsda(data$X, Y = unmap(data$Y),
                                study = data$study, ncomp = 1),
               "'Y' should be a factor or a class vector.", fixed = TRUE)
  expect_error(mint.block.plsda(data$X, Y = Y.with.na,
                                study = data$study, ncomp = 1),
               "Unmapped Y contains samples with no associated class", fixed = TRUE)
  expect_error(mint.block.plsda(data$X, study = data$study),
               "Either 'Y' or 'indY' is needed", fixed = TRUE)
})

test_that("(mint.block.plsda:edge.case): missing 'study'", {
  data <- .mint_block_plsda_data(nvar = 20)
  res <- mint.block.plsda(data$X, data$Y, ncomp = 1)

  expect_true("mint.block.plsda" %in% class(res))
  expect_equal(table(res$study)[["1"]], nrow(data$X$mrna))
})

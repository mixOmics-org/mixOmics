context("mint.block.splsda")

.mint_block_splsda_data <- function(nvar = NULL) {
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

test_that("(mint.block.splsda:basic): breast.TCGA", {
  data <- .mint_block_splsda_data()
  res <- mint.block.splsda(data$X, data$Y, study = data$study, design = "full",
                           keepX = list(mrna = c(2, 2), mirna = c(3, 3)))

  expect_true("mint.block.splsda" %in% class(res))
  .expect_numerically_close(res$variates$mrna[1, 1], -1.48)
  .expect_numerically_close(res$variates$Y[220, 2], -0.5758)
})

test_that("(mint.block.splsda:data): Y supplied through 'indY'", {
  data <- .mint_block_splsda_data(nvar = 20)
  res <- mint.block.splsda(c(data$X, list(Y = data$Y)), indY = 3,
                           study = data$study, ncomp = 1,
                           keepX = list(mrna = 2, mirna = 3))

  expect_true("mint.block.splsda" %in% class(res))
  expect_equal(res$Y, data$Y)
})

test_that("(mint.block.splsda:parameter): partially specified 'keepX'", {
  data <- .mint_block_splsda_data(nvar = 20)
  res <- mint.block.splsda(data$X, data$Y, study = data$study,
                           ncomp = 2, keepX = list(mrna = 2))

  # unspecified components and blocks default to keeping all variables,
  # keepY defaults to the number of outcome levels
  expect_equal(as.numeric(res$keepX$mrna), c(2, ncol(data$X$mrna)))
  expect_equal(as.numeric(res$keepX$mirna), rep(ncol(data$X$mirna), 2))
  expect_equal(as.numeric(res$keepY), rep(nlevels(data$Y), 2))
})

test_that("(mint.block.splsda:output): keepX/keepY, DA structure and weights", {
  data <- .mint_block_splsda_data()
  res <- mint.block.splsda(data$X, data$Y, study = data$study, ncomp = 2,
                           keepX = list(mrna = c(10, 10), mirna = c(20, 20)))

  expect_true(all(c("mint.block.splsda", "mint.block.spls", "block.splsda", "block.spls")
                  %in% class(res)))

  # keepX returned as supplied, keepY defaults to the number of outcome levels
  expect_equal(res$keepX, list(mrna = c(10, 10), mirna = c(20, 20)))
  expect_equal(as.numeric(res$keepY), c(3, 3))
  expect_equal(unname(colSums(res$loadings$mrna != 0)), c(10, 10))
  expect_equal(unname(colSums(res$loadings$mirna != 0)), c(20, 20))

  # X excludes the outcome block; Y is the original factor
  expect_equal(names(res$X), c("mrna", "mirna"))
  expect_true(is.factor(res$Y))
  expect_equal(dim(res$ind.mat), c(220, 3))
  expect_equal(res$indY, 3)
  expect_equal(rownames(res$weights), c("mrna", "mirna"))

  # selectVar excludes the outcome block
  sv <- selectVar(res, comp = 1)
  expect_equal(names(sv), c("mrna", "mirna", "comp"))
})

test_that("(mint.block.splsda:predict): predict with weighted outputs", {
  data <- .mint_block_splsda_data()
  res <- mint.block.splsda(data$X, data$Y, study = data$study, ncomp = 2,
                           keepX = list(mrna = c(10, 10), mirna = c(20, 20)))
  pred <- predict(res, newdata = data$X, study.test = data$study)

  expect_true(all(c("predict", "class", "MajorityVote", "WeightedVote",
                    "AveragedPredict", "WeightedPredict") %in% names(pred)))
  acc <- mean(pred$class$max.dist[["mrna"]][, 2] == as.character(data$Y))
  expect_gt(acc, 0.7)
})

test_that("(mint.block.splsda:error): invalid inputs", {
  data <- .mint_block_splsda_data(nvar = 20)

  expect_error(mint.block.splsda(data$X, Y = unmap(data$Y),
                                 study = data$study, ncomp = 1),
               "'Y' should be a factor or a class vector.", fixed = TRUE)
  expect_error(mint.block.splsda(data$X, Y = rep("Basal", length(data$Y)),
                                 study = data$study, ncomp = 1),
               "'Y' should be a factor with more than one level", fixed = TRUE)
  expect_error(mint.block.splsda(data$X, study = data$study),
               "Either 'Y' or 'indY' is needed", fixed = TRUE)
  expect_error(mint.block.splsda(data$X, data$Y, study = data$study,
                                 keepX = list(bad = c(1, 1)), ncomp = 2),
               "Each entry of 'keepX' must have a unique name", fixed = TRUE)
})

test_that("(mint.block.splsda:edge.case): single component", {
  data <- .mint_block_splsda_data(nvar = 20)
  res <- mint.block.splsda(data$X, data$Y, study = data$study,
                           ncomp = 1, keepX = list(mrna = 2, mirna = 3))

  expect_true("mint.block.splsda" %in% class(res))
  expect_equal(ncol(res$variates$mrna), 1L)
})

context("mint.block.pls")

# Test
test_that("mint.block.pls works and returns the block-object structure", {
  data(breast.TCGA)
  study = c(rep("study1", 150), rep("study2", 70))
  mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  Y = mrna[, 1, drop = FALSE]
  mrna = mrna[, -1]
  data = list(mrna = mrna, mirna = mirna)

  res = mint.block.pls(data, Y, study = study, ncomp = 2)

  expect_true(all(c("mint.block.pls", "block.pls") %in% class(res)))
  expect_equal(res$mode, "regression")
  .expect_numerically_close(res$loadings$mrna[1, 1], -0.0598)
  .expect_numerically_close(res$variates$mrna[1, 1], -7.6191)

  # X contains all blocks including the outcome, located by indY
  expect_equal(names(res$X), c("mrna", "mirna", "Y"))
  expect_equal(res$indY, 3)

  # design and weights are returned
  expect_equal(dim(res$design), c(3, 3))
  expect_equal(rownames(res$design), c("mrna", "mirna", "Y"))
  expect_equal(rownames(res$weights), c("mrna", "mirna"))
  .expect_numerically_close(res$weights[1, 1], 0.5989)

  # predict works on the fitted model
  pred = predict(res, newdata = data, study.test = study)
  expect_true(all(c("predict", "variates", "B.hat") %in% names(pred)))
  expect_equal(unname(dim(pred$predict$mrna)), c(220, 1, 2))
  .expect_numerically_close(pred$predict$mrna[1, 1, 2], 2.7174)
})

test_that("mint.block.pls errors informatively with a vector Y", {
  data(breast.TCGA)
  study = c(rep("study1", 150), rep("study2", 70))
  mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  Y = mrna[, 1]
  mrna = mrna[, -1]
  data = list(mrna = mrna, mirna = mirna)

  expect_error(mint.block.pls(data, Y, study = study, ncomp = 2),
               "'Y' must be a numeric matrix")
})

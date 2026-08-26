context("mint.block.spls")

# Test
test_that("mint.block.spls works and returns correct keepX and keepY", {
  data(breast.TCGA)
  study = c(rep("study1", 150), rep("study2", 70))
  mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  Y = mrna[, 1, drop = FALSE]
  mrna = mrna[, -1]
  data = list(mrna = mrna, mirna = mirna)

  res = mint.block.spls(data, Y, study = study, ncomp = 2,
                        keepX = list(mrna = c(10, 10), mirna = c(20, 20)))

  expect_true(all(c("mint.block.spls", "block.spls") %in% class(res)))
  expect_equal(res$mode, "regression")
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

  # predict works on the fitted model
  pred = predict(res, newdata = data, study.test = study)
  expect_true(all(c("predict", "variates", "B.hat") %in% names(pred)))
  expect_equal(unname(dim(pred$predict$mrna)), c(220, 1, 2))
})

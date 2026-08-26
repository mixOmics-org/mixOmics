

###############################################################################
### ================================ BASIC ================================ ###
###############################################################################

test_that("mint.block.plsda and mint.block.splsda works", {
  
  data(breast.TCGA)
  mrna <- rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna <- rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  X <- list(mrna = mrna, mirna = mirna)
  Y <- c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)
  
  study <- c(rep("study1",150), rep("study2",70))
  
  res.mint.block.plsda <- mint.block.plsda(X, Y, study=study, design="full")
  
  expect_true("mint.block.plsda" %in% class(res.mint.block.plsda))
  .expect_numerically_close(res.mint.block.plsda$variates$mrna[1,1], -7.503)
  .expect_numerically_close(res.mint.block.plsda$variates$Y[220,2], -0.796)
  
  .expect_numerically_close(res.mint.block.plsda$loadings$mrna[1,1], 0.100)
  .expect_numerically_close(res.mint.block.plsda$loadings$mirna[100,2], -0.019)
  
  
  res.mint.block.splsda <- mint.block.splsda(X, Y, study=study, design="full", keepX = list(mrna=c(2,2), mirna=c(3,3)))
  
  expect_true("mint.block.splsda" %in% class(res.mint.block.splsda))
  .expect_numerically_close(res.mint.block.splsda$variates$mrna[1,1], -1.48)
  .expect_numerically_close(res.mint.block.splsda$variates$Y[220,2], -0.5758)
  
})


###############################################################################
### ================================ DATA ================================= ###
###############################################################################





###############################################################################
### ============================== PARAMETER ============================== ###
###############################################################################

test_that("mint.block.plsda returns the DA block-object structure and predicts", {

  data(breast.TCGA)
  mrna <- rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna <- rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  X <- list(mrna = mrna, mirna = mirna)
  Y <- c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)

  study <- c(rep("study1",150), rep("study2",70))

  res <- mint.block.plsda(X, Y, study = study, ncomp = 2)

  expect_true(all(c("mint.block.plsda", "mint.block.pls", "block.plsda", "block.pls")
                  %in% class(res)))

  # X excludes the outcome block; Y is the original factor
  expect_equal(names(res$X), c("mrna", "mirna"))
  expect_true(is.factor(res$Y))
  expect_equal(dim(res$ind.mat), c(220, 3))
  expect_equal(res$indY, 3)
  expect_equal(rownames(res$weights), c("mrna", "mirna"))

  pred <- predict(res, newdata = X, study.test = study)
  expect_true(all(c("predict", "class", "MajorityVote", "WeightedVote") %in% names(pred)))
})

test_that("mint.block.splsda returns correct keepX/keepY and predicts", {

  data(breast.TCGA)
  mrna <- rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
  mirna <- rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
  X <- list(mrna = mrna, mirna = mirna)
  Y <- c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)

  study <- c(rep("study1",150), rep("study2",70))

  res <- mint.block.splsda(X, Y, study = study, ncomp = 2,
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

  # predict works, including the weighted outputs
  pred <- predict(res, newdata = X, study.test = study)
  expect_true(all(c("predict", "class", "MajorityVote", "WeightedVote",
                    "AveragedPredict", "WeightedPredict") %in% names(pred)))
  acc <- mean(pred$class$max.dist[["mrna"]][, 2] == as.character(Y))
  expect_gt(acc, 0.7)
})





###############################################################################
### ================================ ERROR ================================ ###
###############################################################################





###############################################################################
### ============================== EDGE CASES ============================= ###
###############################################################################





###############################################################################

# if the method is a graphical function, include the following:
#dev.off()
#unlink(list.files(pattern = "*.pdf"))
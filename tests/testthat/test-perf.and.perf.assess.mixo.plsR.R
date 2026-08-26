context("perf.assess.mixo.pls")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.mixo_pls() and perf.assess.mixo_pls()

test_that("perf.assess.mixo_pls works same as perf for same components in serial and in parallel for regression mode", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  pls.obg <- pls(Y, X, ncomp = 2)
  
  # Execution using old perf() function in serial
  out.1 <- perf(pls.obg, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SerialParam(RNGseed = 1000), seed = 12)

  # Execution using old perf() function in parallel
  out.2 <- perf(pls.obg, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SnowParam(workers = 2), seed = 12)
  
  # New perf.assess in series
  out.3 <- perf.assess(pls.obg, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SerialParam(RNGseed = 1000), seed = 12)
  
  # New perf.assess in parallel
  out.4 <- perf.assess(pls.obg, validation = "Mfold", folds = 2, nrepeat = 1, BPPARAM = SnowParam(workers = 2), seed = 12)
  
  # check output structure
  expect_equal(names(out.1), c("call", "measures", "features"))
  expect_equal(names(out.2), c("call", "measures", "features"))
  expect_equal(names(out.3), c("call", "measures", "features"))
  expect_equal(names(out.4), c("call", "measures", "features"))
  
  # Expect the same errors
  expect_equal(c(0.041, -0.177), round(out.1$measures$Q2.total$summary[, "mean"], 3))
  expect_equal(out.1$measures$RSS.upred$summary[1,3], out.2$measures$RSS.upred$summary[1,3])
  expect_equal(out.1$measures$RSS.upred$summary[2,3], out.3$measures$RSS.upred$summary[1,3])
  expect_equal(out.1$measures$RSS.upred$summary[2,3], out.4$measures$RSS.upred$summary[1,3])
  
  # should have no feature stability values as pls object
  expect_equal(out.1$features, NULL)
  expect_equal(out.2$features, NULL)
  expect_equal(out.3$features, NULL)
  expect_equal(out.4$features, NULL)

})
  
test_that("perf() and perf.assess() works on pls with nzv features (all modes)", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:1000]
  Y <- liver.toxicity$clinic
  
  # to reproduce error, we need to induce some features to have near zero variance
  X[, c(1, 23, 62, 234, 789)] <- 0
  
  modes <- c("classic", "regression", "canonical")
  trueVals <- list(c(0.021,  -0.014),
                   c(-0.071, -0.260),
                   c(0.089, -0.415))
  
  for (m in 1:3) {
    suppressWarnings(pls.obg <- pls(Y, X, ncomp = 2, mode = modes[m]))
    suppressWarnings(pls.perf.obj <- perf(pls.obg, validation = "Mfold", folds = 2, 
                                          progressBar = F, 
                                          nrepeat = 1, BPPARAM = SerialParam(),
                                          seed = 2))
    
    testVals <- round(pls.perf.obj$measures$Q2.total$summary[, "mean"], 3)
    expect_equal(trueVals[[m]], testVals)
    
    suppressWarnings(pls.perf.assess.obj <- perf.assess(pls.obg, validation = "Mfold", folds = 2, 
                                          progressBar = F, 
                                          nrepeat = 1, BPPARAM = SerialParam(),
                                          seed = 2))
    testVals <- round(pls.perf.assess.obj$measures$Q2.total$summary[, "mean"], 3)
    expect_equal(trueVals[[m]][2], testVals)
  }
  
})

## ------------------------------------------------------------------------ ##
## Test perf.mixo_spls() and perf.assess.mixo_spls()

test_that("perf() works on spls object in serial and parallel", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  model.spls = spls(X, Y, ncomp = 2, mode = 'regression',
                    keepX = c(10, 10), keepY = c(4, 4))
  
  # run perf in series
  out.1 <- perf(model.spls, validation = "Mfold", folds = 2, BPPARAM = SerialParam(), seed = 12)
  
  # run perf in parallel
  out.2 <- perf(model.spls, validation = "Mfold", folds = 2, BPPARAM = SnowParam(workers = 2), seed = 12)
  
  # run perf.assess in series
  out.3 <- perf.assess(model.spls, validation = "Mfold", folds = 2, BPPARAM = SerialParam(), seed = 12)
  
  # run perf.assess in parallel
  out.4 <- perf.assess(model.spls, validation = "Mfold", folds = 2, BPPARAM = SerialParam(), seed = 12)
  
  # check output structure
  expect_equal(names(out.1), c("call", "measures", "features"))
  expect_equal(names(out.2), c("call", "measures", "features"))
  expect_equal(names(out.3), c("call", "measures", "features"))
  expect_equal(names(out.4), c("call", "measures", "features"))
  
  # check error values are expected
  trueVals <- c(0.096, -0.190)
  expect_equal(trueVals, round(out.1$measures$Q2.total$summary[, "mean"], 3))
  expect_equal(trueVals, round(out.2$measures$Q2.total$summary[, "mean"], 3))
  expect_equal(trueVals[2], round(out.3$measures$Q2.total$summary[, "mean"], 3))
  expect_equal(trueVals[2], round(out.4$measures$Q2.total$summary[, "mean"], 3))
  expect_equal(out.1$measures$RSS.upred$summary, out.2$measures$RSS.upred$summary)
  
  # check feature stability
  expect_equal(names(out.1$features$stability.X$comp2[5]), "A_43_P10885")
  expect_equal(names(out.2$features$stability.X$comp2[5]), "A_43_P10885")
  expect_equal(names(out.3$features$stability.X$comp2[5]), "A_43_P10885")
  expect_equal(names(out.4$features$stability.X$comp2[5]), "A_43_P10885")
  
})

## ------------------------------------------------------------------------ ##
## Check get error when try to pass pls or spls object made with invariant mode to perf or perf.assess

test_that("perf gives error for invariant mode", code = {
  
  # set up data
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:500]
  Y <- liver.toxicity$clinic
  pls.obg <- pls(Y, X, ncomp = 2, mode = "invariant")
  spls.obj <- spls(X, Y, ncomp = 2, mode = 'invariant', keepX = c(10, 10), keepY = c(4, 4))
  
  # try running perf and perf.assess on pls
  expect_error(perf(pls.obg), "BiocParallel errors.*mode 'invariant'", fixed = FALSE)
  expect_error(perf.assess(pls.obg), "BiocParallel errors.*mode 'invariant'", fixed = FALSE)
  
  # try running perf and perf.assess on spls
  expect_error(perf(spls.obj), "BiocParallel errors.*mode 'invariant'", fixed = FALSE)
  expect_error(perf.assess(spls.obj), "BiocParallel errors.*mode 'invariant'", fixed = FALSE)

  })

## ------------------------------------------------------------------------ ##
## Regression guard: ggplot2 size -> linewidth

test_that("(plot.perf.pls:edge.case): no ggplot2 size-aesthetic deprecation warning (Q2 criterion)", {
  data("liver.toxicity")
  X <- liver.toxicity$gene[1:100]
  Y <- liver.toxicity$clinic
  pls.mod  <- pls(X, Y, ncomp = 2)
  perf.obj <- perf(pls.mod, validation = "Mfold", folds = 2, nrepeat = 1,
                   progressBar = FALSE)
  expect_no_warning(plot(perf.obj, criterion = "Q2.total"))
})

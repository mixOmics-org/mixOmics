context("perf.assess.diablo")
library(BiocParallel)

## ------------------------------------------------------------------------ ##
## Test perf.sgccda()

test_that("perf.diablo works with with auroc", {
  # set up data and model
    data(nutrimouse)
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    nutrimouse.sgccda <- block.splsda(X = data,
                                      Y = nutrimouse$diet,
                                      design = design,
                                      keepX = list(gene=c(10,10), lipid=c(15,15)),
                                      ncomp = 2)
    
    # run perf model - set seed as 100, ignores RNGseed
    perf.res <- perf(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                      BPPARAM = SerialParam(RNGseed = 100000), seed = 100, progressBar = FALSE)
    true_aucs <- c(0.95, 0.62, 0.68, 0.54, 0.77)
    aucs <- round(unname(perf.res$auc$comp1[,1]), 2)
    expect_equal(aucs, true_aucs)
    
    # run in parallel
    perf.res.parallel = perf(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                      BPPARAM = SnowParam(workers = 2), seed = 100, progressBar = FALSE)
    aucs <- round(unname(perf.res.parallel$auc$comp1[,1]), 2)
    expect_equal(aucs, true_aucs)
    expect_equal(perf.res$weights, perf.res.parallel$weights)
    
    # run perf.assess
    perf.assess.res <- perf.assess(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                            BPPARAM = SerialParam(RNGseed = 100000), seed = 100, progressBar = FALSE)
    expect_equal(perf.res$weights$comp2, perf.assess.res$weights$comp2)
    
    # run perf.assess in parallel
    perf.assess.res.parallel <- perf.assess(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
                                   BPPARAM = SnowParam(workers = 2), seed = 100, progressBar = FALSE)
    expect_equal(perf.res$weights$comp2, perf.assess.res.parallel$weights$comp2)
    
})

## ------------------------------------------------------------------------ ##
## Test  perf.assess.sgccda() give informative error message when one sample in one class

test_that("perf.assess.sgccda error when one sample in one class", code = {
  
  # set up data and model
  data(nutrimouse)
  data = list(gene = nutrimouse$gene[1:10, ], lipid = nutrimouse$lipid[1:10, ])
  design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
  nutrimouse.sgccda <- block.splsda(X = data,
                                    Y = nutrimouse$diet[1:10],
                                    design = design,
                                    keepX = list(gene=c(10,10), lipid=c(15,15)),
                                    ncomp = 2)
  
  test_run <- function() {
    perf.assess(nutrimouse.sgccda, folds = 2, nrepeat = 1, auc = TRUE, 
         BPPARAM = SerialParam(RNGseed = 100000), seed = 100, progressBar = FALSE)}
  
  expect_error(test_run(), 
               "Cannot evaluate performance when a class level ('ref') has only a single associated sample.",
               fixed = TRUE)
})

## ------------------------------------------------------------------------ ##
## Test perf.sgccda() and perf.assess.sgccda() with folds supplied as a list

# model shared by the tests below: 40 samples, 5 classes of 8 samples
.diablo.folds.setup <- function() {
  data(nutrimouse, envir = environment())
  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
  design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
  model <- block.splsda(X = data, Y = nutrimouse$diet, design = design,
                        keepX = list(gene = c(10,10), lipid = c(15,15)),
                        ncomp = 2)
  n <- nrow(data$gene)
  # pretend the samples come from 5 groups (e.g. sites or batches)
  group <- rep(1:5, length.out = n)
  list(model = model, n = n, folds = split(seq_len(n), group))
}

test_that("(perf.sgccda:parameter): folds supplied as a list", {
  setup <- .diablo.folds.setup()
  model <- setup$model; n <- setup$n; my.folds <- setup$folds
  
  res <- perf(model, validation = "Mfold", folds = my.folds, nrepeat = 1,
              BPPARAM = SerialParam(), progressBar = FALSE)
  
  # one fold per element of the list, and a single repeat
  expect_equal(sort(unique(res$weights$fold)), seq_along(my.folds))
  expect_equal(unique(res$weights$rep), 1)
  
  # every sample is predicted exactly once, in the order of the supplied folds
  pred <- res$predict$nrep1$gene$comp1
  expect_equal(nrow(pred), n)
  expect_false(any(duplicated(rownames(pred))))
  expect_setequal(rownames(pred), rownames(model$X$gene))
  expect_equal(rownames(pred), rownames(model$X$gene)[unlist(my.folds)])
  
  # the split is not random any more: the same folds give the same result
  res2 <- perf(model, validation = "Mfold", folds = my.folds, nrepeat = 1,
               BPPARAM = SerialParam(), progressBar = FALSE)
  expect_equal(res$error.rate, res2$error.rate)
  expect_equal(res$MajorityVote.error.rate, res2$MajorityVote.error.rate)
})

test_that("(perf.assess.sgccda:parameter): folds supplied as a list", {
  setup <- .diablo.folds.setup()
  model <- setup$model; n <- setup$n; my.folds <- setup$folds
  
  res <- perf.assess(model, validation = "Mfold", folds = my.folds, nrepeat = 1,
                     BPPARAM = SerialParam(), progressBar = FALSE)
  
  expect_equal(sort(unique(res$weights$fold)), seq_along(my.folds))
  expect_equal(unique(res$weights$rep), 1)
  
  pred <- res$predict$nrep1$gene[[1]]
  expect_equal(nrow(pred), n)
  expect_false(any(duplicated(rownames(pred))))
  expect_setequal(rownames(pred), rownames(model$X$gene))
})

test_that("(perf.sgccda:edge.case): a list of folds with nrepeat > 1 warns and runs once", {
  setup <- .diablo.folds.setup()
  
  expect_warning(
    res <- perf(setup$model, validation = "Mfold", folds = setup$folds, nrepeat = 3,
                BPPARAM = SerialParam(), progressBar = FALSE),
    "'nrepeat' is set to '1'", fixed = TRUE)
  expect_equal(unique(res$weights$rep), 1)
  expect_null(res$error.rate.sd)
  
  expect_warning(
    res.assess <- perf.assess(setup$model, validation = "Mfold", folds = setup$folds,
                              nrepeat = 3, BPPARAM = SerialParam(), progressBar = FALSE),
    "'nrepeat' is set to '1'", fixed = TRUE)
  expect_equal(unique(res.assess$weights$rep), 1)
})

test_that("(perf.sgccda:error): invalid lists of folds give a clear error", {
  setup <- .diablo.folds.setup()
  model <- setup$model; n <- setup$n; my.folds <- setup$folds
  run <- function(folds) perf(model, validation = "Mfold", folds = folds, nrepeat = 1,
                              BPPARAM = SerialParam(), progressBar = FALSE)
  
  # a sample appears in two folds
  dup <- my.folds; dup[[1]][1] <- dup[[2]][1]
  expect_error(run(dup), "Repeated samples in folds", fixed = TRUE)
  # the error is raised before the repeats are dispatched to BiocParallel,
  # so it is not wrapped in a 'BiocParallel errors' message
  expect_error(run(dup), "^Invalid folds")
  
  # a sample is missing
  short <- my.folds; short[[1]] <- short[[1]][-1]
  expect_error(run(short), paste0("total number of samples in folds must be equal to ", n),
               fixed = TRUE)
  
  # an index larger than n
  large <- my.folds; large[[1]][1] <- n + 1L
  expect_error(run(large), paste0("integer indices between 1 and ", n), fixed = TRUE)
  
  # a non-integer index
  frac <- my.folds; frac[[1]] <- frac[[1]] + 0.5
  expect_error(run(frac), paste0("integer indices between 1 and ", n), fixed = TRUE)
  
  # a single fold
  expect_error(run(list(seq_len(n))), "Invalid number of folds", fixed = TRUE)
  
  # an empty fold
  expect_error(run(c(my.folds, list(integer(0)))), "at least one sample", fixed = TRUE)

  # a list of lists
  expect_error(run(list(my.folds, my.folds)), "one vector per fold", fixed = TRUE)
  
  # a class absent from the training set of a fold (one fold per class)
  by.class <- split(seq_len(n), model$Y)
  expect_error(run(by.class), "not represented in the training set of fold 1", fixed = TRUE)
  
  # perf.assess uses the same checks
  expect_error(
    perf.assess(model, validation = "Mfold", folds = dup, nrepeat = 1,
                BPPARAM = SerialParam(), progressBar = FALSE),
    "Repeated samples in folds", fixed = TRUE)
})

context("tune.rcc")
library(BiocParallel)

# Reference workflow: prepare each training fold explicitly before fitting.
.rcc_cv_reference <- function(X, Y, folds, lambda1, lambda2, scale) {
  scores <- lapply(folds, function(omit) {
    train.X <- base::scale(X[-omit, , drop = FALSE], scale = scale)
    train.Y <- base::scale(Y[-omit, , drop = FALSE], scale = scale)
    test.X <- base::scale(X[omit, , drop = FALSE],
                         center = attr(train.X, "scaled:center"),
                         scale = if (scale) attr(train.X, "scaled:scale") else FALSE)
    test.Y <- base::scale(Y[omit, , drop = FALSE],
                         center = attr(train.Y, "scaled:center"),
                         scale = if (scale) attr(train.Y, "scaled:scale") else FALSE)
    # Preprocessing is explicit here; disable automatic standardisation.
    fit <- rcc(train.X, train.Y, ncomp = 1,
               lambda1 = lambda1, lambda2 = lambda2, scale = FALSE)
    test.X[is.na(test.X)] <- 0
    test.Y[is.na(test.Y)] <- 0
    cbind(test.X %*% fit$loadings$X, test.Y %*% fit$loadings$Y)
  })
  scores <- do.call(rbind, scores)
  cor(scores[, 1], scores[, 2])
}

test_that("tune.rcc works with Mfold method", code = {
  # set up data
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  
  # run
  tune.rcc.res <- tune.rcc(X, Y, validation = "Mfold", seed = 20)
  
  # check outputs
  expect_equal(class(tune.rcc.res), "tune.rcc")
  best <- cbind(match(tune.rcc.res$opt.lambda1, tune.rcc.res$grid1),
                match(tune.rcc.res$opt.lambda2, tune.rcc.res$grid2))
  expect_equal(unname(tune.rcc.res$mat[best]),
               rep(max(tune.rcc.res$mat), nrow(best)))
  expect_equal(tune.rcc.res$grid1, c(0.00100, 0.25075, 0.50050, 0.75025, 1.00000))
  
  # check can plot
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(tune.rcc.res))
})

test_that("tune.rcc works with loo method", code = {
  # set up data
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  
  # run
  tune.rcc.res <- tune.rcc(X, Y, validation = "loo", seed = 20)
  
  # check outputs
  expect_equal(class(tune.rcc.res), "tune.rcc")
  best <- cbind(match(tune.rcc.res$opt.lambda1, tune.rcc.res$grid1),
                match(tune.rcc.res$opt.lambda2, tune.rcc.res$grid2))
  expect_equal(unname(tune.rcc.res$mat[best]),
               rep(max(tune.rcc.res$mat), nrow(best)))
  expect_equal(tune.rcc.res$grid1, c(0.00100, 0.25075, 0.50050, 0.75025, 1.00000))
  
  # check can plot
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(tune.rcc.res))
})
  
test_that("tune.rcc works in parallel same as in series", code = {
  # set up data
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  
  # run in series
  tune.rcc.res <- tune.rcc(X, Y, validation = "Mfold",
                           BPPARAM = SerialParam(RNGseed = NULL), seed = 12)
  # run in parallel
  tune.rcc.res.parallel <- tune.rcc(X, Y, validation = "Mfold",
                           BPPARAM = SnowParam(workers = 2, RNGseed = NULL), seed = 12)
  
  # check outputs
  expect_equal(class(tune.rcc.res), "tune.rcc")
  expect_equal(tune.rcc.res$opt.lambda1, tune.rcc.res.parallel$opt.lambda1)
  expect_equal(tune.rcc.res$opt.lambda2, tune.rcc.res.parallel$opt.lambda2)
  expect_equal(tune.rcc.res$opt.score, tune.rcc.res.parallel$opt.score)
  expect_equal(tune.rcc.res$mat, tune.rcc.res.parallel$mat)
  expect_equal(tune.rcc.res$grid2, tune.rcc.res.parallel$grid2)
  
  # check can plot
  pdf(NULL)
  on.exit(dev.off())
  expect_silent(plot(tune.rcc.res))
  expect_silent(plot(tune.rcc.res.parallel))
})

test_that("tune.rcc validates the scale option", {
  data(linnerud)
  for (bad in list(NULL, NA, logical(), c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(tune.rcc(linnerud$exercise, linnerud$physiological, scale = bad),
                 "'scale' must be")
  }
})

test_that("Mfold reports columns that cannot be scaled within a training fold", {
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  # A variable can be nonconstant globally but constant in one training fold.
  X[, 1] <- 1
  X[1:5, 1] <- 2:6
  expect_error(mixOmics:::Mfold(X, Y, 1, 1, list(1:5, 6:20), scale = TRUE),
               "zero variance in 'X'")
})

test_that("Mfold uses training statistics and preserves missing training values", {
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  X[1, ] <- X[1, ] * 10
  Y[2, ] <- Y[2, ] * 3
  folds <- split(seq_len(nrow(X)), rep(1:4, length.out = nrow(X)))
  expect_identical(formals(mixOmics:::Mfold)$scale, TRUE)
  for (scaling in c(FALSE, TRUE)) {
    expected <- .rcc_cv_reference(X, Y, folds, 1, 1, scale = scaling)
    actual <- mixOmics:::Mfold(X, Y, 1, 1, folds, scale = scaling)
    expect_equal(actual, expected, tolerance = 1e-8)
  }
  expect_equal(mixOmics:::Mfold(X, Y, 1, 1, folds), actual)
  global <- mixOmics:::Mfold(base::scale(X), base::scale(Y), 1, 1, folds,
                            scale = FALSE)
  expect_gt(abs(actual - global), 1e-4)

  # These omitted values become training values in subsequent folds.
  X[1, 2] <- NA
  Y[2, 3] <- NA
  for (scaling in c(FALSE, TRUE)) {
    expected <- .rcc_cv_reference(X, Y, folds, 10, 10, scale = scaling)
    actual <- mixOmics:::Mfold(X, Y, 10, 10, folds, scale = scaling)
    expect_equal(actual, expected, tolerance = 1e-8)
    expect_true(is.finite(actual))
  }
})

test_that("Mfold is invariant to shifting the origin of each variable", {
  data(nutrimouse)
  X <- base::scale(as.matrix(nutrimouse$lipid))
  Y <- base::scale(as.matrix(nutrimouse$gene[, 1:10]))
  set.seed(7)
  folds <- split(sample(seq_len(nrow(X))), rep(1:5, length.out = nrow(X)))
  shifted.X <- sweep(X, 2, seq_len(ncol(X)) * 10, "+")
  shifted.Y <- sweep(Y, 2, seq_len(ncol(Y)) * 10, "+")
  for (scaling in c(FALSE, TRUE)) {
    a <- mixOmics:::Mfold(X, Y, 0.1, 0.1, folds, scale = scaling)
    b <- mixOmics:::Mfold(shifted.X, shifted.Y, 0.1, 0.1, folds, scale = scaling)
    expect_equal(a, b, tolerance = 1e-8)
  }
})

test_that("both validation modes pass scaling through and match manual fold scores", {
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  grid <- expand.grid(c(0.1, 1), c(0.2, 1))
  for (validation in c("Mfold", "loo")) {
    set.seed(399)
    folds <- if (validation == "loo") {
      as.list(seq_len(nrow(X)))
    } else {
      split(sample(seq_len(nrow(X))), rep(1:4, length.out = nrow(X)))
    }
    for (scaling in c(FALSE, TRUE)) {
      expected <- vapply(seq_len(nrow(grid)), function(i) {
        .rcc_cv_reference(X, Y, folds, grid[i, 1], grid[i, 2], scale = scaling)
      }, numeric(1))
      fit <- tune.rcc(X, Y, grid1 = c(0.1, 1), grid2 = c(0.2, 1),
                      validation = validation, folds = 4, seed = 399, scale = scaling)
      expect_equal(fit$mat, matrix(expected, 2, 2), tolerance = 1e-8)
      expect_equal(fit$opt.score, max(expected), tolerance = 1e-8)
    }
  }
})

test_that("tune.rcc accepts scale before cross-validation options", {
  expect_identical(names(formals(tune.rcc)),
                   c("X", "Y", "grid1", "grid2", "scale", "validation",
                     "folds", "BPPARAM", "seed"))
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  for (validation in c("Mfold", "loo")) {
    for (scaling in c(FALSE, TRUE)) {
      named <- tune.rcc(X, Y, grid1 = c(0.1, 1), grid2 = c(0.2, 1),
                        scale = scaling, validation = validation, folds = 4,
                        BPPARAM = BiocParallel::SerialParam(), seed = 399)
      positional <- tune.rcc(X, Y, c(0.1, 1), c(0.2, 1), scaling,
                             validation, 4, BiocParallel::SerialParam(), 399)
      .almost_identical(positional, named)

      mixed <- tune.rcc(X, Y, c(0.1, 1), c(0.2, 1), scaling,
                        validation = validation, folds = 4, seed = 399)
      .almost_identical(mixed, named)
    }
  }
})

test_that("rcc tuning defaults to scaling and honours explicit choices in both entry points", {
  expect_identical(formals(tune.rcc)$scale, TRUE)
  expect_identical(formals(tune)$scale, TRUE)
  data(linnerud)
  args <- list(X = linnerud$exercise, Y = linnerud$physiological,
               grid1 = c(0.1, 1), grid2 = c(0.2, 1),
               validation = "Mfold", folds = 4, seed = 399)
  for (validation in c("Mfold", "loo")) {
    args$validation <- validation
    default <- do.call(tune.rcc, args)
    wrapped.default <- suppressMessages(do.call(tune, c(list(method = "rcc"), args)))
    scaled <- do.call(tune.rcc, c(args, list(scale = TRUE)))
    .almost_identical(default, scaled)
    .almost_identical(default, wrapped.default)
    for (scaling in c(FALSE, TRUE)) {
      direct <- do.call(tune.rcc, c(args, list(scale = scaling)))
      wrapped <- suppressMessages(do.call(tune, c(list(method = "rcc"), args,
                                                 list(scale = scaling))))
      .almost_identical(direct, wrapped)
      if (!scaling)
        expect_gt(max(abs(default$mat - direct$mat)), 1e-4)
    }
  }
})

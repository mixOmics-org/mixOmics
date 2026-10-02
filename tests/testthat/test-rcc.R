context("rcc")

test_that("rcc accepts scale before verbose.call", {
  expect_identical(names(formals(rcc)),
                   c("X", "Y", "ncomp", "method", "lambda1", "lambda2",
                     "scale", "verbose.call"))
  expect_identical(formals(rcc)$scale, TRUE)
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  for (scaling in c(FALSE, TRUE)) {
    named <- rcc(X, Y, ncomp = 2, method = "ridge",
                 lambda1 = 0.1, lambda2 = 0.2, scale = scaling)
    positional <- rcc(X, Y, 2, "ridge", 0.1, 0.2, scaling)
    .almost_identical(positional, named)
    expect_true(is.call(positional$call))

    verbose <- rcc(X, Y, 2, "ridge", 0.1, 0.2, scaling, TRUE)
    .almost_identical(verbose, named)
    expect_true(verbose$call$verbose.call)
    expect_identical(verbose$call$scale, scaling)
  }

  # Name both options to retain unscaled fits and verbose output explicitly.
  unscaled.verbose <- rcc(X, Y, 2, "ridge", 0.1, 0.2,
                        scale = FALSE, verbose.call = TRUE)
  expect_identical(unscaled.verbose$scale, list(X = FALSE, Y = FALSE))
  expect_identical(unscaled.verbose$call$scale, FALSE)
  expect_true(unscaled.verbose$call$verbose.call)
})

test_that("unscaled ridge solves the covariance-based regularised CCA problem", {
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  fit <- rcc(X, Y, lambda1 = 0.1, lambda2 = 0.2, scale = FALSE)
  expect_s3_class(fit, "rcc")
  expect_identical(fit$X, X)
  expect_identical(fit$Y, Y)
  expect_equal(fit$center, list(X = colMeans(X), Y = colMeans(Y)))
  expect_identical(fit$scale, list(X = FALSE, Y = FALSE))

  # Check the defining covariance constraints, independently of the
  # Cholesky/SVD implementation and arbitrary component signs.
  Cxx <- cov(X) + diag(0.1, ncol(X))
  Cyy <- cov(Y) + diag(0.2, ncol(Y))
  Cxy <- cov(X, Y)
  A <- fit$loadings$X
  B <- fit$loadings$Y
  expect_equal(unname(crossprod(A, Cxx %*% A)), diag(2))
  expect_equal(unname(crossprod(B, Cyy %*% B)), diag(2))
  expect_equal(unname(crossprod(A, Cxy %*% B)), diag(unname(fit$cor[1:2])))

  # The squared canonical correlations are the eigenvalues of this operator.
  eigenvalues <- eigen(solve(Cxx, Cxy) %*% solve(Cyy, t(Cxy)),
                       only.values = TRUE)$values
  expect_equal(unname(fit$cor)^2, sort(eigenvalues, decreasing = TRUE))
  expect_equal(fit$variates$X, base::scale(X, scale = FALSE) %*% A)
  expect_equal(fit$variates$Y, base::scale(Y, scale = FALSE) %*% B)

  # Explained variance is the fraction of input sum of squares captured
  # by projection onto each component.
  for (block in c("X", "Y")) {
    input <- fit[[block]]
    scores <- fit$variates[[block]]
    expected <- colSums(crossprod(input, scores)^2) /
      (colSums(scores^2) * sum(input^2))
    expect_equal(unname(fit$prop_expl_var[[block]]), unname(expected))
  }
})

test_that("ridge defaults to standardisation and stores its training statistics", {
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  fit <- rcc(X, Y, lambda1 = 0.1, lambda2 = 0.2)
  explicit <- rcc(X, Y, lambda1 = 0.1, lambda2 = 0.2, scale = TRUE)
  .almost_identical(fit, explicit)
  manual <- rcc(base::scale(X), base::scale(Y),
                lambda1 = 0.1, lambda2 = 0.2, scale = FALSE)
  expect_equal(fit$cor, manual$cor)
  expect_equal(abs(diag(cor(fit$variates$X, manual$variates$X))), rep(1, 2))
  expect_equal(abs(diag(cor(fit$variates$Y, manual$variates$Y))), rep(1, 2))
  expect_equal(fit$X, base::scale(X))
  expect_equal(fit$Y, base::scale(Y))
  expect_equal(fit$center, list(X = colMeans(X), Y = colMeans(Y)))
  expect_equal(fit$scale, list(X = apply(X, 2, sd), Y = apply(Y, 2, sd)))
  expect_equal(fit$X %*% fit$loadings$X, fit$variates$X)
  expect_equal(fit$Y %*% fit$loadings$Y, fit$variates$Y)
  expect_equal(fit$prop_expl_var, manual$prop_expl_var)

  # Already scaled input must not cause the old input attributes to be reused.
  expect_equal(manual$center$X, colMeans(base::scale(X)))
  expect_identical(manual$scale, list(X = FALSE, Y = FALSE))
})

test_that("scaled ridge is invariant to variable units with and without penalties", {
  data(nutrimouse)
  X <- as.matrix(nutrimouse$lipid)
  Y <- as.matrix(nutrimouse$gene[, 1:10])
  rescaled.X <- sweep(X, 2, 10^seq(-2, 3, length.out = ncol(X)), "*")
  rescaled.Y <- sweep(Y, 2, 10^seq(2, -2, length.out = ncol(Y)), "*")
  for (lambda in c(0, 0.1, 1)) {
    for (swap in c(FALSE, TRUE)) {
      a <- if (swap) list(Y, X) else list(X, Y)
      b <- if (swap) list(rescaled.Y, rescaled.X) else list(rescaled.X, rescaled.Y)
      fit <- rcc(a[[1]], a[[2]], lambda1 = lambda, lambda2 = lambda, scale = TRUE)
      rescaled <- rcc(b[[1]], b[[2]], lambda1 = lambda, lambda2 = lambda, scale = TRUE)
      expect_equal(fit$cor, rescaled$cor, tolerance = 1e-8)
      expect_equal(abs(diag(cor(fit$variates$X, rescaled$variates$X))), rep(1, 2))
      expect_equal(abs(diag(cor(fit$variates$Y, rescaled$variates$Y))), rep(1, 2))
    }
  }
})

test_that("the scale option leaves shrinkage calculations unchanged", {
  data(linnerud)
  a <- rcc(linnerud$exercise, linnerud$physiological, method = "shrinkage")
  b <- rcc(linnerud$exercise, linnerud$physiological, method = "shrinkage", scale = TRUE)
  unscaled <- rcc(linnerud$exercise, linnerud$physiological,
                  method = "shrinkage", scale = FALSE)
  .almost_identical(a, b)
  .almost_identical(a, unscaled)
  expect_identical(b$scale, list(X = FALSE, Y = FALSE))
})

test_that("rcc validates scaling and reports unscalable columns", {
  data(linnerud)
  X <- as.matrix(linnerud$exercise)
  Y <- as.matrix(linnerud$physiological)
  for (bad in list(NULL, NA, logical(), c(TRUE, FALSE), 1, "TRUE")) {
    expect_error(rcc(X, Y, scale = bad), "'scale' must be")
  }
  constant <- X
  constant[, 1] <- 1
  expect_error(rcc(constant, Y, scale = TRUE), "zero variance in 'X'")
  expect_error(rcc(X, constant, scale = TRUE), "zero variance in 'Y'")
  expect_error(rcc(constant, Y), "zero variance in 'X'")
  expect_silent(rcc(constant, Y, ncomp = 1, lambda1 = 1, lambda2 = 1,
                    scale = FALSE))
  missing <- X
  missing[, 1] <- NA
  expect_error(rcc(missing, Y, scale = TRUE), "at least two non-missing")
  missing[1, 1] <- 1
  expect_error(rcc(X, missing, scale = TRUE), "at least two non-missing")
})

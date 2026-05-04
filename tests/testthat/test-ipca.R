context("ipca")

# Test that ipca works
test_that("ipca works", {
  data(liver.toxicity)
  ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode="deflation")
  expect_equal(ipca.res$loadings$X[1:1], 0.0001637122)
})

# Test 1: Test invalid 'X' input (non-numeric matrix)
test_that("ipca throws error with invalid X input", {
  expect_error(ipca("invalid_input"), "'X' must be a numeric matrix.")
})

# Test 2: Test 'ncomp' validation (greater than the number of variables)
test_that("ipca throws error if ncomp is greater than the number of variables", {
  data(liver.toxicity)
  expect_error(ipca(liver.toxicity$gene, ncomp = 100), "use smaller 'ncomp'")
})

# Test 3: Test 'mode' parameter validation (unsupported mode)
test_that("ipca throws error with unsupported mode", {
  data(liver.toxicity)
  expect_error(ipca(liver.toxicity$gene, mode = "invalid_mode"), "'mode' should be one of 'deflation' or 'parallel'.")
})

# Test 4: Test function with valid 'X' input and default parameters
test_that("ipca runs with valid input", {
  data(liver.toxicity)
  ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode = "deflation")
  expect_s3_class(ipca.res, "ipca")
  expect_equal(ncol(ipca.res$x), 3)  # Check if the number of components is correct
})

# Test 5: Test 'w.init' parameter with custom matrix
test_that("ipca runs with custom w.init matrix", {
  data(liver.toxicity)
  w.init <- matrix(1, 3, 3)  # Custom w.init
  ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode = "deflation", w.init = w.init)
  expect_s3_class(ipca.res, "ipca")
})

# Test 6: Test 'fun' parameter with exp function
test_that("ipca runs with exp function in deflation mode", {
  data(liver.toxicity)
  ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode = "deflation", fun = "exp", max.iter = 20)
  expect_s3_class(ipca.res, "ipca")
  expect_equal(ncol(ipca.res$x), 3)
})

# Test 7: Test 'mode' parameter with parallel algorithm
test_that("ipca runs in parallel mode with logcosh", {
  data(liver.toxicity)
  ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode = "parallel", fun = "logcosh", max.iter = 20)
  expect_s3_class(ipca.res, "ipca")
  expect_equal(ncol(ipca.res$x), 3)
})

# Test 8: Test parallel algorithm with exp function
test_that("ipca runs in parallel mode with exp", {
  data(liver.toxicity)
  ipca.res <- ipca(liver.toxicity$gene, ncomp = 3, mode = "parallel", fun = "exp", max.iter = 20)
  expect_s3_class(ipca.res, "ipca")
  expect_equal(ncol(ipca.res$x), 3)
})

# Test 9: Test for missing values in 'X'
test_that("ipca throws error if missing values in X", {
  data(liver.toxicity)
  liver.toxicity$gene[1, 1] <- NA  # Inject missing value
  expect_error(ipca(liver.toxicity$gene), "missing values in 'X'.")
})

context("sipca")

# Test that sipca works
test_that("sipca works", {
  data(liver.toxicity)
  sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode="deflation", keepX=c(50,50,50))
  expect_equal(sipca.res$loadings$X[1], 0)
})

# Test 1: Test invalid 'X' input (non-numeric matrix)
test_that("sipca throws error with invalid X input", {
  expect_error(sipca("invalid_input"), "data must be in a matrix form")
})

# Test 2: What happens if ncomp larger than variables
test_that("sipca resets ncomp if it's too large", {
  data(liver.toxicity)
  res <- sipca(liver.toxicity$gene, ncomp = 100)
  expect_equal(res$ncomp, min(nrow(liver.toxicity$gene), ncol(liver.toxicity$gene)))  # Check if ncomp is reset
})

# # Test 3: Test invalid 'mode' argument
# test_that("sipca throws error with unsupported mode", {
#   expect_error(sipca(liver.toxicity$gene, mode = "invalid_mode"), 
#                "'arg' should be one of “deflation”, “parallel”")  # Match curly quotes
# })

# Test 4: Test 'w.init' parameter validation (invalid size)
test_that("sipca throws error if w.init is the wrong size", {
  data(liver.toxicity)
  expect_error(sipca(liver.toxicity$gene, w.init = matrix(1, 2, 3)), "w.init is not a matrix or is the wrong size")
})

# Test 5: Test for missing values in 'X'
test_that("sipca throws error if missing values in X", {
  data(liver.toxicity)
  liver.toxicity$gene[1, 1] <- NA  # Inject missing value
  expect_error(sipca(liver.toxicity$gene), "infinite or missing values in 'x'")
})

# Test 6: Test function with valid 'X' input and default parameters
test_that("sipca runs with valid input", {
  data(liver.toxicity)
  sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode = "deflation")
  expect_s3_class(sipca.res, "sipca")
  expect_equal(ncol(sipca.res$x), 3)  # Check if the number of components is correct
})

# Test 7: Test custom 'w.init' matrix
test_that("sipca runs with custom w.init matrix", {
  data(liver.toxicity)
  w.init <- matrix(1, 3, 3)  # Custom w.init
  sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode = "deflation", w.init = w.init)
  expect_s3_class(sipca.res, "sipca")
})

# Test 8: Test 'fun' parameter with exp function
test_that("sipca runs with exp function in deflation mode", {
  data(liver.toxicity)
  sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode = "deflation", fun = "exp")
  expect_s3_class(sipca.res, "sipca")
  expect_equal(ncol(sipca.res$x), 3)
})

# Test 9: Test 'mode' parameter with parallel algorithm
test_that("sipca runs in parallel mode with logcosh", {
  data(liver.toxicity)
  sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode = "parallel", fun = "logcosh")
  expect_s3_class(sipca.res, "sipca")
  expect_equal(ncol(sipca.res$x), 3)
})

# Test 10: Test parallel algorithm with exp function
test_that("sipca runs in parallel mode with exp", {
  data(liver.toxicity)
  sipca.res <- sipca(liver.toxicity$gene, ncomp = 3, mode = "parallel", fun = "exp")
  expect_s3_class(sipca.res, "sipca")
  expect_equal(ncol(sipca.res$x), 3)
})

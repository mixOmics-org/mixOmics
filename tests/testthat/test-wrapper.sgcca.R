context("wrapper.sgcca")

test_that("(wrapper.sgcca:basic): nutrimouse blocks", {
  data(nutrimouse)
  X <- list(gene = nutrimouse$gene[, 1:20],
            lipid = nutrimouse$lipid[, 1:10])
  design <- matrix(c(0, 1, 1, 0), nrow = 2)

  res <- wrapper.sgcca(
    X, design = design, ncomp = 2,
    keepX = list(gene = c(5, 5), lipid = c(5, 5))
  )
  selected <- sapply(res$loadings, function(x) colSums(x != 0))

  expect_s3_class(res, "sgcca")
  expect_equal(names(res$X), names(X))
  expect_equal(lapply(res$X, dim), lapply(X, dim))
  expect_equal(res$ncomp, c(gene = 2, lipid = 2))
  expect_equal(selected, matrix(5, nrow = 2, ncol = 2,
                                dimnames = list(c("comp1", "comp2"), names(X))))
})

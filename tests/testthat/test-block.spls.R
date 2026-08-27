context("block.spls")

# Regression test: keepX/keepY used to be extracted from the internal per-component
# 'keepA' structure, returning a malformed keepX and a NULL keepY
test_that("block.spls returns keepX and keepY as supplied", {
  data(breast.TCGA)
  mrna = breast.TCGA$data.train$mrna
  mirna = breast.TCGA$data.train$mirna
  Y = mrna[, 1, drop = FALSE]
  data = list(mrna = mrna[, -1], mirna = mirna)

  res = block.spls(data, Y, ncomp = 2, keepX = list(mrna = c(10, 10), mirna = c(20, 20)))

  expect_equal(res$keepX, list(mrna = c(10, 10), mirna = c(20, 20)))
  expect_equal(as.numeric(res$keepY), c(1, 1))
  expect_equal(unname(colSums(res$loadings$mrna != 0)), c(10, 10))
  expect_equal(unname(colSums(res$loadings$mirna != 0)), c(20, 20))

  # explicit keepY is honoured
  res2 = block.spls(data, Y, ncomp = 2,
                    keepX = list(mrna = c(5, 5), mirna = c(7, 7)), keepY = c(1, 1))
  expect_equal(res2$keepX, list(mrna = c(5, 5), mirna = c(7, 7)))
  expect_equal(as.numeric(res2$keepY), c(1, 1))
})

context("diabolo")

test_that("block.splsda works", {
  data(nutrimouse)
  Y = nutrimouse$diet

  data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
  design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE,
                  dimnames = list(c("gene", "lipid", "Y"),
                                  c("gene", "lipid", "Y")))


  nutrimouse.sgccda <- block.splsda(X = data,
                                    Y = Y,
                                    design = design,
                                    keepX = list(gene = c(10,10),
                                                 lipid = c(15,15)),
                                    ncomp = 2,
                                    tol = 1e-30)
  expect_length(nutrimouse.sgccda, 24L)

  expect_equal(names(nutrimouse.sgccda),
               c("call", "X", "Y", "ind.mat", "ncomp", "mode", "keepX", "keepY",
                 "variates", "loadings", "crit", "AVE", "names", "init", "tol",
                 "iter", "max.iter", "nzv", "scale", "design", "scheme", "indY",
                 "weights", "prop_expl_var"))

  expect_is(nutrimouse.sgccda$X, "list")
  expect_is(nutrimouse.sgccda$design, "matrix")
  expect_is(nutrimouse.sgccda$X$gene, "matrix")
  expect_is(nutrimouse.sgccda$X$lipid, "matrix")

  expect_equal(nutrimouse.sgccda$design, design)
  expect_equal(dim(nutrimouse.sgccda$X$gene), dim(data$gene))
  expect_equal(dim(nutrimouse.sgccda$X$lipid), dim(data$lipid))
  expect_equal(nutrimouse.sgccda$X$gene[1, 1], 0.524732755905559)
  expect_equal(nutrimouse.sgccda$X$lipid[1, 1], -0.528375020663953)
  expect_is(nutrimouse.sgccda$ind.mat, "matrix")

  expect_equal(dim(nutrimouse.sgccda$ind.mat), c(40L, 5L))

  expect_equal(nutrimouse.sgccda$Y, Y)

  expect_equal(nutrimouse.sgccda$ncomp, c("gene" = 2L, "lipid" = 2L, "Y" = 2L))

  expect_equal(nutrimouse.sgccda$mode, "regression")

  expect_equal(dim(nutrimouse.sgccda$loadings$gene), c(120L, 2L))
  expect_equal(dim(nutrimouse.sgccda$loadings$lipid), c(21L, 2L))
  expect_equal(dim(nutrimouse.sgccda$loadings$Y), c(5L, 2L))

  expect_null(nutrimouse.sgccda$nzv)

  expect_true(nutrimouse.sgccda$scale)

  expect_equal(nutrimouse.sgccda$indY, 3L)

  expect_equal(rowMeans(nutrimouse.sgccda$weights),
               c(gene = 0.694506104274723, lipid = 0.915845972615744))

  expect_length(nutrimouse.sgccda$prop_expl_var, 3L)
  expect_is(nutrimouse.sgccda$prop_expl_var, "list")
  expect_equal(names(nutrimouse.sgccda$prop_expl_var), colnames(design))

  expect_length(nutrimouse.sgccda$AVE, 3L)
  expect_equal(names(nutrimouse.sgccda$AVE),
               c("AVE_X", "AVE_outer", "AVE_inner"))
  expect_equal(nutrimouse.sgccda$AVE$AVE_outer[1], 0.217938372815004)
  expect_equal(nutrimouse.sgccda$AVE$AVE_inner[1], 0.663209598406049)
  expect_equal(nutrimouse.sgccda$AVE$AVE_X$Y[1], c(`comp1` = 0.25))

  expect_length(nutrimouse.sgccda$variates, 3L)
  expect_equal(names(nutrimouse.sgccda$variates), c("gene", "lipid", "Y"))
  expect_is(nutrimouse.sgccda$variates$gene, "matrix")
  expect_is(nutrimouse.sgccda$variates$lipid, "matrix")
  expect_is(nutrimouse.sgccda$variates$Y, "matrix")

  expect_equal(nutrimouse.sgccda$variates$gene[1, 1], 2.9424296984024)
  expect_equal(nutrimouse.sgccda$variates$lipid[1, 1], 2.73351593820324)
  expect_equal(nutrimouse.sgccda$variates$Y[1, 1], 0.639567998302767)
})

test_that("Check.entry.wrapper.mint.block handles zero-variance predictors", {
  samples = paste0("s", seq_len(8))
  X = list(gene = cbind(g1 = seq_len(8), zero = 1, g2 = 8:1),
           lipid = cbind(l1 = 2:9, l2 = 9:2))
  Y = cbind(class1 = rep(c(1, 0), 4), class2 = rep(c(0, 1), 4))
  rownames(X$gene) = rownames(X$lipid) = rownames(Y) = samples

  nzv.check = Check.entry.wrapper.mint.block(
    X = X, Y = Y, ncomp = 1, keepX = list(gene = 3, lipid = 2),
    DA = TRUE, init = "svd.single", scheme = "horst", scale = TRUE,
    near.zero.var = TRUE, mode = "regression", tol = 1e-06, max.iter = 10
  )
  expect_equal(ncol(nzv.check$A$gene), 2)
  expect_equal(nzv.check$keepA$gene, 2)

  expect_error(
    Check.entry.wrapper.mint.block(
      X = X, Y = Y, ncomp = 1, DA = TRUE, init = "svd.single",
      scheme = "horst", scale = TRUE, near.zero.var = FALSE,
      mode = "regression", tol = 1e-06, max.iter = 10
    ),
    "zero variance in block 'gene'", fixed = TRUE
  )
})

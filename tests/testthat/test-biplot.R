context("biplot")


## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2"
library(vdiffr)
set.seed(100)

## pca model
data("nutrimouse")
pca.lipid <- pca(nutrimouse$lipid, ncomp = 3, scale = TRUE)

test_that("biplot works for pca objects", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # cant run because not reproducible
  # simple plot showing sample names
  # set.seed(100)
  # invisible(capture.output(
  #   expect_doppelganger(
  #     title = "biplot plot pca model", 
  #     fig = biplot(pca.lipid, group = nutrimouse$diet))
  # ))
  
  # cant run because not reproducible
  # samples coloured by primary groups, sample names, default colours
  # set.seed(100)
  # invisible(capture.output(
  #   expect_doppelganger(
  #     title = "biplot plot pca with customised colours",
  #     fig = biplot(pca.lipid, group = nutrimouse$genotype,
  #                  col = c("red", "blue")))
  # ))
})

## pls model
data("nutrimouse")
pls.nutrimouse <- pls(X = nutrimouse$gene, Y = nutrimouse$lipid, ncomp = 2)

test_that("biplot works for pls objects", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # cant run because not reproducible
  # # simple plot showing sample names
  # set.seed(100)
  # invisible(capture.output(
  #   expect_doppelganger(
  #     title = "biplot plot pls model", 
  #     fig = biplot(pls.nutrimouse, group = nutrimouse$genotype, block = 'X',
  #                  legend.title = 'Genotype', cutoff = 0.878) )
  # ))
})


## 'plsda' model

data(breast.tumors)
X <- breast.tumors$gene.exp
colnames(X) <- paste0('GENE_', colnames(X))
rownames(X) <- paste0('SAMPLE_', rownames(X))
Y <- breast.tumors$sample$treatment
nrow(X)
grouping <- factor(c(rep("A", 20), rep("B", 20), rep("C", 7)))
plsda.breast <- plsda(X, Y, ncomp = 2)

test_that("biplot works for plsda objects", {
  skip_on_ci() # only run the vdiffr tests locally

  # simple plot showing sample names
  set.seed(100)
  invisible(capture.output(
    expect_doppelganger(
      title = "biplot plot plsda model control pch",
      fig = biplot(plsda.breast, cutoff = 0.72, ind.names = FALSE, pch = 2, pch.size = 5))
  ))

  # cant run because not reproducible
  # samples coloured by primary groups, sample names, default colours
  # set.seed(100)
  # invisible(capture.output(
  #   expect_doppelganger(
  #     title = "biplot plot plsda model with customised pch",
  #     fig = biplot(plsda.breast, cutoff = 0.72, pch = grouping, pch.size = 5))
  # ))
})

## ------------------------------------------------------------------------ ##
## Regression guard: ggplot2 size -> linewidth

test_that("(biplot:edge.case): no ggplot2 size-aesthetic deprecation warning", {
  # vline/hline default to FALSE; enable them to exercise the linewidth code paths
  expect_no_warning(biplot(plsda.breast, vline = TRUE, hline = TRUE))
})
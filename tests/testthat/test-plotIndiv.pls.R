context("plotIndiv.pls")

# Check every sample and displayed component, independently of model scale
# and arbitrary component signs.
.expect_plotIndiv_coordinates <- function(plot, expected, display.names = TRUE) {
  coordinates <- c("x", "y", "z")[seq_len(ncol(expected))]
  expect_equal(unname(as.matrix(plot$df[, coordinates, drop = FALSE])),
               unname(expected))
  if (display.names) {
    expect_equal(as.character(plot$df$names), rownames(expected))
  } else {
    expect_null(plot$df$names)
  }
  if (inherits(plot$graph, "ggplot")) {
    expect_equal(unname(as.matrix(plot$graph$data[, coordinates, drop = FALSE])),
                 unname(expected))
  }
}

## rcc model using the default standardisation
data(nutrimouse)
X <- nutrimouse$lipid
Y <- nutrimouse$gene
rcc.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)

## spls model
data(liver.toxicity)
X <- liver.toxicity$gene
Y <- liver.toxicity$clinic
spls.res <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                      keepY = c(10, 10, 10))

## splsda models
data(breast.tumors)
X <- breast.tumors$gene.exp
Y <- breast.tumors$sample$treatment
splsda.breast <- splsda(X, Y,keepX=c(10,10), ncomp=2)

data(srbct)
X <- srbct$gene
Y <- srbct$class
srbct.splsda <- splsda(X, Y, ncomp = 10)

## ------------------------------------------------------------------------ ##
## Test that outputs are correct when running default style = "ggplot2"

test_that("plotIndiv works for rcc", {
  
  pl.res <- plotIndiv(rcc.res)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               rbind(rcc.res$variates$X[, 1:2],
                                     rcc.res$variates$Y[, 1:2]))
  expect_equal(as.character(pl.res$df$Block),
               rep(c("Block: X", "Block: Y"), each = nrow(rcc.res$X)))
  
  pl.res <- plotIndiv(rcc.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
                      legend = TRUE)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               (rcc.res$variates$X[, 1:2] +
                                  rcc.res$variates$Y[, 1:2]) / 2)
  expect_equal(as.character(pl.res$df$group), as.character(nutrimouse$genotype))
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv maps selected rcc components in each representation space", {
  comp <- c(3, 1)
  X <- rcc.res$variates$X[, comp]
  Y <- rcc.res$variates$Y[, comp]
  expected <- list("X-variate" = X, "Y-variate" = Y,
                   "XY-variate" = (X + Y) / 2, "multi" = rbind(X, Y))
  for (space in names(expected)) {
    pl.res <- plotIndiv(rcc.res, comp = comp, rep.space = space)
    .expect_plotIndiv_coordinates(pl.res, expected[[space]])
  }
})

test_that("plotIndiv also represents an explicitly unscaled rcc model", {
  data(nutrimouse)
  fit <- rcc(nutrimouse$lipid, nutrimouse$gene, ncomp = 3,
              lambda1 = 0.064, lambda2 = 0.008, scale = FALSE)
  pl.res <- plotIndiv(fit, comp = c(2, 3), rep.space = "XY-variate")
  .expect_plotIndiv_coordinates(pl.res,
                               (fit$variates$X[, c(2, 3)] +
                                  fit$variates$Y[, c(2, 3)]) / 2)
})

test_that("plotIndiv works for (s)pls", {
  
  pl.res <- plotIndiv(spls.res, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col = c("red", "blue", "green", "black"),
                      pch = factor(liver.toxicity$treatment$Dose.Group),
                      legend.title.pch = 'Dose',
                      legend = TRUE)
  
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], 4.146771)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # check right number of samples
  expect_equal(dim(spls.res$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
  
})

test_that("plotIndiv works for (s)plsda", {
  
  pl.res <- plotIndiv(splsda.breast)
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -1.075222)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

test_that("plotIndiv works for (s)plsda and ellipses", {
  
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  pl.res <- plotIndiv(srbct.splsda , comp = 1:2, col = c("red", "blue", "green", "black"),
            group = groups, ind.names = FALSE,  # colour points by class
            ellipse = TRUE, # include 95% confidence ellipse for each class
            legend = TRUE, title = '(a) PLSDA with confidence ellipses')
  
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -6.83832)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "RMS"]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "NB"]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "EWS"]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BL"]), "black")
  # check right number of samples
  expect_equal(dim(srbct.splsda$X)[1], dim(pl.res$df)[1])
  # check ellipses
  expect_false(is.null(pl.res$df.ellipse))
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
})

test_that("plotIndiv works for (s)plsda and backgrounds", {
  data(srbct)
  X <- srbct$gene
  Y <- srbct$class
  srbct.splsda <- splsda(X, Y, ncomp = 10)
  groups <- factor(srbct$class, levels = c("RMS", "NB", "EWS", "BL"))
  background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")
  pl.res <- plotIndiv(srbct.splsda, comp = 1:2,
            group = srbct$class, ind.names = FALSE, # colour points by class
            background = background, # include prediction background for each class
            legend = TRUE, title = " (b) PLSDA with prediction background")
  
  # check coordinates
  .expect_numerically_close(pl.res$graph$data$x[1], -6.83832)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check right number of samples
  expect_equal(dim(srbct.splsda$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(groups)))
})

## ------------------------------------------------------------------------ ##

test_that("plotIndiv.rcc works without ind.names", code = {
    data(nutrimouse)
    X <- nutrimouse$lipid
    Y <- nutrimouse$gene
    nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
    
    plotIndiv.res <- plotIndiv(nutri.res, group = nutrimouse$genotype, ind.names = FALSE, legend = TRUE)
    
    expect_is(plotIndiv.res$graph, "ggplot")
})

## ------------------------------------------------------------------------ ##
## Plotting with 'lattice' style

test_that("plotIndiv works for rcc (lattice style)", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  pl.res <- plotIndiv(nutri.res, style = "lattice")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               rbind(nutri.res$variates$X[, 1:2],
                                     nutri.res$variates$Y[, 1:2]))
  
  pl.res <- plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
                      legend = TRUE, style = "lattice")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               (nutri.res$variates$X[, 1:2] +
                                  nutri.res$variates$Y[, 1:2]) / 2)
  expect_equal(as.character(pl.res$df$group), as.character(nutrimouse$genotype))
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pls (lattice style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col = c("red", "blue", "green", "black"),
                      pch = factor(liver.toxicity$treatment$Dose.Group),
                      legend.title.pch = 'Dose',
                      legend = TRUE, style = "lattice")
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 4.146771)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
  
})

test_that("plotIndiv works for (s)plsda (lattice style)", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y,keepX=c(10,10), ncomp=2)
  
  pl.res <- plotIndiv(splsda.breast, style = "lattice")
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], -1.075222)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

## ------------------------------------------------------------------------ ##
## Plotting with 'graphics' style

test_that("plotIndiv works for rcc (graphics style)", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  pl.res <- plotIndiv(nutri.res, style = "graphics")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               rbind(nutri.res$variates$X[, 1:2],
                                     nutri.res$variates$Y[, 1:2]))
  
  pl.res <- plotIndiv(nutri.res, rep.space= 'XY-variate', group = nutrimouse$genotype,
                      legend = TRUE, style = "graphics")
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               (nutri.res$variates$X[, 1:2] +
                                  nutri.res$variates$Y[, 1:2]) / 2)
  expect_equal(as.character(pl.res$df$group), as.character(nutrimouse$genotype))
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pls (graphics style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  pl.res <- plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                      group = factor(liver.toxicity$treatment$Time.Group),
                      legend.title = 'Time',
                      col = c("red", "blue", "green", "black"),
                      pch = factor(liver.toxicity$treatment$Dose.Group),
                      legend.title.pch = 'Dose',
                      legend = TRUE, style = "graphics")
  
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], 4.146771)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
  
})

test_that("plotIndiv works for (s)plsda (graphics style)", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y,keepX=c(10,10), ncomp=2)
  
  pl.res <- plotIndiv(splsda.breast, style = "graphics")
  # check coordinates
  .expect_numerically_close(pl.res$df[1,1], -1.075222)
  # check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

## ------------------------------------------------------------------------ ##
## Plotting with '3d' style

library(rgl)

test_that("plotIndiv works for rcc (3d style)", {
  data(nutrimouse)
  X <- nutrimouse$lipid
  Y <- nutrimouse$gene
  nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
  
  # Clear any existing rgl plots and disable plot rendering for non-interactive environments
  options(rgl.useNULL = TRUE)
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(nutri.res, style = "3d", pch = "sphere")))
  
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               rbind(nutri.res$variates$X[, 1:3],
                                     nutri.res$variates$Y[, 1:3]),
                               display.names = FALSE)
  
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(nutri.res, rep.space = 'XY-variate', group = nutrimouse$genotype,
                                                        legend = TRUE, style = "3d", pch = "sphere")))
  
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check coordinates
  .expect_plotIndiv_coordinates(pl.res,
                               (nutri.res$variates$X[, 1:3] +
                                  nutri.res$variates$Y[, 1:3]) / 2,
                               display.names = FALSE)
  expect_equal(as.character(pl.res$df$group), as.character(nutrimouse$genotype))
  # Check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(nutrimouse$genotype)))
})

test_that("plotIndiv works for (s)pls (3d style)", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  
  # Expect error when passing numbers as pch argument
  expect_error(suppressWarnings(suppressMessages(plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                                                           group = factor(liver.toxicity$treatment$Time.Group),
                                                           legend.title = 'Time',
                                                           col = c("red", "blue", "green", "black"),
                                                           pch = factor(liver.toxicity$treatment$Dose.Group),
                                                           legend.title.pch = 'Dose',
                                                           legend = TRUE, style = "3d"))),
               "pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.",
               fixed = TRUE)
  
  # Plot with correct pch values
  pchs <- factor(liver.toxicity$treatment$Dose.Group)
  levels(pchs) <- c("sphere", "tetra", "octa", "icosa")
  
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(toxicity.spls, rep.space="X-variate", ind.name = FALSE,
                                                        group = factor(liver.toxicity$treatment$Time.Group),
                                                        legend.title = 'Time',
                                                        col = c("red", "blue", "green", "black"),
                                                        pch = pchs,
                                                        legend.title.pch = 'Dose',
                                                        legend = TRUE, style = "3d")))
  
  # Check coordinates
  .expect_numerically_close(pl.res$df[1,1], 4.146771)
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check colour assignments are correct
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 6]), "red")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 18]), "blue")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 24]), "green")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == 48]), "black")
  # Check right number of samples
  expect_equal(dim(toxicity.spls$X)[1], dim(pl.res$df)[1])
  # Check groups
  expect_true(!is.null(pl.res$df$group))
  expect_equal(length(unique(pl.res$df$group)), length(unique(liver.toxicity$treatment$Time.Group)))
})

test_that("plotIndiv works for (s)plsda (3d style)", {
  data(breast.tumors)
  X <- breast.tumors$gene.exp
  Y <- breast.tumors$sample$treatment
  splsda.breast <- splsda(X, Y, keepX=c(10,10), ncomp=3)
  
  clear3d()
  pl.res <- suppressWarnings(suppressMessages(plotIndiv(splsda.breast, style = "3d", pch = "cube")))
  
  # Check coordinates
  .expect_numerically_close(pl.res$df[1,1], -1.075222)
  # Check correct output structure
  expect_equal(names(pl.res), c("df", "df.ellipse", "graph"))
  # Check colours
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "BE"]), "#F68B33")
  expect_equal(unique(pl.res$df$col[pl.res$df$group == "AF"]), "#388ECC")
  # Check right number of samples
  expect_equal(dim(splsda.breast$X)[1], dim(pl.res$df)[1])
})

# Clear the rgl device
if (rgl::rgl.cur() > 0) {
  rgl::close3d()
}

## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2" with rcc object
library(vdiffr)

test_that("plotIndiv works for rcc with sample names (default)", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot showing sample names
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot sample names", 
    fig = plotIndiv(rcc.res))
  ))
  # samples coloured by primary groups, sample names, default colours
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot sample names coloured by primary groups",
    fig = plotIndiv(rcc.res, group = nutrimouse$diet))
  ))
  # samples coloured by primary groups, sample names, user-defined colours
  invisible(capture.output(
  expect_doppelganger(
    title = "PCA plot sample names coloured by primary groups custom cols",
    fig = plotIndiv(rcc.res, group = nutrimouse$diet, col = c("red", "blue", "black", "pink", "grey"), legend = TRUE))
  ))
})

test_that("plotIndiv works for rcc without sample names", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot", 
    fig = plotIndiv(rcc.res, ind.names = FALSE))
  # samples coloured by primary groups, default colours, by default shapes also match primary groups
  ))
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot coloured by primary groups", 
    fig = plotIndiv(rcc.res, ind.names = FALSE,
                    group = nutrimouse$diet, legend = TRUE))
  ))
  # samples coloured by primary groups, user-defined colours, groups reordered, by default shapes also match primary groups
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot coloured by primary groups custom cols reordered groups", 
    fig = plotIndiv(rcc.res, ind.names = FALSE,
                    group = factor(nutrimouse$diet, levels = c("coc", "lin", "ref", "sun", "fish")), 
                    col = c("red", "blue", "black", "pink", "grey"), legend = TRUE))
  ))
})

test_that("plotIndiv works for rcc with ellipse on groups", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot with ellipse coloured by primary groups", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, group = nutrimouse$diet, ellipse = TRUE, legend = TRUE))
  ))
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse confidence
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot with ellipse coloured by primary groups, ellipse level 0.5", 
    fig = plotIndiv(rcc.res, ind.names = TRUE, group = nutrimouse$diet, ellipse = TRUE, 
                    legend = TRUE, ellipse.level = 0.5))
  ))
})

test_that("plotIndiv works for rcc with centroids on groups", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot with centroids coloured by primary groups", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, group = nutrimouse$diet, centroid = TRUE, 
                    legend = TRUE))
  ))
})

test_that("plotIndiv works for rcc controlling colours and pch", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # force pch to be the same for all samples, so groups only coloured differently
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot coloured by primary groups custom cols with set pch (circle) for all samples", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, group = nutrimouse$diet, 
                    legend = TRUE, pch = 1))
  ))
  # control the pch of each of the primary groups by setting pch with length = levels(primary_groups)
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot coloured by primary groups with set pch for each group", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, group = nutrimouse$diet, 
                    legend = TRUE, pch = c(2, 4, 6, 3, 1)))
  ))
  # use pch to show a whole new grouping (secondary groups)
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot coloured by primary groups with pch for secondary groups", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, group = nutrimouse$diet, 
                    legend = TRUE, pch = as.factor(c(rep("A", 20), rep("B", 20)))))
  ))
  # only pch grouping
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot with pch for primary groups, col consistent", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, col = "purple",
                    legend = TRUE, pch = as.factor(c(rep("A", 20), rep("B", 20))), legend.title.pch = "Groups"))
  ))
})

test_that("plotIndiv works for rcc with different rep space", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples projected on XY variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot on XY variate space", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, rep.space = "XY-variate"))
  ))
  # samples projected on X variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot on X variate space", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, rep.space = "X-variate"))
  ))
  # samples projected on Y variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot on Y variate space", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, rep.space = "Y-variate"))
  ))
  # samples projected on multivariate space
  invisible(capture.output(
  expect_doppelganger(
    title = "rCCA plot on multi variate space", 
    fig = plotIndiv(rcc.res, ind.names = FALSE, rep.space = "multi"))
  ))
})

## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2" with pls object

test_that("plotIndiv works for sPLS with sample names (default)", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot showing sample names
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot sample names", 
    fig = plotIndiv(spls.res))
  ))
  # samples coloured by primary groups, sample names, default colours
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot sample names coloured by primary groups",
    fig = plotIndiv(spls.res, group = liver.toxicity$treatment$Treatment.Group))
  ))
})

test_that("plotIndiv works for sPLS without sample names", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot", 
    fig = plotIndiv(spls.res, ind.names = FALSE))
  ))
  # samples coloured by primary groups, default colours, by default shapes also match primary groups
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot coloured by primary groups", 
    fig = plotIndiv(spls.res, ind.names = FALSE,
                    group = liver.toxicity$treatment$Treatment.Group, legend = TRUE))
  ))
})

test_that("plotIndiv works for sPLS with ellipse on groups", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot with ellipse coloured by primary groups", 
    fig = plotIndiv(spls.res, ind.names = FALSE, group = liver.toxicity$treatment$Treatment.Group, ellipse = TRUE, legend = TRUE))
  ))
})

test_that("plotIndiv works for spls with centroids on groups", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot with centroids and stars coloured by primary groups", 
    fig = plotIndiv(spls.res, ind.names = FALSE, group = liver.toxicity$treatment$Treatment.Group, centroid = TRUE, 
                    legend = TRUE, star = TRUE))
  ))
})

test_that("plotIndiv works for sPLS with different rep space", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples projected on XY variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot on XY variate space", 
    fig = plotIndiv(spls.res, ind.names = FALSE, rep.space = "XY-variate"))
  ))
  # samples projected on X variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot on X variate space", 
    fig = plotIndiv(spls.res, ind.names = FALSE, rep.space = "X-variate"))
  ))
  # samples projected on Y variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot on Y variate space", 
    fig = plotIndiv(spls.res, ind.names = FALSE, rep.space = "Y-variate"))
  ))
  # samples projected on multivariate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLS plot on multi variate space", 
    fig = plotIndiv(spls.res, ind.names = FALSE, rep.space = "multi"))
  ))
})

## ------------------------------------------------------------------------ ##
## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2" with plsda object

test_that("plotIndiv works for sPLSDA with background", {
  skip_on_ci() # only run the vdiffr tests locally
  
  background = background.predict(srbct.splsda, comp.predicted=2, dist = "max.dist")
  invisible(capture.output(
  expect_doppelganger(
    title = "mint sPLSDA plot with background", 
    fig = plotIndiv(srbct.splsda, background = background, legend = TRUE))
  ))
})

test_that("plotIndiv works for sPLSDA with sample names (default)", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot showing sample names, by default coloured by groups
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot sample names", 
    fig = plotIndiv(srbct.splsda))
  ))
  # samples coloured by custom groups, sample names
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot sample names coloured by custom groups",
    fig = plotIndiv(srbct.splsda, group = as.factor(c(rep("A", 30), rep("B", 33)))))
  ))
})

test_that("plotIndiv works for sPLSDA without sample names", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # coloured and shapes by default by extracted groups
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot coloured by extracted groups", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, legend = TRUE))
  ))
})

test_that("plotIndiv works for sPLSDA with ellipse from predictions", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, ellipse
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot with ellipse", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, ellipse = TRUE, legend = TRUE))
  ))
  # samples coloured by primary groups and ellipse on custom groups
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot with ellipse on custom groups", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, group = as.factor(c(rep("A", 30), rep("B", 33))),
                    ellipse = TRUE, legend = TRUE))
  ))
})

test_that("plotIndiv works for sPLSDA with centroids on groups", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot with centroids and stars", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, centroid = TRUE, 
                    legend = TRUE, star = TRUE))
  ))
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid, custom cols
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot with centroids and stars, custom cols", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, centroid = TRUE, 
                    legend = TRUE, star = TRUE,
                    col = c("red", "purple", "orange", "yellow")))
  ))
  # samples coloured by primary groups, default colours, by default shapes also match primary groups, centroid, custom cols, pch on secondary groups
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot with centroids and stars, custom cols, pch on second grouping", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, centroid = TRUE, 
                    legend = TRUE, star = TRUE,
                    col = c("red", "purple", "orange", "yellow"),
                    pch = as.factor(c(rep("A", 30), rep("B", 33)))))
  ))
})

test_that("plotIndiv works for sPLSDA with different rep space", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # samples projected on XY variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot on XY variate space", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, rep.space = "XY-variate"))
  ))
  # samples projected on X variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot on X variate space", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, rep.space = "X-variate"))
  ))
  # samples projected on Y variate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot on Y variate space", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, rep.space = "Y-variate"))
  ))
  # samples projected on multivariate space
  invisible(capture.output(
  expect_doppelganger(
    title = "sPLSDA plot on multi variate space", 
    fig = plotIndiv(srbct.splsda, ind.names = FALSE, rep.space = "multi"))
  ))
})

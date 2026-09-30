context("plotVar")
## ------------------------------------------------------------------------ ##
test_that("plotVar works for pls with var.names", {
  data(nutrimouse)
  x <- nutrimouse$gene
  y <- nutrimouse$lipid
  ## custom var.names
  var.names <- list(x = seq_along(x), y = seq_along(y))
  
  pls.res <- pls(x, y) 
  df <- plotVar(pls.res , var.names = var.names, plot = FALSE)
  
  var.names.char.vec <- unname(unlist(lapply(var.names, as.character)))
  
  expect_true(all(df$names == var.names.char.vec))
  
  ## ------------- spls
  spls.res <- spls(x, y , keepX = c(10, 10))
  df <- plotVar(spls.res , var.names = var.names, plot = FALSE)
  expect_true(is(df, 'data.frame'))
  expect_true(all(df$names %in% as.character(unlist(var.names))))
  
  ## ------------- spca
  var.names = list(seq_along(x))
  spca.res <- spca(x, keepX = c(10, 10))
  df <- plotVar(spca.res, var.names = var.names, plot = FALSE)
  expect_true(all(df$names %in% as.character(unlist(var.names))))
  
})

test_that("plotVar works in block.(s)PLS1 cases", {
  
  data(breast.TCGA)
  X <- list(miRNA = breast.TCGA$data.train$mirna[,1:10],
            mRNA = breast.TCGA$data.train$mrna[,1:10])
  
  Y <- matrix(breast.TCGA$data.train$protein[,4], ncol=1)
  rownames(Y) <- rownames(X$miRNA)
  colnames(Y) <- "response"
  
  block.pls.result <- block.spls(X, Y, design = "full",
                                 keepX = list(miRNA=c(3,3),
                                              mRNA=c(3,3)))
  
  plotVar.result <- plotVar(block.pls.result, plot = FALSE)
  
  expect_equal(as.character(unique(plotVar.result$Block)), c("miRNA", "mRNA"))
})
test_that("plotVar leaves the coordinates unchanged when the components are orthogonal", {
  data(nutrimouse)
  pls.res <- pls(nutrimouse$gene, nutrimouse$lipid, ncomp = 3)
  expect_no_message(df <- plotVar(pls.res, comp = c(2, 3), plot = FALSE))
  expected <- rbind(cor(pls.res$X, pls.res$variates$X[, 2:3]),
                    cor(pls.res$Y, pls.res$variates$X[, 2:3]))
  expect_equal(cbind(df$x, df$y), unname(expected), tolerance = 1e-12)
  expect_equal(df$names, rownames(expected))
})

test_that("plotVar rejects the same component twice", {
  data(nutrimouse)
  pls.res <- pls(nutrimouse$gene, nutrimouse$lipid, ncomp = 2)
  expect_error(plotVar(pls.res, comp = c(1, 1), plot = FALSE), "distinct")
})

test_that("plotVar orthogonalises correlated components so the variables stay inside the circle", {
  data(liver.toxicity)
  X <- as.matrix(liver.toxicity$gene)
  Y <- as.matrix(liver.toxicity$clinic)
  res <- suppressMessages(block.pls(list(gene = X), Y, ncomp = 3))
  u <- res$variates$Y[, 2:3]
  # the Y block is deflated with the gene components, so its own components are correlated
  expect_gt(abs(cor(u)[1, 2]), 0.3)
  expect_message(df <- plotVar(res, comp = c(2, 3), blocks = "Y", plot = FALSE), "correlated")
  # inside the circle, with the squared radius equal to the variance of each variable
  # explained by the two components
  radius2 <- df$x^2 + df$y^2
  expect_true(all(radius2 <= 1))
  r2 <- apply(res$X$Y, 2, function(y) summary(lm(y ~ u))$r.squared)
  expect_equal(radius2, unname(r2), tolerance = 1e-10)
  # reversing comp only swaps the coordinates
  swapped <- suppressMessages(plotVar(res, comp = c(3, 2), blocks = "Y", plot = FALSE))
  expect_equal(swapped$x, df$y)
  expect_equal(swapped$y, df$x)
  # symmetric: both components are rotated by the same angle, half of the angle
  # by which they miss being perpendicular
  axes <- .cor_orthogonalised(u, u)
  expect_equal(axes[1, 1], axes[2, 2])
  expect_equal(axes[1, 2], axes[2, 1])
  expect_equal(acos(axes[1, 1]), abs(acos(cor(u)[1, 2]) - pi / 2) / 2)
  # the axis label says so
  pdf(NULL)
  on.exit(dev.off())
  suppressMessages(plotVar(res, comp = c(2, 3), blocks = "Y"))
  p <- ggplot2::last_plot()
  labs <- if ("get_labs" %in% getNamespaceExports("ggplot2")) ggplot2::get_labs(p) else p$labels
  expect_equal(labs$x, "Component 2 (orthogonalised)")
  expect_equal(labs$y, "Component 3 (orthogonalised)")
})

test_that("plotVar with missing values correlates each variable on its observed samples", {
  data(nutrimouse)
  set.seed(1)
  X <- as.matrix(nutrimouse$lipid)
  X[sample(length(X), 0.3 * length(X))] <- NA
  res <- suppressMessages(pca(X, ncomp = 2))
  # pairwise correlations with the components as they are can leave the circle
  plain <- cor(res$X, res$variates$X, use = "pairwise")
  expect_true(any(rowSums(plain^2) > 1))
  expect_message(df <- plotVar(res, plot = FALSE), "missing values")
  expect_true(all(df$x^2 + df$y^2 <= 1))
  # the squared radius of each variable is the variance explained by the two
  # components on the samples where that variable is observed
  r2 <- apply(res$X, 2, function(x) summary(lm(x ~ res$variates$X[, 1:2]))$r.squared)
  expect_equal(df$x^2 + df$y^2, unname(r2), tolerance = 1e-10)
})

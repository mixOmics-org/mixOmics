context("network")

# rCCA similarities are inner products of variable correlations with the
# canonical bisectors; the cutoff removes weak similarities.
.rcc_network_reference <- function(fit, comp, cutoff) {
  bisectors <- fit$variates$X[, comp, drop = FALSE] +
    fit$variates$Y[, comp, drop = FALSE]
  similarities <- cor(fit$X, bisectors, use = "pairwise") %*%
    t(cor(fit$Y, bisectors, use = "pairwise"))
  similarities[abs(similarities) < cutoff] <- 0
  similarities
}

test_that("network works for the default standardised rcc model", {
  data(nutrimouse)
  nutri.res <- rcc(nutrimouse$lipid, nutrimouse$gene, ncomp = 3,
                   lambda1 = 0.064, lambda2 = 0.008)

  # network appends the extension to name.save.
  tmp.file <- tempfile("network")
  on.exit(unlink(paste0(tmp.file, ".jpeg")))
  network.rcc.res <- network(nutri.res, comp = 1:3, cutoff = 0.6, save = "jpeg",
                             name.save = tmp.file)
  expected <- .rcc_network_reference(nutri.res, comp = 1:3, cutoff = 0.6)
  expect_equal(names(network.rcc.res), c("gR", "M", "cutoff"))
  expect_equal(network.rcc.res$M, expected)
  expect_equal(network.rcc.res$cutoff, 0.6)
  expect_s3_class(network.rcc.res$gR, "igraph")
  expect_equal(igraph::ecount(network.rcc.res$gR), sum(expected != 0))
  expect_equal(sort(igraph::edge_attr(network.rcc.res$gR, "weight")),
               sort(expected[expected != 0]))
  expect_true(file.exists(paste0(tmp.file, ".jpeg")))
  expect_gt(file.info(paste0(tmp.file, ".jpeg"))$size, 0)
})

test_that("network also represents an explicitly unscaled rcc model", {
  data(nutrimouse)
  nutri.res <- rcc(nutrimouse$lipid, nutrimouse$gene, ncomp = 3,
                   lambda1 = 0.064, lambda2 = 0.008, scale = FALSE)
  pdf(NULL)
  on.exit(dev.off())
  net <- network(nutri.res, comp = 1:3, cutoff = 0.6, plot.graph = FALSE)
  expected <- .rcc_network_reference(nutri.res, comp = 1:3, cutoff = 0.6)
  expect_equal(net$M, expected)
  expect_equal(net$cutoff, 0.6)
  expect_equal(igraph::ecount(net$gR), sum(expected != 0))
  expect_equal(sort(igraph::edge_attr(net$gR, "weight")),
               sort(expected[expected != 0]))
})

test_that("network works for spls", {
  data(liver.toxicity)
  X <- liver.toxicity$gene
  Y <- liver.toxicity$clinic
  toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
                        keepY = c(10, 10, 10))
  ## create a tmp file
  tmp.file <- tempfile("network", fileext = ".jpeg")
  network.spls.res <- network(toxicity.spls, comp = 1:3, cutoff = 0.8,
          color.node = c("mistyrose", "lightcyan"),
          shape.node = c("rectangle", "circle"),
          color.edge = color.spectral(100),
          lty.edge = "solid", lwd.edge =  1,
          show.edge.labels = FALSE, interactive = FALSE, save = "jpeg", 
          name.save = tmp.file)
  
  expect_equal(names(network.spls.res), c("gR", "M", "cutoff"))
  .expect_numerically_close(sum(network.spls.res$M), 45.6061, digits = 3)
  unlink(tmp.file)
})

test_that("network works when `shape.node` == 'none'", {
  data(breast.TCGA)
  X <- list(mRNA = breast.TCGA$data.train$mrna[,1:100],
            proteomics = breast.TCGA$data.train$protein[,1:100])
  Y <- breast.TCGA$data.train$subtype
  
  model.pls <- pls(X$mRNA, X$proteomics)
  model.block <- block.splsda(X, Y)
  
  tmp.file <- tempfile("network", fileext = ".jpeg")
  
  net <- network(model.pls, cutoff = 0.5,
                 shape.node = c("none", "none"), 
                 save = "jpeg", 
                 name.save = tmp.file)
  expect_equal(names(net), c("gR", "M", "cutoff"))
  .expect_numerically_close(sum(net$M), 11.17439, digits = 3)
  unlink(tmp.file)
  
  net <- network(model.block, cutoff = 0.5,
                 shape.node = c("none", "none"), 
                 save = "jpeg", 
                 name.save = tmp.file)
  expect_equal(names(net), c("gR", "M_mRNA_proteomics", "cutoff"))
  .expect_numerically_close(sum(net$M), 14.36216, digits = 3)
  unlink(tmp.file)
})

test_that("catches invalid `graph.scale` values", {
  data(breast.TCGA)
  X <- list(mRNA = breast.TCGA$data.train$mrna[,1:100],
            proteomics = breast.TCGA$data.train$protein[,1:100])
  Y <- breast.TCGA$data.train$subtype
  
  model <- block.splsda(X, Y)
  
  expect_error(network(model, cutoff = 0.5,
                       graph.scale = -1),
               "graph.scale")
  expect_error(network(model, cutoff = 0.5,
                       graph.scale = 2),
               "graph.scale")
})

test_that("catches invalid `size.node` values", {
  data(breast.TCGA)
  X <- list(mRNA = breast.TCGA$data.train$mrna[,1:100],
            proteomics = breast.TCGA$data.train$protein[,1:100])
  Y <- breast.TCGA$data.train$subtype
  
  model <- block.splsda(X, Y)
  
  expect_error(network(model, cutoff = 0.5,
                       size.node = -1),
               "size.node")
  expect_error(network(model, cutoff = 0.5,
                       size.node = 2),
               "size.node")
})

unlink(list.files(pattern = "*.pdf"))


test_that("network plot.graph parameter does not affect numerical output", {
  data("nutrimouse")
  X <- nutrimouse$gene
  Y <- nutrimouse$lipid
  
  pls.obj <- pls(X, Y)
  
  network.obj.F <- network(pls.obj, plot.graph = F)
  network.obj.T <- network(pls.obj, plot.graph = T)
  
  expect_equal(network.obj.F$M, network.obj.T$M)
})

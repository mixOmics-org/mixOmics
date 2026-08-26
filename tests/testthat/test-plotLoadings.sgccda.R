context("plotLoadings.sgccda")

# Create a sgccda object for testing
data(nutrimouse)
Y = nutrimouse$diet
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
diablo.obj = wrapper.sgccda(X = data,
                            Y = Y,
                            design = design,
                            keepX = list(gene = c(10,10), lipid = c(15,15)),
                            ncomp = 2)

# Unit test 1: Test default behavior with graphics style
test_that("Test default behavior with graphics style", {
  skip_on_cran()
  skip_on_ci()

  png(tempfile(), width = 1200, height = 1000, res = 150)
  old_par <- par(no.readonly = TRUE)  # Save current par settings
  par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
  result <- plotLoadings(diablo.obj, comp = 2, style = "graphics")
  par(old_par)  # Restore original par settings
  dev.off()
  expect_equal(length(result), 2)  # Expect 2 blocks (gene and lipid)
  expect_s3_class(result[[1]], "data.frame")
  expect_equal(ncol(result[[1]]), 1)  # Expect 1 column
})

# Unit test 2: Test block-specific functionality
test_that("Test block-specific plotting", {
  skip_on_cran()
  skip_on_ci()

  png(tempfile(), width = 1200, height = 1000, res = 150)
  old_par <- par(no.readonly = TRUE)  # Save current par settings
  par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
  result <- plotLoadings(diablo.obj, block = "lipid", style = "graphics")
  par(old_par)  # Restore original par settings
  dev.off()
  expect_equal(length(result), 1)  # Only one block
  expect_equal(names(result), "lipid")
})

# Unit test 3: Test contrib parameter
test_that("Test contrib parameter functionality", {
  skip_on_cran()
  skip_on_ci()

  png(tempfile(), width = 1200, height = 1000, res = 150)
  old_par <- par(no.readonly = TRUE)  # Save current par settings
  par(mar = c(8, 4, 4, 2))  # Increase bottom margin to fit long names
  result <- plotLoadings(diablo.obj, contrib = "max", style = "graphics")
  par(old_par)  # Restore original par settings
  dev.off()
  expect_equal(length(result), 2)  # Both blocks
  expect_true(all(sapply(result, function(x) "color" %in% colnames(x))))  # Check for color column
})

test_that("plotLoadings supports vector xlim with sgccda contrib", {
  skip_on_cran()
  skip_on_ci()

  png(tempfile(), width = 1200, height = 1000, res = 150)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    dev.off()
  })
  result <- plotLoadings(diablo.obj, contrib = "max", xlim = c(-0.5, 0.5))
  expect_equal(length(result), 2)
  expect_true(all(sapply(result, function(x) "color" %in% colnames(x))))
})

# Error: plotLoadings encountered margin errors. Ensure feature names are not too long and the "Plots" pane is enlarged.
# # Unit test 4: Test method parameter
# test_that("Test method parameter functionality", {
#   result <- plotLoadings(diablo.obj, contrib = "max", method = "median", style = "graphics")
#   expect_equal(length(result), 2)
#   expect_equal(nrow(result[[1]]), 10)
# })

# Error: plotLoadings encountered margin errors. Ensure feature names are not too long and the "Plots" pane is enlarged.
# # Unit test 5: Test show.ties parameter
# test_that("Test show.ties parameter functionality", {
#   result <- plotLoadings(diablo.obj, contrib = "max", show.ties = FALSE, style = "graphics")
#   expect_equal(length(result), 2)
#   expect_equal(ncol(result[[1]]), 14)
# })

# Error: plotLoadings encountered margin errors. Ensure feature names are not too long and the "Plots" pane is enlarged.
# # Unit test 7: Test legend parameters
# test_that("Test legend parameters functionality", {
#   result <- plotLoadings(diablo.obj, contrib = "max", legend = TRUE, 
#                         legend.color = c("red", "blue", "green", "yellow", "black"), 
#                         legend.title = "Test Legend", 
#                         style = "graphics")
#   expect_equal(length(result), 2)
#   expect_equal(nrow(result[[1]]), 10)
# })

## ------------------------------------------------------------------------ ##
## vdiffr testing - "graphics"
library(vdiffr)

test_that("plotLoadings works for sgccda graphics", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo simple - graphics", 
      fig = plotLoadings(diablo.obj))
  ))
  
  # block-specific plot with contrib
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo block specific with contrib - graphics",
      fig = plotLoadings(diablo.obj, block = "lipid", contrib = "max"))
  ))
  
  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo change labels and label sizes - graphics",
      fig = plotLoadings(diablo.obj, X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.axis = 1))
  ))
  
  # change legend
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo change legend - graphics",
      fig = plotLoadings(diablo.obj, contrib = "max", 
                        legend.color = c("red", "blue", "green", "yellow", "black"), 
                        legend.title = "Test Legend"))
  ))
})

## ------------------------------------------------------------------------ ##
## vdiffr testing - "ggplot2"

test_that("plotLoadings works for sgccda ggplot2", {
  skip_on_ci() # only run the vdiffr tests locally
  
  # simple plot
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo simple - ggplot2", 
      fig = plotLoadings(diablo.obj, style = "ggplot2"))
  ))
  
  # block-specific plot with contrib
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo block specific with contrib - ggplot2",
      fig = plotLoadings(diablo.obj, style = "ggplot2", block = "lipid", contrib = "max"))
  ))
  
  # change labels and label sizes
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo change labels and label sizes - ggplot2",
      fig = plotLoadings(diablo.obj, style = "ggplot2", X.label = "X", Y.label = "Y", title = "Title Test",
                         size.labs = 2, size.title = 3, size.axis = 1))
  ))
  
  # change legend
  invisible(capture.output(
    expect_doppelganger(
      title = "Loadings plot diablo change legend - ggplot2",
      fig = plotLoadings(diablo.obj, style = "ggplot2", contrib = "max", 
                        legend.color = c("red", "blue", "green", "yellow", "black"), 
                        legend.title = "Test Legend"))
  ))
})

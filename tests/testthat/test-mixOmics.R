context("mixOmics")

test_that("Test that PLS works with SE", {
    data(liver.toxicity)
    X <- liver.toxicity$gene
    Y <- liver.toxicity$clinic
    Y.factor <- as.factor(liver.toxicity$treatment[, 4])
    #
    out_ref <- mixOmics(X, Y, ncomp = 2)
    # SummarizedExperiment
    library(SummarizedExperiment)
    assay <- t(X)
    coldata <- DataFrame(Y)
    coldata[["outcome"]] <- Y.factor
    se <- SummarizedExperiment(
        assays = SimpleList(abundance = assay),
        colData = coldata
    )
    columns <- colnames(colData(se))
    columns <- columns[ !columns %in% c("outcome") ]
    out = mixOmics(se, assay.type = "abundance", col.var = columns, ncomp = 2)
    #
    expect_equal(out$Y, out_ref$Y)
})

test_that("Test that block sPLS-DA works with MAE", {
    data(nutrimouse)
    data <- list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    #
    out_ref <- mixOmics(X = data, Y = nutrimouse$diet, ncomp = 3,
                        keepX = list(gene = c(10,10), lipid = c(15,15)))
    # MultiAssayExperiment
    library(SummarizedExperiment)
    library(MultiAssayExperiment)
    data(nutrimouse)
    Y = unmap(nutrimouse[["diet"]])
    se1 <- SummarizedExperiment(
        assays = SimpleList(abundance = t(nutrimouse[["gene"]])))
    se2 <- SummarizedExperiment(
        assays = SimpleList(abundance = t(nutrimouse[["lipid"]])))
    coldata <- DataFrame(diet = nutrimouse[["diet"]])
    rownames(coldata) <- colnames(se1)
    mae <- MultiAssayExperiment(
        experiments = ExperimentList(gene = se1, lipid = se2),
        colData = coldata
    )
    columns <- colnames(colData(mae))
    out <- mixOmics(
        mae, experiments = c("gene", "lipid"),
        assay.type = c("abundance", "abundance"),
        col.var = "diet", ncomp = 3,
        keepX = list(gene = c(10,10), lipid = c(15,15))
    )
    #
    expect_equal(out$loadings, out_ref$loadings)
})
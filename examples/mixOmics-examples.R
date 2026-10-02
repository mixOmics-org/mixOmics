## -- directed towards PLS framework because X is a matrix and the study argument is missing
# ----------------------------------------------------
data(liver.toxicity)
X = liver.toxicity$gene
Y = liver.toxicity$clinic
Y.factor = as.factor(liver.toxicity$treatment[, 4])

# directed towards PLS
out = mixOmics(X, Y, ncomp = 2)

# directed towards sPLS because of keepX and/or keepY
out = mixOmics(X, Y, ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))

# directed towards PLS-DA because Y is a factor
out = mixOmics(X, Y.factor, ncomp = 2)

# directed towards sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(X, Y.factor, ncomp = 2, keepX = c(20, 20))


\dontrun{
## -- directed towards block.pls framework because X is a list
# ----------------------------------------------------
data(nutrimouse)
Y = unmap(nutrimouse$diet)
data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)

# directed towards block PLS
out = mixOmics(X = data, Y = Y,ncomp = 3)

# directed towards block sPLS because of keepX and/or keepY
out = mixOmics(X = data, Y = Y,ncomp = 3,
keepX = list(gene = c(10,10), lipid = c(15,15)))

# directed towards block PLS-DA because Y is a factor
out = mixOmics(X = data, Y = nutrimouse$diet, ncomp = 3)

# directed towards block sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(X = data, Y = nutrimouse$diet, ncomp = 3,
keepX = list(gene = c(10,10), lipid = c(15,15)))


## -- directed towards mint.pls framework because of the study factor
# ----------------------------------------------------
data(stemcells)
# directed towards PLS
out = mixOmics(X = stemcells$gene, Y = unmap(stemcells$celltype), ncomp = 2)

# directed towards mint.PLS
out = mixOmics(X = stemcells$gene, Y = unmap(stemcells$celltype),
ncomp = 2, study = stemcells$study)

# directed towards mint.sPLS because of keepX and/or keepY
out = mixOmics(X = stemcells$gene, Y = unmap(stemcells$celltype),
ncomp = 2, study = stemcells$study, keepX = c(10, 5, 15))

# directed towards mint.PLS-DA because Y is a factor
out = mixOmics(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2,
study = stemcells$study)

# directed towards mint.sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(X = stemcells$gene, Y = stemcells$celltype, ncomp = 2,
study = stemcells$study, keepX = c(10, 5, 15))

## -- Use SummarizedExperiment as an input
# ----------------------------------------------------

data(liver.toxicity)
X = liver.toxicity[["gene"]]
Y = liver.toxicity[["clinic"]]
Y.factor = as.factor(liver.toxicity[["treatment"]][, 4])

# Create SummarizedExperiment object
library(SummarizedExperiment)
assay <- t(X)
coldata <- DataFrame(Y)
coldata[["outcome"]] <- Y.factor
se <- SummarizedExperiment(
    assays = SimpleList(abundance = assay),
    colData = coldata
)

# directed towards PLS
columns <- colnames(colData(se))
columns <- columns[ !columns %in% c("outcome") ]
out = mixOmics(se, assay.type = "abundance", col.var = columns, ncomp = 2)

# directed towards sPLS because of keepX and/or keepY
out = mixOmics(
    se, assay.type = "abundance", col.var = columns,
    ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))

# directed towards PLS-DA because Y is a factor
out = mixOmics(se, assay.type = "abundance", col.var = "outcome", ncomp = 2)

# directed towards sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(
    se, assay.type = "abundance", col.var = "outcome",
    ncomp = 2, keepX = c(20, 20))

## -- Use MultiAssayExperiment as an input
# ----------------------------------------------------

# Create MultiAssayExperiment
data(nutrimouse)
Y = unmap(nutrimouse[["diet"]])
se1 <- SummarizedExperiment(
    assays = SimpleList(abundance = t(nutrimouse[["gene"]])))
se2 <- SummarizedExperiment(
    assays = SimpleList(abundance = t(nutrimouse[["lipid"]])))
coldata <- DataFrame(Y)
colnames(coldata) <- levels(Y)
coldata[["diet"]] <- nutrimouse[["diet"]]
rownames(coldata) <- colnames(se1)
mae <- MultiAssayExperiment(
    experiments = ExperimentList(gene = se1, lipid = se2),
    colData = coldata
)

## -- directed towards block.pls framework because X is a MultiAssayExperiment
# ----------------------------------------------------

# directed towards block PLS
columns <- colnames(colData(mae))
columns <- columns[ !columns %in% c("diet") ]
out = mixOmics(
    mae, experiments = c("gene", "lipid"),
    assay.type = c("abundance", "abundance"),
    col.var = columns, ncomp = 3)

# directed towards block sPLS because of keepX and/or keepY
out = mixOmics(
    mae, experiments = c("gene", "lipid"),
    assay.type = c("abundance", "abundance"),
    col.var = columns, ncomp = 3,
    keepX = list(gene = c(10,10), lipid = c(15,15))
)

# directed towards block PLS-DA because Y is a factor
out <- mixOmics(
    mae, experiments = c("gene", "lipid"),
    assay.type = c("abundance", "abundance"),
    col.var = "diet", ncomp = 3)

# directed towards block sPLS-DA because Y is a factor and there is a keepX
out <- mixOmics(
    mae, experiments = c("gene", "lipid"),
    assay.type = c("abundance", "abundance"),
    col.var = "diet", ncomp = 3,
    keepX = list(gene = c(10,10), lipid = c(15,15))
)

# Create MultiAssayExperiment
data(nutrimouse)
assay <- t(stemcells[["gene"]])
coldata <- DataFrame(unmap(stemcells[["celltype"]]))
colnames(coldata) <- levels(unmap(stemcells[["celltype"]]))
rownames(coldata) <- colnames(assay)
se <- SummarizedExperiment(
    assays = SimpleList(abundance = assay),
    colData = coldata
)
se[["celltype"]] <- stemcells[["celltype"]]
# Split into MultiAssayExperiment
if( !require("mia") ){
    BiocManager::install("mia")
}
se_list <- splitOn(se, by = "columns", f = stemcells[["study"]])
names(se_list) <- paste0("study_", names(se_list))
mae <- MultiAssayExperiment(experiments = se_list)

## -- directed towards mint.pls framework because of the specified MINT argument
# ----------------------------------------------------

# Get columns from the first experiment/study
columns <- colnames(colData(mae[[1]]))
columns <- columns[ !columns %in% c("celltype") ]

# directed towards PLS (note that we are using SummarizedExperiment)
out = mixOmics(se, assay.type = "abundance", col.var = columns, ncomp = 2)

# directed towards mint.PLS
out = mixOmics(
    mae, experiments = c("study_1", "study_3"),
    assay.type = c("abundance", "abundance"),
    col.var = columns, ncomp = 2, MINT = TRUE)

# directed towards mint.sPLS because of keepX and/or keepY
out = mixOmics(
    mae, experiments = c("study_1", "study_3"),
    assay.type = c("abundance", "abundance"),
    col.var = columns, ncomp = 3, MINT = TRUE,
    keepX = c(10, 5, 15)
)

# directed towards mint.PLS-DA because we specify class
out = mixOmics(
    mae, experiments = c("study_1", "study_3"),
    assay.type = c("abundance", "abundance"),
    col.var = "celltype", ncomp = 2, MINT = TRUE
)

# directed towards mint.sPLS-DA because Y is a factor and there is a keepX
out = mixOmics(
    mae, experiments = c("study_1", "study_3"),
    assay.type = c("abundance", "abundance"),
    col.var = "celltype", ncomp = 3, MINT = TRUE, keepX = c(10, 5, 15)
)

}
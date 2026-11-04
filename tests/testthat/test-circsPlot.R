test_that("circosPlot works", code = {
    data(nutrimouse)
    Y = nutrimouse$diet
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
    
    # wrapper.sgccda is just exported version of block.splsda()
    # nutrimouse.sgccda <- wrapper.sgccda(X=data,
    #                                     Y = Y,
    #                                     design = design,
    #                                     keepX = list(gene=c(8,8), lipid=c(4,4)),
    #                                     ncomp = 2,
    #                                     scheme = "horst")
    nutrimouse.sgccda <- block.splsda(X=data, Y = Y, design = design,
                                        keepX = list(gene=c(8,8), lipid=c(4,4)),
                                        ncomp = 2)
    
    
    cp_res <- circosPlot(nutrimouse.sgccda, cutoff = 0.7, ncol.legend = 2, size.legend = 1.1,
                        color.Y = 1:5, color.blocks = c("green","brown"), color.cor = c("magenta", "purple"))
    expect_is(cp_res, "matrix")
    
})

test_that("circosPlot works with similar feature names in different blocks", code = {

    create_similar_feature_names <- function(data_list)
    {
        lapply(data_list, function(x){
            colnames(x) <- paste0('feature_', seq_len(ncol(x)))
            x
        })
    }
    
    
    
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, 
                mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    
    data <- create_similar_feature_names(data)
    list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2), protein = rep(10, 2))
    TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype, 
                                     ncomp = 2, keepX = list.keepX, design = 'full')
    cp_res <- .quiet(circosPlot(TCGA.block.splsda, cutoff = 0.7))
    expect_is(cp_res, "matrix")
})

test_that("circosPlot works when using the indY parameter", code = {
    
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna, 
                mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)
    
    list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2), protein = rep(10, 2))
    TCGA.block.spls = block.spls(X = data, indY = 3, 
                                     ncomp = 2, keepX = list.keepX, design = 'full')
    cp_res <- circosPlot(TCGA.block.spls, cutoff = 0.7, group = breast.TCGA$data.train$subtype)
    
    expect_is(cp_res, "matrix")
})


test_that("circosPlot works with comp = 1 when a block has a single selected variable", code = {
    # Check that selecting a single variable preserves its name and avoids indexing errors.
    data("breast.TCGA")
    data = list(mrna = breast.TCGA$data.train$mrna,
                mirna = breast.TCGA$data.train$mirna,
                protein = breast.TCGA$data.train$protein)

    list.keepX = list(mrna = c(10, 10), mirna = c(1, 10), protein = c(10, 10))
    TCGA.block.splsda = block.splsda(X = data, Y = breast.TCGA$data.train$subtype,
                                     ncomp = 2, keepX = list.keepX, design = 'full')

    cp_res <- circosPlot(TCGA.block.splsda, comp = 1, cutoff = 0.5)

    expect_is(cp_res, "matrix")
    # one row/column per variable selected on comp 1 across all blocks
    expect_equal(ncol(cp_res), sum(sapply(list.keepX, `[`, 1)))
    # preserve every selected feature name in block and input-variable order
    expected.names <- unlist(lapply(TCGA.block.splsda$loadings[names(data)],
                                   function(x) rownames(x)[x[, 1] != 0]),
                             use.names = FALSE)
    expect_identical(rownames(cp_res), expected.names)
    expect_identical(colnames(cp_res), expected.names)
})

test_that("circosPlot keeps the links of every block listed in blocks.link", code = {
    # three blocks of two features each, all pairs similar
    chr <- data.frame(chrom = rep(c("chra", "chrb", "chrc"), each = 2),
                      chromStart = rep(c(0, 1), 3), chromEnd = rep(c(1, 2), 3),
                      name = paste0("f", 1:6), color = "black")
    simMat <- matrix(0.9, 6, 6, dimnames = list(chr$name, chr$name))
    
    all.links <- genLinks(chr, simMat, threshold = 0.5)
    abc.links <- genLinks(chr, simMat, threshold = 0.5, blocks.link = c("a", "b", "c"))
    expect_equal(nrow(abc.links), nrow(all.links))
    
    ab.links <- genLinks(chr, simMat, threshold = 0.5, blocks.link = c("a", "b"))
    expect_setequal(as.character(ab.links$feat2), paste0("f", 1:4))
})

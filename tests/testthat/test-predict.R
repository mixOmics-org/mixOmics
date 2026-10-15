test_that("predict.mixo_pls works", code = {
    data("linnerud")
    X <- linnerud$exercise
    Y <- linnerud$physiological
    linn.pls <- pls(X, Y, ncomp = 2, mode = "classic")
    
    indiv1 <- as.integer(colMeans(X)) + 3
    indiv2 <- as.integer(colMeans(X)) - 3
    newdata <- rbind(indiv1, indiv2)
    colnames(newdata) <- colnames(X)
    
    # yhat_approx <- colMeans(Y)
    
    pred <- predict(linn.pls, newdata)
    expect_equal(round(pred$predict[1]), round(175.8865 ))
    # pred$predict
    # yhat_approx
    
})


test_that("predict.mixo_plsda works", code = {
    data("liver.toxicity")
    X <- liver.toxicity$gene
    Y <- as.factor(liver.toxicity$treatment[, 4])
    
    
    ## if training is perfomed on 4/5th of the original data
    samp <- sample(1:5, nrow(X), replace = TRUE)
    test <- which(samp == 1)   # testing on the first fold
    train <- setdiff(1:nrow(X), test)
    
    plsda.train <- plsda(X[train, ], Y[train], ncomp = 2)
    test.predict <- predict(plsda.train, X[test, ], dist = "max.dist")
    expect_is(test.predict, "predict")
})

test_that("predict.block.splsda", code = {
    # example with block.splsda=diablo=sgccda and a missing block
    data(nutrimouse)
    # need to unmap Y for an unsupervised analysis, where Y is included as a data block in data
    Y.mat = unmap(nutrimouse$diet)
    data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y.mat)
    # with this design, all blocks are connected
    design = matrix(c(0,1,1,1,0,1,1,1,0), ncol = 3, nrow = 3,
                    byrow = TRUE, dimnames = list(names(data), names(data)))
    
    # train on 75
    ind.train=NULL
    for(i in 1:nlevels(nutrimouse$diet))
        ind.train=c(ind.train,which(nutrimouse$diet==levels(nutrimouse$diet)[i])[1:6])
    
    #training set
    gene.train=nutrimouse$gene[ind.train,]
    lipid.train=nutrimouse$lipid[ind.train,]
    Y.mat.train=Y.mat[ind.train,]
    Y.train=nutrimouse$diet[ind.train]
    data.train=list(gene=gene.train,lipid=lipid.train,Y=Y.mat.train)
    
    #test set
    gene.test=nutrimouse$gene[-ind.train,]
    lipid.test=nutrimouse$lipid[-ind.train,]
    Y.mat.test=Y.mat[-ind.train,]
    Y.test=nutrimouse$diet[-ind.train]
    data.test=list(gene=gene.test,lipid=lipid.test)
    
    # example with block.splsda=diablo=sgccda and a missing block
    res.train = block.splsda(X=list(gene=gene.train,lipid=lipid.train),Y=Y.train,
                             ncomp=3,keepX=list(gene=c(10,10,10),lipid=c(5,5,5)))
    
    expect_warning((test.predict <- predict(res.train, newdata=data.test[2], method = "max.dist")))
    expect_is(test.predict$WeightedPredict, 'array')
})


test_that("predict.mint.splsda works", code = {
    ## example with mint.splsda
    data(stemcells)
    
    #training set
    ind.test = which(stemcells$study == "3")
    gene.train = stemcells$gene[-ind.test,]
    Y.train = stemcells$celltype[-ind.test]
    study.train = factor(stemcells$study[-ind.test])
    
    #test set
    gene.test = stemcells$gene[ind.test,]
    Y.test = stemcells$celltype[ind.test]
    study.test = factor(stemcells$study[ind.test])
    
    res = mint.splsda(X = gene.train, Y = Y.train, ncomp = 3, keepX = c(10, 5, 15),
                      study = study.train)
    
    pred = predict(res, newdata = gene.test, study.test = study.test)
    # saveRDS(pred, file = 'inst/testdata/predict.mint.splsda.rda')
    pred_ref <- readRDS(system.file("testdata", "predict.mint.splsda.rda", package = 'mixOmics'))
    expect_equal(pred_ref, pred)
})

test_that("(predict:error): catches when colnames(newdata) differs from colnames(X)", {
    
    data(breast.TCGA) # load in the data

    # extract data
    X.train = list(mirna = breast.TCGA$data.train$mirna,
                   mrna = breast.TCGA$data.train$mrna)
    
    X.test = list(mirna = breast.TCGA$data.test$mirna,
                  mrna = breast.TCGA$data.test$mrna)
    
    Y.train = breast.TCGA$data.train$subtype
    
    # use optimal values from the case study on mixOmics.org
    optimal.ncomp = 2
    optimal.keepX = list(mirna = c(10,5),
                         mrna = c(26, 16))
    
    # set design matrix
    design = matrix(0.1, ncol = length(X.train), nrow = length(X.train),
                    dimnames = list(names(X.train), names(X.train)))
    diag(design) = 0
    
    # generate model
    final.diablo.model = block.splsda(X = X.train, Y = Y.train, ncomp = optimal.ncomp, # set the optimised DIABLO model
                                      keepX = optimal.keepX, design = design)
    
    
    # create new test data with one dataframe being reordered
    new.var.order = sample(1:dim(X.test$mrna)[2])
    X.test.reorder <- X.test
    X.test.reorder$mrna <- X.test.reorder$mrna[, new.var.order]
    
    
    X.test.new.feat <- X.test
    colnames(X.test.new.feat$mrna)[1] <- "random.feature.name"
    
    # ---------------------------------------------------------------------------- #
    
    # should raise error about mismatching ORDERS of features
    expect_error(predict(final.diablo.model, newdata = X.test.reorder),
                 "The order of features in 'object$X$mrna' is different to 'newdata$mrna'. Please ensure that you adjust these to the same order",
                 fixed = T)

    # should raise error about mismatching SETS of features
    expect_error(predict(final.diablo.model, newdata = X.test.new.feat),
                 "The set of features in 'object$X$mrna' is different to 'newdata$mrna'. Please ensure they have the same features.",
                 fixed = T)
})


# test_that("predict.block.splsda works on reordered test data", code = {
#     data(breast.TCGA) # load in the data
#     
#     X.train = list(mirna = breast.TCGA$data.train$mirna,
#                    mrna = breast.TCGA$data.train$mrna)
#     
#     X.test = list(mirna = breast.TCGA$data.test$mirna,
#                   mrna = breast.TCGA$data.test$mrna)
#     
#     Y.train = breast.TCGA$data.train$subtype
#     Y.test = breast.TCGA$data.test$subtype
#     
#     optimal.ncomp = 2
#     optimal.keepX = list(mirna = c(10,5),
#                          mrna = c(26, 16))
#     
#     design = matrix(0.1, ncol = length(X.train), nrow = length(X.train), 
#                     dimnames = list(names(X.train), names(X.train)))
#     diag(design) = 0 
#     
#     final.diablo.model = block.splsda(X = X.train, Y = Y.train, ncomp = optimal.ncomp, # set the optimised DIABLO model
#                                       keepX = optimal.keepX, design = design)
#     
#     
#     # create new test data with one dataframe being reordered
#     new.var.order = sample(1:dim(X.test$mirna)[2])
#     
#     X.test.dup <- X.test
#     X.test.dup$mirna <- X.test.dup$mirna[, new.var.order]
#     
#     predict.diablo = predict(final.diablo.model, newdata = X.test)
#     predict.diablo.reordered = predict(final.diablo.model, newdata = X.test.dup)
#     
#     homogenity <- matrix(NA, nrow = 2, ncol = 3)
#     colnames(homogenity) <- c("max.dist", "centroids.dist", "mahalanobis.dist")
#     rownames(homogenity) <- c("mirna", "mrna")
#     
#     for (dist in colnames(homogenity)) {
#         for (block in rownames(homogenity)) {
#             homogenity[block, dist] = all(predict.diablo$class[[dist]][[block]] == 
#                                               predict.diablo.reordered$class[[dist]][[block]])
#         }
#     }
#     
#     expect_true(all(homogenity))
# })

test_that("predict.block.(s)pls(da) works", code = {
    data("nutrimouse")
    X <- list(gene1 = nutrimouse$gene, gene2 = nutrimouse$gene + 0.5, lipid = nutrimouse$lipid)
    indY <- 3
    labels <- nutrimouse$diet
    nutri.pls <- block.pls(X, indY = indY, ncomp = 2, design = "full", scale = TRUE)
    nutri.spls <- block.spls(X, indY = indY, ncomp = 2, design = "full", scale = TRUE, keepX = list(gene1 = c(5, 5), gene2 = c(5, 5)))
    nutri.plsda <- block.plsda(X, Y = labels, ncomp = 2, design = "full", scale = TRUE)
    nutri.splsda <- block.splsda(X, Y = labels, ncomp = 2, design = "full", scale = TRUE, keepX = list(gene1 = c(5, 5), gene2 = c(5, 5), lipid = c(5, 5)))
    
    newdata <- lapply(X, function (data) {
        indiv1 <- colMeans(data) + 0.3
        indiv2 <- colMeans(data) - 0.3
        ret <- rbind(indiv1, indiv2)
        colnames(ret) <- colnames(data)
        ret
    })
    
    pred <- list(pls = predict(nutri.pls, newdata[-indY])$WeightedPredict,
                 spls = predict(nutri.spls, newdata[-indY])$WeightedPredict,
                 plsda = predict(nutri.plsda, newdata)$WeightedPredict,
                 splsda = predict(nutri.splsda, newdata)$WeightedPredict)
    # saveRDS(pred, file = "inst/testdata/predict.block.pls.rda")
    pred_ref <- readRDS(system.file("testdata", "predict.block.pls.rda", package = 'mixOmics'))
    expect_equal(pred, pred_ref)
})


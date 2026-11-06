################################################################################
# Author :
#   Florian Rohart,
#
# created: 22-04-2015
# last modified: 04-10-2017
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################


# perform the mint.pls on a subset of variables on one only dimension,
# deflate the intial matrices X and Y (already center by study)

# mean centering with attach and without modify.na, need to look at how to
# remove some of means/sigma when nearZerVar is used
# we can have a list of studies for Discriminant Analyses,
# not for pls/spls as they would be overlapping batch effects

# ==============================================================================
# internal_wrapper.mint.block: this function is a particular setting of
#   internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint.block
# ==============================================================================
# used in (mint).block approaches

internal_wrapper.mint.block <-
    function(X,
             Y,
             indY,
             study,
             ncomp,
             keepX,
             keepY,
             DA=NULL,
             test.keepX=NULL,
             test.keepY=NULL,
             design,
             scheme,
             mode,
             scale = TRUE,
             init ,
             tol = 1e-06,
             max.iter = 100,
             near.zero.var = FALSE,
             misdata = NULL, is.na.A = NULL, ind.NA = NULL, ind.NA.col = NULL,
             all.outputs=TRUE
    )
    {
        if (missing(scheme))
            scheme= "horst"
        
        if (missing(mode))
            mode="regression"
        
        
        # checks (near.zero.var is done there)
        check <- Check.entry.wrapper.mint.block(X = X, Y = Y, indY = indY,
                                             ncomp = ncomp, keepX = keepX, keepY = keepY, DA = DA,
                                             study = study, design = design, init = init, scheme = scheme, scale = scale,
                                             near.zero.var = near.zero.var, mode = mode, tol = tol,
                                             max.iter = max.iter)
        
        # get some values after checks
        A = check$A
        indY = check$indY
        study = check$study
        design = check$design
        ncomp = check$ncomp
        keepA = check$keepA
        keepA.save = keepA
        init = check$init
        nzv.A = check$nzv.A
        
        #--------------------------------------------------------------------------#
        #-- keepA ----------------------------------------------------#
        
        # shaping keepA, will need to be done somewhere before eventually
        
        if(!is.null(test.keepX) & !is.null(test.keepY))
        {
            test.keepA = lapply(c(test.keepX, Y=test.keepY),sort)
            #sort test.keepX so as to be sure to chose the smallest in case of
            # several minimum
        } else {test.keepA=NULL}
        
        keepAA = vector("list", length = max(ncomp)) # one keepA per comp
        names(keepAA) = paste0("comp",1:max(ncomp))
        for(comp in 1:max(ncomp)) # keepA[[block]] [1:ncomp]
            keepAA[[comp]] = lapply(keepA, function(x) x[comp])
        
        if(!is.null(test.keepA))
            keepAA[[max(ncomp)]] = test.keepA
        
        keepA = lapply(keepAA, expand.grid)
        
        #print(keepA)
        # keepA[[comp]] is a matrix where each row is all the keepX the test over
        # the block (each block is a column)
        
        #-- keepA ----------------------------------------------------#
        #--------------------------------------------------------------------------#
        
        
        # A: list of matrices
        # indY: integer, pointer to one of the matrices of A
        # design: design matrix, links between matrices. Diagonal must be 0
        # ncomp: vector of ncomp, per matrix
        # scheme: a function "g", refer to the article (thanks Benoit)
        # scale: do you want to scale ? mean is done by default
        # init: one of "svd" or "random", initialisation of the algorithm
        # tol: nobody cares about this
        # mode: canonical, classic, invariant, regression
        # max.iter: nobody cares about this
        # study: factor for each matrix of A, must be a vector
        # keepA: keepX of spls for each matrix of A. must be a list.
        #   Each entry must be of the same length (max ncomp)
        # near.zero.var: do you want to remove variables with very small variance
        
        result <- internal_mint.block(A = A,
                                   indY = indY,
                                   design = design,
                                   ncomp = ncomp,
                                   scheme = scheme,
                                   scale = scale,
                                   init = init,
                                   tol = tol,
                                   tau = NULL,
                                   mode = mode,
                                   max.iter = max.iter,
                                   study = study,
                                   keepA = keepA,
                                   misdata = misdata, is.na.A = is.na.A, ind.NA = ind.NA,
                                   ind.NA.col = ind.NA.col, all.outputs= all.outputs)
        
        if(near.zero.var)
            result$nzv=nzv.A
        
        result$keepX = keepA.save
        
        class(result) = c("sparse.mint.block")
        return(invisible(result))
        
        
    }




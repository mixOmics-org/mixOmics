################################################################################
# Author :
#   Florian Rohart,
#
# created: 22-04-2015
# last modified: 05-10-2017
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################


# ==============================================================================
# internal_wrapper.mint: perform a vertical PLS on a combination of experiments,
#   input as a matrix in X
# this function is a particular setting of internal_mint.block, the formatting
#   of the input is checked in Check.entry.pls
# internal function. Do not export in NAMESPACE.
# ==============================================================================
# used in (mint).(s)pls(da)

# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# study: grouping factor indicating which samples are from the same study
# keepX: number of \eqn{X} variables kept in the model on the last components.
# keepY: number of \eqn{Y} variables kept in the model on the last components.
# mode: input mode, one of "canonical", "classic", "invariant" or "regression".
#   Default to "regression"
# scale: Logical. If scale = TRUE, each block is standardized to zero means and
#   unit variances (default: TRUE).
# near.zero.var: Logical, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# max.iter: integer, the maximum number of iterations.
# tol: Convergence stopping value.
# logratio: one of "none", "CLR"
# DA: indicate whether it's a DA analysis, only used for the multilvel approach
#   with withinVariation
# multilevel: multilevel is passed to multilevel(design=) in withinVariation.
#   Y is omitted and should be included in multilevel design

internal_wrapper.mint <- 
  function(X,
           Y,
           study,
           ncomp = 2,
           keepX,
           keepY,
           test.keepX=NULL,
           test.keepY=NULL,
           mode,
           scale = FALSE,
           near.zero.var = FALSE,
           max.iter = 100,
           tol = 1e-06,
           logratio = "none", 
           DA = FALSE,
           multilevel = NULL,
           misdata = NULL, is.na.A = NULL, ind.NA = NULL, ind.NA.col = NULL,
           all.outputs=FALSE,
           remove.object=NULL
  )
  {
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || length(ncomp)>1)
      stop("invalid number of variates, 'ncomp'.")
    
    #-- validation des arguments --#
    
    check = Check.entry.pls(X, Y, ncomp, keepX, keepY, mode=mode, scale=scale,
                            near.zero.var=near.zero.var, max.iter=max.iter ,tol=tol ,logratio=logratio,
                            DA=DA, multilevel=multilevel)
    X = check$X
    input.X = X # save the checked X, before logratio/multileve/scale
    Y = check$Y
    ncomp = check$ncomp
    mode = check$mode
    keepX = check$keepX
    keepY = check$keepY
    nzv.A = check$nzv.A
    multilevel <- check$multilevel
    rm(check) # free memory
    #remove `X' from the previous environment
    if(!is.null(remove.object))
      rm(list=remove.object, envir=parent.frame()) # free memory
    
    #test.keepX and test.keepY must be checked before (in tune)
    
    #set the default study factor
    if (missing(study))
    {
      study = factor(rep(1,nrow(X)))
    } else {
      study = as.factor(study)
      #if(nlevels(study) == 1)
      #stop("'study' has a single level, no need to use the MINT approach")
    }
    if (length(study) != nrow(X))
      stop(paste0("'study' must be a factor of length ",nrow(X),"."))
    
    if (any(table(study) <= 1))
      stop("At least one study has only one sample, please consider removing
    before calling the function again")
    if (any(table(study) < 5))
      warning("At least one study has less than 5 samples, mean centering might
    not do as expected")
    
    design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
    
    
    #-----------------------------#
    #-- logratio transformation --#
    
    X = logratio.transfo(X=X, logratio=logratio)
    
    #as X may have changed
    if (ncomp > min(ncol(X), nrow(X)))
      stop("'ncomp' should be smaller than ", min(ncol(X), nrow(X)),
           call. = FALSE)
    
    #-- logratio transformation --#
    #-----------------------------#
    
    
    #--------------------------------------------------------------------------#
    #-- multilevel approach ---------------------------------------------------#
    
    if (!is.null(multilevel))
    {
      if (!DA)
      {
        Xw = withinVariation(X, design = multilevel)
        Yw = withinVariation(Y, design = multilevel)
        X = Xw
        Y = Yw
      } else {
        Xw = withinVariation(X, design = multilevel)
        X = Xw
        
        #-- Need to set Y variable for 1 or 2 factors
        Y = multilevel[, -1,drop=FALSE]
        if (ncol(Y)>0)
          Y = apply(Y, 1, paste, collapse = ".")
        #  paste is to combine in the case we have 2 levels
        
        Y = as.factor(Y)
        Y.factor = Y
        Y = unmap(Y)
        colnames(Y) = levels(Y)
        rownames(Y) = rownames(X)
        # if DA keepY should be all the levels
        # (which is not happening in the check because of multilevel
        keepY = rep(ncol(Y),ncomp)
      }
    }
    #-- multilevel approach ---------------------------------------------------#
    #--------------------------------------------------------------------------#
    
    
    #--------------------------------------------------------------------------#
    #-- keepA ----------------------------------------------------#
    
    # shaping keepA, contains all the keepX/keepY models to be constructed
    
    if(!is.null(test.keepX) & !is.null(test.keepY))
    {
      test.keepA = lapply(list(X=test.keepX, Y=test.keepY),sort)
      #sort test.keepX so as to be sure to chose the smallest in case of
      # several minimum
    } else {test.keepA=NULL}
    
    keepA = vector("list", length = ncomp) # one keepA per comp
    names(keepA) = paste0("comp",1:ncomp)
    for(comp in 1:length(keepX)) # keepA[[block]] [1:ncomp]
      keepA[[comp]] = lapply(list(X=keepX, Y=keepY), function(x) x[comp])
    
    if(!is.null(test.keepA))
      keepA[[ncomp]] = test.keepA
    
    keepA = lapply(keepA, expand.grid)
    
    # keepA[[comp]] is a matrix where each row is all the keepX the test over
    # the block (each block is a column)
    
    #-- keepA ----------------------------------------------------#
    #--------------------------------------------------------------------------#
    
    
    #--------------------------------------------------------------------------#
    #-- pls approach ----------------------------------------------------#
    result = internal_mint.block(A = list(X = X, Y = Y), indY = 2, mode = mode,
                                 ncomp = c(ncomp, ncomp), tol = tol, max.iter = max.iter,
                                 design = design, keepA = keepA, scale = scale, scheme = "horst",init="svd",
                                 study = study, misdata = misdata, is.na.A = is.na.A, ind.NA = ind.NA,
                                 ind.NA.col = ind.NA.col, all.outputs= all.outputs, remove.object=c("X"))
    
    #-- pls approach ----------------------------------------------------#
    #--------------------------------------------------------------------------#
    
    # result contains all loadings and variates of the test.keepX and
    # test.keepY (if not null)
    # if no test.keepX and test.keepY, then it's classical outputs
    
    result$keepX = keepX
    result$keepY = keepY
    result$ncomp = ncomp
    if(near.zero.var)
      result$nzv = nzv.A
    
    if(!is.null(multilevel) & DA)
      result$Y.factor = Y.factor
    
    result$input.X = input.X
    
    class(result) = c("mint.spls.hybrid")
    return(invisible(result))
    
  }




################################################################################
# Authors:
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

# --------------------------------------
# check and construct keepA
# --------------------------------------
get.keepA <- function(X,keepX,ncomp)
{
    # X:data
    # keepA
    
    keepA = list()
    if (missing(keepX) || length(keepX)==0)
    {
        #if keepX is missing, pls-like: keepX=ncol(X)
        for (q in 1:length(X))
            keepA[[q]]=rep(ncol(X[[q]]),max(ncomp)) #keepX
        
        names(keepA)=names(X)
    } else {
        
        if (!is.list(keepX))
            stop("'keepX' must be a list")
        
        if (length(keepX) > length(X))
            stop(paste0("length(keepX) is higher than the number of blocks in X,
        which is ", length(X), "."))
        
        # error if no names on keepX or not matching the names in X
        if(length(unique(names(keepX)))!=length(keepX) |
           sum(is.na(match(names(keepX),names(X)))) > 0)
            stop("Each entry of 'keepX' must have a unique name corresponding to
        a block of 'X'")
        
        # I want to match keepX to X by names
        ind.match = match(names(X), names(keepX))
        
        for (q in 1:length(X))
        {
            
            if (!is.na(ind.match[q]))
                # means there is a keepX with the same name as X[q]
                #(q <= length(keepX))
            {
                #checking entries of keepX
                if (is.list(keepX[[ind.match[q]]]))
                    stop(paste0("keepX[[",ind.match[q],"]]' must be a vector"))
                
                if (any(keepX[[ind.match[q]]] > ncol(X[[q]])))
                    stop(paste0("each component of 'keepX[[",ind.match[q],"]]'
                must be lower or equal to ncol(X[[",q,"]])=",ncol(X[[q]]),"."))
                
                if (any(keepX[[ind.match[q]]] < 0))
                    stop(paste0("each component of 'keepX[[",ind.match[q],"]]'
                must be non negative."))
                
                if (length(keepX[[ind.match[q]]]) > ncomp[q])
                    stop(paste0("length of 'keepX[[",ind.match[q],"]]'
                must be lower or equal to ncomp[",q,"]=",ncomp[q], "."))
                
                keepA[[q]] = keepX[[ind.match[q]]]
                if (length(keepA[[q]]) < max(ncomp))
                    keepA[[q]] = c(keepA[[q]], rep(ncol(X[[q]]),
                                                   max(ncomp) - length(keepA[[q]])))
                #complete the keepX already provided
                
            }else{
                keepA[[q]]=rep(ncol(X[[q]]),max(ncomp))
            }
        }
        
    }
    
    #print(keepA)
    
    names(keepA) = names(X)
    return(list(keepA=keepA))
}


# --------------------------------------
# Check.entry.single: check input parameter for a single matrix that comes
#   from a list of matrix
# --------------------------------------
# X: a single numeric matrix
# ncomp: the number of components to include in the model
# q: position of X in a list of matrix, as used in horizontal analysis
#   (e.g. block.pls)

Check.entry.single = function(X,  ncomp, q)
{
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
        stop(paste0("'X[[", q, "]]' must be a numeric matrix."))
    
    if(! any(class(X) %in% "matrix"))
        X = as.matrix(X)
    
    if (!is.numeric(X))
        stop(paste0("'X[[", q, "]]'  must be a numeric matrix."))
    
    N = nrow(X)
    P = ncol(X)
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
        stop(paste0(
            "invalid number of variates 'ncomp' for matrix 'X[[", q, "]]'."))
    
    ncomp = round(ncomp)
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = 1:N
        rownames(X)  = ind.names
    }
    
    if (length(unique(rownames(X))) != nrow(X))
        stop("samples should have a unique identifier/rowname")
    if (length(unique(X.names)) != P)
        stop("Unique indentifier is needed for the columns of X")
    
    return(list(X=X, ncomp=ncomp, X.names=X.names, ind.names=ind.names))
}



# --------------------------------------
# Check.entry.pls: check the entry associated with (s)PLS.
#   Used in 'internal_wrapper.mint.R'
# --------------------------------------
# X: numeric matrix of predictors
# Y: numeric vector or matrix of responses
# ncomp: the number of components to include in the model. Default to 2.
# keepX: number of \eqn{X} variables kept in the model on the last components.
# keepY: number of \eqn{Y} variables kept in the model on the last components.
# mode: deflation mode
# scale: Logical. If scale = TRUE, each block is standardized to zero means and
#   unit variances (default: TRUE).
# near.zero.var: Logical, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# max.iter: integer, the maximum number of iterations.
# tol: Convergence stopping value.
# logratio: whether a log transformation will be performed, and which one,
# DA: whether the input is used for a Discrimant Analysis
# multilevel: whether a multilevel analysis will be performed
#   (repeated measurements)


Check.entry.pls = function(X, Y, ncomp, keepX, keepY, test.keepX, test.keepY,
                           mode, scale, near.zero.var, max.iter, tol, logratio, DA, multilevel)
{
    
    if (missing(mode))
        mode = "regression"
    
    if (length(mode)>1)
        mode = mode[1]
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
        stop("Choose one of the four following modes: canonical, invariant,
    classic or regression")
    
    
    #-- validation des arguments --#
    if (length(dim(X)) != 2)
        stop("'X' must be a numeric matrix.")
    
    if(! any(class(X) %in% "matrix"))
        X = as.matrix(X)
    
    if (!(logratio %in% c("none", "CLR")))
        stop("Choose one of the two following logratio transformation: none or CLR")
    
    # if DA and the unmapped Y has rows without associated class, only run this check if Y is not null (ie multilevel)
    if (DA && !is.null(Y)) {
      if (length(which(rowSums(Y)==0)) != 0) {
        stop("Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector")
      }
    }
    
    if(!is.null(multilevel))
    {
        #multilevel analysis: withinVariation and then pls-like
        # if it's DA analysis, Y and 'multilevel' are combined
        if(DA)
        {
            Y = multilevel
        }else{
            multilevel = data.frame(multilevel)
            if ((nrow(X) != nrow(multilevel)))
                stop("unequal number of rows in 'X' and 'multilevel'.")
            
            Y = as.matrix(Y)
            if (!is.numeric(X) || !is.numeric(Y))
                stop("'X' and/or 'Y' must be a numeric matrix.")
        }
    }else{
        Y = as.matrix(Y)
        if (!is.numeric(X) || !is.numeric(Y))
            stop("'X' and/or 'Y' must be a numeric matrix.")
    }
    N = nrow(X)
    Q = ncol(Y)
    P= ncol(X)
    
    if ((N != nrow(Y)))
        stop("Unequal number of rows in 'X' and 'Y'.")
    
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0 || length(ncomp)>1)
        stop("invalid number of variates, 'ncomp'.")
    
    if(mode == "canonical" & ncomp>ncol(Y))
        stop("For `canonical mode', 'ncomp' needs to be lower than ncol(Y)= ",
             ncol(Y))
    
    
    ncomp = round(ncomp)
    if(ncomp > P)
    {
        warning("Reset maximum number of variates 'ncomp' to ncol(X) = ", P,
                ".")
        ncomp = P
    }
    
    if (!is.numeric(tol) | tol<=0)
        stop("tol must be non negative")
    
    if (!is.numeric(max.iter) | max.iter<=0)
        stop("max.iter must be non negative")
    
    
    # add colnames and rownames if missing
    X.names = dimnames(X)[[2]]
    if (is.null(X.names))
    {
        X.names = paste("X", 1:P, sep = "")
        dimnames(X)[[2]] = X.names
    }
    
    
    
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names))
    {
        ind.names = dimnames(Y)[[1]]
        if (is.null(ind.names))
        {
            ind.names = 1:N
            rownames(X) = rownames(Y) = ind.names
            
        } else {
            rownames(X) = ind.names
        }
        
        #rownames(X) = ind.names
    } else {
        rownames(Y) = ind.names
    }
    
    #rownames(X) = rownames(Y) = ind.names
    
    
    #if (dim(Y)[2] == 1) Y.names = "Y"
    Y.names = dimnames(Y)[[2]]
    if (is.null(Y.names))
    {
        if (dim(Y)[2] == 1)
        {
            Y.names = "Y"
        } else {
            Y.names = paste("Y", 1:Q, sep = "")
        }
        
        dimnames(Y)[[2]]=Y.names
    }
    
    if (length(unique(X.names)) != P)
        stop("Unique indentifier is needed for the columns of X")
    
    if (length(unique(Y.names)) != Q)
        stop("Unique indentifier is needed for the columns of Y")
    
    
    # check keepX
    if (missing(keepX))
    {
        keepX = rep(P, ncomp)
    } else {
        len <- length(keepX)
        if (len<ncomp)
            keepX = c(keepX, rep(P, ncomp - len))
        else if (len>ncomp) {
            cat("keepX is of length", len, "while ncomp is", ncomp, "\n")
            cat("trimming keepX to ")
            cat(.colour.txt(paste0('c(', paste0(keepX[seq_len(ncomp)], collapse = ','), ')'), sQuote = FALSE))
            
        }
           
    }
    
    # check keepY
    if (missing(keepY))
    {
        keepY = rep(Q, ncomp)
    } else {
        if (length(keepY) < ncomp)
            keepY = c(keepY, rep(Q, ncomp - length(keepY)))
        #complete the keepY already provided
    }
    
    
    
    
    if (any(keepX<0))
        stop("each component of 'keepX' must be non negative ")
    if (any(keepY<0))
        stop("each component of 'keepY' must be non negative ")
    
    if (any(keepX > ncol(X)))
        stop("each component of 'keepX' must be no greater than ", P, ".")
    if (any(keepY > ncol(Y)))
        stop("each component of 'keepY' must be no greater than ", Q, ".")
    
    
    if (!is.logical(scale))
        stop("'scale' must be either TRUE or FALSE")
    
    if (!is.logical(near.zero.var))
        stop("'near.zero.var' must be either TRUE or FALSE")
    
    
    ### near.zero.var, remove the variables with very small variances
    if (near.zero.var == TRUE)
    {
        nzv.A = nearZeroVar(X)
        
        if (length(nzv.A$Position) > 0)
        {
            names.remove.X = colnames(X)[nzv.A$Position]
            X = X[, -nzv.A$Position, drop=FALSE]
            warning("Zero- or near-zero variance predictors.\n Reset
            predictors matrix to not near-zero variance predictors.\n
            See $nzv for problematic predictors.")
            if (ncol(X) == 0)
                stop("No more variables in X")
            
            #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
            if (any(keepX > ncol(X)))
            {
                ind = which(keepX > ncol(X))
                keepX[ind] = ncol(X)
            }
        }
        
    }else{nzv.A=NULL}
    
    
    return(list(X=X, Y=Y, ncomp=ncomp, X.names=X.names, Y.names=Y.names,
                ind.names=ind.names, mode=mode, multilevel=multilevel, keepX=keepX, keepY=keepY, nzv.A=nzv.A))
}


# --------------------------------------
# Check.entry.wrapper.mint.block: used in 'internal_wrapper.mint.block.R'
# --------------------------------------
# X: a list of data sets (called 'blocks') matching on the same samples.
#   Data in the list should be arranged in samples x variables, with samples
#   order matching in all data sets. \code{NA}s are not allowed.
# Y: a factor or a class vector for the discrete outcome.
# indY: to supply if Y is missing, indicate the position of the outcome in X.
# ncomp: numeric vector of length the number of blocks in \code{X}.
#   The number of components to include in the model for each block
#   (does not necessarily need to take the same value for each block).
#   By default set to 2 per block.
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of
#   X[[i]]-variables kept in the model on the last components.
# keepY: Only used if Y is provided. Each entry keepY[i] is the number of
#   Y-variables kept in the model on the last components.
# study: grouping factor indicating which samples are from the same study
# design: the input design.
# init: intialisation of the algorithm, one of "svd" or "svd.single".
#   Default to "svd"
# scheme: the input scheme, one of "horst", "factorial" or ""centroid".
#   Default to "centroid"
# scale: Logical. If scale = TRUE, each block is standardized to zero means
#   and unit variances (default: TRUE).
# near.zero.var: Logical, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# mode: input mode, one of "canonical", "classic", "invariant" or "regression".
#   Default to "regression"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.

Check.entry.wrapper.mint.block = function(X,
                                          Y,
                                          indY, #only use if Y not provided
                                          ncomp,
                                          keepX,
                                          keepY,
                                          DA=NULL,
                                          study, #mint
                                          design, #block
                                          init,
                                          scheme, #block
                                          scale,
                                          near.zero.var,
                                          mode,
                                          tol,
                                          max.iter)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if (!is.list(X))
        stop("X must be a list")
    
    # check names on X are unique
    if (length(unique(names(X))) != length(X))
        stop("Each block of 'X' must have a unique name.")
    #names(X)=paste0("block", 1:length(X)) #add names to the blocks if no names
    # or not unique name for each block
    
    if (length(unique(unlist(lapply(X, nrow)))) != 1)
        stop("Unequal number of rows among the blocks of 'X'")
    
    #check for numeric vector instead of matrix
    if (length(unlist(lapply(X, nrow))) != length(X))
        #means there is at least one NULL
    {
        nrow.all = unlist(lapply(X, nrow)) # could be NULL
        ind.null = which(is.na(match(names(X), names(nrow.all))))
        #gives the vectors
        
        # we need matrices with rownames and colnames
        stop("At least one block of 'X' is not a matrix")
    }
    
    if ((missing(indY) & missing(Y)))
        stop("Either 'Y' or 'indY' is needed")
        
    # if DA and the unmapped Y has rows without associated class
    if (DA) {
      Y.tmp <- NULL
      if (missing(Y)) {
        Y.tmp <- X[[indY]]
      } else {
        Y.tmp <- Y
      }
      
      if (length(which(rowSums(Y.tmp)==0)) != 0) {
        stop("Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector")
      }
      rm(Y.tmp)
    }
    
    # if DA and the unmapped Y has rows without associated class
    if (DA) {
        Y.tmp <- NULL
        if (missing(Y)) {
            Y.tmp <- X[[indY]]
        } else {
            Y.tmp <- Y
        }
        
        if (length(which(rowSums(Y.tmp)==0)) != 0) {
            stop("Unmapped Y contains samples with no associated class. May be caused by NAs in input Y vector")
        }
        rm(Y.tmp)
    }
    
    if (missing(ncomp))
        ncomp = 1
    
    #check length(ncomp)=1
    if (length(ncomp) != 1)
        stop("'ncomp' must be a single value")
    
    # transform ncomp to length(X)
    ncomp = rep(ncomp, length(X))
    
    
    #check dimnames and ncomp per block of A
    for (q in 1:length(X))
    {
        check = Check.entry.single(X[[q]], ncomp[q], q=q)
        X[[q]] = check$X
        ncomp[q] = check$ncomp
    }
    
    #check ncomp[q] < ncol(X[[q]])
    for (q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[", q, "]'
            to ncol(X[[", q, "]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    # construction of keepA for X, we include Y later on if needed (i.e. if Y is
    #   provided, otherwise Y is in X[[indY]])
    check = get.keepA(X=X, keepX=keepX, ncomp=ncomp)
    keepA = check$keepA
    
    
    if (!missing(Y))# Y is not missing, so we don't care about indY
    {
        if (!missing(indY))
            warning("'Y' and 'indY' are provided, 'Y' is used.")
        
        if (is.list(Y) | length(dim(Y)) != 2 | !is.numeric(Y))
            stop("'Y' must be a numeric matrix.")
        
        #check dimnames and ncomp per block of A
        check = Check.entry.single(Y, max(ncomp), q=1)
        Y = check$X
        
        if (nrow(Y)!=nrow(X[[1]]))
            stop("Different number of samples between the blocks and Y")
        
        # if not missing, we transform keepY in list to input them in get.keepA
        if (!missing(keepY))
            keepY = list(Y = keepY)
        
        check.temp.Y = get.keepA(X = list(Y = Y), keepX = keepY, ncomp = ncomp)
        keepY.temp = check.temp.Y$keepA
        
        keepA[[length(X)+1]] = keepY.temp[[1]] #add keepY in keepA
        
        # check design matrix before adding Y in
        A = X
        ### Start check design matrix
        if (missing(design))
        {
            design = 1 - diag(length(A)+1)
            rownames(design) = colnames(design) = c(names(A), "Y")
        } else if (is.null(dim(design)))
        {
            design <- .create_design(X = X, design = design, indY = indY)
        } else if (ncol(design) != nrow(design) || ncol(design) < length(X) ||
                   ncol(design) > (length(X) + 1) ){#|| any(!design %in% c(0,1))) {
            stop(paste0("'design' must be a square matrix with ", length(X),
                        "columns."))
        } else if (ncol(design) == length(X)) {
            message("Design matrix has changed to include Y; each block will be
            linked to Y.")
            design = rbind(cbind(design, 1), 1)
            diag(design) = 0
            rownames(design) = colnames(design) = c(names(A), "Y")
        }
        rownames(design) = colnames(design) = c(names(A), "Y")
        ### End check design matrix
        
        # build the list A by adding Y, and creating indY
        A[[length(A)+1]] = Y
        names(A)[length(A)] = names(keepA)[length(A)] = "Y"
        indY = length(A)
        
        if (mode == "canonical")
            ncomp = c(ncomp, min(ncomp, ncol(Y) - 1)) #adjust ncomp for Y
        
        if (mode == "regression")
            ncomp = c(ncomp, max(ncomp)) #adjust ncomp for Y
        
    } else {        #missing(Y) but indY not missing
        
        if (!missing(keepY))
            warning("indY is provided: 'keepY' is ignored and only
            'keepX' is used.")
        
        A = X #list provided as input
        
        ### Start check design matrix
        if (missing(design))
        {
            design = 1 - diag(length(A))
        } else if (is.null(dim(design)))
        {
            design <- .create_design(X = X, design = design, indY = indY)
        }  else if (ncol(design) != nrow(design) || ncol(design) < length(A) ||
                   ncol(design) > (length(A) + 1) )#|| any( !design %in% c(0,1)))
        {
            stop(paste0("'design' must be a square matrix with ", length(A),
                        "columns."))
        }
        rownames(design) = colnames(design) = names(A)
        ### End check design matrix
        
        # check indY
        if (!is.numeric(indY) | indY > length(A) | length(indY) > 1)
            stop(paste0("'indY' must be a numeric value lower or equal to ",
                        length(A), ", the number of blocks in A."))
        
    }
    
    # -----------------------------------------------------------------------
    # at this stage, we have A, indY, keepA, ncomp verified
    # -----------------------------------------------------------------------
    
    # force all rownames to be the same
    ind.names=lapply(A,rownames)
    check = sapply(1:length(A),function(j){identical(ind.names[[1]],
                                                     ind.names[[j]])})
    
    if(sum(check) != length(check))
        stop("Please check the rownames of the data, there seems to be some
    discrepancies")
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
        stop("Choose one of the four following modes: canonical, invariant,
    classic or regression")
    
    #set the default study factor
    if (missing(study))
    {
        study = factor(rep(1, nrow(A[[1]])))
    } else {
        study = as.factor(study)
    }
    if (length(study) != nrow(A[[1]]))
        stop(paste0("'study' must be a factor of length ", nrow(A[[1]]), "."))
    
    if (any(table(study) <= 1))
        stop("At least one study has only one sample, please consider removing
        before calling the function again")
    if (any(table(study) < 5))
        warning("At least one study has less than 5 samples, mean centering
        might not do as expected")
    
    if (missing(init))
        init = "svd.single"
    
    if (!init%in%c("svd","svd.single"))
        stop("init should be one of 'svd' or 'svd.single'")
    
    if (!abs(indY - round(indY) < 1e-25))
        stop ("indY must be an integer")
    if (indY > length(A))
        stop ("indY must point to a block of A")
    
    # =====================================================
    # with or without tau (RGGCA or mint.block.spls algo)
    # =====================================================
    x = unlist(lapply(A,nrow))
    if (!isTRUE(all.equal(max(x), min(x))))
        stop("The samplesize must be the same for all blocks")
    
    #check scheme
    if (!(scheme %in% c("horst", "factorial", "centroid")))
    {
        stop("Choose one of the three following schemes: horst, centroid or
            factorial")
    }
    
    if (!is.numeric(tol) | tol <= 0)
        stop("tol must be non negative")
    if (!is.numeric(max.iter) | max.iter <= 0)
        stop("max.iter must be non negative")
    
    if (!is.logical(scale))
        stop("scale must be either TRUE or FALSE")
    if (!is.logical(near.zero.var))
        stop("near.zero.var must be either TRUE or FALSE")
    
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A, nearZeroVar)
        for(q in 1:length(A))
        {
          if (length(nzv.A[[q]]$Position) <= 0) { next }
          if (DA && q == indY) { next }
          
          names.remove.X = colnames(A[[q]])[nzv.A[[q]]$Position]
          A[[q]] = A[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
          #if (verbose)
          #warning("Zero- or near-zero variance predictors.\n
          #Reset predictors matrix to not near-zero variance predictors.\n
          # See $nzv for problematic predictors.")
          if (ncol(A[[q]]) == 0)
              stop(paste0("No more variables in",A[[q]]))
          
          #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
          if (any(keepA[[q]] > ncol(A[[q]])))
          {
              ind = which(keepA[[q]] > ncol(A[[q]]))
              keepA[[q]][ind] = ncol(A[[q]])
          }
            
        }
    } else {
        nzv.A=NULL
    }
    
    for(q in 1:length(A))
    {
      vars <- apply(A[[q]], 2, sd)^2
      if (length(which(vars==0)) >0) {
        stop(sprintf("There are features with zero variance in block '%s'. If nearZeroVar() function or 'near.zero.var' parameter hasn't been used,  please use it. If you have used one of these, you may need to manually filter out these features.", names(A)[q]), call.=FALSE)
      }
    }
    return(list(A=A, ncomp=ncomp, study=study, keepA=keepA,
                indY=indY, design=design, init=init, nzv.A=nzv.A))
}


# --------------------------------------
# Check.entry.sgcca: used in 'wrapper.sgcca.R'
# --------------------------------------
# X: a list of data sets (called 'blocks') matching on the same samples.
#   Data in the list should be arranged in samples x variables, with samples
#   order matching in all data sets. \code{NA}s are not allowed.
# design: the input design.
# ncomp: numeric vector of length the number of blocks in \code{X}.
#   The number of components to include in the model for each block
#   (does not necessarily need to take the same value for each block).
#   By default set to 2 per block.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid".
#   Default to "centroid"
# mode: input mode, one of "canonical", "classic", "invariant" or "regression".
#   Default to "regression"
# scale: Logical. If scale = TRUE, each block is standardized to zero means and
#   unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single".
#   Default to "svd"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: Logical, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of
#   X[[i]]-variables kept in the model on the last components.

Check.entry.sgcca = function(X,
                             design,
                             ncomp,
                             scheme,
                             mode,
                             scale,
                             init,
                             tol,
                             max.iter,
                             near.zero.var,
                             keepX)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if (!is.list(X))
        stop("X must be a list of at least two matrices")
    
    if (length(X)<2)
        stop("X must be a list of at least two matrices")
    
    if (length(unique(unlist(lapply(X, nrow)))) != 1)
        stop("Unequal number of rows among the blocks of 'X'")
    
    # check names on X are unique
    if (length(unique(names(X))) != length(X))
        stop("Each block of 'X' must have a unique name.")
    #names(X)=paste0("block", 1:length(X)) #add names to the blocks if no
    #names or not unique name for each block
    
    if (missing(ncomp))
        ncomp = 1
    
    #check length(ncomp)=1
    if (length(ncomp) != 1)
        stop("'ncomp' must be a single value")
    
    # transform ncomp to length(X)
    ncomp = rep(ncomp, length(X))
    
    #check dimnames and ncomp per block of A
    for (q in 1:length(X))
    {
        check = Check.entry.single(X[[q]], ncomp[q], q=q)
        X[[q]] = check$X
        ncomp[q] = check$ncomp
    }
    
    #check ncomp[q]<ncol(X[[q]])
    for (q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[", q,
                           "]' to ncol(X[[", q, "]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    A = X#input
    
    if (missing(init))
        init="svd"
    
    if (!init%in%c("svd", "svd.single"))
        stop("init should be one of 'svd' or 'svd.single'")
    
    if (!(mode %in% c("canonical", "invariant", "classic", "regression")))
        stop("Choose one of the four following modes: canonical, invariant, classic
    or regression")
    
    
    #check scheme
    if (missing(scheme))
        scheme= "horst"
    
    if (!(scheme %in% c("horst", "factorial","centroid")))
    {
        stop("Choose one of the three following schemes: horst, centroid or
        factorial")
    }
    
    if(missing(design))
        design = 1 - diag(length(A))
    
    # check design matrix
    if (nrow(design) != ncol(design))
        stop(paste0("'design' must be a square matrix."))
    
    if (nrow(design) != length(A))
        stop(paste0("'design' must be a square matrix with", length(A), "columns."))
    
    if (tol <= 0)
        stop("tol must be non negative")
    
    if (max.iter <= 0)
        stop("max.iter must be non negative")
    
    if (!is.logical(scale))
        stop("scale must be either TRUE or FALSE")
    if (!is.logical(near.zero.var))
        stop("near.zero.var must be either TRUE or FALSE")
    
    # construction of keepA
    check = get.keepA(X=A, keepX=keepX, ncomp=ncomp)
    keepA = check$keepA
    
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A,nearZeroVar)
        for (q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X = colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
                #if (verbose)
                #warning("Zero- or near-zero variance predictors.\n
                #Reset predictors matrix to not near-zero variance predictors.\n
                #See $nzv for problematic predictors.")
                if (ncol(A[[q]]) == 0)
                    stop(paste0("No more variables in", A[[q]]))
                
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if (any(keepA[[q]]>ncol(A[[q]])))
                {
                    ind = which(keepA[[q]] > ncol(A[[q]]))
                    keepA[[q]][ind] = ncol(A[[q]])
                }
            }
            
        }
    } else {
        nzv.A=NULL
    }
    
    
    return(list(A=A, ncomp=ncomp, design=design, init=init, scheme=scheme,
                nzv.A=nzv.A, keepA=keepA))
    
}



# --------------------------------------
# Check.entry.rgcca: used in 'wrapper.rgcca.R'
# --------------------------------------
# X: a list of data sets (called 'blocks') matching on the same samples.
#   Data in the list should be arranged in samples x variables, with samples
#   order matching in all data sets. \code{NA}s are not allowed.
# design: the input design.
# tau:
# ncomp: numeric vector of length the number of blocks in \code{X}.
#   The number of components to include in the model for each block
#   (does not necessarily need to take the same value for each block).
#   By default set to 2 per block.
# scheme: the input scheme, one of "horst", "factorial" or ""centroid".
#   Default to "centroid"
# scale: Logical. If scale = TRUE, each block is standardized to zero means and
#    unit variances (default: TRUE).
# init: intialisation of the algorithm, one of "svd" or "svd.single".
#   Default to "svd"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: Logical, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# keepX: A vector of same length as X.  Each entry keepX[i] is the number of
#   X[[i]]-variables kept in the model on the last components.
Check.entry.rgcca = function(X,
                             design,
                             tau,
                             ncomp,
                             scheme,
                             scale,
                             init,
                             tol,
                             max.iter,
                             near.zero.var,
                             keepX)
{
    #need to give the default values of mint.block.spls to mixOmics
    
    if (!is.list(X))
        stop("X must be a list of at list two matrices")
    
    if (length(X)<2)
        stop("X must be a list of at list two matrices")
    
    if (length(unique(unlist(lapply(X, nrow)))) != 1)
        stop("Unequal number of rows among the blocks of 'X'")
    
    # check names on X are unique
    if (length(unique(names(X))) != length(X))
        stop("Each block of 'X' must have a unique name.")
    #names(X)=paste0("block", 1:length(X)) #add names to the blocks if no names
    #or not unique name for each block
    
    if (is.null(tau))
        stop("'tau' is needed")
    
    if (missing(ncomp))
        ncomp = 1
    
    #check length(ncomp)=1
    if (length(ncomp) != 1)
        stop("'ncomp' must be a single value")
    
    # transform ncomp to length(X)
    ncomp = rep(ncomp, length(X))
    
    #check dimnames and ncomp per block of A
    for (q in 1:length(X))
    {
        check = Check.entry.single(X[[q]], ncomp[q], q=q)
        X[[q]] = check$X
        ncomp[q] = check$ncomp
    }
    
    
    #check ncomp[q]<ncol(X[[q]])
    for(q in 1:length(X))
    {
        ncomp[q] = round(ncomp[q])
        if(ncomp[q] > ncol(X[[q]]))
        {
            warning(paste0("Reset maximum number of variates 'ncomp[", q,
                           "]' to ncol(X[[", q, "]])= ", ncol(X[[q]]), "."))
            ncomp[q] = ncol(X[[q]])
        }
    }
    
    A = X#input
    
    if (is.numeric(tau))
    {
        if (any(tau < 0) | any(tau > 1))
            stop("'tau' contains either values between 0 and 1, or 'optimal'.")
        
        if (is.vector(tau))
        {
            if (length(tau) != length(A))
                stop(paste0("'tau' must be of length ", length(A), "."))
            
            tau = matrix(rep(tau, max(ncomp)), nrow = max(ncomp),
                         ncol = length(tau), byrow = TRUE)
        }
    } else {
        if (tau != "optimal")
            stop("'tau' contains either values between 0 and 1, or 'optimal'.")
    }
    
    if (missing(init))
        init = "svd.single"
    
    if (init != "svd.single")
        stop("init should be 'svd.single'.")
    
    # check scheme
    if(missing(scheme)) scheme = "horst"
    if (!(scheme %in% c("horst", "factorial", "centroid")))
    {
        stop("Choose one of the three following schemes: horst, centroid or
        factorial")
    }
    
    
    if (missing(design))
        design = 1 - diag(length(A))
    
    #check design matrix
    if (nrow(design) != ncol(design))
        stop(paste0("'design' must be a square matrix."))
    if (nrow(design) != length(A))
        stop(paste0("'design' must be a square matrix with", length(A), "columns."))
    
    
    if (missing(near.zero.var))
        near.zero.var = FALSE
    
    if (tol <= 0)
        stop("tol must be non negative")
    
    if (max.iter <= 0)
        stop("max.iter must be non negative")
    
    if (!is.logical(scale))
        stop("scale must be either TRUE or FALSE")
    if (!is.logical(near.zero.var))
        stop("near.zero.var must be either TRUE or FALSE")
    
    # construction of keepA
    check = get.keepA(X=A, keepX=keepX, ncomp=ncomp)
    keepA = check$keepA
    
    ### near.zero.var, remove the variables with very small variances
    if(near.zero.var == TRUE)
    {
        nzv.A = lapply(A,nearZeroVar)
        for (q in 1:length(A))
        {
            if (length(nzv.A[[q]]$Position) > 0)
            {
                names.remove.X = colnames(A[[q]])[nzv.A[[q]]$Position]
                A[[q]] = A[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
                #if (verbose)
                #warning("Zero- or near-zero variance predictors.\n
                #Reset predictors matrix to not near-zero variance predictors.\n
                #See $nzv for problematic predictors.")
                if (ncol(A[[q]]) == 0)
                    stop(paste0("No more variables in", A[[q]]))
                
                #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
                if (any(keepA[[q]] > ncol(A[[q]])))
                {
                    ind = which(keepA[[q]] > ncol(A[[q]]))
                    keepA[[q]][ind] = ncol(A[[q]])
                }
            }
            
        }
    } else {
        nzv.A = NULL
    }
    
    
    return(list(A=A, ncomp=ncomp, design=design, init=init, scheme=scheme,
                nzv.A=nzv.A, keepA=keepA))
}





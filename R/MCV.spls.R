#############################################################################################################
# Authors:
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 24-08-2016
# last modified: 05-10-2017
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
#############################################################################################################


# ========================================================================================================
# tune.splsda: chose the optimal number of parameters per component on a splsda method
# ========================================================================================================

# X: numeric matrix of predictors
# Y: a factor or a class vector for the discrete outcome
# multilevel: repeated measurement only
# validation: Mfold or loo cross validation
# folds: if validation = Mfold, how many folds?
# nrepeat: number of replication of the Mfold process
# ncomp: the number of components to include in the model. Default to 1.
# choice.keepX: a vector giving keepX on the components that were already tuned
# test.keepX: grid of keepX among which to chose the optimal one
# measure: one of c("overall","BER"). Accuracy measure used in the cross validation processs
# dist: distance to classify samples. see predict
# auc:
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: Logical, see the internal \code{\link{nearZeroVar}} function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations
# progressBar: show progress,
# class.object
# cl: if parallel, the clusters
# scale: Logical. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE).
# misdata: optional. any missing values in the data? list, misdata[[q]] for each data set
# is.na.A: optional. where are the missing values? list, is.na.A[[q]] for each data set (if misdata[[q]] == TRUE)
# ind.NA: optional. which rows have missing values? list, ind.NA[[q]] for each data set.
# ind.NA.col: optional. which col have missing values? list, ind.NA.col[[q]] for each data set.
# parallel: logical.

#' @importFrom matrixStats colVars
MCVfold.spls <- function(
    X,
    Y,
    multilevel = NULL, # repeated measurement only
    validation,
    folds,
    nrepeat = 1,
    ncomp,
    choice.keepX = NULL, # keepX chosen on the first components
    test.keepX, # a vector of value(keepX) to test on the last component. There needs to be names(test.keepX)
    test.keepY, # a vector of value(keepX) to test on the last component. There needs to be names(test.keepX)
    measure = c("overall"), # one of c("overall","BER")
    dist = "max.dist",
    auc = FALSE,
    tol = 1e-06,
    max.iter = 100,
    near.zero.var = FALSE,
    progressBar = FALSE,
    class.object,
    scale,
    misdata,
    is.na.A,
    #ind.NA,
    #ind.NA.col,
    BPPARAM
)
{   

  
  #-- Define function to process one fold --------------------------------------#
  #---------------------------------------------------------------------------#
  # this is called inside the for-loop over repeats below
  
    function.j.folds <- function(j)#for (j in 1:M)
    {
      if (progressBar ==  TRUE){setTxtProgressBar(pb, (M*(nrep-1)+j-1)/(M*nrepeat))} # update progress bar
      
      # identify the test (omitted) samples
      omit = which(repeated.measure %in% folds[[j]] == TRUE) # finds the row indices that belong to the current fold
      
      # split the data into training and test sets
      X.train = X[-omit, ]
      X.test = X[omit, , drop = FALSE]
      
      if(any(class.object == "DA")){
        Y.train = Y[-omit]
        Y.train.mat = unmap(Y.train)
        q = ncol(Y.train.mat)
        colnames(Y.train.mat) = levels(Y.train)
        Y.test = Y[omit]
      } else {
        Y.train.mat = Y[-omit, , drop = FALSE]
        q = ncol(Y.train.mat)
        Y.test = Y[omit, , drop = FALSE]
      }
      
      # remove variables with no or near-zero variance
      remove <- NULL
      var.train <- colVars(X.train, na.rm <- TRUE)
      ind.var <- which(var.train == 0)
      if (length(ind.var) > 0){
        remove <- c(remove, colnames(X.train)[ind.var])
        X.train = X.train[, -c(ind.var),drop = FALSE]
        X.test = X.test[, -c(ind.var),drop = FALSE]
        # reduce choice.keepX and test.keepX if needed
        if (any(choice.keepX > ncol(X.train)))
          choice.keepX[which(choice.keepX>ncol(X.train))] = ncol(X.train)
        # reduce test.keepX if needed
        if (any(test.keepX > ncol(X.train)))
          test.keepX[which(test.keepX>ncol(X.train))] = ncol(X.train)
      }
      if(near.zero.var == TRUE){
        remove.zero = nearZeroVar(X.train)$Position
        if (length(remove.zero) > 0){
          remove = c(remove, colnames(X.train)[remove.zero])
          X.train = X.train[, -c(remove.zero),drop = FALSE]
          X.test = X.test[, -c(remove.zero),drop = FALSE]
          # reduce choice.keepX and test.keepX if needed
          if (any(choice.keepX > ncol(X.train)))
            choice.keepX[which(choice.keepX>ncol(X.train))] = ncol(X.train)
          # reduce test.keepX if needed
          if (any(test.keepX > ncol(X.train)))
            test.keepX[which(test.keepX>ncol(X.train))] = ncol(X.train)
        }
      }
      
      # handle missing data - splits missing data between training and test sets
      if(any(misdata)){
        if(any(class.object == "DA")){
          if(length(remove)>0){
            ind.remove = which(colnames(X) %in% remove)
            is.na.A.train = is.na.A[-omit, -ind.remove, drop=FALSE]
            is.na.A.test = list(X=is.na.A[omit, -ind.remove, drop=FALSE])
          } else {
            is.na.A.train = is.na.A[-omit,, drop=FALSE]
            is.na.A.test = list(X=is.na.A[omit,, drop=FALSE])
          }
          temp = which(is.na.A.train, arr.ind=TRUE)
          ind.NA.train = unique(temp[,1])
          ind.NA.col.train = unique(temp[,2])
          is.na.A.train = list(X=is.na.A.train, Y=NULL)
          ind.NA.train = list(X=ind.NA.train, Y=NULL)
          ind.NA.col.train = list(X=ind.NA.col.train, Y=NULL)
        } else{
          if(length(remove)>0){
            ind.remove = which(colnames(X) %in% remove)
            is.na.A.train = list(X=is.na.A[[1]][omit, -ind.remove, drop=FALSE], Y=is.na.A[[1]][omit,, drop=FALSE])
            is.na.A.test = list(X=is.na.A[[1]][omit, -ind.remove, drop=FALSE]) #only for X
          }else {
            is.na.A.train = lapply(is.na.A, function(x){x[-omit,, drop=FALSE]})
            is.na.A.test = list(X=is.na.A[[1]][omit,, drop=FALSE]) #only for X
          }
          temp = lapply(is.na.A.train,function(x){which(x,arr.ind=TRUE)})
          ind.NA.train = lapply(temp,function(x){unique(x[,1])})
          ind.NA.col.train = lapply(temp,function(x){unique(x[,2])})
        }
        
      } else {
        is.na.A.train = is.na.A.test =NULL
        ind.NA.train = NULL
        ind.NA.col.train = NULL
      }
      
      # prepare arrays to store predictions on the fold
      class.comp.j = list()
      if(any(class.object == "DA")){
        prediction.comp.j = array(0, c(length(omit), nlevels(Y), length(test.keepX)), dimnames = list(rownames(X.test), levels(Y), names(test.keepX)))
        for(ijk in dist)
          class.comp.j[[ijk]] = matrix(0, nrow = length(omit), ncol = length(test.keepX))# prediction of all samples for each test.keepX and  nrep at comp fixed
      } else {
        prediction.comp.j = array(0, c(length(omit), ncol(Y), length(test.keepX)), dimnames = list(rownames(X.test), colnames(Y), test.keepX))
      }
      
      # fit the model on the training data - `result' returns loadings and variates for all test.keepX on the ncomp component
      result = suppressWarnings(internal_wrapper.mint(X=X.train, Y=Y.train.mat, study=factor(rep(1,nrow(X.train))), ncomp=ncomp,
                                                      keepX=choice.keepX, keepY=rep(ncol(Y.train.mat), ncomp-1), test.keepX=test.keepX, test.keepY=test.keepY,
                                                      mode="regression", scale=scale, near.zero.var=near.zero.var,
                                                      max.iter=max.iter, logratio="none", DA=TRUE, multilevel=NULL,
                                                      misdata = misdata, is.na.A = is.na.A.train, ind.NA = ind.NA.train,
                                                      ind.NA.col = ind.NA.col.train, all.outputs=FALSE))
      
      # scale the test data in the same way as the training set has been done above
      if (!is.null(attr(result$A[[1]], "scaled:center"))){X.test = sweep(X.test, 2, STATS = attr(result$A[[1]], "scaled:center"))}
      if (scale){X.test = sweep(X.test, 2, FUN = "/", STATS = attr(result$A[[1]], "scaled:scale"))}
      means.Y = matrix(attr(result$A[[2]], "scaled:center"),nrow=nrow(X.test),ncol=q,byrow=TRUE);
      if (scale)
      {sigma.Y = matrix(attr(result$A[[2]], "scaled:scale"),nrow=nrow(X.test),ncol=q,byrow=TRUE)}else{sigma.Y=matrix(1,nrow=nrow(X.test),ncol=q)}
      
      # record prediction results for each test.keepX
      keepA = result$keepA
      test.keepA = keepA[[ncomp]]
      
      # save variates and loadings for all test.keepA
      result.temp = list(variates = result$variates, loadings = result$loadings)
      
      # creates temporary splsda object to use the predict function
      result$X = result$A$X
      
      # add the "splsda" or "spls" class
      if(any(class.object == "DA")){
        class(result) = c("mixo_splsda","mixo_spls","DA")
        result$ind.mat = result$A$Y
        result$Y = factor(Y.train)
      } else{
        class(result) = c("mixo_spls")
        result$Y = result$A$Y
      }
      result$A = NULL
      
      # make predictions for each candidate model on the test set
      for(i in seq_len(nrow(test.keepA))){ # loop through each candidate i.e. row of test.keepA
        
        # pick the appropriate loadings and variates based on the candidate's keepX value
        names.to.pick = NULL
        if(ncomp>1) {
          names.to.pick = unlist(lapply(seq_len(ncomp-1), function(x){paste(paste0("comp",x),apply(keepA[[x]],1,function(x) paste(x,collapse="_")), sep=":")}))
        }
        names.to.pick.ncomp = paste(paste0("comp",ncomp),paste(as.numeric(keepA[[ncomp]][i,]),collapse="_"), sep=":")
        names.to.pick = c(names.to.pick, names.to.pick.ncomp)
        
        #change variates and loadings for each test.keepA
        result$variates = lapply(result.temp$variates, function(x){if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick, drop=FALSE]}else{x}})
        result$loadings = lapply(result.temp$loadings, function(x){if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick, drop=FALSE]}else{x}})
        
        # record selected features
        if (any(class.object == "mixo_splsda") & length(test.keepX) ==  1) # only done if splsda and if only one test.keepX as not used if more so far
          # note: if plsda, 'features' includes everything: to optimise computational time, we don't evaluate for plsda object
          features.j = selectVar(result, comp = ncomp)$name
        
        # predict responses for the test set using the scaled X.test and specified distance metric -> stored in prediction.comp.j
        test.predict.sw <- predict.mixo_spls(result, newdata.scale = X.test, dist = dist, misdata.all=misdata[1], is.na.X = is.na.A.train, is.na.newdata = is.na.A.test)
        prediction.comp.j[, , i] =  test.predict.sw$predict[, , ncomp]
        if(any(class.object == "DA")){
          for(ijk in dist)
            class.comp.j[[ijk]][, i] =  test.predict.sw$class[[ijk]][, ncomp] # if DA predicted labels stored in class.comp.j
        }
      } # end loop
      
      return(list(class.comp.j = class.comp.j, prediction.comp.j = prediction.comp.j, features = features.j, omit = omit))
      # class predictions for each candidate (if DA), numeric predictions for each candidate, any selected features, the indices of the test samples for this fold
      
    }
  
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#

    # progress bar
    if (progressBar ==  TRUE)
    {pb <- txtProgressBar(style = 3)
    nBar <- 1} else {pb = FALSE}

    # set AUC
    if(any(measure == "AUC")) {auc <- TRUE}

    # set design matric
    design <- matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
    #      [,1] [,2]
    # [1,]    0    1
    # [2,]    1    0
    
    # decide on the number of response variables (or dummy variables) to keep for each component beyond the first
    if(ncomp>1 &  any(class.object == "DA")) {keepY <- rep(nlevels(Y), ncomp-1)} else {keepY <- rep(ncol(Y),ncomp-1)}

    # save row name information
    rownames.X <- rownames(X)

    # number of cross-validation folds
    M <- length(folds)
    folds.input <- folds
    
    #-- initialise outputs ----------------------------------------#
    #---------------------------------------------------------------------------#

    # initialise outputs
    features <- features.j <- NULL
    auc.all <- prediction.comp <- class.comp <- list()
    
    # initialise the output matrices based on whether object is DA or not
    if(any(class.object == "DA")){
        for(ijk in dist)
            class.comp[[ijk]] = array(0, c(nrow(X), nrepeat, length(test.keepX)))# prediction of all samples for each test.keepX and  nrep at comp fixed
    } else {
        prediction.keepX = vector("list", length=length(test.keepX))
        for(i in seq_len(length(test.keepX)))
            prediction.keepX[[i]] = array(0,c(nrow(X), ncol(Y), nrepeat), dimnames = list(rownames(X), colnames(Y), paste0("nrep.", seq_len(nrepeat))))
    }
    
    #-- loop over the number of repeats ----------------------------------------#
    #---------------------------------------------------------------------------#
    
    for(nrep in seq_len(nrepeat))
    {
        ##### Initialise empty arrays to store predictions #####
      
        n = nrow(X) # number of samples
        
        if(any(class.object == "DA")){
            # create empty `prediction.comp` array of dimensions n x nlevels(Y) x length(test.keepX)
            prediction.comp[[nrep]] = array(0, c(n, nlevels(Y), length(test.keepX)), dimnames = list(rownames.X, levels(Y), names(test.keepX)))
            colnames(prediction.comp[[nrep]]) = levels(Y) # column names are levels of Y
            # create empty `auc.all` array of dimensions nlevels(Y) or just 1 if have 2 levels x 2 (AUC and p-value) x length(test.keepX)
            if(nlevels(Y)>2)
            { auc.all[[nrep]] = array(0, c(nlevels(Y),2, length(test.keepX)), dimnames = list(paste(levels(Y), "vs Other(s)"), c("AUC","p-value"), names(test.keepX)))
            }else{auc.all[[nrep]] = array(0, c(1,2, length(test.keepX)), dimnames = list(paste(levels(Y)[1], levels(Y)[2], sep = " vs "), c("AUC","p-value"), names(test.keepX)))}
            
        } else {
            # create empty `prediction.comp` array of dimensions n x ncol(Y) x length(test.keepX)
            prediction.comp[[nrep]] = array(0, c(nrow(X), ncol(Y), length(test.keepX)), dimnames = list(rownames(X), colnames(Y), test.keepX))
            colnames(prediction.comp[[nrep]]) = colnames(Y)
        }
        
        rownames(prediction.comp[[nrep]]) = rownames.X # rownames of `prediction.comp` match rownames of X
        
        #####  Deal with multilevel data  #####
        
        repeated.measure = seq_len(n)
        if (!is.null(multilevel)) # assigns identifier to each sample grouping rows by repeated measures
        {repeated.measure = multilevel[,1]
        n = length(unique(repeated.measure))}
        
        #####  Define the folds  #####

        # Mfold validation
        if (validation ==  "Mfold"){
          
          if (nrep > 1) {folds = folds.input} # reset the folds if this is not the first repetition
          if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n){stop("Invalid number of folds.")} # checks fold value is valid
          
          else {M = round(folds)
              if (is.null(multilevel)){
                # for DA object not multilevel
                    if(any(class.object == "DA")){
                        temp = stratified.subsampling(Y, folds = M) # splits samples into folds while trying to balance classes
                        folds = temp$SAMPLE
                        if(temp$stop > 0 & nrep == 1){ # show this warning only once
                          warning("At least one class is not represented in one fold, which may unbalance the error rate.\n  Consider a number of folds lower than the minimum in table(Y): ", min(table(Y)))
                          rm(temp)}
                    } else {
                        # for non-DA object not multilevel
                        folds = suppressWarnings(split(sample(seq_len(n)), rep(seq_len(n), length = M)))
                    }
                # for non multilevel object    
                } else {
                    folds = split(sample(seq_len(n)), rep(seq_len(M), length = n)) # needs to have all repeated samples in the same fold
                }
          }
          
        # LOO validation
        } else if (validation ==  "loo") {
            folds = split(seq_len(n), rep(seq_len(n), length = n))
            M = n
        }
        
        # creates object `folds`, list with `M` elements each containing indices (or group IDs if multilevel) for the test set for that fold
        M = length(folds)
        
        #####  Prepare prediction error matrix  #####
        
        error.sw = matrix(0, nrow <- M, ncol <- length(test.keepX)) # dim number of folds (`M`) x length(test.keepX)
        rownames(error.sw) = paste0("fold",seq_len(M))
        colnames(error.sw) = test.keepX
        
        #####  Run function on each fold #####

        if (!is.null(BPPARAM)) {result.all = bplapply(seq_len(M), function.j.folds, BPPARAM = BPPARAM)
        } else {result.all = lapply(seq_len(M), function.j.folds)}
        
        #####  Combine the results #####
        
        for(j in seq_len(M)){
            omit = result.all[[j]]$omit
            prediction.comp.j = result.all[[j]]$prediction.comp.j
            prediction.comp[[nrep]][omit, , ] = prediction.comp.j
            
            if(any(class.object == "DA")){
                class.comp.j = result.all[[j]]$class.comp.j
                for(ijk in dist)
                    class.comp[[ijk]][omit,nrep, ] = class.comp.j[[ijk]]
                if (length(test.keepX) ==  1) # only done if splsda and if only one test.keepX as not used if more so far
                    features = c(features, result.all[[j]]$features)
            }}
        if(!any(class.object == "DA")){
            # create a array, for each keepX: n*q*nrep
            for(i in seq_len(length(test.keepX)))
                prediction.keepX[[i]][,,nrep] = prediction.comp[[nrep]][,,i, drop=FALSE]
        }
        
        #####  Progress bar #####
        
        if (progressBar ==  TRUE){setTxtProgressBar(pb, (M*nrep)/(M*nrepeat))}
        
        ##### Compute AUC #####
        
        if(auc){
            data=list()
            for (i in seq_len(length(test.keepX))){
                data$outcome = Y
                data$data = prediction.comp[[nrep]][, , i]
                auc.all[[nrep]][, , i] = statauc(data)[[1]] #as.matrix(statauc(data)[[1]]) # [[1]] because [[2]] is the graph
            }}
        
    } #end nrep 1:nrepeat
    #### Summary of objects created by this loop:
    # `prediction.comp[[nrep]]` - n_samples x n_classes (for DA) or n_columns of Y (for regression) x length(test.keepX)
    # `class.comp[[ijk]]` - for each distance metric, array of predictions for the test samples (test data subsets x number of repeats)
    # `auc.all[[nrep]]` - if AUC is computed, array stores AUC and p-value for each candidate model
    # `error.sw` - temporary matrix (per repetition) that holds error rates for each fold and candidate
    # `result.all` - list of length M (number of folds) where each element is a list returned by inner function for that fold
    
    # rename slots of output arrays
    if(auc) {names(auc.all) <- paste0("nrep.", seq_len(nrepeat))}
    names(prediction.comp) <- paste0("nrep.", seq_len(nrepeat))
    
    #-------------------------------------------------------------#
    #----- average AUC over the nrepeat, for each test.keepX -----#
    
    if(auc){
        # matrix of AUC per keepX, per nrepeat, averaged over classes
        auc.mean.sd.over.class =  matrix(0, nrow= length(test.keepX), ncol=nrepeat, dimnames = list(names(test.keepX), paste0("nrep.", seq_len(nrepeat))))
        
        # matrix of AUC per keepX, per class, averaged over nrepeat
        if(nlevels(Y)>2){
            auc.mean.sd.over.nrepeat =  array(0, c(nlevels(Y),2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }else{
            auc.mean.sd.over.nrepeat =  array(0, c(1,2, length(test.keepX)), dimnames = list(rownames(auc.all[[1]]), c("AUC.mean","AUC.sd"), names(test.keepX)))
        }
        
        for(i in seq_len(length(test.keepX))){
            temp = NULL
            for(nrep in seq_len(nrepeat)){
                temp = cbind(temp, auc.all[[nrep]][, 1, i])
                auc.mean.sd.over.class[i,nrep] = mean (auc.all[[nrep]][,1,i])}
            auc.mean.sd.over.nrepeat[, 1, i] = apply(temp,1,mean)
            auc.mean.sd.over.nrepeat[, 2, i] = apply(temp,1,sd)
        }
    } else {
        auc.mean.sd.over.nrepeat = auc.all  = NULL
        
    }
    
    #-------------------------------------------------------------#
    #----- initialise output lists -----#
    
    result = list()
    error.mean = error.sd = error.per.class.keepX.opt.comp = keepX.opt = test.keepX.out = test.keepY.out = mat.error.final = choice.keepX.out = choice.keepY.out = list()
    
    #-------------------------------------------------------------#
    #----- `overall` measure -----#
    
    if (any(measure == "overall")){
      
        # for each distance in distance measures
        for(ijk in dist){
            rownames(class.comp[[ijk]]) = rownames.X
            colnames(class.comp[[ijk]]) = paste0("nrep.", seq_len(nrepeat))
            dimnames(class.comp[[ijk]])[[3]] = paste0("test.keepX.",names(test.keepX))
            
            # compute the misclassification error
            error = apply(class.comp[[ijk]],c(3,2),function(x){sum(as.character(Y) != x)})
            rownames(error) = names(test.keepX)
            colnames(error) = paste0("nrep.",seq_len(nrepeat))
            
            # average error over repeats 
            error.mean[[ijk]] = apply(error,1,mean)/length(Y)
            if (!nrepeat ==  1){error.sd[[ijk]] = apply(error,1,sd)/length(Y)}
            mat.error.final[[ijk]] = error/length(Y)  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
            
            # select optimal candidate `keepX.opt`
            keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
            
            # compute per-class error (confusion matrices)
            error.per.class.keepX.opt.comp[[ijk]] <- apply(class.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x){
                conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", seq_len(nrepeat))
            
            # store chosen candidate and build final output
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
            result$"overall"$error.rate.mean = error.mean
            if (!nrepeat ==  1){result$"overall"$error.rate.sd = error.sd}
            result$"overall"$confusion = error.per.class.keepX.opt.comp
            result$"overall"$mat.error.rate = mat.error.final
            result$"overall"$keepX.opt = test.keepX.out
        }
    }
    
    if (any(measure ==  "BER"))
    {
        for(ijk in dist)
        {
            rownames(class.comp[[ijk]]) = rownames.X
            colnames(class.comp[[ijk]]) = paste0("nrep.", seq_len(nrepeat))
            dimnames(class.comp[[ijk]])[[3]] = paste0("test.keepX.",names(test.keepX))
            
            error = apply(class.comp[[ijk]],c(3,2),function(x)
            {
                conf = get.confusion_matrix(truth = factor(Y),predicted = x)
                get.BER(conf)
            })
            rownames(error) = names(test.keepX)
            colnames(error) = paste0("nrep.",seq_len(nrepeat))
            
            # average BER over the nrepeat
            error.mean[[ijk]] = apply(error,1,mean)
            if (!nrepeat ==  1)
                error.sd[[ijk]] = apply(error,1,sd)
            
            mat.error.final[[ijk]] = error  # BER for each test.keepX (rows) and each nrepeat (columns)
            
            keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1]
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] = apply(class.comp[[ijk]][, , keepX.opt[[ijk]], drop = FALSE], 2, function(x)
            {
                conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
            })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            colnames(error.per.class.keepX.opt.comp[[ijk]]) = paste0("nrep.", seq_len(nrepeat))
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
            
            result$"BER"$error.rate.mean = error.mean
            if (!nrepeat ==  1)
                result$"BER"$error.rate.sd = error.sd
            
            result$"BER"$confusion = error.per.class.keepX.opt.comp
            result$"BER"$mat.error.rate = mat.error.final
            result$"BER"$keepX.opt = test.keepX.out
            
        }
        
        
    }
    
    if (any(measure ==  "AUC"))
    {
        # no loop in ijk because AUC does not use any distances
        # still need to put [[1]] to be similar to other distances
        
        #save(list=ls(),file="temp.Rdata")
        # average the AUCs over the class (AUC per keepX)
        error.mean[[1]] = apply(auc.mean.sd.over.class, 1, mean)
        
        #choose highest AUC
        keepX.opt[[1]] = which(error.mean[[1]] ==  max(error.mean[[1]]))[1]
        
        # match to the best keepX
        test.keepX.out[[1]]= test.keepX[keepX.opt[[1]]]
        choice.keepX.out[[1]] = c(choice.keepX, test.keepX.out)
        
        result$"AUC"$error.rate.mean = error.mean
        
        # we are calculating a mean over nrepeat from means over the groups
        # we do not outputs SD so far (which one should we use?)
        
        #if (!nrepeat ==  1)
        #result$"AUC"$error.rate.sd = list(list())
        #result$"AUC"$confusion = list(list())
        
        result$"AUC"$mat.error.rate = list(auc.mean.sd.over.class)
        result$"AUC"$keepX.opt = test.keepX.out
    }
    
    
    
    if (any(measure == "MSE"))
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)
        
        # MSE error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
        {
            temp = (Y-x)^2
            temp2=t(t(temp)/varY)
            temp3=apply(temp2,2, sum)
            temp3
        })})
        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
        
        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y)
        if (!nrepeat ==  1)
            error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
        
        result$"MSE"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
            result$"MSE"$error.rate.sd = error.sd
        
        result$"MSE"$mat.error.rate = mat.error.final
        result$"MSE"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "MAE")) # MAE (Mean Absolute Error: MSE without the square)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)
        
        # MAE error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
        {
            temp = abs(Y-x)
            temp2=t(t(temp)/varY)
            temp3=apply(temp2,2, sum)
            temp3
        })})
        
        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
        
        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y)
        if (!nrepeat ==  1)
            error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
        
        result$"MAE"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
            result$"MAE"$error.rate.sd = error.sd
        
        result$"MAE"$mat.error.rate = mat.error.final
        result$"MAE"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "Bias")) # Bias (average of the differences)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)
        
        # Bias error for each nrep and each test.keepX: summing over all samples
        varY=apply(Y,2,var)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
        {
            temp = (Y-x)
            temp2=t(t(temp)/varY)
            temp3=abs(apply(temp2,2, sum)) # absolute value of the bias
            temp3
        })})
        
        mat.error.final[[ijk]] = lapply(error, function(x){x/nrow(Y)})  # percentage of misclassification error for each test.keepX (rows) and each nrepeat (columns)
        
        # we want to average the error per keepX over nrepeat and choose the minimum error
        error.mean[[ijk]] = sapply(error, mean)/nrow(Y)
        if (!nrepeat ==  1)
            error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)/nrow(Y)}), mean)
        
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
        
        result$"Bias"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
            result$"Bias"$error.rate.sd = error.sd
        
        result$"Bias"$mat.error.rate = mat.error.final
        result$"Bias"$keepX.opt = test.keepX.out
    }
    
    if (any(measure == "R2")) # R2 (square of the correlation of the truth and the predicted values, averaged over the columns of Y)
    {
        ijk=1 # in case more measure later on
        names(prediction.keepX) = paste0("test.keepX.",test.keepX)
        
        # R2 error for each test.keepX (list),each Y (rows) and each nrepeat (columns)
        error = lapply(prediction.keepX, function(z){apply(z,c(3),function(x)
        {
            temp = diag(cor(Y,x))^2
            temp
        })}) # list for each test.keepX of nrow=ncol(Y) and ncol=nrepeat
        
        mat.error.final[[ijk]] = error
        
        # we want to average the error over the nrepeat and the columns of Y
        error.mean[[ijk]] = sapply(error, mean)
        if (!nrepeat ==  1)
            error.sd[[ijk]] = sapply(lapply(error, function(x){apply(matrix(x,nrow=ncol(Y)), 1, sd)}), mean)
        
        keepX.opt[[ijk]] = which(error.mean[[ijk]] ==  max(error.mean[[ijk]]))[1] # chose the lowest keepX if several minimum
        
        
        test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
        
        choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
        
        result$"R2"$error.rate.mean = error.mean
        if (!nrepeat ==  1)
            result$"R2"$error.rate.sd = error.sd
        
        result$"R2"$mat.error.rate = mat.error.final
        result$"R2"$keepX.opt = test.keepX.out
    }
    
    
    result$prediction.comp = prediction.comp
    if(any(class.object == "DA")){
        if(auc){
            result$auc = auc.mean.sd.over.nrepeat
            result$auc.all = auc.all
        }
        result$class.comp = class.comp
        if (length(test.keepX) ==  1)
            result$features$stable = sort(table(as.factor(features))/M/nrepeat, decreasing = TRUE)
    }
    return(result)
}

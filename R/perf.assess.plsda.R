#############################################################################################################
# Authors:
#   Amrit Singh, University of British Columbia, Vancouver.
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 01-04-2015
# last modified: 27-05-2016
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
#############################################################################################################


# ----------------------------------------------------------------------------------------------------------
# perf.assess.plsda - Function to evaluate the performance of fitted PLS-DA (cross-validation)
# ----------------------------------------------------------------------------------------------------------

## ------------------------------- (s)PLSDA ------------------------------- ##
#' @rdname perf.assess
#' @method perf.assess mixo_plsda
#' @export
perf.assess.mixo_plsda <- function(object,
                            validation = c("Mfold", "loo"),
                            folds,
                            nrepeat = 1,
                            dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
                            auc = FALSE,
                            progressBar = FALSE,
                            BPPARAM = SerialParam(),
                            seed = NULL,
                            ...)
{
    
    #-- initialising arguments --#
    set.seed(seed)
    BPPARAM$RNGseed <- seed
    # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
    X = object$input.X
    level.Y = object$names$colnames$Y  #to make sure the levels are ordered
    Y = object$Y
    ncomp = object$ncomp
    n = nrow(X)
    
    logratio = object$logratio
    if (is.null(logratio))
        logratio = "none"
    
    multilevel = object$multilevel # repeated measurement and Y
    near.zero.var = !is.null(object$nzv) # if near.zero.var was used, we set it to TRUE. if not used, object$nzv is NULL
    
    #-- tells which variables are selected in X and in Y --#
    
    if (is(object, "mixo_splsda"))
    {
        keepX = object$keepX
    } else {
        keepX = rep(ncol(X), ncomp)
    }
    
    tol = object$tol
    max.iter = object$max.iter
    scale = object$scale
    
    # check input arguments
    if (hasArg(method.predict))
        stop("'method.predict' argument has been replaced by 'dist' to match the 'tune' function")
    method.predict = NULL # to pass R CMD check
    
    dist = match.arg(dist, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    if (any(dist == "all"))
    {
        nmthdd = 3
        dist = c("max.dist", "centroids.dist", "mahalanobis.dist")
    } else {
        nmthdd = length(dist)
    }
    
    if (length(validation) > 1 )
        validation = validation [1]
    if (!(validation %in% c("Mfold", "loo")))
        stop("Choose 'validation' among the two following possibilities: 'Mfold' or 'loo'")
    
    if (validation == "loo")
    {
        if (nrepeat != 1)
            warning("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
        nrepeat = 1
    }

    if (any(table(object$Y) <= 1)) {
        stop(paste("Cannot evaluate performance when a class level ('", 
                   names(table(object$Y))[which(table(object$Y) == 1)],
                   "') has only a single associated sample.", sep = ""))
    }
    
    if (!is.logical(progressBar))
        stop("'progressBar' must be either TRUE or FALSE")
    
    measure = c("overall","BER") # one of c("overall","BER")
    
    
    if (!(logratio %in% c("none", "CLR")))
        stop("Choose one of the two following logratio transformation: 'none' or 'CLR'")
    #fold is checked in 'MCVfold'
    
    #---------------------------------------------------------------------------#
    #-- logration + multilevel approach ----------------------------------------#

    # we can do logratio and multilevel on the whole data as these transformation are done per sample
    X = logratio.transfo(X = X, logratio = logratio)
    if (!is.null(multilevel))
    {
        Xw = withinVariation(X, design = multilevel)
        X = Xw
    }
    
    # -------------------------------------
    # first check for near zero var on the whole data set
    if (near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]
            
            if (ncol(X)==0)
                stop("No more predictors after Near Zero Var has been applied!")
            
            if (any(keepX > ncol(X)))
                keepX = ncol(X)
            
        }
    }
    # and then we start from the X data set with the nzv removed
    
    #---------------------------------------------------------------------------#
    #-- NA calculation      ----------------------------------------------------#
    
    misdata = c(X=anyNA(X), Y=FALSE) # Detection of missing data. we assume no missing values in the factor Y
    
    if (any(misdata))
    {
        is.na.A = is.na(X)
        
        #ind.NA = which(apply(is.na.A, 1, sum) > 0) # calculated only once
        #ind.NA.col = which(apply(is.na.A, 2, sum) > 0) # calculated only once
    } else {
        is.na.A = NULL
        #ind.NA = ind.NA.col = NULL
    }
    
    #-- Set up empty lists for outputs to go      ------------------------------#
    #---------------------------------------------------------------------------#
    list.features = list()
    
    mat.error.rate = mat.sd.error = mat.mean.error = error.per.class.keepX.opt = error.per.class.keepX.opt.mean = list()
    error.per.class = list()
    final=list()
    
    for (measure_i in measure)
    {
        mat.sd.error[[measure_i]] = matrix(0, nrow = 1, ncol = length(dist),
                                           dimnames = list(c(paste0('comp', ncomp)), dist))
        mat.mean.error[[measure_i]] = matrix(0,nrow = 1, ncol = length(dist),
                                             dimnames = list(c(paste0('comp', ncomp)), dist))
        error.per.class.keepX.opt[[measure_i]] = list()
        error.per.class.keepX.opt.mean[[measure_i]] = list()
        mat.error.rate[[measure_i]]=list()
        for(ijk in dist)
        {
            mat.error.rate[[measure_i]][[ijk]] = matrix(0, nrow = 1, ncol = nrepeat,
                                                        dimnames = list(c(paste0('comp', ncomp)), c(paste0('nrep', 1 : nrepeat))))
            
            error.per.class.keepX.opt[[measure_i]][[ijk]] <- array(0, c(nlevels(Y), nrepeat, 1), 
                                                                dimnames = list(levels(Y), paste0("nrep", 1:nrepeat), paste0('comp', ncomp)))
            
            error.per.class.keepX.opt.mean[[measure_i]][[ijk]] = matrix(nrow = nlevels(Y), ncol = 1,
                                                                        dimnames = list(c(levels(Y)), c(paste0('comp', ncomp))))
        }
    }
    
    if(auc == TRUE)
    {
        auc.mean=list()
        auc.all=list()
    }
    
    prediction.all = class.all = auc.mean = auc.all = list()
    for(ijk in dist)
    {
        class.all[[ijk]] = array(0, c(nrow(X),  nrepeat, 1),
                                 dimnames = list(rownames(X),c(paste0('nrep', 1 : nrepeat)),c(paste0('comp', ncomp))))
    }
    
    class.object=class(object)


    #-- Calculate error values - just on detected ncomp      -------------------#
    #---------------------------------------------------------------------------#

        if (progressBar == TRUE)
            cat("\ncomp", ncomp, "\n")
        
        
        if(ncomp > 1)
        {
            choice.keepX = keepX[1 : (ncomp - 1)]
        } else {
            choice.keepX = NULL
        }
        test.keepX = keepX[ncomp]
        names(test.keepX) = test.keepX
        #test.keepX is a value
        
        # estimate performance of the model
        result = MCVfold.spls (X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = ncomp,
                               choice.keepX = choice.keepX, test.keepX = test.keepX, test.keepY = nlevels(Y),
                               measure = measure, dist = dist, scale=scale,
                               near.zero.var = near.zero.var,
                               auc = auc, progressBar = progressBar, class.object = class.object,
                               BPPARAM = BPPARAM,
                               misdata = misdata, is.na.A = is.na.A)#, ind.NA = ind.NA, ind.NA.col = ind.NA.col)
        
        # ---- extract stability of features ----- # NEW
        if (is(object, "mixo_splsda"))
            list.features[[1]] = result$features$stable
        
        for (ijk in dist)
        {
            for (measure_i in measure)
            {
                mat.error.rate[[measure_i]][[ijk]] = result[[measure_i]]$mat.error.rate[[ijk]][1,]
                mat.mean.error[[measure_i]][, ijk]=result[[measure_i]]$error.rate.mean[[ijk]]
                if (!is.null(result[[measure_i]]$error.rate.sd))
                {
                    mat.sd.error[[measure_i]][, ijk]=result[[measure_i]]$error.rate.sd[[ijk]]
                } else {
                    mat.sd.error= NULL
                }
                # confusion matrix for keepX.opt, for each nrep
                error.per.class.keepX.opt[[measure_i]][[ijk]] = result[[measure_i]]$confusion[[ijk]]
                
                # confusion matrix for keepX.opt, averaged over all nrep
                error.per.class.keepX.opt.mean[[measure_i]][[ijk]] = apply(result[[measure_i]]$confusion[[ijk]],1 , mean)
            }
            
            #prediction of each samples for each fold and each repeat, on each comp
            class.all[[ijk]] = result$class.comp[[ijk]][,,1]
        }
        prediction.all[[1]] = array(unlist(result$prediction.comp),c(nrow(result$prediction.comp[[1]]), ncol(result$prediction.comp[[1]]), nrepeat),
                                       dimnames = c(dimnames(result$prediction.comp[[1]])[1:2], list(paste0("nrep",1:nrepeat))))#[[1]][, , 1] #take only one component [[1]] and one of test.keepX [,,1]
        
        if(auc == TRUE)
        {
            auc.all[[1]] = lapply(result$auc.all, function(x) x[,,1])
            auc.mean[[1]] = result$auc[, , 1]
        }
    
    names(prediction.all) = paste0("comp", ncomp)

    #-- Put together outputs     -----------------------------------------------#
    #---------------------------------------------------------------------------#
    
    result = list(error.rate = mat.mean.error,
                  error.rate.sd = mat.sd.error,
                  error.rate.all = mat.error.rate,
                  error.rate.class = error.per.class.keepX.opt.mean[[1]],
                  error.rate.class.all = error.per.class.keepX.opt[[1]],
                  predict = prediction.all,
                  class = class.all)
    
    if(auc)
    {
        names(auc.mean) = c(paste0("comp", ncomp))
        result$auc = auc.mean
        
        names(auc.all) = c(paste0("comp", ncomp))
        result$auc.all =auc.all
    }
    
    if (is(object, "mixo_splsda"))
    {
        names(list.features) = paste0("comp", ncomp)
        result$features$stable = list.features
    }
    
    if (progressBar == TRUE)
        cat('\n')
    
    # added
    if (near.zero.var == TRUE)
        result$nzvX = nzv$Position
    
    if (is(object, "mixo_splsda"))
    {
        method = "splsda.mthd"
    } else if (is(object, "mixo_plsda")) {
        method = "plsda.mthd"
    } else {
        warning("Something that should not happen happened. Please contact us.")
    }
    # edit class to 'perf' and remove 'perf.plsda.mthd' class so that plotting functionality is stopped
    class(result) = "perf"
    result$call = match.call()
    
    
    #updated outputs
    return(invisible(result))
}

#' @rdname perf.assess
#' @method perf.assess mixo_splsda
#' @export
perf.assess.mixo_splsda <- perf.assess.mixo_plsda

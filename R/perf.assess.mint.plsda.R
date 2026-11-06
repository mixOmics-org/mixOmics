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
# perf.assess.mint.plsda - Function to evaluate the performance of the fitted PLS (cross-validation)
#   inputs: object - object obtain from running mint.plsda
#           dist - to evaluate the classification performance
#           validation - type of validation
#           folds - number of folds if validation = "Mfold"
# ----------------------------------------------------------------------------------------------------------
#' ## -------------------------- perf.mint(s)plsda --------------------------- ##

#' @rdname perf.assess
#' @method perf.assess mint.plsda
#' @export
perf.assess.mint.plsda <- function (object,
                             dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
                             auc = FALSE,
                             progressBar = FALSE,
                             ...
)
{   #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #------------------#
    #-- check entries --#
    X = object$X
    Y = object$Y
    study = object$study
    ncomp = object$ncomp
    scale = object$scale
    
    keepX = apply(object$loadings$X, 2, function(x){sum(x!=0)})
    
    tol = object$tol
    max.iter = object$max.iter
    
    dist = match.arg(dist, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    if (any(dist == "all"))
    {
        nmthdd = 3
        dist = c("max.dist", "centroids.dist", "mahalanobis.dist")
    } else {
        nmthdd = length(dist)
    }
    
    if (!is.logical(progressBar))
        stop("'progressBar' must be either TRUE or FALSE")

    near.zero.var = !is.null(object$nzv) # if near.zero.var was used, we set it to TRUE. if not used, object$nzv is NULL

    if (any(table(object$Y) <= 1)) {
      stop(paste("Cannot evaluate performance when a class level ('", 
                names(table(object$Y))[which(table(object$Y) == 1)],
                "') has only a single associated sample.", sep = ""))
    }
    
    #-- end checking --#
    #------------------#
    
    # -------------------------------------
    # added: first check for near zero var on the whole data set
    if (near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]
            
            if (ncol(X)==0)
                stop("No more predictors after Near Zero Var has been applied!")
        }
    }
    # and then we start from the X data set with the nzv removed
    
    # Init data structures to store results
    prediction.all = class.all = list()
    for(ijk in dist)
    {
        class.all[[ijk]] = matrix(nrow = nrow(X), ncol = 1,
                                  dimnames = list(rownames(X), c(paste0('comp', ncomp))))
    }
    
    if(auc)
    {
        auc.mean=list()
        auc.mean.study=list()
    }
    
    study.specific = global = list()
    for (study_i in 1:nlevels(study)) #LOO on the study factor
    {
        study.specific[[study_i]] =list()
        study.specific[[study_i]]$BER = global$BER = matrix(0, nrow = 1, ncol = length(dist),
                                                            dimnames = list(c(paste0('comp', ncomp)), dist))
        
        study.specific[[study_i]]$overall = global$overall = matrix(0, nrow = 1, ncol = length(dist),
                                                                    dimnames = list(c(paste0('comp', ncomp)), dist))
        
        study.specific[[study_i]]$error.rate.class = list()
        for(ijk in dist)
            study.specific[[study_i]]$error.rate.class[[ijk]] = global$error.rate.class[[ijk]] = matrix(0, nrow = nlevels(Y), ncol = 1,
                                                                                                        dimnames = list(levels(Y), c(paste0('comp', ncomp))))
        
    }
    names(study.specific) = levels(study)
    
    # successively tune the components until ncomp: comp1, then comp2, ...
    for(comp in ncomp : ncomp)
    {
        
        already.tested.X = keepX[1:comp]
        
        if (progressBar == TRUE)
            cat("\ncomp",comp, "\n")
        
        #-- set up a progress bar --#
        if (progressBar ==  TRUE)
        {
            pb = txtProgressBar(style = 3)
        } else {
            pb = FALSE
        }
        
        M = nlevels(study)
        names.study = levels(study)
        features = NULL
        prediction.comp = matrix(0, nrow = nrow(X), ncol = nlevels(Y), dimnames = list(rownames(X), levels(Y)))
        
        class.comp = list()
        for(ijk in dist)
            class.comp[[ijk]] = matrix(0, nrow = nrow(X), ncol = 1)# prediction of all samples for each test.keepX and  nrep at comp fixed
        
        if(auc)
            auc.mean.study[[comp]] = list()
        
        for (study_i in 1:M) #LOO on the study factor
        {
            if (progressBar ==  TRUE)
                setTxtProgressBar(pb, (study_i-1)/M)
            
            omit = which(study %in% names.study[study_i])
            X.train = X[-omit,]
            Y.train = factor(Y[-omit])
            study.learn.CV = factor(as.character(study[-omit]))
            
            X.test = X[omit, , drop = FALSE] #note: drop is useless as there should always be more than a single sample in a study
            Y.test = Y[omit]
            study.test.CV = factor(as.character(study[omit]))
            
            #---------------------------------------#
            #-- near.zero.var ----------------------#
            if(near.zero.var == TRUE)
            {
                remove.zero = nearZeroVar(X.train)$Position
                
                if (length(remove.zero) > 0)
                {
                    X.train = X.train[, -c(remove.zero),drop = FALSE]
                    X.test = X.test[, -c(remove.zero),drop = FALSE]
                }
            }
            #-- near.zero.var ----------------------#
            #---------------------------------------#
            
            if (progressBar ==  TRUE)
                setTxtProgressBar(pb, (study_i-1)/M)
            
            object.res = mint.splsda(X.train, Y.train, study = study.learn.CV, ncomp = comp,
                                     keepX = already.tested.X,
                                     scale = scale)
            
            test.predict.sw <- predict.mixo_spls(object.res, newdata = X.test, dist = dist, study.test = study.test.CV)
            prediction.comp[omit, match(levels(Y.train),levels(Y))] =  test.predict.sw$predict[, , comp]
            
            for(ijk in dist)
                class.comp[[ijk]][omit,1] =  test.predict.sw$class[[ijk]][, comp] #levels(Y)[test.predict.sw$class[[ijk]][, ncomp]]
            
            
            if (progressBar ==  TRUE)
                setTxtProgressBar(pb, (study_i)/M)
            
            # result per study
            #BER
            study.specific[[study_i]]$BER[1,] = sapply(test.predict.sw$class, function(x){
            conf = get.confusion_matrix(truth = Y[omit], all.levels = levels(Y), predicted = x[,ncomp])
            get.BER(conf)
            })

            #overall
            study.specific[[study_i]]$overall[1,] = sapply(test.predict.sw$class, function(x){
            conf = get.confusion_matrix(truth = Y[omit], all.levels = levels(Y), predicted = x[,ncomp])
            out = sum(apply(conf, 1, sum) - diag(conf)) / length(Y[omit])
            })

            #classification for each level of Y
            temp = lapply(test.predict.sw$class, function(x){
            conf = get.confusion_matrix(truth = Y[omit], all.levels = levels(Y), predicted = x[,ncomp])
            out = (apply(conf, 1, sum) - diag(conf)) / summary(Y[omit])
            })
            for (ijk in dist)
            study.specific[[study_i]]$error.rate.class[[ijk]][,1] = temp[[ijk]]

            #AUC per study
            if(auc)
            {
            data = list()
            data$outcome = Y[omit]
            data$data = prediction.comp[omit, ]
            auc.mean.study[[comp]][[study_i]] = statauc(data)
            }

            # average of ER and BER across studies, weighted by study sample size
            global$BER[1,] <- global$BER[1,] + study.specific[[study_i]]$BER[1, ] * table(study)[study_i]/length(study)
            global$overall[1,] <- global$overall[1,] + study.specific[[study_i]]$overall[1, ] * table(study)[study_i]/length(study)
            
        } # end study_i 1:M (M folds)
        
        for (ijk in dist)
        {
            #prediction of each samples for each fold and each repeat, on each comp
            class.all[[ijk]][,1] = class.comp[[ijk]][,1]
        }
        prediction.all[[1]] = prediction.comp
        
        #classification for each level of Y
        temp = lapply(class.comp, function(x){
            conf = get.confusion_matrix(truth = factor(Y), predicted = x)
            out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
        })
        for (ijk in dist)
            global$error.rate.class[[ijk]][,1] = temp[[ijk]]
        
        #AUC global
        if(auc)
        {
            names(auc.mean.study[[comp]]) = names.study
            
            data = list()
            data$outcome = Y
            data$data = prediction.comp
            auc.mean[[comp]] = statauc(data)
        }
        
        
    } # end comp
    names(prediction.all) = paste0('comp', ncomp)
    
    result = list(study.specific.error = study.specific,
                  global.error = global,
                  predict = prediction.all,
                  class = class.all)
    
    if(auc)
    {
        names(auc.mean) = names(auc.mean.study) = paste0('comp', ncomp)
        result$auc = auc.mean
        result$auc.study = auc.mean.study
    }
    
    if (progressBar == TRUE)
        cat('\n')
    
    # added
    if (near.zero.var == TRUE)
        result$nzvX = nzv$Position
    
    class(result) = c("perf")
    result$call = match.call()
    
    return(invisible(result))
}

#' @rdname perf.assess
#' @method perf.assess mint.splsda
#' @export
perf.assess.mint.splsda <- perf.assess.mint.plsda
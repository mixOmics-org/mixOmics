################################################################################
# Author :
#   Florian Rohart,
#
# created: 04-07-2015
# last modified: 22-04-2016
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################


# ==============================================================================
# perform a Leave-one-out cross validation on the study to tune the number of
# variables (keepX) to keep in a mint.splsda analysis
# ==============================================================================


LOGOCV <- function(X,
                  Y,
                  ncomp,
                  study,
                  choice.keepX = NULL,
                  test.keepX = c(5, 10, 15),
                  dist = "max.dist",
                  measure = c("BER"), # one of c("overall","BER")
                  auc = auc,
                  max.iter = 100,
                  progressBar = FALSE,
                  near.zero.var = FALSE,
                  scale)
{
    # X input
    # Y factor input
    # ncomp: which component we are tuning
    # study: study effect
    # choice.keepX: how many variables are kept on the first ncomp-1 components
    # test.keepX: grid of keepX that is to be tested in the CV
    # dist= which distance should be used to classify the samples?
    # showProgress=TRUE, show the progress of the iteration
    
    
    if(missing(X))
        stop("missing X")
    
    if(missing(Y))
        stop("missing Y")
    
    if(missing(study))
        stop("missing study")
    
    if(missing(ncomp))
        ncomp = 1
    
    if(missing(scale))
        scale = FALSE
    
    if(missing(dist))
        dist = "max.dist"
    
    if(length(Y) != nrow(X))
        stop("X and Y have to be of same length")
    
    #-- set up a progress bar --#
    if (progressBar ==  TRUE)
    {
        pb = txtProgressBar(style = 3)
        nBar = 1
    } else {
        pb = FALSE
    }
    
    M = nlevels(study)
    names.study = levels(study)
    features = NULL
    prediction.comp = array(0, c(nrow(X), nlevels(Y), length(test.keepX)),
                            dimnames = list(rownames(X), levels(Y), test.keepX))
    
    class.comp = list()
    for(ijk in dist)
        class.comp[[ijk]] = matrix(0, nrow = nrow(X), ncol = length(test.keepX))
    # prediction of all samples for each test.keepX and  nrep at comp fixed
    
    PRED=INDICE = matrix(0, nrow = length(test.keepX), ncol = M)
    rownames(PRED) = rownames(INDICE) = test.keepX
    nbr.temp = matrix(0, nrow = length(test.keepX), ncol = nlevels(Y))
    
    for (study_i in 1:M) #LOO on the study factor
    {
        if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (study_i-1)/M)
        
        omit = which(study %in% names.study[study_i])
        X.train = X[-omit,]
        Y.train = factor(Y[-omit])
        study.learn.CV = factor(as.character(study[-omit]))
        
        X.test = X[omit, , drop = FALSE]
        #note: drop is useless as there should always be more than a single
        # sample in a study
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
        
        for (i in 1:length(test.keepX))
        {
            if (progressBar ==  TRUE)
                setTxtProgressBar(pb, (study_i-1)/M + (i-1)/length(test.keepX)/M)
            
            object.res = mint.splsda(X.train, Y.train,
                                                      study = study.learn.CV, ncomp = ncomp, keepX =
                                                          c(choice.keepX, test.keepX[i]),
                                                      scale = scale, max.iter = max.iter)
            
            # record selected features
            if (length(test.keepX) ==  1)
                # only done if only one test.keepX as not used if more so far
                features = c(features, selectVar(object.res, comp = ncomp)$name)
            
            test.predict.sw <- predict.mixo_spls(object.res, newdata = X.test,
                                                 method = dist, study.test = study.test.CV)
            # Y.train can be missing factors, so the prediction
            # 'test.predict.sw' might be missing factors compared to the
            # full prediction.comp
            prediction.comp[omit, match(levels(Y.train),levels(Y)) , i] =
                test.predict.sw$predict[, , ncomp]
            
            for(ijk in dist)
                class.comp[[ijk]][omit,i] =  test.predict.sw$class[[ijk]][, ncomp]
            #levels(Y)[test.predict.sw$class[[ijk]][, ncomp]]
        }#end test.keepX
        
        if (progressBar ==  TRUE)
            setTxtProgressBar(pb, (study_i)/M)
        
    } # end study_i 1:M (M folds)
    
    result = list()
    
    auc.mean = error.mean = error.sd = error.per.class.keepX.opt.comp =
        keepX.opt = test.keepX.out = choice.keepX.out =
        error.per.study.keepX.opt = list()
    
    if(auc)
    {
        data=list()
        for (i in 1:length(test.keepX))
        {
            
            data$outcome=Y
            data$data=prediction.comp[, , i]
            
            auc.mean[[i]]=statauc(data)
            
        }
        names(auc.mean)=test.keepX
    }
    
    if (any(measure == "overall"))
    {
        for(ijk in dist)
        {
            rownames(class.comp[[ijk]]) = rownames(X)
            colnames(class.comp[[ijk]]) = paste0("test.keepX.",test.keepX)
            
            #finding the best keepX depending on the error measure:
            #   overall or BER
            # classification error for each nrep and each test.keepX:
            #   summing over all samples
            error = apply(class.comp[[ijk]], 2, function(x)
            {
                sum(as.character(Y) != x)
            })
            
            # we divide the error by the number of samples and choose the
            #   minimum error
            error.mean[[ijk]] = error/length(Y)
            keepX.opt[[ijk]] = which(error.mean[[ijk]] ==
                                         min(error.mean[[ijk]]))[1]
            # chose the lowest keepX if several minimum
            
            # overall error per study
            temp = matrix(0, ncol = length(test.keepX), nrow = nlevels(study))
            for (study_i in 1:M) #LOO on the study factor
            {
                omit = which(study %in% names.study[study_i])
                temp[study_i,] = apply(class.comp[[ijk]][omit,], 2, function(x)
                {
                    sum(as.character(Y)[omit] != x)/length(omit)
                })
            }
            error.per.study.keepX.opt[[ijk]] = temp[,keepX.opt[[ijk]]]
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] =
                apply(class.comp[[ijk]][, keepX.opt[[ijk]], drop = FALSE], 2,
                      function(x)
                      {
                          conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                          out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
                      })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            
            choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
            
            result$"overall"$error.rate.mean = error.mean
            result$"overall"$error.per.study.keepX.opt =
                error.per.study.keepX.opt
            result$"overall"$confusion = error.per.class.keepX.opt.comp
            result$"overall"$keepX.opt = test.keepX.out
        }
    }
    
    if (any(measure ==  "BER"))
    {
        for(ijk in dist)
        {
            rownames(class.comp[[ijk]]) = rownames(X)
            colnames(class.comp[[ijk]]) = paste0("test.keepX.",test.keepX)
            
            # we calculate a BER for each study and each test.keepX,
            # then average over the study factor
            error = matrix(0, ncol = length(test.keepX), nrow = nlevels(study))
            for (study_i in 1:M) #LOO on the study factor
            {
                omit = which(study %in% names.study[study_i])
                error[study_i,] = apply(class.comp[[ijk]][omit,], 2, function(x)
                {
                    conf = get.confusion_matrix(truth = factor(Y[omit]),
                                                all.levels = levels(factor(Y)), predicted = x)
                    get.BER(conf)
                })
            }
            
            # weighted average BER over the studies
            error.mean[[ijk]] = apply(error, 2, function(x) {
                sum(x * table(study)/length(study))
            })
            
            
            keepX.opt[[ijk]] =
                which(error.mean[[ijk]] ==  min(error.mean[[ijk]]))[1]
            
            # error per study
            error.per.study.keepX.opt[[ijk]] = error[,keepX.opt[[ijk]]]
            
            # confusion matrix for keepX.opt
            error.per.class.keepX.opt.comp[[ijk]] =
                apply(class.comp[[ijk]][, keepX.opt[[ijk]], drop = FALSE], 2,
                      function(x)
                      {
                          conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                          out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
                      })
            
            rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
            
            
            test.keepX.out[[ijk]] = test.keepX[keepX.opt[[ijk]]]
            choice.keepX.out[[ijk]] = c(choice.keepX, test.keepX.out)
            
            result$"BER"$error.rate.mean = error.mean
            result$"BER"$error.per.study.keepX.opt = error.per.study.keepX.opt
            result$"BER"$confusion = error.per.class.keepX.opt.comp
            result$"BER"$keepX.opt = test.keepX.out
        }
        
    }
    
    result$prediction.comp = prediction.comp
    result$class.comp = class.comp
    result$auc = auc.mean
    result$features$stable =
        sort(table(as.factor(features))/M, decreasing = TRUE)
    return(result)
    
}#end function




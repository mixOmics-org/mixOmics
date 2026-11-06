################################################################################
# Authors:
#   Florian Rohart,
#   Kim-Anh Le Cao,
#   Francois Bartolo,
#
# created: 18-09-2017
# last modified: 05-10-2017
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################


# ==============================================================================
# tune.block.splsda: chose the optimal number of parameters per component on a
# block.splsda method
# ==============================================================================

# X: a list of data sets (called 'blocks') matching on the same samples.
#   Data in the list should be arranged in samples x variables, with
#   samples order matching in all data sets. \code{NA}s are not allowed.
# Y: a factor or a class vector for the discrete outcome.
# validation: Mfold or loo cross validation
# folds: if validation = Mfold, how many folds?
# nrepeat: number of replication of the Mfold process
# ncomp: the number of components to include in the model. Default to 1.
# choice.keepX: a list, each choice.keepX[[i]] is a vector giving keepX on
#   the components that were already tuned for block i
# test.keepX: a list of value(keepX) to test on the last component.
#   There needs to be names(test.keepX)
# measure: one of c("overall","BER"). Accuracy measure used in the cross
#   validation processs
# weighted: optimise the weighted or not-weighted prediction
# dist: distance to classify samples. see predict
# scheme: the input scheme, one of "horst", "factorial" or ""centroid".
#   Default to "centroid"
# design: the input design.
# init: intialisation of the algorithm, one of "svd" or "svd.single".
#   Default to "svd"
# tol: Convergence stopping value.
# max.iter: integer, the maximum number of iterations.
# near.zero.var: Logical, see the internal \code{\link{nearZeroVar}} function
#   (should be set to TRUE in particular for data with many zero values).
# progressBar: show progress,
# scale: Logical. If scale = TRUE, each block is standardized to zero means and
#   unit variances (default: TRUE).
# misdata: optional. any missing values in the data? list,
#   misdata[[q]] for each data set
# is.na.A: optional. where are the missing values? list,
#   is.na.A[[q]] for each data set (if misdata[[q]] == TRUE)
# ind.NA: optional. which rows have missing values? list,
#   ind.NA[[q]] for each data set.
# ind.NA.col: optional. which col have missing values? list,
#   ind.NA.col[[q]] for each data set.
# BPPARAM A \linkS4class{BiocParallelParam} object indicating the type
# of parallelisation.
#' @importFrom matrixStats colVars
MCVfold.block.splsda <- 
    function(
        X,
        Y,
        validation,
        folds,
        nrepeat = 1,
        ncomp,
        choice.keepX = NULL,
        test.keepX,
        measure = c("overall"),
        weighted = TRUE,
        dist = "max.dist",
        scheme,
        design,
        init,
        tol,
        max.iter = 100,
        near.zero.var = FALSE,
        progressBar = FALSE,
        scale,
        misdata,
        is.na.A,
        BPPARAM = SerialParam())
    {
        ## evaluate all args -- essential for SnowParam()
        mget(names(formals()), sys.frame(sys.nframe()))
        #-- checking general input parameters ------------------------------------#
        #--------------------------------------------------------------------------#
        #-- set up a progress bar --#
        if (progressBar ==  TRUE)
        {
            pb = txtProgressBar(style = 3)
        } else {
            pb = FALSE
        }
        test.keepA <- NULL
        #design = matrix(c(0,1,1,0), ncol = 2, nrow = 2, byrow = TRUE)
        
        if(ncomp>1)
        {
            keepY = rep(nlevels(Y), ncomp-1)
        } else {keepY = NULL}
        # prediction of all samples for each test.keepX and  nrep at comp fixed
        folds.input = folds
        
        n = nrow(X[[1]])
        repeated.measure = 1:n
        
        #-- define the folds --#
        if (validation ==  "Mfold")
        {
            
            if (is.null(folds) || !is.numeric(folds)) {
                stop("'folds' need to be non-NULL and numeric")
            } 
            else if (folds < 2) {
                stop("'folds' needs to be at least 2")  
            } 
            else if (folds > n) {
                stop("'folds' cannot be greater than the number of input samples") 
            } 
            else {
                M = round(folds)
            }
        } else if (validation ==  "loo") {
            M = n
            if(nrepeat != 1) { stop("nrepeat should be set to 1 with validation='loo'\n") }
        }
        
        all_folds <- lapply(seq_len(nrepeat), function(nrep) {
            #-- define the folds --#
            if (validation ==  "Mfold")
            {
                temp = stratified.subsampling(Y, folds = M)
                folds = temp$SAMPLE
                if(temp$stop > 0 & nrep == 1) # to show only once
                    warning("\nAt least one class is not represented in one fold,
                    which may unbalance the error rate.\n  Consider a number of
                    folds lower than the minimum in table(Y): ", min(table(Y)))
            } else if (validation ==  "loo") {
                folds = split(1:n, rep(1:n, length = n))
            }
            return(folds)
        })
       
        repeat_cv <- function(nrep, all_folds = all_folds, M = M)
        {
            folds = all_folds[[nrep]]
            class.comp.rep <- list()
            for(ijk in dist)
                class.comp.rep[[ijk]] = array(0, c(nrow(X[[1]]),
                                                   nrow(expand.grid(test.keepX))))
            # we don't record all the predictions for all fold and all blocks,
            #   too much data
            
            error.sw = matrix(0,nrow = M, ncol = length(test.keepX))
            rownames(error.sw) = paste0("fold",1:M)
            colnames(error.sw) = names(test.keepX)
            # for the last keepX (i) tested, prediction combined for all M folds so
            # as to extract the error rate per class
            # prediction.all = vector(length = nrow(X))
            # in case the test set only includes one sample, it is better to advise
            #the user to perform loocv
            # stop.user = FALSE
            
            #result.all=list()
            fonction.j.folds =function(j)#for(j in 1:M)
            {
                if (progressBar ==  TRUE)
                    setTxtProgressBar(pb, (M*(nrep-1)+j-1)/(M*nrepeat))
                
                #print(j)
                #set up leave out samples.
                omit = which(repeated.measure %in% folds[[j]] == TRUE)
                
                # get training and test set
                X.train = lapply(X, function(x){x[-omit, ]})
                Y.train = Y[-omit]
                Y.train.mat = unmap(Y.train)
                Q = ncol(Y.train.mat)
                colnames(Y.train.mat) = levels(Y.train)
                rownames(Y.train.mat) = rownames(X.train[[1]])
                X.test = lapply(X, function(x){x[omit, , drop = FALSE]})
                #matrix(X[omit, ], nrow = length(omit))
                #removed to keep the colnames in X.test
                Y.test = Y[omit]
                
                
                #---------------------------------------#
                #-- near.zero.var ----------------------#
                remove = vector("list",length=length(X))
                # first remove variables with no variance inside each X.train/X.test
                var.train = lapply(X.train, function(x){colVars(x, na.rm=TRUE)})
                for(q in 1:length(X))
                {
                    ind.var = which(var.train[[q]] == 0)
                    if (length(ind.var) > 0)
                    {
                        remove[[q]] =
                            c(remove[[q]], colnames(X.train[[q]])[ind.var])
                        
                        X.train[[q]] = X.train[[q]][, -c(ind.var),drop = FALSE]
                        X.test[[q]] = X.test[[q]][, -c(ind.var),drop = FALSE]
                        
                        # reduce choice.keepX and test.keepX if needed
                        if (any(choice.keepX[[q]] > ncol(X.train[[q]])))
                            choice.keepX[[q]][which(choice.keepX[[q]]>
                                                        ncol(X.train[[q]]))] = ncol(X.train[[q]])
                        
                        # reduce test.keepX if needed
                        if (any(test.keepX[[q]] > ncol(X.train[[q]])))
                            test.keepX[[q]][which(test.keepX[[q]]>
                                                      ncol(X.train[[q]]))] = ncol(X.train[[q]])
                        
                    }
                }
                
                # near zero var on X.train
                if(near.zero.var == TRUE)
                {
                    nzv.A = lapply(X.train, nearZeroVar)
                    for(q in 1:length(X.train))
                    {
                        if (length(nzv.A[[q]]$Position) > 0)
                        {
                            names.remove.X =
                                colnames(X.train[[q]])[nzv.A[[q]]$Position]
                            remove[[q]] = c(remove[[q]], names.remove.X)
                            
                            X.train[[q]] =
                                X.train[[q]][, -nzv.A[[q]]$Position, drop=FALSE]
                            X.test[[q]] =
                                X.test[[q]][, -nzv.A[[q]]$Position,drop = FALSE]
                            
                            if (ncol(X.train[[q]]) == 0)
                                stop(paste0("No more variables in",X.train[[q]]))
                            
                            #need to check that the keepA[[q]] is now not higher
                            # than ncol(A[[q]])
                            if (any(test.keepX[[q]] > ncol(X.train[[q]])))
                                test.keepX[[q]][which(test.keepX[[q]]>
                                                          ncol(X.train[[q]]))] = ncol(X.train[[q]])
                        }
                    }
                }
                
                #-- near.zero.var ----------------------#
                #---------------------------------------#
                
                
                #------------------------------------------#
                # split the NA in training and testing
                if(any(misdata))
                {
                    
                    is.na.A.train = is.na.A.test = ind.NA.train = ind.NA.col.train =
                        vector("list", length = length(X))
                    
                    for(q in 1:length(X))
                    {
                        if(misdata[q])
                        {
                            if(length(remove[[q]])>0){
                                ind.remove =
                                    which(colnames(X[[q]]) %in% remove[[q]])
                                is.na.A.train[[q]] =
                                    is.na.A[[q]][-omit, -ind.remove, drop=FALSE]
                                is.na.A.test[[q]] =
                                    is.na.A[[q]][omit, -ind.remove, drop=FALSE]
                            } else {
                                is.na.A.train[[q]] =
                                    is.na.A[[q]][-omit, , drop=FALSE]
                                is.na.A.test[[q]] = is.na.A[[q]][omit, , drop=FALSE]
                            }
                            temp = which(is.na.A.train[[q]], arr.ind=TRUE)
                            ind.NA.train[[q]] = unique(temp[,1])
                            ind.NA.col.train[[q]] = unique(temp[,2])
                        }
                    }
                    names(is.na.A.train) = names(is.na.A.test) =
                        names(ind.NA.train) = names(ind.NA.col.train) = names(is.na.A)
                    
                    
                    if(FALSE){
                        is.na.A.train = ind.NA.train = ind.NA.col.train =
                            vector("list", length = length(X))
                        
                        is.na.A.train= lapply(is.na.A, function(x)
                        {x[-omit,, drop=FALSE]})
                        is.na.A.test = lapply(is.na.A, function(x)
                        {x[omit,,drop=FALSE]})
                        for(q in 1:length(X))
                        {
                            if(misdata[q])
                            {
                                
                                temp = which(is.na.A.train[[q]], arr.ind=TRUE)
                                ind.NA.train[[q]] = unique(temp[,1])
                                ind.NA.col.train[[q]] = unique(temp[,2])
                            }
                        }
                    }
                } else {
                    is.na.A.train = is.na.A.test = NULL
                    ind.NA.train = NULL
                    ind.NA.col.train = NULL
                }
                
                # split the NA in training and testing
                #------------------------------------------#
                
                is.na.A.temp = ind.NA.temp = ind.NA.col.temp =
                    vector("list", length = length(X)+1)
                # inside wrapper.mint.block, X and Y are combined,
                # so the ind.NA need  length(X)+1
                is.na.A.temp[1:length(X)] = is.na.A.train
                ind.NA.temp[1:length(X)] = ind.NA.train
                ind.NA.col.temp[1:length(X)] = ind.NA.col.train
                
                # shape input for `internal_mint.block' (keepA, test.keepA, etc)
                result <- suppressMessages(
                    internal_wrapper.mint.block(
                        X=X.train, Y=Y.train.mat, study=factor(rep(1,length(Y.train))),
                        ncomp=ncomp, keepX=choice.keepX,
                        keepY=rep(ncol(Y.train.mat), ncomp-1),
                        test.keepX=test.keepX, test.keepY=ncol(Y.train.mat),
                        mode="regression", scale=scale, near.zero.var=near.zero.var,
                        design=design, max.iter=max.iter, scheme =scheme, init=init,
                        tol=tol, DA = TRUE,
                        misdata = misdata, is.na.A = is.na.A.temp, ind.NA = ind.NA.temp,
                        ind.NA.col = ind.NA.col.temp, all.outputs=FALSE))
                
                # `result' returns loadings and variates for all test.keepX on
                # the ncomp component
                
                # need to find the best keepX/keepY among all the tested models
                #save(list=ls(),file="temp.Rdata")
                
                # we prep the test set for the successive prediction: scale and
                # is.na.newdata
                
                # scale X.test
                ind.match = 1:length(X.train)# for missing blocks in predict.R
                if (!is.null(attr(result$A[[1]], "scaled:center")))
                    X.test[which(!is.na(ind.match))] = lapply(which(!is.na(ind.match)),
                                                              function(x){sweep(X.test[[x]], 2, STATS =
                                                                                    attr(result$A[[x]], "scaled:center"))})
                if (scale)
                    X.test[which(!is.na(ind.match))] = lapply(which(!is.na(ind.match)),
                                                              function(x){sweep(X.test[[x]], 2, FUN = "/", STATS =
                                                                                    attr(result$A[[x]], "scaled:scale"))})
                
                means.Y = matrix(attr(result$A[[result$indY]], "scaled:center"),
                                 nrow=nrow(X.test[[1]]),ncol=Q,byrow=TRUE);
                if (scale)
                {sigma.Y = matrix(attr(result$A[[result$indY]], "scaled:scale"),
                                  nrow=nrow(X.test[[1]]),ncol=Q,byrow=TRUE)
                
                }else{
                    sigma.Y=matrix(1,nrow=nrow(X.test[[1]]),ncol=Q)}
                
                names(X.test)=names(X.train)
                
                # record prediction results for each test.keepX
                keepA = result$keepA
                test.keepA = keepA[[ncomp]]
                
                class.comp.j = list()
                for(ijk in dist)
                    class.comp.j[[ijk]] =
                    matrix(0, nrow = length(omit), ncol = nrow(test.keepA))
                #prediction of all samples for each test.keepX and  nrep
                # at comp fixed
                
                
                # creates temporary block.splsda object to use the predict function
                class(result) =
                    c("block.splsda", "block.spls", "sgccda", "sgcca", "DA")
                
                result$X = result$A[-result$indY]
                result$ind.mat = result$A[result$indY][[1]]
                result$Y = factor(Y.train)
                
                #save variates and loadings for all test.keepA
                result.temp =
                    list(variates = result$variates, loadings = result$loadings)
                
                #time2 = proc.time()
                for(i in 1:nrow(test.keepA))
                {
                    #print(i)
                    
                    #only pick the loadings and variates relevant to that test.keepX
                    
                    names.to.pick = NULL
                    if(ncomp>1)
                        names.to.pick = unlist(lapply(1:(ncomp-1), function(x){
                            paste(paste0("comp",x),apply(keepA[[x]],1,function(x)
                                paste(x,collapse="_")), sep=":")
                            
                        }))
                    
                    names.to.pick.ncomp = paste(paste0("comp",ncomp),
                                                paste(as.numeric(keepA[[ncomp]][i,]),collapse="_"), sep=":")
                    names.to.pick = c(names.to.pick, names.to.pick.ncomp)
                    
                    
                    result$variates = lapply(result.temp$variates, function(x)
                    {if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick,
                                           drop=FALSE]}else{x}})
                    result$loadings = lapply(result.temp$loadings, function(x)
                    {if(ncol(x)!=ncomp) {x[,colnames(x)%in%names.to.pick,
                                           drop=FALSE]}else{x}})
                    
                    result$weights = get.weights(result$variates, indY=result$indY)
                    
                    # do the prediction, we are passing to the function some
                    # invisible parameters:the scaled newdata and the missing values
                    
                    test.predict.sw <- predict.block.spls(result,
                                                          newdata.scale = X.test, dist = dist, misdata.all=misdata,
                                                          is.na.X = is.na.A.train, is.na.newdata = is.na.A.test,
                                                          noAveragePredict=FALSE)
                    
                    if(weighted ==TRUE) #WeightedVote
                    {
                        for(ijk in dist)
                            class.comp.j[[ijk]][, i] =
                                test.predict.sw$WeightedVote[[ijk]][, ncomp]
                    } else {#MajorityVote
                        for(ijk in dist)
                            class.comp.j[[ijk]][, i] =
                                test.predict.sw$MajorityVote[[ijk]][, ncomp]
                    }
                } # end i
                
                return(list(
                    class.comp.j = class.comp.j, omit = omit, keepA = keepA))
            } # end fonction.j.folds
            ## if number of CPUs greater than nrepeat, also parallel on folds on non-Windows
            # excess_cpus <- ifelse(is(BPPARAM, 'MulticoreParam') && (BPPARAM$workers - nrepeat) > 0,
            #                       BPPARAM$workers - nrepeat,
            #                       0)
            # if (excess_cpus > 0 & do_excess) {
            #     result.all = bplapply(1:M, fonction.j.folds, BPPARAM = MulticoreParam(workers = excess_cpus))
            # } else {
            #TODO function.j.folds and repeat_cv should be modified to output tidy data.frames
            result.all = lapply(1:M, fonction.j.folds) ## too much overhead if we also parallelise folds
            # }
            
            keepA <- result.all[[1]]$keepA[[ncomp]]
            # combine the results
            for(j in 1:M)
            {
                omit = result.all[[j]]$omit
                #prediction.comp.j = result[[j]]$prediction.comp.j
                class.comp.j = result.all[[j]]$class.comp.j
                
                #prediction.comp[[nrep]][omit, , ] = prediction.comp.j
                for(ijk in dist) {
                    class.comp.rep[[ijk]][omit, ] = class.comp.j[[ijk]]
                }
                
                
            }
            
            if (progressBar ==  TRUE)
                setTxtProgressBar(pb, nrep/nrepeat)
            
            return(list(class.comp.rep=class.comp.rep, keepA=keepA))
        } #end nrep 1:nrepeat
        
        class.comp.reps <- bplapply(seq_len(nrepeat), repeat_cv, all_folds = all_folds, M = M, BPPARAM = BPPARAM)
 
        list2array <- function(cc) {
            ## function to make an array of results of all repeats ino the former form
            ## before nrepeat loop becomes a function
            dims <- dim(cc[[1]][[1]][[1]])
            dist.array <- list()
            for (dist in names(cc[[1]][[1]])) {
                dist.array[[dist]] <- array(0, dim = c(dims[1], length(cc), dims[2]))
                for (rep in seq_len(length(cc))) {
                    dist.array[[dist]][,rep,] <- cc[[rep]][[1]][[dist]]
                }
            }
            return(dist.array)
        }
        
        class.comp <- list2array(class.comp.reps)
        #names(prediction.comp) =
        # class.comp[[ijk]] is a matrix containing all prediction for test.keepX,
        # all nrepeat and all distance, at comp fixed
        test.keepA <- class.comp.reps[[1]][["keepA"]]
        keepA.names = apply(test.keepA[,seq_along(X), drop=FALSE],1,function(x)
            paste(x,collapse="_"))#, sep=":")
        
        result = list()
        error.mean = error.sd = error.per.class.keepX.opt.comp = keepX.opt =
            test.keepX.out = mat.error.final = list()
        #save(list=ls(), file="temp22.Rdata")
        
        if (any(measure == "overall"))
        {
            for(ijk in dist)
            {
                rownames(class.comp[[ijk]]) = rownames(X)
                colnames(class.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
                dimnames(class.comp[[ijk]])[[3]] = keepA.names
                
                #finding the best keepX depending on the error measure:
                # overall or BER
                # classification error for each nrep and each test.keepX:
                # summing over all samples
                error = apply(class.comp[[ijk]],c(3,2),function(x)
                {
                    length(Y) - sum(as.character(Y) == x, na.rm=TRUE)
                })
                
                # we want to average the error per keepX over nrepeat and choose
                # the minimum error
                error.mean[[ijk]] = apply(error,1,mean)/length(Y)
                if (!nrepeat ==  1)
                    error.sd[[ijk]] = apply(error,1,sd)/length(Y)
                
                mat.error.final[[ijk]] = error/length(Y)
                # percentage of misclassification error for each test.keepX (rows)
                # and each nrepeat (columns)
                
                
                min.error = min(error.mean[[ijk]])
                min.keepX = rownames(error)[which(error.mean[[ijk]] == min.error)]
                # vector of all keepX combination that gives the minimum error
                
                a = lapply(min.keepX, function(x)
                {as.numeric(strsplit(x, "_")[[1]][1:length(X)])})
                
                #transform keepX in percentage of variable per dataset, so we choose
                # the minimal overall
                p = sapply(X,ncol)
                percent = sapply(a, function(x) sum(x/p))
                ind.opt = which.min(percent) # we take only one
                a = a[[ind.opt]]# vector of each optimal keepX for all block on
                # component comp.real[comp]
                
                # best keepX
                opt.keepX.comp = as.list(a)
                names(opt.keepX.comp) = names(X)
                
                choice.keepX = lapply(1:length(X), function(x)
                {c(choice.keepX[[x]],opt.keepX.comp[[x]])})
                names(choice.keepX) = names(X)
                
                keepX.opt[[ijk]] = which(error.mean[[ijk]] == min.error)[ind.opt]
                
                
                # confusion matrix for keepX.opt
                error.per.class.keepX.opt.comp[[ijk]] =apply(class.comp[[ijk]][, ,
                                                                               keepX.opt[[ijk]], drop = FALSE], 2, function(x)
                                                                               {
                                                                                   conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                                                                                   out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
                                                                               })
                
                test.keepX.out[[ijk]] = keepA.names[keepX.opt[[ijk]]]
                #strsplit(keepA.names[keepX.opt[[ijk]]],":")[[1]][2]
                # single entry of the keepX for each block
                
                result$"overall"$error.rate.mean = error.mean
                if (!nrepeat ==  1)
                    result$"overall"$error.rate.sd = error.sd
                
                result$"overall"$confusion = error.per.class.keepX.opt.comp
                result$"overall"$mat.error.rate = mat.error.final
                result$"overall"$ind.keepX.opt = keepX.opt
                result$"overall"$keepX.opt = test.keepX.out
                result$"overall"$choice.keepX = choice.keepX
                
            }
        }
        
        if (any(measure ==  "BER"))
        {
            for(ijk in dist)
            {
                rownames(class.comp[[ijk]]) = rownames(X[[1]])
                colnames(class.comp[[ijk]]) = paste0("nrep.", 1:nrepeat)
                dimnames(class.comp[[ijk]])[[3]] = keepA.names
                
                error = apply(class.comp[[ijk]],c(3,2),function(x)
                {
                    conf = get.confusion_matrix(truth = factor(Y),predicted = x)
                    get.BER(conf)
                })
                rownames(error) = keepA.names
                colnames(error) = paste0("nrep.",1:nrepeat)
                
                # average BER over the nrepeat
                error.mean[[ijk]] = apply(error,1,mean)
                if (!nrepeat ==  1)
                    error.sd[[ijk]] = apply(error,1,sd)
                
                mat.error.final[[ijk]] = error
                # BER for each test.keepX (rows) and each nrepeat (columns)
                
                
                min.error = min(error.mean[[ijk]])
                min.keepX = rownames(error)[which(error.mean[[ijk]] == min.error)]
                # vector of all keepX combination that gives the minimum error
                
                a = lapply(min.keepX, function(x)
                {as.numeric(strsplit(x, "_")[[1]][1:length(X)])})
                
                #transform keepX in percentage of variable per dataset, so we
                # choose the minimal overall
                p = sapply(X,ncol)
                percent = sapply(a, function(x) sum(x/p))
                ind.opt = which.min(percent) # we take only one
                a = a[[ind.opt]]# vector of each optimal keepX for all block on
                # component comp.real[comp]
                
                # best keepX
                opt.keepX.comp = as.list(a)
                names(opt.keepX.comp) = names(X)
                
                choice.keepX = lapply(1:length(X), function(x)
                {c(choice.keepX[[x]],opt.keepX.comp[[x]])})
                names(choice.keepX) = names(X)
                
                keepX.opt[[ijk]] = which(error.mean[[ijk]] == min.error)[ind.opt]
                
                # confusion matrix for keepX.opt
                error.per.class.keepX.opt.comp[[ijk]] = apply(class.comp[[ijk]][, ,
                                                                                keepX.opt[[ijk]], drop = FALSE], 2, function(x)
                                                                                {
                                                                                    conf = get.confusion_matrix(truth = factor(Y), predicted = x)
                                                                                    out = (apply(conf, 1, sum) - diag(conf)) / summary(Y)
                                                                                })
                
                rownames(error.per.class.keepX.opt.comp[[ijk]]) = levels(Y)
                colnames(error.per.class.keepX.opt.comp[[ijk]]) =
                    paste0("nrep.", 1:nrepeat)
                
                test.keepX.out[[ijk]] = keepA.names[keepX.opt[[ijk]]]
                
                result$"BER"$error.rate.mean = error.mean
                if (!nrepeat ==  1)
                    result$"BER"$error.rate.sd = error.sd
                
                result$"BER"$confusion = error.per.class.keepX.opt.comp
                result$"BER"$mat.error.rate = mat.error.final
                result$"BER"$ind.keepX.opt = keepX.opt
                result$"BER"$keepX.opt = test.keepX.out
                result$"BER"$choice.keepX = choice.keepX
            }
            
            
        }
        
        #result$prediction.comp = prediction.comp
        result$class.comp = class.comp
        return(result)
    }

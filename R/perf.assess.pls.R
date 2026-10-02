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
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# ----------------------------------------------------------------------------------------------------------
# perf.assess.pls - Function to evaluate the performance of fitted PLS (cross-validation)
# ----------------------------------------------------------------------------------------------------------

## -------------------------------- (s)PLS -------------------------------- ##
#' @rdname perf.assess
#' @method perf.assess mixo_pls
#' @export
perf.assess.mixo_pls <- function(object,
                          validation = c("Mfold", "loo"),
                          folds,
                          nrepeat = 1,
                          progressBar = FALSE,
                          BPPARAM = SerialParam(),
                          seed = NULL,
                          ...)
{

    # checking args and initialize params
    ncomp = object$ncomp
    spls.model <- is(object, 'mixo_spls')
    progressBar <- .check_logical(progressBar)
    BPPARAM$RNGseed <- seed
    
    # run CV in parallel depending on BPPARAM
    repeat_names <- .name_list(char = seq_len(nrepeat))
    result <- bplapply(X = repeat_names, FUN = function(nrep) {
        ## progress bar
        if (progressBar == TRUE){.progressBar(nrep/nrepeat)}
        ## CV
        .perf.assess.mixo_pls_cv(object, validation = validation, folds = folds, nrep = nrep)
    }, BPPARAM = BPPARAM)
    
    # extract measures
    measures <- lapply(result, function(x){x$measures})
    measures <- Reduce(rbind, measures)
    measures <- as.data.frame(measures)

    # Add this line to remove rows with NAs that correspond to components < ncomp
    measures <- dplyr::filter(measures, .data$comp == ncomp)
    
    ## R CMD check stuff
    measure <- feature <- comp <- block <- stability <- value <- NULL
    lower <- upper <- keepX <- keepY <- NULL
    
    # extract measures
    measure.names <- .name_list(unique(measures$measure))
    measures <- lapply(measure.names, function(meas) {
        ## ------ value of measures across repeats
        df <- measures %>% 
            filter(measure == meas) %>% 
            mutate(measure = NULL) %>% 
            as.data.frame()
        
        ## ------ summary of measures across repeats
        df.summ <- df %>%  
            group_by(feature, comp) %>% 
            summarise(mean = mean(value, na.rm = TRUE), 
                      sd = sd(value, na.rm = TRUE)) %>% 
            as.data.frame()
        
        list(values = df, summary = df.summ)
    })

    ## ------ feature stability
    if (spls.model)
    {
        features <- lapply(result, function(x){
            x$features
        })
        
        features <- Reduce(rbind, features) %>% 
            group_by(feature, comp, block) %>% 
            summarise(stability = mean(stability, na.rm = TRUE))
        
        features <- as.data.frame(features)
        features <- lapply(list(stability.X = 'X', stability.Y = 'Y'), function(z){
            lapply(.name_list(unique(features$comp)), function(n.comp){
                
                df <- features
                df <- filter(df, block == z & comp == n.comp)
                df <- df[,c('feature', 'stability')]
                vec <- df$stability
                names(vec) <- df$feature
                sort(vec, decreasing = TRUE)
            })
        })
    } else
    {
        features <- NULL
    }
    
    result <- list(measures = measures,
                   features = features)
    mc <- mget(names(formals())[-1], sys.frame(sys.nframe()))
    ## replace function, object with unevaluated call
    mc <- as.call(c(as.list(match.call())[1:2], mc))
    result <- c(list(call = mc), result)
    
    # change this so cant plot the output
    class(result) <- "perf"
    
    return(result)
    
}

#' @rdname perf.assess
#' @method perf.assess mixo_spls
#' @export
perf.assess.mixo_spls  <- perf.assess.mixo_pls

#' @noRd
#' @keywords Internal
.perf.assess.mixo_pls_cv <- function(object,
                              validation = c("Mfold", "loo"),
                              folds,
                              nrep = 1,
                              ...)

{
# changes to bypass the loop for the Q2

    ## R CMD check stuff
    measure <- feature <- comp <- block <- stability <- value <- NULL
    lower <- upper <- keepX <- keepY <- NULL
    
    ## -------- checks -------- ##
    if (object$mode == 'invariant')
        stop("'perf' is only available for (s)pls with modes: 'regression', 'canonical' or 'classic'.  Object has mode 'invariant'", call. = FALSE)
    
    validation = match.arg(validation)
    
    ## ---------- CV ---------- ##
    ## ------------- initialise
    # these are the centered and scaled matrices output from pls, we remove $nzv if needed
    if (length(object$nzv$Position)>0)
    {
        X = object$X[, -object$nzv$Position]
    } else {
        X = object$X
    }
    Y = object$Y
    
    scale = object$scale
    tol = object$tol
    max.iter = object$max.iter
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    
    if (any(is.na(X)) || any(is.na(Y)))
        stop("missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.", call. = FALSE)
    
    
    #-- tells which variables are selected in X and in Y --#
    if (is(object, "mixo_spls"))
    {
        keepX = object$keepX
        keepY = object$keepY
    } else {
        keepX = rep(ncol(X), ncomp)
        keepY = rep(ncol(Y), ncomp)
    }
    
    #-- define the folds --#
    if (validation == "Mfold")
    {
        if (is.list(folds))
        {
            # folds supplied as a list: validated by the shared helper
            folds = .check_folds(folds, n = n)
            M = length(folds)
        } else {
            if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n)
            {
                stop("Invalid number of folds.", call. = FALSE)
            } else {
                M = round(folds)
                folds = split(sample(1:n), rep(1:M, length = n))
            }
        }
    } else {
        folds = split(1:n, rep(1:n, length = n))
        M = n
    }
    
    #-- initialize new objects --#
    if (mode == 'canonical'){
        RSS = rbind(rep(n - 1, p), matrix(nrow = ncomp, ncol = p))
        # RSS.indiv is the reconstructed matrix X
        #RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = p)})
        #RSS.indiv[[1]] = X
        press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = p)})
        PRESS.inside = Q2 = matrix(nrow = ncomp, ncol = p)
    }else{
        RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
        # RSS.indiv is the reconstructed matrix Y
        #RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
        #RSS.indiv[[1]] = Y # KA changed
        press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = q)})
        PRESS.inside = Q2 = matrix(nrow = ncomp, ncol = q)
    }
    
    MSEP.mat = Ypred = array(0, c(n, q, ncomp))
    MSEP = R2 = matrix(nrow = ncomp, ncol = q)
    
    # to store the predicted components
    t.pred.cv = matrix(nrow = nrow(X), ncol = ncomp)
    u.pred.cv = matrix(nrow = nrow(X), ncol = ncomp)
    
    # to record feature stability, a list of form
    # list(X = list(comp1 = c(feature1 = 0.99, ...), 
    #               comp2 = c(feature2 = 0.98, ...)), 
    #      Y = ...)
    features <-
        lapply(list(X = X, Y = Y), function(Z){
            features <- vector(mode = 'numeric', length = ncol(Z))
            names(features) <- colnames(Z)
            features <- lapply(seq_len(ncomp), function(x) features)
            names(features) <- paste0('comp', seq_len(ncomp))
            
            return(features)
        })
    
    
    # ====  loop on h = ncomp is only for the calculation of Q2 on each component
    ## loop adds data to new row in RSS 

    for (h in 1:ncomp) # this loop needs to run on all components because Q2 is calculated using the RSS of ncomp-1
    {
        #-- initialising arguments --#
        tt = object$variates$X[, h]
        u = object$variates$Y[, h]
        b = object$loadings$Y[, h]
        #nx = p - keepX[h]
        #ny = q - keepY[h]
        
        # only used for matrices deflation across dimensions
        c = crossprod(X, tt)/drop(crossprod(tt))  #object$mat.c[, h]
        d = crossprod(Y, tt)/drop(crossprod(tt))  #object$mat.d[, h]
        e = crossprod(Y, u)/drop(crossprod(u))    
        
        # deflate matrices
        X = X - tt %*% t(c)
        
        #-- mode classic
        if (mode == "classic")
            Y = Y - tt %*% t(b)
        #-- mode regression
        if (mode == "regression")
            Y = Y - tt %*% t(d)
        #-- mode canonical 
        if (mode == "canonical")
            Y = Y - u %*% t(e)
        #-- mode invariant: Y is unchanged
        
        # update RSS for X/Y deflated
        if(mode == 'canonical'){  # based on X
            RSS[h + 1, ] =  colSums((X)^2)   # ==  colSums((X - tt %*% t(c))^2) if we had not deflated
        }else{ # regression, invariant, classic
            RSS[h + 1, ] = colSums((Y)^2)  # 
        }
        
    } # end h to calculate RSS   
    
    
    
    # ======== loop on i for cross validation ===================#
    for (i in 1:M) # M is number of folds
    {
        # initialise the train / test datasets
        omit = folds[[i]]
        X.train = object$X[-omit, , drop = FALSE]
        Y.train = object$Y[-omit, , drop = FALSE]
        X.test = object$X[omit, , drop = FALSE]
        Y.test = object$Y[omit, , drop = FALSE]
        
        # New loop to calculate prediction - theoretically should be able to just calculate ncomp directly but this gives a different result for whatever reason
        # So just keeping the loop but exacting rows corresponding to ncomp at the end
        for (h in 1:ncomp)
        { 
            #-- for MSEP and R2 criteria, no loop on the component as we do a spls with ncomp
            ##if (h == 1)
            #{
            #nzv = (apply(X.train, 2, var) > .Machine$double.eps) # removed in v6.0.0 so that MSEP, R2 and Q2 are obtained with the same data
            # re-added in >6.1.3 to remove constant variables
            nzv.X = (apply(X.train, 2, var) > .Machine$double.eps)
            nzv.Y = (apply(Y.train, 2, var) > .Machine$double.eps)
            
            # creating a keepX/Y.temp that can change for each fold, depending on nzv.X/Y
            keepX.temp = keepX
            keepY.temp = keepY
            if(any(keepX.temp > sum(nzv.X)))
                keepX.temp[which(keepX.temp>sum(nzv.X))] = sum(nzv.X)
            if(any(keepY.temp > sum(nzv.Y)))
                keepY.temp[which(keepY.temp>sum(nzv.Y))] = sum(nzv.Y)
            # TODO clarify the iterative nzv process in docs -- give it a better name (these are actually !nzv)
            # here h = 1 because we deflate at each step then extract the vectors for each h comp
            spls.res = spls(X.train[, nzv.X, drop = FALSE], Y.train[, nzv.Y, drop = FALSE], ncomp = 1, mode = mode, max.iter = max.iter, tol = tol, 
                            keepX = keepX.temp[h], keepY = keepY.temp[h], near.zero.var = FALSE, scale = scale)
            Y.hat = predict.mixo_spls(spls.res, X.test[, nzv.X, drop = FALSE])$predict
            
            # added the stop msg
            if(sum(is.na(Y.hat))>0) stop('Predicted Y values include NA')  
            
            # replaced h by 1; Y.hat is the prediction of the test samples for all q variable in comp h = 1
            Ypred[omit, nzv.Y, h] = Y.hat[, , 1]
            MSEP.mat[omit, nzv.Y, h] = (Y.test[, nzv.Y] - Y.hat[, , 1])^2
            
            
            # Q2 criterion: buidling directly from spls object
            u.cv = spls.res$variates$Y[, 1]
            t.cv = spls.res$variates$X[, 1]
            a.cv = spls.res$loadings$X[, 1]
            b.cv = spls.res$loadings$Y[, 1, drop = FALSE]
            
            # reg coefficients:
            c.cv = crossprod(X.train, u.cv) / drop(crossprod(u.cv)) 
            d.cv = crossprod(Y.train, t.cv) / drop(crossprod(t.cv)) # d.cv \neq to b.cv as d.cv is normed wrt to t.cv
            e.cv = crossprod(Y.train, u.cv) / drop(crossprod(u.cv)) 
            
            # calculate predicted components and store
            t.pred = c(X.test %*% a.cv)
            t.pred.cv[omit,h] = t.pred    # needed for tuning
            b.pred = crossprod(Y.test, t.pred)
            b.pred.cv = b.pred/ drop(sqrt(crossprod(b.pred)))
            u.pred.cv[omit,h] = Y.test[, nzv.Y] %*% b.cv  # needed for tuning, changed instead of b.pred.cv
            
            # predicted reg coeff, could be removed
            e.pred.cv = crossprod(as.matrix(Y.test), Y.test %*% b.pred.cv) / drop(crossprod(Y.test %*% b.pred))
            d.pred.cv = crossprod(as.matrix(Y.test), t.pred) / drop(crossprod(t.pred)) 
            
            # deflate matrices X
            X.train = X.train - t.cv %*% t(c.cv)
            X.test = X.test - t.pred %*% t(c.cv)
            # deflate matrices X      
            #-- mode classic
            if (mode == "classic"){
                Y.train[, nzv.Y] = Y.train[, nzv.Y] - t.cv %*% t(b.cv)  # could be pred on b
                Y.test[, nzv.Y] = Y.test[, nzv.Y] - t.pred %*% t(b.cv)
            }
            #-- mode regression
            if (mode == "regression"){
                Y.train = Y.train - t.cv %*% t(d.cv) # could be pred d.pred.cv? does not decrease enough
                Y.test[, nzv.Y] = Y.test[, nzv.Y] - Y.hat[, , 1]   # == Y.test - t.pred %*% t(d.cv) 
            }
            
            #-- mode canonical  ## KA added
            if (mode == "canonical"){
                Y.train = Y.train - u.cv %*% t(e.cv)  # could be pred on e
                Y.test = Y.test - (Y.test[, nzv.Y] %*% b.cv) %*% t(e.cv)  # here u.pred = Y.test %*% b.cv (b.pred.cv gives similar results)
            }
            #-- mode invariant: Y is unchanged
            
            # calculate predicted matrix X.hat or Y.hat based on X.test
            if(mode == 'canonical'){  # Xa c' = t c'
                #X.hat.cv = t.pred %*% t(c.cv), calculated earlier
                press.mat[[h]][omit, ] = X.test        # == X.test - X.hat.cv
            }else{ #  if(mode == 'regression'){  # Xa d' = t d'
                #Y.hat.cv = t.pred %*% t(d.cv), calculated earlier
                press.mat[[h]][omit, ] = Y.test        # == Y.test - Y.hat.cv
            }  
            
            # Record selected features in each set
            if (is(object,"mixo_spls"))
            {
                X.feature <- as.numeric(names(features$X[[h]]) %in% selectVar(spls.res, comp = 1)$X$name)
                Y.feature <- as.numeric(names(features$Y[[h]]) %in% selectVar(spls.res, comp = 1)$Y$name)
                # TODO using comp = 1 after deflation: this is problematic if, say, folds = 3, keepX = c(2, 100) (max 4 features (2 folds x 2 features) should be output for comp2 before calculating stability)
                features$X[[h]] <- features$X[[h]] + X.feature / length(folds)
                features$Y[[h]] <- features$Y[[h]] + Y.feature / length(folds)
            }
            
        } #  end loop on h ncomp
    } # end i (cross validation)
    
    
    # store results for each comp
    for (h in ncomp:ncomp){
        #-- compute the Q2 criterion --#
        # norm is equivalent to summing here the squared press values:
        PRESS.inside[h, ] = apply(press.mat[[h]], 2, function(x){norm(x, type = "2")^2})
        
        if(mode != 'canonical'){
            Q2[h, ] = 1 - PRESS.inside[h, ] / RSS[h, ] # note that RSS has an extra row so here Q is being calculated with PRESS of h and RSS of h-1
            MSEP[h, ] = apply(as.matrix(MSEP.mat[, , h]), 2, mean)
            R2[h, ] = (diag(cor(object$Y, Ypred[, , h])))^2
        } # if mode == canonical, do not output

    }
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]),
                      nrow = 1, ncol = ncomp,
                      dimnames = list("Q2.total", paste0("comp", seq_len(ncomp))))
    
    # set up dimnames and outputs
    result = list()
    
    if(mode != 'canonical'){
        rownames(MSEP) = rownames(R2) = rownames(Q2) = paste0("comp", seq_len(ncomp))
        colnames(MSEP) = colnames(R2) = colnames(Q2) = object$names$colnames$Y
        
        result$MSEP = t(MSEP)
        result$RMSEP = sqrt(t(MSEP))
        #result$MSEP.mat = MSEP.mat  
        result$R2 = t(R2)
        result$Q2 = t(Q2)  # remove this output when canonical mode?
    }
    
    result$Q2.total =  Q2.total
    RSS <- t(RSS) ## bc all others are transposed
    PRESS = t(PRESS.inside)
    result$RSS <- RSS[,-1, drop = FALSE] ## drop q/p
    result$PRESS <- PRESS
    if (ncol(object$Y) > 1)
    {
        # TODO ensure these are in fact no more necessary
        #result$d.cv = d.cv  # KA added  
        #result$b.cv = b.cv  # KA added 
        #result$c.cv = c.cv  # KA added 
        #result$u.cv = u.cv  # KA added 
        #result$a.cv = a.cv  # KA added 
        #result$t.pred.cv = t.pred.cv  # needed for tuning
        #result$u.pred.cv = u.pred.cv  # needed for tuning
        
        # extract the predicted components per dimension, take abs value
        result$cor.tpred = diag(abs(cor(t.pred.cv, object$variates$X)))
        result$cor.tpred = t(data.matrix(result$cor.tpred, rownames.force = TRUE))
        result$cor.upred = diag(abs(cor(u.pred.cv, object$variates$Y)))
        result$cor.upred = t(data.matrix(result$cor.upred, rownames.force = TRUE))
        
        # RSS: no abs values here
        result$RSS.tpred = apply((t.pred.cv - object$variates$X)^2, 2, sum)/(nrow(X) -1)
        result$RSS.tpred  = t(data.matrix(result$RSS.tpred, rownames.force = TRUE))
        result$RSS.upred = apply((u.pred.cv - object$variates$Y)^2, 2, sum)/(nrow(X) -1)
        result$RSS.upred  = t(data.matrix(result$RSS.upred, rownames.force = TRUE))
    }
    result <- mapply(result, names(result), FUN = function(arr, measure) {
        arr <- data.matrix(arr)
        col.names <- seq_len(ncomp)
        if (ncol(arr) == ncomp)
            colnames(arr) <- col.names
        else
            stop("unexpected dimension for entry in perf measures: ", measure)
        if (nrow(arr) == 1)
            rownames(arr) <- measure
        arr
    }, SIMPLIFY = FALSE)
    
    ## melt by comp
    result <- lapply(result, FUN = function(arr, nrep) {
        arr <- melt(arr)
        colnames(arr) <- c('feature', 'comp', 'value')
        if (nlevels(arr$feature) == 1) ## for Y-level measures (ass opossed to Y_feature level) such as Q2.total
            arr$feature <- factor('Y')
        arr$nrep <- nrep
        arr
    }, nrep = nrep)
    col.names <- names(result[[1]])
    #' @importFrom reshape2 melt
    result <- melt(result, id.vars = col.names)
    colnames(result) <- c(col.names, 'measure')
    
    result <- list(measures = result)

    #---- stability of features -----#
    if (is(object, "mixo_spls"))
    {
        features <- lapply(features, function(x)
        {
            x <- lapply(x, function(stab) round(stab, 2))
            df <- data.frame(x)
            df <- data.frame(feature = rownames(df), df)
            df
        })
        features <- melt(features, id.vars = 'feature', value.name = 'stability', variable.name = 'comp')
        ## add block name column instead of default 'L1'
        colnames(features) <- c(rev(rev(colnames(features))[-1]), 'block')
        features$nrep <- nrep

        # filter features to only include those based on ncomp
        features <- dplyr::filter(features, comp == paste0("comp", ncomp))
        
        result$features <- features
    }
    return(invisible(result))
}
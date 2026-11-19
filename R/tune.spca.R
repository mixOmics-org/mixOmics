# This file is part of mixOmics.
#
# mixOmics is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, version 3.
#
# mixOmics is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero
# General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with mixOmics. If not, see <https://www.gnu.org/licenses/>.
#
# SPDX-License-Identifier: AGPL-3.0-only

#' Tune number of selected variables for spca
#'
#' This function performs sparse pca and optimises the number of variables to
#' keep on each component using repeated cross-validation.
#'
#' Essentially, for the first component, and for a grid of the number of
#' variables to select (\code{keepX}), a number of repeats and folds, data are
#' split to train and test and the extracted components are compared against
#' those from a spca model with all the data to ascertain the optimal
#' \code{keepX}. In order to keep at least 3 samples in each test set for
#' reliable scaling of the test data for comparison, \code{folds} must be <=
#' \code{floor(nrow(X)/3)}
#'
#' The number of selected variables for the following components will then be
#' sequentially optimised. If the number of observations are small (e.g. < 30),
#' it is recommended to use Leave-One-Out Cross-Validation which can be
#' achieved by setting \code{folds = nrow(X)}.
#' @inheritParams spca
#' @param folds Number of folds in 'Mfold' cross-validation. See details.
#' @param nrepeat Number of times the Cross-Validation process is repeated.
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set.
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating the type
#'   of parallelisation. See examples.
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' @importFrom BiocParallel SerialParam bplapply
#' @return A \code{tune.spca} object containing: \describe{ 
#' \item{call}{ The
#'   function call}
#' \item{choice.keepX}{The selected number of components on
#'   each component} 
#' \item{cor.comp}{The correlations between the components
#'   from the cross-validated studies and those from the study which used all of
#'   the data in training.} }
#' @export
#'
#' @example ./examples/tune.spca-examples.R
tune.spca <- function(X, 
                      ncomp = 2, 
                      nrepeat = 1, 
                      folds, 
                      test.keepX, 
                      center = TRUE, 
                      scale = TRUE, 
                      BPPARAM = SerialParam(),
                      seed = NULL)
{
    ## evaluate all args
    BPPARAM$RNGseed <- seed
    mget(names(formals()), sys.frame(sys.nframe()))
    X <- data.matrix(X, rownames.force = TRUE)
    X <- scale(X, center = center, scale = scale)
    ncomp <-      .check_ncomp(ncomp = ncomp, X = X)
    test.keepX <- .check_test.keepX(test.keepX = test.keepX, X = X)
    ## check cv args
    if (!(is.numeric(folds) && (folds > 2 & folds <= floor(nrow(X)/3))))
        stop("'folds' must be an integer > 2 and <= floor(nrow(X)/3)=", floor(nrow(X)/3))
    
    if (!(is.numeric(nrepeat) && (nrepeat > 0)))
        stop("'nrepeat' must be a positive integer")
    
    ## optimal keepX for all components
    keepX.opt <- NULL
    ## a list of cor.df for each component
    cor.df.list <- rep(list(NA), ncomp)
    names(cor.df.list) <- paste0('comp', seq_len(ncomp))
    
    all.keepX <- test.keepX
    names(all.keepX) <- paste0('keepX_', all.keepX)
    
    if (any(is.na(X))) {
      X[which(is.na(X))] <- 0
      warning("There were NAs present in the input dataframe. These were converted to 0 values. If you don't want these as 0, handle missing values prior to tuning.", call. = FALSE)
    }
    
    ## ------ component loop
    for(ncomp in seq_len(ncomp)) {
        iter_keepX <- function(keepX.value) {
            # full data
            spca.full = spca(X, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = FALSE)
            ## ------ repeated cv
            repeat_cv_j <- function(j) {
                repeat.j.folds <- suppressWarnings(split(sample(seq_len(nrow(X))),seq_len(folds)))
                ## ------ mean cor for CV
                cor.pred = sapply(repeat.j.folds, function(test.fold.inds){
                    t.comp.actual <- spca.full$variates$X[,ncomp][test.fold.inds]
                    ## split data to train/test
                    X.train = X[-test.fold.inds,,drop=FALSE]
                    X.test = X[test.fold.inds,,drop=FALSE]
                    # ---- run sPCA 
                    ## train
                    spca.train = spca(X.train, ncomp = ncomp, keepX = c(keepX.opt, keepX.value), center = center, scale = FALSE)
                    # ---- deflation on X with only the fold left out
                    for(k in seq_len(ncomp)){ 
                        # loop to calculate deflated matrix and predicted comp
                        # calculate the predicted comp on the fold left out
                        # calculate reg coeff, then deflate
                        if(k != 1){
                          # calculate deflation beyond comp 1
                          # recalculate the loading vector (here c.sub) on the test set 
                          # (perhaps we could do this instead on the training set by extracting from spca.train$loadings$X[,k]?)
                          c.sub = crossprod(X.test, t.comp.pred) / drop(crossprod(t.comp.pred)) 
                          X.test = X.test - t.comp.pred %*% t(c.sub) 
                          # update predicted comp based on deflated matrix
                        } 
                      t.comp.pred = X.test %*% spca.train$loadings$X[,k]
                    }
                    # calculate predicted component and compare with component from sPCA on full data on the left out set
                    # cor with the component on the full data, abs value
                    cor.comp = abs(cor(t.comp.pred, t.comp.actual, use = 'pairwise.complete.obs'))
                    # need to flip the sign in one of the comp to calculate the RSS
                    #RSS.comp = sum((c(sign(cor.comp))*(t.comp.pred) - spca.full$variates$X[,comp][i])^2)
                    # also look at correlation between loading vectors
                    # cor.loading = abs(cor(spca.full$loadings$X[,comp], spca.train$loadings$X[,comp]))
                    # get the cor(pred component, comp) per fold for a given comp
                    return(cor.comp)
                })  # end CV, get the cor(pred component, comp) per fold for a given comp
                
                # output is correlation between reconstructed components cor.comp on the left out samples and if we had those samples in PCA,
                # average correlations across folds for the repeat
                return(mean(cor.pred, na.rm = TRUE))
            }
            out <-  lapply(seq_len(nrepeat), FUN=repeat_cv_j)
            return(unlist(out))
        }
        test.keepX.cors <- bplapply(X=all.keepX, FUN=iter_keepX, BPPARAM = BPPARAM)
        test.keepX.cors <- data.frame(test.keepX.cors)
        test.keepX.cors <- t(test.keepX.cors)
        colnames(test.keepX.cors) <- paste0('repeat_', seq_len(nrepeat))
        cor.df.list[[ncomp]] <- test.keepX.cors
        ## use a one-sided t.test using repeat correlations to assess if addition of keepX improved the correlation
        ## and get the index of optimum keepX
        t.test.df <- data.frame(t(cor.df.list[[ncomp]]))
        if (nrepeat > 2) 
            keepX.opt.comp.ind <-  t.test.process(t.test.df, alpha = 0.1, alternative = 'less') 
        else 
            keepX.opt.comp.ind <-  which.max(colMeans(t.test.df))
        
        keepX.opt.comp <- all.keepX[keepX.opt.comp.ind]
        ## update keepX.optimum for next comp
        keepX.opt <- c(keepX.opt, keepX.opt.comp)
    }
    if (nrepeat > 2)
    {
        choice.keepX <- keepX.opt
        names(choice.keepX) <- paste0('comp', seq_len(ncomp))
        
    } else
    {
        choice.keepX <- "not computable with 'nrepeat < 3'"
    }
    ## get the mean and sd values
    cor.comp = mapply(df=cor.df.list, opt = keepX.opt, FUN = function(df, opt){
        out <- data.frame(keepX = all.keepX, cor.mean = apply(df, 1, mean), cor.sd = apply(df, 1, sd))
        out$opt.keepX <- NA
        out[which(out$keepX == opt),]$opt.keepX <- ifelse(nrepeat > 2 ,TRUE, NA)
        return(out)
    }, SIMPLIFY = FALSE)
    # evaluate all for output except X to save memory
    mc <- mget(names(formals())[-1], sys.frame(sys.nframe()))
    mc <- as.call(c(as.list(match.call())[1:2], mc))
    
    result <- list(
        call = mc,
        choice.keepX = choice.keepX,
        cor.comp = cor.comp)
    class(result) <- 'tune.spca'
    return(result)
}

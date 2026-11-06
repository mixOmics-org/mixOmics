## tune.pls1
#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2013
# last modified: 05-10-2017
#
# Copyright (C) 2013
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
#############################################################################################################


# ========================================================================================================
# tune.spls: chose the optimal number of parameters per component on a spls method, based on MSE
# ========================================================================================================

# I start with no selection on Y. Otherwise we need to be able to tell which submodel of Y is better. I'm afraid a sum(MSE_Yi) is gonna lead to very sparse model (only 1Y), to test.

#' Tuning functions for sPLS method
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
#' grid to determine optimal values for the sparsity parameters in \code{spls}.
#' 
#' 
#' This tuning function should be used to tune the parameters in the
#' \code{spls} function (number of components and the number of variables in
#' \code{keepX} to select).
#' 
#' If \code{validation = "loo"}, leave-one-out cross-validation is performed.
#' By default \code{folds} is set to the number of unique individuals. If
#' \code{validation = "Mfold"}, M-fold cross-validation is performed. How many
#' folds to generate is selected by specifying the number of folds in
#' \code{folds}.
#' 
#' Four measures of accuracy are available: Mean Absolute Error (\code{MAE}),
#' Mean Square Error(\code{MSE}), \code{Bias} and \code{R2}. Both MAE and MSE
#' average the model prediction error. MAE measures the average magnitude of
#' the errors without considering their direction. It is the average over the
#' fold test samples of the absolute differences between the Y predictions and
#' the actual Y observations. The MSE also measures the average magnitude of
#' the error. Since the errors are squared before they are averaged, the MSE
#' tends to give a relatively high weight to large errors. The Bias is the
#' average of the differences between the Y predictions and the actual Y
#' observations and the R2 is the correlation between the predictions and the
#' observations. All those measures are averaged across all Y variables in the
#' PLS2 case. We are still improving the function to tune an sPLS2 model,
#' contact us for more details and examples.
#' 
#' The function outputs the optimal number of components that achieve the best
#' performance based on the chosen measure of accuracy. The assessment is
#' data-driven and similar to the process detailed in (Rohart et al., 2016),
#' where one-sided t-tests assess whether there is a gain in performance when
#' adding a component to the model.
#' 
#' See also \code{?perf} for more details.
#' 
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y \code{if(method = 'spls')} numeric vector or matrix of continuous
#' responses (for multi-response models) \code{NA}s are allowed.
#' @param ncomp the number of components to include in the model.
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set
#' @param already.tested.X Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{X} data set on
#' the firsts components.
#' @param validation character.  What kind of (internal) validation to use,
#' matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is
#' \code{"Mfold"}.
#' @param folds the folds in the Mfold cross-validation. See Details.
#' @param measure One of \code{MSE} (Mean Squared Error), 
#' \code{MAE} (Mean Absolute Error: MSE without the square), 
#' \code{Bias} (average of the differences), 
#' \code{MAPE} (average of the absolute errors,
#'  as a percentage of the actual values) or \code{R2}. 
#'  Default to \code{MSE}. See details.
#' @param scale boleean. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var boolean, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default value is FALSE
#' @param nrepeat Number of times the Cross-Validation process is repeated.
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements) that indicates the repeated measures on each individual, i.e.
#' the individuals ID. See Details.
#' @param light.output if set to FALSE, the prediction/classification of each
#' sample for each of \code{test.keepX} and each comp is returned.
#' @template arg/BPPARAM
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' @return A list that contains: \item{error.rate}{returns the prediction error
#' for each \code{test.keepX} on each component, averaged across all repeats
#' and subsampling folds. Standard deviation is also output. All error rates
#' are also available as a list.} \item{choice.keepX}{returns the number of
#' variables selected (optimal keepX) on each component.}
#' \item{choice.ncomp}{returns the optimal number of components for the model
#' fitted with \code{$choice.keepX} and \code{$choice.keepY} }
#' \item{measure}{reminds which criterion was used} \item{predict}{Prediction
#' values for each sample, each \code{test.keepX,test.keepY}, each comp and
#' each repeat. Only if light.output=FALSE}
#' @author Kim-Anh Lê Cao, Benoit Gautier, Francois Bartolo, Florian Rohart,
#' Al J Abadi
#' @seealso \code{\link{splsda}}, \code{\link{predict.splsda}} and
#' http://www.mixOmics.org for more details.
#' @references mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' 
#' PLS and PLS citeria for PLS regression: Tenenhaus, M. (1998). La regression
#' PLS: theorie et pratique. Paris: Editions Technic.
#' 
#' Chavent, Marie and Patouille, Brigitte (2003). Calcul des coefficients de
#' regression et du PRESS en regression PLS1. Modulad n, 30 1-11. (this is the
#' formula we use to calculate the Q2 in perf.pls and perf.spls)
#' 
#' Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction
#' (MSEP) Estimates for Principal Component Regression (PCR) and Partial Least
#' Squares Regression (PLSR). Journal of Chemometrics 18(9), 422-429.
#' 
#' sparse PLS regression mode:
#' 
#' Lê Cao, K. A., Rossouw D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. Statistical
#' Applications in Genetics and Molecular Biology 7, article 35.
#' 
#' One-sided t-tests (suppl material):
#' 
#' Rohart F, Mason EA, Matigian N, Mosbergen R, Korn O, Chen T, Butcher S,
#' Patel J, Atkinson K, Khosrotehrani K, Fisk NM, Lê Cao K-A&, Wells CA&
#' (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ
#' 4:e1845.
#' @keywords regression multivariate
#' @examples
#' 
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' 
#' \dontrun{
#' tune = tune.spls(X, Y, ncomp=4, test.keepX = c(5,10,15), measure = "MSE",
#' nrepeat=3, progressBar = TRUE)
#' 
#' tune$choice.ncomp
#' tune$choice.keepX
#' 
#' # plot the results
#' plot(tune)
#' }
#' @noRd
tune.spls1 <- #TODO naming for tune.spls1 is wrong. it should hint the block where variable selection is performed
    function (X, 
              Y,
              ncomp = 1,
              test.keepX = c(5, 10, 15),
              already.tested.X,
              validation = "Mfold",
              folds = 10,
              measure = "MSE", # can do a R2 per Y (correlation, linear regression R2), need to call MSEP (see perf.spls).
              scale = TRUE,
              progressBar = FALSE,
              tol = 1e-06,
              max.iter = 100,
              near.zero.var = FALSE,
              nrepeat = 1,
              multilevel = NULL,
              light.output = TRUE,
              BPPARAM = SerialParam(),
              seed = NULL
    )
    {    #-- checking general input parameters --------------------------------------#
        #---------------------------------------------------------------------------#
        
        #------------------#
      BPPARAM$RNGseed <- seed
      set.seed(seed)
      
        if(!is(X, "matrix"))
            X= as.matrix(X)
        
        if(!is(Y, "matrix"))
            Y= as.matrix(Y)
        
        #-- check entries --#
        if (length(dim(X)) != 2 || !is.numeric(X))
            stop("'X' must be a numeric matrix.")
        
        if (length(dim(Y)) != 2 || !is.numeric(Y))
            stop("'Y' must be a numeric matrix.")
        
        #-- progressBar
        if (!is.logical(progressBar))
            stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
        
        #-- measure
        choices = c("MSE", "MAE","Bias","R2")
        measure = choices[pmatch(measure, choices)]
        if (is.na(measure))
            stop("'measure' must be one of 'MSE', 'MAE', 'Bias' or 'R2' ")
        
        if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
            stop("invalid number of variates, 'ncomp'.")
        
        
        #-- validation
        choices = c("Mfold", "loo")
        validation = choices[pmatch(validation, choices)]
        if (is.na(validation))
            stop("'validation' must be either 'Mfold' or 'loo'")
        
        if (validation == 'loo')
        {
            if (nrepeat != 1)
                warning("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
            nrepeat = 1
        }
        
        
        #-- already.tested.X
        if (missing(already.tested.X))
        {
            already.tested.X = NULL
        } else {
            if(is.null(already.tested.X) | length(already.tested.X)==0)
                stop("''already.tested.X' must be a vector of keepX values")
            
            if(is.list(already.tested.X))
                stop("''already.tested.X' must be a vector of keepX values")
            
            message(paste("Number of variables selected on the first", length(already.tested.X), "component(s):", paste(already.tested.X,collapse = " ")))
        }
        
        if(length(already.tested.X) >= ncomp)
            stop("'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",length(already.tested.X) , call. = FALSE)
        
        if (any(is.na(validation)) || length(validation) > 1)
            stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
        
        #-- test.keepX
        if (is.null(test.keepX) | length(test.keepX) == 1 | !is.numeric(test.keepX))
            stop("'test.keepX' must be a numeric vector with more than two entries", call. = FALSE)
        
        # add colnames and rownames if missing
        X.names = dimnames(X)[[2]]
        if (is.null(X.names))
        {
            X.names = paste("X", seq_len(ncol(X)), sep = "")
            dimnames(X)[[2]] = X.names
        }
        
        ind.names = dimnames(X)[[1]]
        if (is.null(ind.names))
        {
            ind.names = seq_len(nrow(X))
            rownames(X)  = ind.names
        }
        
        if (length(unique(rownames(X))) != nrow(X))
            stop("samples should have a unique identifier/rowname")
        if (length(unique(X.names)) != ncol(X))
            stop("Unique indentifier is needed for the columns of X")
        
        
        #-- end checking --#
        #------------------#
        
        
        
        #---------------------------------------------------------------------------#
        #-- NA calculation      ----------------------------------------------------#
        
        misdata = c(X=anyNA(X), Y=anyNA(Y)) # Detection of missing data. we assume no missing values in the factor Y
        
        if (any(misdata))
        {
            is.na.A.X = is.na(X)
            is.na.A.Y = is.na(Y)
            is.na.A = list(X=is.na.A.X, Y=is.na.A.Y)
        } else {
            is.na.A = NULL
        }
        #-- NA calculation      ----------------------------------------------------#
        #---------------------------------------------------------------------------#
        
        
        
        test.keepX = sort(unique(test.keepX)) #sort test.keepX so as to be sure to chose the smallest in case of several minimum
        names(test.keepX) = test.keepX
        
        # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
        if ((!is.null(already.tested.X)))
        {
            comp.real = (length(already.tested.X) + 1):ncomp
        } else {
            comp.real = seq_len(ncomp)
        }
        
        
        mat.error.rate = list()
        error.per.class = list()
        
        mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
                              dimnames = list(test.keepX, c(paste0('comp', comp.real))))
        mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
                                dimnames = list(test.keepX, c(paste0('comp', comp.real))))
        
        # first: near zero var on the whole data set
        if(near.zero.var == TRUE)
        {
            nzv = nearZeroVar(X)
            if (length(nzv$Position > 0))
            {
                warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
                X = X[, -nzv$Position, drop=TRUE]
                
                if(ncol(X)==0)
                    stop("No more predictors after Near Zero Var has been applied!")
                
            }
        }
        
        if (is.null(multilevel) | (!is.null(multilevel) && ncol(multilevel) == 2))
        {
            if(light.output == FALSE)
                prediction.all = list()
            
            
            class.object="mixo_spls"
            
            # successively tune the components until ncomp: comp1, then comp2, ...
            for(comp in seq_len(length(comp.real)))
            {
                
                if (progressBar == TRUE)
                    cat("\ncomp",comp.real[comp], "\n")
                
                result = MCVfold.spls(X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = 1 + length(already.tested.X),
                                      choice.keepX = already.tested.X,
                                      test.keepX = test.keepX, test.keepY = ncol(Y), measure = measure, dist = NULL, auc=FALSE, scale=scale,
                                      near.zero.var = near.zero.var, progressBar = progressBar, tol = tol, max.iter = max.iter,
                                      BPPARAM,
                                      misdata = misdata, is.na.A = is.na.A, class.object=class.object)
                
                #save(list=ls(),file="temp.Rdata")
                # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.spls' can work with multiple distances
                mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
                mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
                if (!is.null(result[[measure]]$error.rate.sd[[1]]))
                    mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
                
                
                # best keepX
                already.tested.X = c(already.tested.X, result[[measure]]$keepX.opt[[1]])
                
                if(light.output == FALSE)
                {
                    #prediction of each samples for each fold and each repeat, on each comp
                    prediction.all[[comp]] = result$prediction.comp
                }
                
                
                
            } # end comp
            
            names(mat.error.rate) = c(paste0('comp', comp.real))
            names(already.tested.X) = c(paste0('comp', seq_len(ncomp)))
            
            if (progressBar == TRUE)
                cat('\n')
            
            # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
            if(nrepeat > 2 & length(comp.real) >1)
            {
                keepX = already.tested.X
                error.keepX = NULL
                for(comp in seq_len(length(comp.real)))
                {
                    ind.row = match(keepX[[comp.real[comp]]],test.keepX)
                    error.keepX = cbind(error.keepX, apply(matrix(mat.error.rate[[comp]][[ind.row]],ncol=nrepeat),2,mean))
                }
                colnames(error.keepX) = c(paste0('comp', comp.real))
                rownames(error.keepX) = c(paste0('nrep.', seq_len(nrepeat)))
                
                opt = t.test.process(error.keepX)
                
                ncomp_opt = comp.real[opt]
            } else {
                ncomp_opt = error.keepX = NULL
            }
            
            
            result = list(
                error.rate = mat.mean.error,
                error.rate.sd = mat.sd.error,
                error.rate.all = mat.error.rate,
                choice.keepX = already.tested.X,
                choice.ncomp = list(ncomp = ncomp_opt, values = error.keepX))
            
            if(light.output == FALSE)
            {
                names(prediction.all) = c(paste0('comp', comp.real))
                result$predict = prediction.all
            }
            
            result$measure = measure
            result$call = match.call()
            
            class(result) = "tune.spls1"
            
            return(result)
        } else {
            # if multilevel with 2 factors, we can not do as before because withinvariation depends on the factors, we maximase a correlation
            message("Not implemented at the moment")
            
        }
    }





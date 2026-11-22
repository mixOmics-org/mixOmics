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

# ========================================================================================================
# tune.splsda: tuning hyperparameters on a splsda method
# ========================================================================================================

#' Tuning functions for sPLS-DA method
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
#' grid to determine optimal values for the parameters in \code{splsda}.
#'
#' 
#' This tuning function should be used to tune the parameters in the
#' \code{splsda} function (number of components and number of variables in
#' \code{keepX} to select).
#' 
#' For a sPLS-DA, M-fold or LOO cross-validation is performed with stratified
#' subsampling where all classes are represented in each fold.
#' 
#' If \code{validation = "loo"}, leave-one-out cross-validation is performed.
#' By default \code{folds} is set to the number of unique individuals.
#' 
#' The function outputs the optimal number of components that achieve the best
#' performance based on the overall error rate or BER. The assessment is
#' data-driven and similar to the process detailed in (Rohart et al., 2016),
#' where one-sided t-tests assess whether there is a gain in performance when
#' adding a component to the model. Our experience has shown that in most case,
#' the optimal number of components is the number of categories in \code{Y} -
#' 1, but it is worth tuning a few extra components to check (see our website
#' and case studies for more details).
#' 
#' For sPLS-DA multilevel one-factor analysis, M-fold or LOO cross-validation
#' is performed where all repeated measurements of one sample are in the same
#' fold. Note that logratio transform and the multilevel analysis are performed
#' internally and independently on the training and test set.
#' 
#' For a sPLS-DA multilevel two-factor analysis, the correlation between
#' components from the within-subject variation of X and the \code{cond} matrix
#' is computed on the whole data set. The reason why we cannot obtain a
#' cross-validation error rate as for the spls-DA one-factor analysis is
#' because of the difficulty to decompose and predict the within matrices
#' within each fold.
#' 
#' For a sPLS two-factor analysis a sPLS canonical mode is run, and the
#' correlation between components from the within-subject variation of X and Y
#' is computed on the whole data set.
#' 
#' If \code{validation = "Mfold"}, M-fold cross-validation is performed. How
#' many folds to generate is selected by specifying the number of folds in
#' \code{folds}.
#' 
#' If \code{auc = TRUE} and there are more than 2 categories in \code{Y}, the
#' Area Under the Curve is averaged using one-vs-all comparison. Note however
#' that the AUC criteria may not be particularly insightful as the prediction
#' threshold we use in sPLS-DA differs from an AUC threshold (sPLS-DA relies on
#' prediction distances for predictions, see \code{?predict.splsda} for more
#' details) and the supplemental material of the mixOmics article (Rohart et
#' al. 2017). If you want the AUC criterion to be insightful, you should use
#' \code{measure==AUC} as this will output the number of variable that
#' maximises the AUC; in this case there is no prediction threshold from
#' sPLS-DA (\code{dist} is not used). If \code{measure==AUC}, we do not output
#' SD as this measure can be a mean (over \code{nrepeat}) of means (over the
#' categories).
#' 
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017).
#' 
#' If test.keepX is set to NULL, the \code{perf()} function will be run internally,
#' which performs cross-validation to identify optimal number of components and 
#' distance measure. Running tuning initially using \code{test.keepX = NULL} speeds
#' up the parameter tuning workflow, as then a lower ncomp value can be used for 
#' variable selection tuning. 
#' 
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y \code{if(method = 'spls')} numeric vector or matrix of continuous
#' responses (for multi-response models) \code{NA}s are allowed.
#' @param ncomp the number of components to include in the model.
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set. If set to NULL, tuning will be performed on ncomp
#' using all variables in the \eqn{X} data set. 
#' @param already.tested.X Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{X} data set on
#' the firsts components.
#' @param scale Logical. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var Logical, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default value is FALSE
#' @param logratio one of ('none','CLR'). Default to 'none'
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements) that indicates the repeated measures on each individual, i.e.
#' the individuals ID. See Details.
#' @param validation character. What kind of (internal) validation to use,
#'   matching one of \code{"Mfold"} or \code{"loo"} (short for 'leave-one-out').
#'   Default is \code{"Mfold"}.
#' @param folds the folds in the Mfold cross-validation. See Details.
#' @param nrepeat Number of times the Cross-Validation process is repeated.
#' @param dist distance metric to use for \code{splsda} to estimate the
#' classification error rate, should one of \code{"centroids.dist"},
#' \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details). 
#' If \code{test.keepX = NULL} multiple distances can be inputted or "all". 
#' @param measure Three misclassification measure are available: overall
#' misclassification error \code{overall}, the Balanced Error Rate \code{BER}
#' or the Area Under the Curve \code{AUC}. Only used when \code{test.keepX} is not NULL. 
#' @param signif.threshold numeric between 0 and 1 indicating the significance
#' threshold required for improvement in error rate of the components. Default
#' to 0.01.
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' performance of the model based on the optimisation measure \code{measure}.
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param light.output if set to FALSE, the prediction/classification of each
#' sample for each of \code{test.keepX} and each comp is returned.
#' @template arg/BPPARAM
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' @return Depending on the type of analysis performed, a list that contains:
#' \item{error.rate}{returns the prediction error for each \code{test.keepX} on
#' each component, averaged across all repeats and subsampling folds. Standard
#' deviation is also output. All error rates are also available as a list.}
#' \item{choice.keepX}{returns the number of variables selected (optimal keepX)
#' on each component.} \item{choice.ncomp}{returns the optimal number of
#' components for the model fitted with \code{$choice.keepX} }
#' \item{error.rate.class}{returns the error rate for each level of \code{Y}
#' and for each component computed with the optimal keepX}
#' 
#' If test.keepX = FALSE,produces a matrix of classification
#' error rate estimation. The dimensions correspond to the components in the
#' model and to the prediction method used, respectively. Note that error rates
#' reported in any component include the performance of the model in earlier
#' components for the specified \code{keepX} parameters (e.g. error rate
#' reported for component 3 for \code{keepX = 20} already includes the fitted
#' model on components 1 and 2 for \code{keepX = 20}).
#' 
#' \item{predict}{Prediction values for each sample, each \code{test.keepX},
#' each comp and each repeat. Only if light.output=FALSE}
#' \item{class}{Predicted class for each sample, each \code{test.keepX}, each
#' comp and each repeat. Only if light.output=FALSE}
#' 
#' \item{auc}{AUC mean and standard deviation if the number of categories in
#' \code{Y} is greater than 2, see details above. Only if auc = TRUE}
#' 
#' \item{cor.value}{only if multilevel analysis with 2 factors: correlation
#' between latent variables.}
#' @author Kim-Anh Lê Cao, Benoit Gautier, Francois Bartolo, Florian Rohart,
#' Al J Abadi
#' @seealso \code{\link{splsda}}, \code{\link{predict.splsda}} and
#' https://mixomics.org/ for more details.
#' @references mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' @export
#' @example ./examples/tune.splsda-examples.R

tune.splsda <- 
    function (X, Y,
              ncomp = 1,
              # sparsity params
              test.keepX = NULL, # in which case perf() is run internally to only tune ncomp
              already.tested.X,
              # params related to spls model building
              scale = TRUE,
              logratio = c('none','CLR'),
              max.iter = 100,
              tol = 1e-06,
              near.zero.var = FALSE,
              multilevel = NULL,
              # params related to CV
              validation = "Mfold",
              folds = 10,
              nrepeat = 1,
              signif.threshold = 0.01, 
              # params related to PA
              dist = "max.dist",
              measure = "BER", # or 'overall'
              auc = FALSE,
              # params related to running
              progressBar = FALSE,
              light.output = TRUE,
              BPPARAM = SerialParam(),
              seed = NULL
    )
    {   #-- checking general input parameters --------------------------------------#
        #---------------------------------------------------------------------------#
        
        BPPARAM$RNGseed <- seed
        set.seed(seed)
      
        #-- check significance threshold
        signif.threshold <- .check_alpha(signif.threshold)
        
        #------------------#
        #-- check entries --#
        if(!is(X, "matrix"))
            X = as.matrix(X)
        
        if (length(dim(X)) != 2 || !is.numeric(X))
            stop("'X' must be a numeric matrix.")
        
        
        # Testing the input Y
        if (is.null(multilevel))
        {
            if (is.null(Y))
                stop("'Y' has to be something else than NULL.")
            
            if (is.null(dim(Y)))
            {
                Y = factor(Y)
            }  else {
                stop("'Y' should be a factor or a class vector.")
            }
            
            if (nlevels(Y) == 1)
                stop("'Y' should be a factor with more than one level")
            
        } else {
            # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
            multilevel = data.frame(multilevel)
            
            if ((nrow(X) != nrow(multilevel)))
                stop("unequal number of rows in 'X' and 'multilevel'.")
            
            if (ncol(multilevel) != 1)
                stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
            
            if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0,1,2))# multilevel 1 or 2 factors
                stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
            
        }
        
        
        #-- progressBar
        if (!is.logical(progressBar))
            stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
        
        
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
        
        if (nrepeat < 3 & ncomp > 1)
            message("Note that the number of components cannot be reliably tuned with nrepeat < 3 or validaion = 'loo'.")
        
        #-- logratio
        logratio <- match.arg(logratio)

        #--already.tested.X    
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
        
        #-- measure
        choices = c("BER", "overall","AUC")
        measure = choices[pmatch(measure, choices)]
        if (is.na(measure))
            stop("'measure' must be either 'BER', 'overall' or 'AUC' ")
        
        if (any(is.na(validation)) || length(validation) > 1)
            stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)

        #-- test.keepX 
        #-> if test.keepX set to NULL run perf() i.e. parameter tuning of just ncomp
        if (is.null(test.keepX)){
            print("test.keepX set to NULL, tuning only for number of components...")
            splsda_res <- splsda(X, Y, ncomp, scale = scale, tol = tol, max.iter = max.iter, 
            near.zero.var = near.zero.var, logratio = logratio, multilevel = multilevel)
            perf_res <- perf(splsda_res, 
                validation = validation, folds = folds, nrepeat = nrepeat,
                dist = dist, signif.threshold = signif.threshold,
                BPPARAM = BPPARAM, seed = seed, progressBar = progressBar,
                auc = auc)
            return(perf_res)
        #-> if test.keepX is not NULL, run tune function as before
        } else {

            # reformat multilevel design for tuning variables
            if (!is.null(multilevel)){
                multilevel = data.frame(multilevel, Y)
                multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
            }
            
            #-- test.keepX
            if (length(test.keepX) == 1 | !is.numeric(test.keepX))
                stop("'test.keepX' must be a numeric vector with more than two entries", call. = FALSE)
            
            # remove some test.keepX if needed
            if (any(test.keepX > ncol(X))){
                test.keepX = test.keepX[-which(test.keepX>ncol(X))]
                if (length(test.keepX) < 2)
                    stop("Some entries of 'test.keepX' were higher than the number of
            variables in 'X' and were removed, the resulting 'test.keepX' has now
            too few entries (<2)", call. = FALSE)
            }
            
            # add colnames and rownames if missing
            X.names = dimnames(X)[[2]]
            if (is.null(X.names))
            {
                X.names = paste0("X", 1:ncol(X))
                dimnames(X)[[2]] = X.names
            }
            
            ind.names = dimnames(X)[[1]]
            if (is.null(ind.names))
            {
                ind.names = 1:nrow(X)
                rownames(X)  = ind.names
            }
            
            if (length(unique(rownames(X))) != nrow(X))
                stop("samples should have a unique identifier/rowname")
            if (length(unique(X.names)) != ncol(X))
                stop("Unique indentifier is needed for the columns of X")
            
            rm(X.names);rm(ind.names)
            
            #-- end checking --#
            #------------------#
            
            
            #---------------------------------------------------------------------------#
            #-- logration + multilevel approach ----------------------------------------#
            # we can do logratio and multilevel on the whole data as these transformation are done per sample
            X = logratio.transfo(X = X, logratio = logratio)
            
            if (!is.null(multilevel)) # logratio is applied per sample, multilevel as well, so both can be done on the whole data
            {
                
                Xw = withinVariation(X, design = multilevel)
                X = Xw
                
                #-- Need to set Y variable for 1 or 2 factors
                Y = multilevel[, -1, drop=FALSE]
                if (ncol(Y) >= 1)
                    Y = apply(Y, 1, paste, collapse = ".")  #  paste is to combine in the case we have 2 levels
                
                Y = as.factor(Y)
            }
            #-- logration + multilevel approach ----------------------------------------#
            #---------------------------------------------------------------------------#
            
            
            #---------------------------------------------------------------------------#
            #-- NA calculation      ----------------------------------------------------#
            
            misdata = c(X=anyNA(X), Y=FALSE) # Detection of missing data. we assume no missing values in the factor Y
            
            if (any(misdata))
            {
                is.na.A = is.na(X)
                
                #ind.NA = which(apply(is.na.A, 1, sum) > 0) # calculated only once
                #ind.NA.col = which(apply(is.na.A, 2, sum) >0) # indice of the col that have missing values. used in the deflation
            } else {
                is.na.A = NULL
                #ind.NA = ind.NA.col = NULL
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
                comp.real = 1:ncomp
            }
            
            
            choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist")
            dist = match.arg(dist, choices, several.ok = TRUE)
            
            
            mat.error.rate = list()
            error.per.class = list()
            AUC = list()
            
            if(nrepeat>1 & all(measure != "AUC")){
                mat.sd.error = matrix(0,nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
                                    dimnames = list(c(test.keepX), c(paste0('comp', comp.real))))
            } else{
                mat.sd.error=NULL
            }
            mat.mean.error = matrix(nrow = length(test.keepX), ncol = ncomp-length(already.tested.X),
                                    dimnames = list(c(test.keepX), c(paste0('comp', comp.real))))
            
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
                    prediction.all = class.all = list()
                
                if(auc)
                {
                    auc.mean.sd=list()
                    if(light.output == FALSE)
                        auc.all=list()
                }
                
                class.object=c("mixo_splsda","DA")
                error.per.class.keepX.opt = list()

                if(all(measure != "AUC")){
                    error.per.class.keepX.opt.mean = matrix(0, nrow = nlevels(Y), ncol = length(comp.real),
                                                            dimnames = list(c(levels(Y)), c(paste0('comp', comp.real))))
                } else {
                    error.per.class.keepX.opt.mean=NULL
                }
                # successively tune the components until ncomp: comp1, then comp2, ...
                for(comp in 1:length(comp.real))
                {
                    
                    if (progressBar == TRUE & is(BPPARAM, 'SerialParam'))
                        cat("\ncomp",comp.real[comp], "\n")
                    
                    result = MCVfold.spls(X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = 1 + length(already.tested.X),
                                        choice.keepX = already.tested.X,
                                        test.keepX = test.keepX, test.keepY = nlevels(Y), measure = measure, dist = dist, scale=scale,
                                        near.zero.var = near.zero.var, progressBar = progressBar, tol = tol, max.iter = max.iter, auc = auc,
                                        BPPARAM = BPPARAM,
                                        misdata = misdata, is.na.A = is.na.A, class.object=class.object)
                    
                    # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.spls' can work with multiple distances
                    mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
                    mat.mean.error[, comp]=result[[measure]]$error.rate.mean[[1]]
                    if (!is.null(result[[measure]]$error.rate.sd[[1]]))
                        mat.sd.error[, comp]=result[[measure]]$error.rate.sd[[1]]
                    
                    # confusion matrix for keepX.opt
                    if (all(measure!="AUC")){
                        error.per.class.keepX.opt[[comp]]=result[[measure]]$confusion[[1]]
                        error.per.class.keepX.opt.mean[, comp]=apply(result[[measure]]$confusion[[1]], 1, mean)
                    }
                    
                    # best keepX
                    already.tested.X = c(already.tested.X, result[[measure]]$keepX.opt[[1]])
                    
                    if(light.output == FALSE)
                    {
                        #prediction of each samples for each fold and each repeat, on each comp
                        class.all[[comp]] = result$class.comp[[1]]
                        prediction.all[[comp]] = result$prediction.comp
                    }
                    
                    if(auc)
                    {
                        auc.mean.sd[[comp]] = result$auc
                        if(light.output == FALSE)
                            auc.all[[comp]] = result$auc.all
                    }
                    
                } # end comp
                
                names(mat.error.rate) = c(paste0('comp', comp.real))
                if (all(measure!="AUC"))
                    names(error.per.class.keepX.opt) = c(paste0('comp', comp.real))
                
                names(already.tested.X) = c(paste0('comp', 1:ncomp))
                
                if (progressBar == TRUE)
                    cat('\n')
                
                # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
                if(nrepeat > 2 & length(comp.real) >1 & all(measure!="AUC"))
                {
                    keepX = already.tested.X
                    error.keepX = NULL
                    for(comp in 1:length(comp.real))
                    {
                        ind.row = match(keepX[[comp.real[comp]]],test.keepX)
                        error.keepX = cbind(error.keepX, mat.error.rate[[comp]][ind.row,])
                    }
                    colnames(error.keepX) = c(paste0('comp', comp.real))
                    
                    opt = t.test.process(error.keepX)
                    
                    ncomp_opt = comp.real[opt]
                }  else if(nrepeat > 2 & length(comp.real) >1 & all(measure=="AUC")){
                    #hacking t.test.process for AUC
                    keepX = already.tested.X
                    error.keepX = NULL
                    for(comp in 1:length(comp.real))
                    {
                        ind.row = match(keepX[[comp.real[comp]]],test.keepX)
                        error.keepX = cbind(error.keepX, mat.error.rate[[comp]][ind.row,])
                    }
                    colnames(error.keepX) = c(paste0('comp', comp.real))
                    
                    #1- to change the AUC to an 'error' and check for decrease
                    opt = t.test.process(1-error.keepX, alpha = signif.threshold)
                    
                    ncomp_opt = comp.real[opt]
                } else {
                    ncomp_opt = error.keepX = NULL
                }
                
                result = list(
                    error.rate = mat.mean.error,
                    error.rate.sd = mat.sd.error,
                    error.rate.all = mat.error.rate,
                    choice.keepX = already.tested.X,
                    choice.ncomp = list(ncomp = ncomp_opt, values = error.keepX),
                    error.rate.class = error.per.class.keepX.opt.mean,
                    error.rate.class.all = error.per.class.keepX.opt)
                
                if(light.output == FALSE)
                {
                    names(class.all) = names(prediction.all) = c(paste0('comp', comp.real))
                    result$predict = prediction.all
                    result$class = class.all
                }
                if(auc)
                {
                    
                    #we add the AUC outputs only if it was not the measure
                    #otherwise there is twice the same output
                    if(all(measure != "AUC")){
                        names(auc.mean.sd) = c(paste0('comp', comp.real))
                        result$auc = auc.mean.sd
                    }
                    if(light.output == FALSE)
                    {
                        names(auc.all) = c(paste0('comp', comp.real))
                        result$auc.all =auc.all
                    }
                }
                result$measure = measure
                result$call = match.call()
                result$call$nrepeat = nrepeat
                result$nrepeat = nrepeat
                
                class(result) = "tune.splsda"
                
                return(result)
            } else {
                # if multilevel with 2 factors, we can not do as before because withinvariation depends on the factors, we maximase a correlation
                message("For a two-factor analysis, the tuning criterion is based on the maximisation of the correlation between the components on the whole data set")
                
                cor.value = vector(length = length(test.keepX))
                names(cor.value) = test.keepX
                
                for (i in 1:length(test.keepX))
                {
                    spls.train = splsda(X, Y, ncomp = ncomp, keepX = c(already.tested.X, test.keepX[i]), logratio = logratio, near.zero.var = FALSE)
                    
                    # Note: this is performed on the full data set
                    # (could be done with resampling (bootstrap) (option 1) and/or prediction (option 2))
                    cor.value[i] = cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
                    #
                }
                return(list(cor.value = cor.value))
                
            }
        } # end else statement (i.e. test.keepX is not NULL)
    }

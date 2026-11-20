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
# tune.plsda: tuning hyperparameters on a plsda method
# ========================================================================================================

#' Tuning functions for PLS-DA method
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
#' grid to determine optimal values for the parameters in \code{plsda}.
#'
#' 
#' This tuning function should be used to tune the parameters in the
#' \code{plsda} function (number of components and distance metric to select).
#' 
#' For a PLS-DA, M-fold or LOO cross-validation is performed with stratified
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
#' For PLS-DA multilevel one-factor analysis, M-fold or LOO cross-validation
#' is performed where all repeated measurements of one sample are in the same
#' fold. Note that logratio transform and the multilevel analysis are performed
#' internally and independently on the training and test set.
#' 
#' For a PLS-DA multilevel two-factor analysis, the correlation between
#' components from the within-subject variation of X and the \code{cond} matrix
#' is computed on the whole data set. The reason why we cannot obtain a
#' cross-validation error rate as for the pls-DA one-factor analysis is
#' because of the difficulty to decompose and predict the within matrices
#' within each fold.
#' 
#' For a PLS two-factor analysis a PLS canonical mode is run, and the
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
#' threshold we use in PLS-DA differs from an AUC threshold (PLS-DA relies on
#' prediction distances for predictions, see \code{?predict.plsda} for more
#' details) and the supplemental material of the mixOmics article (Rohart et
#' al. 2017).
#' 
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017).
#' 
#' The tune.plsda() function calls older function perf() to perform this cross-validation,
#' for more details see the perf() help pages. 
#' 
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y \code{if(method = 'spls')} numeric vector or matrix of continuous
#' responses (for multi-response models) \code{NA}s are allowed.
#' @param ncomp the number of components to include in the model.
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
#' classification error rate, should be a subset of \code{"centroids.dist"},
#' \code{"mahalanobis.dist"} or \code{"max.dist"} (see Details).
#' @param signif.threshold numeric between 0 and 1 indicating the significance
#' threshold required for improvement in error rate of the components. Default
#' to 0.01.
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param light.output if set to FALSE, the prediction/classification of each
#' sample for each of \code{test.keepX} and each comp is returned.
#' @template arg/BPPARAM
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' @return matrix of classification error rate estimation. 
#' The dimensions correspond to the components in the
#' model and to the prediction method used, respectively.
#' 
#' \item{auc}{Averaged AUC values
#' over the \code{nrepeat}}
#' 
#' \item{cor.value}{only if multilevel analysis with 2 factors: correlation
#' between latent variables.}
#' 
#' @author Kim-Anh Lê Cao, Benoit Gautier, Francois Bartolo, Florian Rohart,
#' Al J Abadi
#' @seealso \code{\link{splsda}}, \code{\link{predict.splsda}} and
#' http://www.mixOmics.org for more details.
#' @references mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' @export
#' @example ./examples/tune.plsda-examples.R

tune.plsda <- 
    function (X, Y,
              ncomp = 1,
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
              dist = "all",
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

        #-- validation
        if (any(is.na(validation)) || length(validation) > 1)
            stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)


        #------------------#
        #-- run perf --#
        plsda_res <- plsda(X, Y, ncomp, scale = scale, tol = tol, max.iter = max.iter, 
            near.zero.var = near.zero.var, logratio = logratio, multilevel = multilevel)
        perf_res <- perf(plsda_res, 
                validation = validation, folds = folds, nrepeat = nrepeat,
                dist = dist, signif.threshold = signif.threshold,
                BPPARAM = BPPARAM, seed = seed, progressBar = progressBar,
                auc = auc)
        return(perf_res)
    }
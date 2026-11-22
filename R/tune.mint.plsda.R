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
# tune.mint.plsda: chose the optimal number of parameters per component on a mint.plsda method
# ========================================================================================================
#' Estimate the parameters of mint.plsda method
#' 
#' Computes Leave-One-Group-Out-Cross-Validation (LOGOCV) scores on a
#' user-input grid to determine optimal values for the parameters in
#' \code{mint.plsda}.
#' 
#' This function performs a Leave-One-Group-Out-Cross-Validation (LOGOCV),
#' where each of \code{study} is left out once.
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
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017).
#' 
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y Outcome. Numeric vector or matrix of responses (for multi-response
#' models)
#' @param ncomp Number of components to include in the model (see Details).
#' Default to 1
#' @param study grouping factor indicating which samples are from the same
#' study
#' @param dist only applies to an object inheriting from \code{"plsda"} or
#' \code{"splsda"} to evaluate the classification performance of the model.
#' Should be a subset of \code{"max.dist"}, \code{"centroids.dist"},
#' \code{"mahalanobis.dist"}. Default is \code{"all"}. See
#' \code{\link{predict}}.
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' performance of the model.
#' @param progressBar by default set to \code{TRUE} to output the progress bar
#' of the computation.
#' @param scale Logical. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var Logical, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Default value is FALSE
#' @param light.output if set to FALSE, the prediction/classification of each
#' sample for each of \code{test.keepX} and each comp is returned.
#' @param signif.threshold numeric between 0 and 1 indicating the significance
#' threshold required for improvement in error rate of the components. Default
#' to 0.01.
#' @return The returned value is a list with components:
#' \item{study.specific.error}{A list that gives BER, overall error rate and
#' error rate per class, for each study} \item{global.error}{A list that gives
#' BER, overall error rate and error rate per class for all samples}
#' \item{predict}{A list of length \code{ncomp} that produces the predicted
#' values of each sample for each class} \item{class}{A list which gives the
#' predicted class of each sample for each \code{dist} and each of the
#' \code{ncomp} components. Directly obtained from the \code{predict} output.}
#' \item{auc}{AUC values} \item{auc.study}{AUC values for each study in mint
#' models}.
#' 
#' @author Florian Rohart, Al J Abadi
#' @seealso \code{\link{mint.plsda}} and https://mixomics.org/ for more
#' details.
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords multivariate dplot
#' @export
#' @example ./examples/tune.mint.plsda-examples.R

tune.mint.plsda <-
    function (X, 
              Y,
              ncomp = 1,
              study,
              # model building params
              scale = TRUE,
              tol = 1e-06,
              max.iter = 100,
              near.zero.var = FALSE,
              # CV params
              signif.threshold = 0.01,
              # PA params
              dist = c("max.dist", "centroids.dist", "mahalanobis.dist"),
              auc = FALSE,
              # running params
              progressBar = FALSE,
              light.output = TRUE # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
    )
    {    
        #-- checking general input parameters --------------------------------------#
        #---------------------------------------------------------------------------#

        ## R CMD check stuff
        BPPARAM <- seed <- NULL
        
        #------------------#
        #-- check entries --#
        if(missing(X))
            stop("'X'is missing", call. = FALSE)
        
        X = as.matrix(X)
        
        if (length(dim(X)) != 2 || !is.numeric(X))
            stop("'X' must be a numeric matrix.", call. = FALSE)
        
        
        # Testing the input Y
        if(missing(Y))
            stop("'Y'is missing", call. = FALSE)
        if (is.null(Y))
            stop("'Y' has to be something else than NULL.", call. = FALSE)
        
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.", call. = FALSE)
        }
        
        if (nlevels(Y) == 1)
            stop("'Y' should be a factor with more than one level", call. = FALSE)
        
        #-- check significance threshold
        signif.threshold <- .check_alpha(signif.threshold)
        
        #-- progressBar
        if (!is.logical(progressBar))
            stop("'progressBar' must be a logical constant (TRUE or FALSE).", call. = FALSE)
        
        
        if (is.null(ncomp) || !is.numeric(ncomp) || ncomp <= 0)
            stop("invalid number of variates, 'ncomp'.")
        
        # -- check using the check of mint.splsda
        Y.mat = unmap(Y)
        colnames(Y.mat) = levels(Y)
        
        check = Check.entry.pls(X, Y = Y.mat, ncomp = ncomp, mode="regression", scale=scale,
                                near.zero.var=near.zero.var, max.iter=max.iter ,tol=tol ,logratio="none" ,DA=TRUE, multilevel=NULL)
        X = check$X
        ncomp = check$ncomp
        
        # -- study
        #set the default study factor
        if (missing(study))
            stop("'study' is missing", call. = FALSE)
        
        if (length(study) != nrow(X))
            stop(paste0("'study' must be a factor of length ",nrow(X),"."))
        
        if (any(table(study) <= 1))
            stop("At least one study has only one sample, please consider removing before calling the function again", call. = FALSE)
        if (any(table(study) < 5))
            warning("At least one study has less than 5 samples, mean centering might not do as expected")
        
        if(sum(apply(table(Y,study)!=0,2,sum)==1) >0)
            stop("At least one study only contains a single level of the multi-levels outcome Y. The MINT algorithm cannot be computed.")
        
        if(sum(apply(table(Y,study)==0,2,sum)>0) >0)
            warning("At least one study does not contain all the levels of the outcome Y. The MINT algorithm might not perform as expected.")
        
        #-- light.output
        if (!is.logical(light.output))
            stop("'light.output' must be either TRUE or FALSE", call. = FALSE)
        
        #-- run perf to tune ncomp --------------------------------------#
        #---------------------------------------------------------------------------#
        mint.plsda_res <- mint.plsda(X, Y, ncomp = ncomp, study = study,
                    scale = scale, tol = tol, max.iter = max.iter, near.zero.var = near.zero.var)
        perf_res <- perf(mint.plsda_res, 
                    dist = dist,
                    BPPARAM = BPPARAM, seed = seed, progressBar = progressBar)
        return(perf_res)
    }
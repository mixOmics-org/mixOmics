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
# tune.block.plsda: Tuning hyperparameters on a block plsda method
# ========================================================================================================
#' Tuning function for block.plsda method (N-integration with Discriminant Analysis)
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores based on a
#' user-input grid to determine the optimal parameters for
#' method \code{block.plsda}.
#' 
#' This tuning function should be used to tune the number of components in the \code{block.plsda} function (N-integration with Discriminant Analysis).
#' 
#' M-fold or LOO cross-validation is performed with stratified subsampling
#' where all classes are represented in each fold.
#' 
#' If \code{validation = "Mfold"}, M-fold cross-validation is performed. The
#' number of folds to generate is to be specified in the argument \code{folds}.
#' 
#' If \code{validation = "loo"}, leave-one-out cross-validation is performed.
#' By default \code{folds} is set to the number of unique individuals.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017). Details
#' about the PLS modes are in \code{?pls}.
#' 
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#' 
#' @inheritParams block.plsda
#' @inheritParams tune
#' @param dist Distance metric. Should be a subset of "max.dist", "centroids.dist", "mahalanobis.dist" or "all".
#' Default is "all"
#' @param weighted tune using either the performance of the Majority vote or
#' the Weighted vote.
#' @param signif.threshold numeric between 0 and 1 indicating the significance
#' threshold required for improvement in error rate of the components. Default
#' to 0.01.
#' @param ... Optional arguments:
#' \itemize{
#'  \item \bold{seed} Integer. Seed number for reproducible parallel code.
#'  Default is \code{NULL}.
#' }
#' run in parallel when repeating the cross-validation, which is usually the
#' most computationally intensive process. If there is excess CPU, the
#' cross-vaidation is also parallelised on *nix-based OS which support
#' \code{mclapply}.
#' Note that the argument 'scheme' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' @return returns:
#' \item{error.rate}{Prediction error rate for each block of \code{object$X}
#' and each \code{dist}} \item{error.rate.per.class}{Prediction error rate for
#' each block of \code{object$X}, each \code{dist} and each class}
#' \item{predict}{Predicted values of each sample for each class, each block
#' and each component} \item{class}{Predicted class of each sample for each
#' block, each \code{dist}, each component and each nrepeat} \item{features}{a
#' list of features selected across the folds (\code{$stable.X} and
#' \code{$stable.Y}) for the \code{keepX} and \code{keepY} parameters from the
#' input object.} \item{AveragedPredict.class}{if more than one block, returns
#' the average predicted class over the blocks (averaged of the \code{Predict}
#' output and prediction using the \code{max.dist} distance)}
#' \item{AveragedPredict.error.rate}{if more than one block, returns the
#' average predicted error rate over the blocks (using the
#' \code{AveragedPredict.class} output)} \item{WeightedPredict.class}{if more
#' than one block, returns the weighted predicted class over the blocks
#' (weighted average of the \code{Predict} output and prediction using the
#' \code{max.dist} distance). See details for more info on weights.}
#' \item{WeightedPredict.error.rate}{if more than one block, returns the
#' weighted average predicted error rate over the blocks (using the
#' \code{WeightedPredict.class} output.)} \item{MajorityVote}{if more than one
#' block, returns the majority class over the blocks. NA for a sample means that
#' there is no consensus on the predicted class for this particular sample over
#' the blocks.} \item{MajorityVote.error.rate}{if more than one block, returns
#' the error rate of the \code{MajorityVote} output}
#' \item{WeightedVote}{if more than one block, returns the weighted majority
#' class over the blocks. NA for a sample means that there is no consensus on
#' the predicted class for this particular sample over the blocks.}
#' \item{WeightedVote.error.rate}{if more than one block, returns the error
#' rate of the \code{WeightedVote} output} \item{weights}{Returns the weights
#' of each block used for the weighted predictions, for each nrepeat and each
#' fold} \item{choice.ncomp}{For supervised models; returns the optimal number
#' of components for the model for each prediction distance using one-sided
#' t-tests that test for a significant difference in the mean error rate (gain
#' in prediction) when components are added to the model. See more details in
#' Rohart et al 2017 Suppl. For more than one block, an optimal ncomp is
#' returned for each prediction framework.}
#' 
#' @author Florian Rohart, Amrit Singh, Kim-Anh Lê Cao, AL J Abadi
#' @seealso \code{\link{block.splsda}} and http://www.mixOmics.org for more
#' details.
#' @references Method:
#' 
#' Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê
#' Cao K.A. (2016). DIABLO: multi omics integration for biomarker discovery.
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' @export
#' @example ./examples/tune.block.plsda-examples.R
tune.block.plsda <- 
  function (
            # basic params
            X,
            Y,
            indY,
            ncomp = 2,
            # model building params
            tol = 1e-06,
            max.iter = 100,
            near.zero.var = FALSE,
            design,
            scale = TRUE,
            # cross validation params
            validation = "Mfold",
            folds = 10,
            nrepeat = 1,
            signif.threshold=0.01,
            # measure of performance params
            dist = "all",
            weighted = TRUE, # optimise the weighted or not-weighted prediction
            # processing params
            progressBar = FALSE,
            light.output = TRUE, # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
            BPPARAM = SerialParam(),
            seed = NULL,
            ...)
  {
    if (hasArg('cpus')) #defunct
    {
    stop("'cpus' has been replaced by BPPARAM. See documentation.")  
    }
    BPPARAM$RNGseed <- seed
    set.seed(seed)

    # hardcode init and scheme
    scheme <- 'horst'
    init <- 'svd.single'

    ## ----------- checks -----------
    
    # check input 'Y' and transformation in a dummy matrix
    if (!missing(Y))
    {
      if (is.null(dim(Y)))
      {
        Y = factor(Y)
      } else {
        stop("'Y' should be a factor or a class vector.")
      }
      
      if (nlevels(Y) == 1)
        stop("'Y' should be a factor with more than one level")
      
    } else if (!missing(indY)) {
      Y = X[[indY]]
      if (is.null(dim(Y)))
      {
        Y = factor(Y)
      } else {
        stop("'Y' should be a factor or a class vector.")
      }
      
      if (nlevels(Y) == 1)
        stop("'X[[indY]]' should be a factor with more than one level")
      
      X = X[-indY] #remove Y from X to pass the arguments simpler to block.splsda
      
    } else if (missing(indY)) {
      stop("Either 'Y' or 'indY' is needed")
      
    }
    ## check using internal #TODO we need to unify the checks
    Y.check <- unmap(Y)
    Y.check <- matrix(Y.check, nrow = nrow(Y.check), dimnames = list(rownames(X[[1]]), NULL))
    Check.entry.wrapper.mint.block(X = X, Y = Y.check, indY = indY,
                                   ncomp = ncomp, DA=TRUE,
                                   design = design, init = init, scheme = scheme, scale = scale,
                                   near.zero.var = near.zero.var, mode = 'regression', tol = tol,
                                   max.iter = max.iter)
    
    ## ensure all X blocks are matrices, keeping dimnames
    X <- lapply(X, function(z){
      zm <- z
      if (!is.matrix(zm)) {
        zm <- as.matrix(zm)
        dimnames(zm) <- dimnames(z)
      }
      return(zm)
    })
    
    #-- progressBar
    if (!is.logical(progressBar))
      stop("'progressBar' must be a logical constant (TRUE or FALSE).",
           call. = FALSE)
    
    #-- ncomp
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
        message("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
      nrepeat = 1
    }
    
    #-- check significance threshold
    signif.threshold <- .check_alpha(signif.threshold)
    
    #-- validation
    if (any(is.na(validation)) || length(validation) > 1)
      stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    ## ----------- Cross-validation -----------

    block.plsda_res <- block.plsda(X, Y, ncomp = ncomp, 
                  scale = scale, tol = tol, max.iter = max.iter, near.zero.var = near.zero.var, design = design)
    perf_res <- perf(block.plsda_res, 
                validation = validation, folds = folds, nrepeat = nrepeat,
                dist = dist,
                BPPARAM = BPPARAM, seed = seed, progressBar = progressBar)
    return(perf_res)

}
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
# tune.block.splsda: Tuning hyperparameters on a block splsda method
# ========================================================================================================
#' Tuning function for block.splsda method (N-integration with sparse
#' Discriminant Analysis)
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores based on a
#' user-input grid to determine the optimal parameters for
#' method \code{block.splsda}.
#' 
#' This tuning function should be used to tune the number of components and the 
#' keepX parameters in the \code{block.splsda} function (N-integration with sparse Discriminant
#' Analysis).
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
#' All combination of test.keepX values are tested. A message informs how many
#' will be fitted on each component for a given test.keepX.
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
#' @inheritParams block.splsda
#' @inheritParams tune
#' @inheritParams tune.spca
#' @param test.keepX A named list with the same length and names as X 
#' (without the outcome Y, if it is provided in X and designated using 
#' \code{indY}).  Each entry of this list is a numeric vector for the different 
#' keepX values to test for that specific block. If set to NULL, ncomp is tuned. 
#' @param already.tested.X Optional, if \code{ncomp > 1} A named list of 
#' numeric vectors each of length \code{n_tested} indicating the number of 
#' variables to select from the \eqn{X} data set on the first \code{n_tested} 
#' components.
#' @param measure only used when \code{test.keepX} is not NULL. Measure used when plotting, 
#' should be 'BER' or 'overall'
#' @param dist distance metric to estimate the classification error rate, should be one of 
#' "centroids.dist", "mahalanobis.dist" or "max.dist" (see Details). If \code{test.keepX = NULL}, 
#' can also input "all" or more than one distance metric
#' @param weighted tune using either the performance of the Majority vote or
#' the Weighted vote.
#' @param signif.threshold numeric between 0 and 1 indicating the significance
#' threshold required for improvement in error rate of the components. Default
#' to 0.01.
#' @return A list that contains: \item{error.rate}{returns the prediction error
#' for each \code{test.keepX} on each component, averaged across all repeats
#' and subsampling folds. Standard deviation is also output. All error rates
#' are also available as a list.} \item{choice.keepX}{returns the number of
#' variables selected (optimal keepX) on each component, for each block.}
#' \item{choice.ncomp}{returns the optimal number of components for the model
#' fitted with \code{$choice.keepX}. } \item{error.rate.class}{returns the
#' error rate for each level of \code{Y} and for each component computed with
#' the optimal keepX}
#' 
#' \item{predict}{Prediction values for each sample, each \code{test.keepX},
#' each comp and each repeat. Only if light.output=FALSE}
#' \item{class}{Predicted class for each sample, each \code{test.keepX}, each
#' comp and each repeat. Only if light.output=FALSE}
#' 
#' \item{cor.value}{compute the correlation between latent variables for
#' two-factor sPLS-DA analysis.}
#' 
#' If \code{test.keepX = NULL}, returns:
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
#' @seealso \code{\link{block.splsda}} and https://mixomics.org/ for more
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
#' @example ./examples/tune.block.splsda-examples.R
tune.block.splsda <- 
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
            # sparse method params
            test.keepX,
            already.tested.X,
            # cross validation params
            validation = "Mfold",
            folds = 10,
            nrepeat = 1,
            signif.threshold=0.01,
            # measure of performance params
            dist = "max.dist",
            measure = "BER", # one of c("overall","BER")
            weighted = TRUE, # optimise the weighted or not-weighted prediction
            # processing params
            progressBar = FALSE,
            light.output = TRUE, # if FALSE, output the prediction and classification of each sample during each folds, on each comp, for each repeat
            BPPARAM = SerialParam(),
            seed = NULL)
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
    
    #-- measure
    measure.input = measure
    if (!measure %in% c("overall", "BER"))
      stop("'measure' must be 'overall' or 'BER'")
    
    #-- check significance threshold
    signif.threshold <- .check_alpha(signif.threshold)
    
    #-- already.tested.X
    if (missing(already.tested.X))
    {
      already.tested.X = NULL
    } else {
      if (is.null(already.tested.X))
        stop("'already.tested.X' must be a vector of keepX values ")
      
      # we require the same number of already tuned components on each block
      if (length(unique(sapply(already.tested.X, length))) > 1)
        stop(
          "The same number of components must be already tuned for each block, in 'already.tested.X'"
        )
      
      if (any(sapply(already.tested.X, function(x)
        is.list(x))) == TRUE)
        stop(" Each entry of 'already.tested.X' must be a vector of keepX values")
      
      if (length(already.tested.X[[1]]) >= ncomp)
        stop(
          "'ncomp' needs to be higher than the number of components already tuned, which is length(already.tested.X)=",
          length(already.tested.X) ,
          call. = FALSE
        )
    }
    
    #-- validation
    if (any(is.na(validation)) || length(validation) > 1)
      stop("'validation' should be one of 'Mfold' or 'loo'.", call. = FALSE)
    
    #-- test.keepX
    if (is.null(test.keepX)){
      print("test.keepX is set to NULL, tuning only for number of components...")
      block.splsda_res <- block.splsda(X, Y, ncomp = ncomp,
                  scale = scale, tol = tol, max.iter = max.iter, near.zero.var = near.zero.var, design = design)
      perf_res <- perf(block.splsda_res, 
                validation = validation, folds = folds, nrepeat = nrepeat,
                dist = dist,
                BPPARAM = BPPARAM, seed = seed, progressBar = progressBar)
      return(perf_res)

    } else {

      #-- dist (for tune can only have one dist, for perf can have multiple)
      dist = match.arg(
      dist,
      choices = c("max.dist", "centroids.dist", "mahalanobis.dist"),
      several.ok = TRUE
    )

      if (missing(test.keepX))
      {
        test.keepX = lapply(X, function(x) {
          max.test.keepX <- min(30, ncol(x))
          if (max.test.keepX > 15)
            return(seq(5, max.test.keepX, 5))
          else
            return(seq(1, max.test.keepX, 2))
        })
        
      } else {
        if (length(test.keepX) != length(X))
          stop(
            paste(
              "test.keepX should be a list of length ",
              length(X),
              ", corresponding to the blocks: ",
              paste(names(X), collapse = ", "),
              sep = ""
            )
          )
        
        #aa = sapply(test.keepX, length)
        #if (any(is.null(aa) | aa == 1 | !is.numeric(aa)))
        #stop("Each entry of 'test.keepX' must be a numeric vector with more than two values", call. = FALSE)
        
      }
      
      l = sapply(test.keepX, length)
      n = names(test.keepX)
      temp = data.frame(l, n)
      
      
      message(
        paste(
          "\nYou have provided a sequence of keepX of length: ",
          paste(apply(temp, 1, function(x)
            paste(x, collapse = " for block ")), collapse = " and "),
          ".\nThis results in ",
          prod(sapply(test.keepX, length)),
          " models being fitted for each component and each nrepeat, this may take some time to run, be patient!",
          sep = ""
        )
      )
      
      if (is (BPPARAM, 'SerialParam'))
      {
        message(paste0(
          "\nYou can look into the 'BPPARAM' argument to speed up computation time."
        ))
        
      } else {
        if (progressBar == TRUE)
          message(paste0(
            "\nAs code is running in parallel, the progressBar is not available."
          ))
      }
      
      ## ----------- END checks -----------#
      
      ## ----------- NA calculation ----------- 
      
      misdata = c(sapply(X, anyNA), Y = FALSE) # Detection of missing data. we assume no missing values in the factor Y
      
      is.na.A = vector("list", length = length(X))
      for (q in seq_along(X))
      {
        if (misdata[q])
        {
          is.na.A[[q]] = is.na(X[[q]])
          #ind.NA[[q]] = which(apply(is.na.A[[q]], 1, sum) > 0) # calculated only once
          #ind.NA.col[[q]] = which(apply(is.na.A[[q]], 2, sum) >0) # indice of the col that have missing values. used in the deflation
        }
      }
      
      ## ----------- END NA calculation ----------- #
      
      
      # if some components have already been tuned (eg comp1 and comp2), we're only tuning the following ones (comp3 comp4 .. ncomp)
      if ((!is.null(already.tested.X)) & length(already.tested.X) > 0)
      {
        comp.real = (length(already.tested.X[[1]]) + 1):ncomp
        #check and match already.tested.X to X
        if (length(already.tested.X[[1]]) > 0)
        {
          if (length(unique(names(already.tested.X))) != length(already.tested.X) |
              sum(is.na(match(names(
                already.tested.X
              ), names(X)))) > 0)
            stop(
              "Each entry of 'already.tested.X' must have a unique name corresponding to a block of 'X'"
            )
          
        }
        
      } else {
        comp.real = seq_len(ncomp)
      }
      
      # near zero var on the whole data sets. It will be performed inside each fold as well
      if (near.zero.var == TRUE)
      {
        nzv.A = lapply(X, nearZeroVar)
        for (q in seq_along(X))
        {
          if (length(nzv.A[[q]]$Position) > 0)
          {
            names.remove.X = colnames(X[[q]])[nzv.A[[q]]$Position]
            X[[q]] = X[[q]][, -nzv.A[[q]]$Position, drop = FALSE]
            warning(
              "Zero- or near-zero variance predictors.\n Reset predictors matrix to not near-zero variance predictors.\n See $nzv for problematic predictors."
            )
            if (ncol(X[[q]]) == 0)
              stop(paste0("No more variables in", X[[q]]))
            
            #need to check that the keepA[[q]] is now not higher than ncol(A[[q]])
            if (any(test.keepX[[q]] > ncol(X[[q]])))
              test.keepX[[q]][which(test.keepX[[q]] > ncol(X[[q]]))] = ncol(X[[q]])
          }
          
        }
      }
      N.test.keepX = nrow(expand.grid(test.keepX))
      
      mat.error.rate = list()
      
      mat.sd.error = matrix(0,
                            nrow = N.test.keepX,
                            ncol = ncomp - length(already.tested.X[[1]]))
      
      mat.mean.error = matrix(nrow = N.test.keepX,
                              ncol = ncomp - length(already.tested.X[[1]]))
      
      
      mat.error.rate = list()
      error.per.class.keepX.opt = list()
      error.per.class.keepX.opt.mean = matrix(
        0,
        nrow = nlevels(Y),
        ncol = length(comp.real),
        dimnames = list(c(levels(Y)), c(paste0('comp', comp.real)))
      )
      
      error.opt.per.comp = matrix(
        nrow = nrepeat,
        ncol = length(comp.real),
        dimnames = list(paste("nrep", seq_len(nrepeat), sep = "."), paste0("comp", comp.real))
      )
      
      if (light.output == FALSE)
        class.all = list()
      
      ## ----------- tune components ----------- 
      
      # successively tune the components until ncomp: comp1, then comp2, ...
      for (comp in seq_along(comp.real))
      {
        tune_comp <- comp.real[comp]
        if (progressBar == TRUE)
          cat(sprintf("\ntuning component %s\n", tune_comp))
        
        result = MCVfold.block.splsda(
          X,
          Y,
          validation = validation,
          folds = folds,
          nrepeat = nrepeat,
          ncomp = tune_comp,
          choice.keepX = already.tested.X,
          scheme = scheme,
          design = design,
          init = init,
          tol = tol,
          test.keepX = test.keepX,
          measure = measure,
          dist = dist,
          scale = scale,
          weighted = weighted,
          near.zero.var = near.zero.var,
          progressBar = progressBar,
          max.iter = max.iter,
          misdata = misdata,
          is.na.A = is.na.A,
          BPPARAM = BPPARAM
        )
        
        
        ## returns error.rate for all test.keepX
        
        # in the following, there is [[1]] because 'tune' is working with only 1 distance and 'MCVfold.block.splsda' can work with multiple distances
        mat.error.rate[[comp]] = result[[measure]]$mat.error.rate[[1]]
        mat.mean.error[, comp] = result[[measure]]$error.rate.mean[[1]]
        if (!is.null(result[[measure]]$error.rate.sd[[1]]))
          mat.sd.error[, comp] = result[[measure]]$error.rate.sd[[1]]
        
        # confusion matrix for keepX.opt
        error.per.class.keepX.opt[[comp]] = result[[measure]]$confusion[[1]]
        error.per.class.keepX.opt.mean[, comp] = apply(result[[measure]]$confusion[[1]], 1, mean)
        
        # error rate for best keepX
        error.opt.per.comp[, comp] = mat.error.rate[[comp]][result[[measure]]$ind.keepX.opt[[1]], ]
        
        # best keepX
        already.tested.X = result[[measure]]$choice.keepX
        
        if (light.output == FALSE)
        {
          #prediction of each samples for each fold and each repeat, on each comp
          class.all[[comp]] = result$class.comp[[1]]
        }
      }
      
      ## ----------- END tune components ----------- #
      
      ## ----------- output ----------- 
      
      rownames(mat.mean.error) = rownames(result[[measure]]$mat.error.rate[[1]])
      colnames(mat.mean.error) = paste0("comp", comp.real)
      names(mat.error.rate) = c(paste0("comp", comp.real))
      names(error.per.class.keepX.opt) = c(paste0("comp", comp.real))
      if (nrepeat > 1)
      {
        rownames(mat.sd.error) = rownames(result[[measure]]$mat.error.rate[[1]])
        colnames(mat.sd.error) = paste0("comp", comp.real)
      }
      
      
      # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
      if (nrepeat > 2 & length(comp.real) > 1)
      {
        error.keepX = error.opt.per.comp
        opt = t.test.process(error.opt.per.comp, alpha = signif.threshold)
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
        error.rate.class = error.per.class.keepX.opt
      )
      
      result$measure = measure.input
      result$call = match.call()
      
      class(result) = "tune.block.splsda"
      
      return(result)
      
    }
  }
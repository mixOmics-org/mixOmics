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
# tune.spls: Tuning hyperparameters on a spls method
# ========================================================================================================
#' 
#' Tuning functions for sPLS method
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
#' grid to determine optimal values for the parameters in \code{spls}.
#'
#' 
#' This tuning function should be used to tune the parameters in the
#' \code{spls} function (number of components and number of variables to select).
#' 
#' 
#' @template section/folds
#' @template section/nrepeat
#' @template section/measure-pls
#' @template section/t-test-process
#' 
#' @section more:
#' See also \code{?perf} for more details.
#' 
#' @param X numeric matrix of predictors with the rows as individual observations.
#' @param Y numeric matrix of response(s) with the rows as individual observations matching \code{X}.
#' @template arg/ncomp
#' @template arg/test.keepX-X.matrix
#' @template arg/test.keepY
#' @template arg/validation
#' @template arg/folds
#' @template arg/nrepeat
#' @param measure The tuning measure to use. Cannot be NULL when applied to sPLS1 object. See details.
#' @templateVar modes \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}
#' @template arg/mode
#' @param scale Logical. If scale = TRUE, each block is standardized to zero means and unit variances (default: TRUE
#' @param tol Positive numeric used as convergence criteria/tolerance during the iterative process. Default to 1e-06.
#' @param max.iter Integer, the maximum number of iterations. Default to 100.
#' @param near.zero.var Logical, see the internal nearZeroVar function (should be set to TRUE in particular for data with many zero values). Setting this argument to FALSE (when appropriate) will speed up the computations. Default value is FALSE.
#' @param logratio Character, one of ('none','CLR') specifies the log ratio transformation to deal with compositional values that may arise from specific normalisation in sequencing data. Default to 'none'. See ?logratio.transfo for details.
#' @param multilevel Numeric, design matrix for repeated measurement analysis, where multilevel decomposition is required. For a one factor decomposition, the repeated measures on each individual, i.e. the individuals ID is input as the first column. For a 2 level factor decomposition then 2nd AND 3rd columns indicate those factors. See examples.
#' @template arg/progressBar
#' @template arg/BPPARAM
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' @param ... Optional parameters passed to \code{\link{spls}}
#' @return 
#' If \code{test.keepX != NULL} and \code{test.keepY != NULL} returns a list that contains: 
#' \item{cor.pred}{The correlation of predicted vs
#'   actual components from X (t) and Y (u) for each
#'   component}\item{RSS.pred}{The Residual Sum of Squares of predicted vs
#'   actual components from X (t) and Y (u) for each component}
#'   \item{choice.keepX}{returns the number of variables selected for X (optimal
#'   keepX) on each component.} \item{choice.keepY}{returns the number of
#'   variables selected for Y (optimal keepY) on each component.}
#'   \item{choice.ncomp}{returns the optimal number of components for the model
#'   fitted with \code{$choice.keepX} and \code{$choice.keepY} } \item{call}{The
#'   functioncal call including the parameteres used.}
#' 
#' If \code{test.keepX = NULL} and \code{test.keepY = NULL} returns a list with the following components for every repeat:
#' \item{MSEP}{Mean Square Error Prediction for each \eqn{Y} variable, only 
#' applies to object inherited from \code{"pls"}, and \code{"spls"}. Only 
#' available when in regression (s)PLS.} 
#' \item{RMSEP}{Root Mean Square Error Prediction for each \eqn{Y} variable, only 
#' applies to object inherited from \code{"pls"}, and \code{"spls"}. Only 
#' available when in regression (s)PLS.} 
#' \item{R2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables for models 
#' with \eqn{1, \ldots ,}\code{ncomp} components, only applies to object
#' inherited from \code{"pls"}, and \code{"spls"}. Only available when in 
#' regression (s)PLS.}
#' \item{Q2}{if \eqn{Y} contains one variable, a vector of \eqn{Q^2} values
#' else a list with a matrix of \eqn{Q^2} values for each \eqn{Y}-variable.
#' Note that in the specific case of an sPLS model, it is better to have a look
#' at the Q2.total criterion, only applies to object inherited from
#' \code{"pls"}, and \code{"spls"}. Only available when in regression (s)PLS.} 
#' \item{Q2.total}{a vector of \eqn{Q^2}-total values for models with \eqn{1, 
#' \ldots ,}\code{ncomp} components, only applies to object inherited from 
#' \code{"pls"}, and \code{"spls"}. Available in both (s)PLS modes.}
#' \item{RSS}{Residual Sum of Squares across all selected features and the 
#' components.}
#' \item{PRESS}{Predicted Residual Error Sum of Squares across all selected 
#' features and the components.}
#' \item{features}{a list of features selected across the 
#' folds (\code{$stable.X} and \code{$stable.Y}) for the \code{keepX} and
#' \code{keepY} parameters from the input object. Note, this will be \code{NULL} 
#' if using standard (non-sparse) PLS.} 
#' \item{cor.tpred, cor.upred}{Correlation between the 
#' predicted and actual components for X (t) and Y (u)} 
#' \item{RSS.tpred, RSS.upred}{Residual Sum of Squares between the
#' predicted and actual components for X (t) and Y (u)} 
#' 
#' @author Kim-Anh Lê Cao, Al J Abadi, Benoit Gautier, Francois Bartolo,
#' Florian Rohart,
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
#' @export
#' @example ./examples/tune.spls-examples.R
#' 
# TODO change this so that it simply wraps perf
tune.spls <- 
  function(X,
           Y,
           test.keepX = NULL,
           test.keepY = NULL,
           ncomp,
           # params related to spls model building
           mode = c('regression', 'canonical', 'classic'),
           scale = TRUE,
           logratio = "none",
           tol = 1e-09,
           max.iter = 100,
           near.zero.var = FALSE,
           multilevel = NULL,
           # params related to CV
           validation = c('Mfold', 'loo'),
           nrepeat = 1,
           folds,
           measure = NULL,
           # params related to running tune
           BPPARAM = SerialParam(),
           seed = NULL,
           progressBar = FALSE,
           ...
  ) {

    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    out = list()
    mode <- match.arg(mode)
    
    BPPARAM$RNGseed <- seed

    # hardcode to streamline
    limQ2 <- 0.0975
    
    X <- .check_numeric_matrix(X, block_name = 'X')
    Y <- .check_numeric_matrix(Y, block_name = 'Y')
    
    check_cv <- .check_cv_args(validation = validation, 
                               nrepeat = nrepeat, folds = folds, 
                               N = nrow(X))
    validation <- check_cv$validation
    nrepeat <- check_cv$nrepeat
    folds <- check_cv$folds
    
    #-- test.keepX 
    #-> if test.keepX set to NULL run perf() i.e. parameter tuning of just ncomp using all features
    if (is.null(test.keepX) && is.null(test.keepY)){
      print("test.keepX and test.keepY are set to NULL, tuning only for number of components...")
      spls_res <- spls(X, Y, ncomp, keepX = ncol(X), keepY = ncol(Y),
                  mode = mode, scale = scale, logratio = logratio, tol = tol, max.iter = max.iter, near.zero.var = near.zero.var, multilevel = multilevel)
      perf_res <- perf(spls_res, 
                validation = validation, folds = folds, nrepeat = nrepeat,
                dist = dist,
                BPPARAM = BPPARAM, seed = seed, progressBar = progressBar)
      return(perf_res)

    } else {

    test.keepX <- .change_if_null(arg = test.keepX, default = ncol(X))
    test.keepY <- .change_if_null(arg = test.keepY, default = ncol(Y))
    
    test.keepX <- unique(test.keepX)
    test.keepY <- unique(test.keepY)
    
    spls.model <- !(length(test.keepX) == 1 & length(test.keepY) == 1)
    
    if (ncol(Y) == 1) {
      res <- tune.spls1(X = X, 
                        Y = Y,
                        ncomp = ncomp,
                        test.keepX = test.keepX,
                        validation = validation,
                        folds = folds,
                        measure = measure, # can do a R2 per Y (correlation, linear regression R2), need to call MSEP (see perf.spls).
                        progressBar = progressBar,
                        nrepeat = nrepeat,
                        BPPARAM = BPPARAM,
                        seed = seed,
                        ...
      )
      ## --- call
      res$call <- NULL
      ## eval all but X and Y
      mc <- mget(names(formals())[-1:-2], sys.frame(sys.nframe()))
      ## replace function, X and Y with unevaluated call
      mc <- as.call(c(as.list(match.call())[1:3], mc))
      res <- c(list(call = mc), res)
      class(res) <- 'tune.spls1'
      return(res)
    }
    
    measure <- match.arg(measure, choices = c('cor', 'RSS'))
    
    choice.keepX = choice.keepY = NULL
    measure.table.cols <- c('comp', 'keepX', 'keepY', 'repeat', 't', 'u', 'cor', 'RSS')
    # if (spls.model)
    # {
    measure.pred <- expand.grid(keepX = test.keepX, 
                                keepY = test.keepY,
                                V = c('u', 't'),
                                measure = c('cor', 'RSS'),
                                comp = seq_len(ncomp),
                                optimum.keepA = FALSE)
    # } else {
    #   # Q2 only
    #   measure.pred <- expand.grid(
    #                               comp = seq_len(ncomp),
    #                               optimum = FALSE)
    # }
    measure.pred <- data.frame(measure.pred)
    measure.pred <- cbind(measure.pred, 
                          value.v = I(rep(list(data.frame(matrix(NA_real_, ncol= 1, nrow = nrepeat))), 
                                          times = nrow(measure.pred))),
                          value.Q2.total = I(rep(list(data.frame(matrix(NA_real_, ncol= 1, nrow = nrepeat))), 
                                                 times = nrow(measure.pred)))
    )
    # measure.pred$value <- lapply(measure.pred$value, as.data.frame)
    
    use_progressBar <- progressBar & (is(BPPARAM, 'SerialParam'))
    n_keepA <- length(test.keepX) * length(test.keepY)
    
    ## initialise optimal keepX/Y
    measure.pred$optimum.keepA <- FALSE 
    measure.pred[
      measure.pred$keepX == test.keepX[1] &
        measure.pred$keepY == test.keepY[1]
      ,]$optimum.keepA <- TRUE
    
    
    # list of all keepX/keepY passed to this via bplapply. calculates the 
    # correlation and RSS of components for given input keepX/Y and returns
    # a subset of measure.pred containing these updated values
    iter_keep.vals <- function(keep.vals) { # keep.vals = list(keepX, keepY)
      
      keepX <- keep.vals$keepX # set keep values
      keepY <- keep.vals$keepY
      
      if (use_progressBar) {
        n_tested <- n_tested + 1
        prog_level <- n_tested / n_keepA
        .progressBar(prog_level, title = 'of features tested')
      }
      # generate model using input parameters
      pls.model <- spls(X = X, Y = Y, 
                        keepX = c(choice.keepX, keepX), 
                        keepY = c(choice.keepY, keepY), 
                        ncomp = comp, mode = mode, 
                        scale = scale, tol = tol, max.iter = max.iter,
                        near.zero.var = near.zero.var, logratio = logratio,
                        multilevel = multilevel,
                        ...)
      
      # run perf in serial to avoid nested parallel processes, but make sure seed is passed into perf
      pls.perf <- perf(pls.model, validation = validation, folds = folds, nrepeat = nrepeat,
                       BPPARAM = SerialParam(), seed = seed)
      
      
      # ensures the only rows that are returned are those manipulated in this
      # iteration. otherwise, later component iterations will undo this work
      adjusted.rows <- c()
      for (measure_i in c('cor', 'RSS')) # for each measure...
      {
        
        for (v in c('u', 't')) # for each set of components ...
        {
          ## extract the relevant row index in measure.pred FOR VALUE.V
          vpred.idx <- rownames(measure.pred[measure.pred$comp == comp & 
                                               measure.pred$keepX == keepX &
                                               measure.pred$keepY == keepY &
                                               measure.pred$V == v &
                                               measure.pred$measure == measure_i,])
          
          # extract relevant values
          measure.vpred <- pls.perf$measures[[sprintf("%s.%spred", measure_i, v)]]$values
          measure.vpred <- measure.vpred[measure.vpred$comp == comp,]
          measure.pred[vpred.idx,]$value.v[[1]] <- measure.vpred$value # adjust measure.pred
          
          ## extract the relevant row index in measure.pred FOR Q2
          Q2.idx <- rownames(measure.pred[measure.pred$comp == comp & 
                                            measure.pred$keepX == keepX &
                                            measure.pred$keepY == keepY &
                                            # TODO the following filtering criteria needs review - not sure if Q2 is per c('u' , 't')
                                            measure.pred$V == v &
                                            measure.pred$measure == measure_i
                                          ,])
          # extract relevant values
          value.Q2.total <- pls.perf$measures$Q2.total$values
          value.Q2.total <- filter(value.Q2.total, comp == comp)$value
          measure.pred[Q2.idx,]$value.Q2.total[[1]] <- value.Q2.total  # adjust measure.pred
          
          adjusted.rows <- c(adjusted.rows, vpred.idx, Q2.idx) # append adjusted rows to list
        }
        
      }
      return(measure.pred[unique(adjusted.rows),])
    }
    
    # takes in the final measure.preds data.frame and adjusts the optimum.keepA 
    # feature to reflect the optimal values depending on the t.test
    DetermineOptimalKeepVals <- function(measure.pred, comp, keep.vals, nrepeat) {
      
      ## output TRUE if value improves opt for cor or RSS, else FALSE
      .check_improvement <- function(opt, value, measure, nrepeat) {
        
        if (nrepeat > 2) {
          t.test.res <- tryCatch(t.test(x = opt, y = value, alternative = ifelse(measure == 'cor', 'less', 'greater')),
                                 error = function(e) e)
          improved <- t.test.res$p.value < 0.05
        } else
        {
          ## compare values
          improved <-  if (measure == 'cor') mean(opt) < mean(value) else mean(opt) > mean(value)
        }
        improved
      }
      
      for (keep.val in keep.vals) {
        trial.keepX <- keep.val$keepX # set keep values
        trial.keepY <- keep.val$keepY
        
        # current best for X components
        optimum.u <- measure.pred[measure.pred$comp == comp & 
                                    measure.pred$optimum.keepA == TRUE &
                                    measure.pred$V == 'u' &
                                    measure.pred$measure == measure
                                  ,]$value.v[[1]]
        # current best for Y components
        optimum.t <- measure.pred[measure.pred$comp == comp &
                                    measure.pred$optimum.keepA == TRUE &
                                    measure.pred$V == 't' &
                                    measure.pred$measure == measure
                                  ,]$value.v[[1]]
        # new keepX/Y pair to try for X components
        value.u <- measure.pred[measure.pred$comp == comp &
                                  measure.pred$keepX == trial.keepX &
                                  measure.pred$keepY == trial.keepY &
                                  measure.pred$V == 'u' &
                                  measure.pred$measure == measure
                                ,]$value.v[[1]]
        # new keepX/Y pair to try for Y components
        value.t <- measure.pred[measure.pred$comp == comp &
                                  measure.pred$keepX == trial.keepX &
                                  measure.pred$keepY == trial.keepY &
                                  measure.pred$V == 't' &
                                  measure.pred$measure == measure
                                ,]$value.v[[1]]
        
        ## workaround for constant values in t.test 
        offset.eps <- seq(1, 2, length.out = length(optimum.u))/1e6
        optimum.t <- optimum.t + offset.eps
        optimum.u <- optimum.u + offset.eps
        value.t <- value.t + offset.eps
        value.u <- value.u + offset.eps
        
        # see if there is significant improvement in measure based on new keepX/Y pair
        u.improved <-.check_improvement(opt = optimum.u, value = value.u, measure = measure, nrepeat = nrepeat)
        t.improved <-.check_improvement(opt = optimum.t, value = value.t, measure = measure, nrepeat = nrepeat)
        improved <- if (mode == 'canonical') u.improved & t.improved else u.improved
        if (improved) # if so ...
        {
          ## set previous to FALSE
          measure.pred[measure.pred$comp == comp
                       ,]$optimum.keepA <- FALSE
          ## update optimum
          measure.pred[measure.pred$comp == comp &
                         measure.pred$keepX == trial.keepX &
                         measure.pred$keepY == trial.keepY
                       ,]$optimum.keepA <- TRUE
        }
      }
      
      return(measure.pred)
    }
    
    # generate list of lists containing all keepX/Y pairs to test
    keep.vals <- list()
    i <- 1
    for(keepX in test.keepX){
      for(keepY in test.keepY){
        keep.vals[[i]] <- list(keepX = keepX, keepY=keepY)
        i<-i+1
      }
    }
    
    for (comp in seq_len(ncomp)){ # for each component to test...
      if (use_progressBar) {
        n_tested <- 0
        cat(sprintf("\ntuning component: %s\n", comp))
      }
      
      # calculate measure values for specific comp, keepX/Y
      m.p <- bplapply(X=keep.vals, FUN = iter_keep.vals, 
                      BPPARAM = BPPARAM) 
      
      for (d in m.p) { # output will be a list, so adjust global measure.pred
                       # to contain the desired values
        measure.pred[rownames(d), ] <- d
      }
      
      # changes measure.pred$optimum.keepA to reflect results of t.tests
      measure.pred <- DetermineOptimalKeepVals(measure.pred, comp, keep.vals, nrepeat)
      
      # extract optimal keepX for this component
      choice.keepX.ncomp <-  measure.pred[measure.pred$comp == comp & 
                                            measure.pred$optimum.keepA == TRUE &
                                            measure.pred$measure == measure & 
                                            measure.pred$V == 't' ## doesn't matter t or u
                                          ,]$keepX
      # extract optimal keepY for this component
      choice.keepY.ncomp <-  measure.pred[measure.pred$comp == comp & 
                                            measure.pred$optimum.keepA == TRUE &
                                            measure.pred$measure == measure &
                                            measure.pred$V == 't'
                                          ,]$keepY
      choice.keepX = c(choice.keepX, choice.keepX.ncomp)
      choice.keepY = c(choice.keepY, choice.keepY.ncomp)
      
    } # end comp
    
    choice.ncomp <- 1
    for (comp in seq_len(ncomp))
    {
      Q2.total <- measure.pred[measure.pred$comp == comp & 
                                 measure.pred$keepX == keepX &
                                 measure.pred$keepY == keepY &
                                 measure.pred$V == 'u'
                               ,]$value.Q2.total[[1]]
      Q2.opt <- measure.pred[measure.pred$comp == choice.ncomp & 
                               measure.pred$keepX == keepX &
                               measure.pred$keepY == keepY &
                               measure.pred$V == 'u'
                             ,]$value.Q2.total[[1]]
      keep.comp <- mean(Q2.total) >= limQ2
      if (keep.comp)
        # we want the first component that drops Q2 below limQ2
        choice.ncomp <- comp + 1
    }
    names(choice.keepX) <- names(choice.keepY) <- paste0('comp', seq_len(ncomp))
    out$choice.keepX = choice.keepX
    out$choice.keepY = choice.keepY
    out$choice.ncomp = choice.ncomp
    
    ## add mean and sd
    val <- unclass(measure.pred$value.v)
    measure.pred$mean <- sapply(measure.pred$value.v, function(x){
      mean(c(as.matrix(x)), na.rm = TRUE)
    })
    
    measure.pred$sd <- sapply(measure.pred$value.v, function(x){
      sd(c(as.matrix(x)), na.rm = TRUE)
    })
    
    out$measure.pred = measure.pred
    ### evaluate all for output except X and Y to save memory
    ## eval all but X and Y
    mc <- mget(names(formals())[-1:-2], sys.frame(sys.nframe()))
    ## replace function, X and Y with unevaluated call
    mc <- as.call(c(as.list(match.call())[1:3], mc))
    out <- c(list(call = mc), out)
    class <- c('tune.pls')
    class(out) <- if (spls.model) class <- c(class, 'tune.spls')
    return(out)
    }
  }

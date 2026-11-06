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
# tune.pls: Tuning hyperparameters on a pls method
# ========================================================================================================
#' 
#' Tuning functions for PLS method
#' 
#' Computes M-fold or Leave-One-Out Cross-Validation scores on a user-input
#' grid to determine optimal values for the parameters in \code{spls}.
#'
#' 
#' This tuning function should be used to tune the number of components to select for spls models.
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
#' @param ... Optional parameters passed to \code{\link{pls}}
#' @return 
#' Returns a list with the following components for every repeat:
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
#' @author Kim-Anh Lê Cao, Al J Abadi, Benoit Gautier, Francois Bartolo and Florian Rohart.
#' @seealso \code{\link{splsda}}, \code{\link{predict.splsda}}, and http://www.mixOmics.org for more details.
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
#' Sparse PLS regression mode:
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
#' 
#' @keywords regression multivariate
#' @export
#' @example ./examples/tune.pls-examples.R
#' 
tune.pls <- 
  function(X,
           Y,
           ncomp,
           # params related to CV
           validation = c('Mfold', 'loo'),
           nrepeat = 1,
           folds,
           measure = NULL,
           # params related to spls model building
           mode = c('regression', 'canonical', 'classic'),
           scale = TRUE,
           logratio = "none",
           tol = 1e-06,
           max.iter = 100,
           near.zero.var = FALSE,
           multilevel = NULL,
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
    
    #-- build model and run perf() --------------------------------------#
    #---------------------------------------------------------------------------#
    pls_res <- pls(X, Y, ncomp,
                  mode = mode, scale = scale, logratio = logratio, tol = tol, max.iter = max.iter, near.zero.var = near.zero.var, multilevel = multilevel)
    perf_res <- perf(pls_res, 
                validation = validation, folds = folds, nrepeat = nrepeat,
                BPPARAM = BPPARAM, seed = seed, progressBar = progressBar)
    return(perf_res)

    } 
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

#' Parallelized Tuning function for multilevel sPLS method using BiocParallel
#' 
#' For a multilevel spls analysis, the tuning criterion is based on the
#' maximisation of the correlation between the components from both data sets
#' 
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param Y \code{if(method = 'spls')} numeric vector or matrix of continuous
#' responses (for multi-response models) \code{NA}s are allowed.
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements) that indicates the repeated measures on each individual, i.e.
#' the individuals ID. See Details.
#' @param ncomp the number of components to include in the model.
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}.
#' @param test.keepX numeric vector for the different number of variables to
#' test from the \eqn{X} data set
#' @param test.keepY numeric vector for the different number of variables to
#' test from the \eqn{Y} data set
#' @param already.tested.X Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{X} data set on
#' the firsts components.
#' @param already.tested.Y Optional, if \code{ncomp > 1} A numeric vector
#' indicating the number of variables to select from the \eqn{Y} data set on
#' the firsts components.
#' @param BPPARAM BiocParallelParam object to manage parallelization
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' @return
#' \item{cor.value}{correlation between latent variables}
#' @export
tune.splslevel <- function (X, Y,
                            multilevel,
                            ncomp = NULL,
                            mode = "regression",
                            test.keepX = rep(ncol(X), ncomp),
                            test.keepY = rep(ncol(Y), ncomp),
                            already.tested.X = NULL,
                            already.tested.Y = NULL,
                            BPPARAM = BiocParallel::SerialParam(),
                            seed = NULL) {
  
  BPPARAM$RNGseed <- seed
  set.seed(seed)

    message("For a multilevel spls analysis, the tuning criterion is based on the maximisation of the correlation between the components from both data sets")
    
    Y <- as.matrix(Y)
    if (length(dim(Y)) != 2 || !is.numeric(Y))
        stop("'Y' must be a numeric matrix.")
    
    if (!is.null(already.tested.X) && is.null(already.tested.Y))
        stop("Input already.tested.Y is missing")
    
    if (!is.null(already.tested.Y) && is.null(already.tested.X))
        stop("Input already.tested.X is missing")
    
    if (length(already.tested.X) != (ncomp - 1))
        stop("The number of already.tested.X parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
    
    if (length(already.tested.Y) != (ncomp - 1))
        stop("The number of already.tested.Y parameters should be ", ncomp - 1, " since you set ncomp = ", ncomp)
    
    if ((!is.null(already.tested.X)) && (!is.numeric(already.tested.X)))
        stop("Expecting a numerical value in already.tested.X", call. = FALSE)
    
    if ((!is.null(already.tested.Y)) && (!is.numeric(already.tested.Y)))
        stop("Expecting a numerical value in already.tested.Y", call. = FALSE)
    
    Xw <- suppressMessages(withinVariation(X = X, design = multilevel))
    Yw <- suppressMessages(withinVariation(X = Y, design = multilevel))
    
    param_grid <- expand.grid(test.keepX = test.keepX, test.keepY = test.keepY, KEEP.OUT.ATTRS = FALSE)

    cor.values <- BiocParallel::bplapply(1:nrow(param_grid), function(idx) {
        i <- param_grid$test.keepX[idx]
        j <- param_grid$test.keepY[idx]
        
        if (ncomp == 1) {
            spls.train <- spls(Xw, Yw, ncomp = ncomp,
                               keepX = i, keepY = j,
                               mode = mode)
        } else {
            spls.train <- spls(Xw, Yw, ncomp = ncomp,
                               keepX = c(already.tested.X, i),
                               keepY = c(already.tested.Y, j),
                               mode = mode)
        }
        
        cor(spls.train$variates$X[, ncomp], spls.train$variates$Y[, ncomp])
    }, BPPARAM = BPPARAM)
    
    cor.value <- matrix(unlist(cor.values), nrow = length(test.keepX), ncol = length(test.keepY),
                        dimnames = list(paste("varX ", test.keepX, sep = ""),
                                        paste("varY ", test.keepY, sep = "")))
    
    return(list(cor.value = cor.value))
}

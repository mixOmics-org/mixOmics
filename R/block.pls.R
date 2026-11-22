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

# =============================================================================
# block.pls: perform a horizontal PLS on a combination of datasets,
# input as a list in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint.block
# =============================================================================

#' N-integration with Projection to Latent Structures models (PLS)
#' 
#' Integration of multiple data sets measured on the same samples or
#' observations, ie. N-integration. The method is partly based on Generalised
#' Canonical Correlation Analysis.
#' 
#' \code{block.spls} function fits a horizontal integration PLS model with a
#' specified number of components per block). An outcome needs to be provided,
#' either by \code{Y} or by its position \code{indY} in the list of blocks
#' \code{X}. Multi (continuous)response are supported. \code{X} and \code{Y}
#' can contain missing values. Missing values are handled by being disregarded
#' during the cross product computations in the algorithm \code{block.pls}
#' without having to delete rows with missing data. Alternatively, missing data
#' can be imputed prior using the \code{\link{impute.nipals}} function.
#' 
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and \code{?pls} for more details). Note that the argument 'scheme' 
#' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' 
#' Note that our method is partly based on Generalised Canonical Correlation
#' Analysis and differs from the MB-PLS approaches proposed by Kowalski et al.,
#' 1989, J Chemom 3(1) and Westerhuis et al., 1998, J Chemom, 12(5).
#' 
#' @inheritParams pls
#' @param X A named list of data sets (called 'blocks') measured on the same 
#' samples. Data in the list should be arranged in matrices, samples x variables, 
#' with samples order matching in all data sets.
#' @param Y Matrix response for a multivariate regression framework. Data
#' should be continuous variables (see \code{?block.plsda} for supervised
#' classification and factor response).
#' @param indY To supply if \code{Y} is missing, indicates the position of 
#' the matrix response in the list \code{X}.
#' @param ncomp the number of components to include in the model. Default to 2.
#' Applies to all blocks.
#' @param design numeric matrix of size (number of blocks in X) x (number of
#'   blocks in X) with values between 0 and 1. Each value indicates the strenght
#'   of the relationship to be modelled between two blocks; a value of 0
#'   indicates no relationship, 1 is the maximum value. Alternatively, one of
#'   c('null', 'full') indicating a disconnected or fully connected design,
#'   respecively, or a numeric between 0 and 1 which will designate all
#'   off-diagonal elements of a fully connected design (see examples in
#'   \code{block.splsda}). If \code{Y} is provided instead of \code{indY}, the
#'   \code{design} matrix is changed to include relationships to \code{Y}.
#' @template arg/verbose.call
#' @return \code{block.pls} returns an object of class \code{'block.pls'}, a
#' list that contains the following components:
#' 
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{variates}{list
#' containing the variates of each block of X.} \item{loadings}{list containing
#' the estimated loadings for the variates.} \item{names}{list containing the
#' names to be used for individuals and variables.} \item{nzv}{list containing
#' the zero- or near-zero predictors information.} \item{iter}{Number of
#' iterations of the algorithm for each component}
#' \item{prop_expl_var}{Percentage of explained variance for each
#' component and each block}
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{plotIndiv}}, \code{\link{plotArrow}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.spls}},
#' \code{\link{block.plsda}} and https://mixomics.org/ for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' 
#' Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical
#' Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' @keywords regression multivariate
#' @example ./examples/block.pls-examples.R
#' @export
block.pls <- function(X,
                      Y,
                      indY,
                      ncomp = 2,
                      design,
                      mode,
                      scale = TRUE,
                      tol = 1e-06,
                      max.iter = 100,
                      near.zero.var = FALSE,
                      all.outputs = TRUE,
                      verbose.call = FALSE)
{
    
    # call to 'internal_wrapper.mint.block'
    result = internal_wrapper.mint.block(X=X, Y=Y, indY=indY, ncomp=ncomp,
                                         design=design, scheme="horst", mode=mode, scale=scale,
                                         init="svd.single", tol=tol, max.iter=max.iter ,near.zero.var=near.zero.var,
                                         all.outputs = all.outputs, DA = FALSE)
    
    # calculate weights for each dataset
    weights = get.weights(result$variates, indY = result$indY)
    
    # choose the desired output from 'result'
    out=list(call = match.call(),
             X = result$A,
             indY = result$indY,
             ncomp = result$ncomp,
             mode = result$mode,
             variates = result$variates,
             loadings = result$loadings,
             crit = result$crit,
             AVE = result$AVE,
             names = result$names,
             init = result$init,
             tol = result$tol,
             iter = result$iter,
             max.iter = result$max.iter,
             nzv = result$nzv,
             scale = result$scale,
             design = result$design,
             scheme = result$scheme,
             weights = weights,
             prop_expl_var = result$prop_expl_var)
    
    if (verbose.call) {
        c <- out$call
        out$call <- mget(names(formals()))
        out$call <- append(c, out$call)
        names(out$call)[1] <- "simple.call"
    }
    
    # give a class
    class(out) = c("block.pls","sgcca")
    
    return(invisible(out))
    
}

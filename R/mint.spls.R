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
# mint.spls: perform a vertical sPLS on a combination of experiments, input as a matrix in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint, which then call 'internal_mint.block'
# ========================================================================================================

#' P-integration with variable selection
#' 
#' Function to integrate and combine multiple independent studies measured on
#' the same variables or predictors (P-integration) using variants of
#' multi-group sparse PLS for variable selection (unsupervised analysis).
#' 
#' \code{mint.spls} fits a vertical sparse PLS-DA models with \code{ncomp}
#' components in which several independent studies measured on the same
#' variables are integrated. The aim is to explain the continuous outcome
#' \code{Y} and selecting correlated features between both data sets \code{X}
#' and \code{Y}. The \code{study} factor indicates the membership of each
#' sample in each study. We advise to only combine studies with more than 3
#' samples as the function performs internal scaling per study.
#' 
#' Multi (continuous)response are supported. \code{X} and \code{Y} can contain
#' missing values. Missing values are handled by being disregarded during the
#' cross product computations in the algorithm \code{mint.spls} without having
#' to delete rows with missing data. Alternatively, missing data can be imputed
#' prior using the \code{nipals} function.
#' 
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and more details in \code{?pls}).
#' 
#' Variable selection is performed on each component for each block of
#' \code{X}, and for \code{Y} if specified, via input parameter \code{keepX}
#' and \code{keepY}.
#' 
#' Useful graphical outputs are available, e.g. \code{\link{plotIndiv}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}.
#' 
#' @inheritParams mint.pls
#' @param Y Matrix or vector response for a multivariate regression framework.
#' Data should be continuous variables (see \code{mint.splsda} for supervised
#' classification and factor response)
#' @param keepX numeric vector indicating the number of variables to select in
#' \code{X} on each component. By default all variables are kept in the model.
#' @param keepY numeric vector indicating the number of variables to select in
#' \code{Y} on each component. By default all variables are kept in the model.
#' @template arg/verbose.call
#' @return \code{mint.spls} returns an object of class
#' \code{"mint.spls","spls"}, a list that contains the following components:
#' 
#' \item{X}{numeric matrix of predictors combining multiple independent studies
#' on the same set of predictors. \code{NA}s are allowed.} \item{Y}{the
#' centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{study}{The study grouping factor} \item{mode}{the algorithm used to
#' fit the model.} \item{keepX}{Number of variables used to build each
#' component of X} \item{keepY}{Number of variables used to build each
#' component of Y} \item{variates}{list containing the variates of X - global
#' variates.} \item{loadings}{list containing the estimated loadings for the
#' variates - global loadings.} \item{variates.partial}{list containing the
#' variates of X relative to each study - partial variates.}
#' \item{loadings.partial}{list containing the estimated loadings for the
#' partial variates - partial loadings.} \item{names}{list containing the names
#' to be used for individuals and variables.} \item{nzv}{list containing the
#' zero- or near-zero predictors information.} \item{iter}{Number of iterations
#' of the algorithm for each component} \item{prop_expl_var}{The amount
#' of the variance explained by each variate / component divided by the total
#' variance in the \code{data} for each study (after removing the possible
#' missing values) using the definition of 'redundancy'. Note that contrary to
#' \code{PCA}, this amount may not decrease in the following components as the
#' aim of the method is not to maximise the variance, but the covariance between
#' data sets (including the dummy matrix representation of the outcome variable
#' in case of the supervised approaches).}
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' @author Florian Rohart, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}},
#' \code{\link{mint.pls}}, \code{\link{mint.plsda}}, \code{\link{mint.splsda}}
#' and http://www.mixOmics.org/mixMINT for more details.
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#' 
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms
#' for multi-group PLS. J. Chemometrics, 28(3), 192-201.
#' @keywords regression multivariate
#' @export
#' @examples
#' 
#' data(stemcells)
#' 
#' # for the purpose of this example, we artificially
#' # create a continuous response Y by taking gene 1.
#' 
#' res = mint.spls(X = stemcells$gene[,-1], Y = stemcells$gene[,1], ncomp = 3,
#' keepX = c(10, 5, 15), study = stemcells$study)
#' 
#' plotIndiv(res)
#' 
#' #plot study-specific outputs for all studies
#' plotIndiv(res, study = "all.partial")
#' 
#' \dontrun{
#' #plot study-specific outputs for study "2"
#' plotIndiv(res, study = "2", col = 1:3, legend = TRUE)
#' }
mint.spls <- function(X,
                      Y,
                      ncomp = 2,
                      mode = c("regression", "canonical", "invariant", "classic"),
                      study,
                      keepX = rep(ncol(X), ncomp),
                      keepY = rep(ncol(Y), ncomp),
                      scale = TRUE,
                      tol = 1e-06,
                      max.iter = 100,
                      near.zero.var = FALSE,
                      all.outputs = TRUE,
                      verbose.call = FALSE)
{
    
    # call to 'internal_wrapper.mint'
    result <- internal_wrapper.mint(
        X = X,
        Y = Y,
        ncomp = ncomp,
        scale = scale,
        near.zero.var = near.zero.var,
        study = study,
        mode = mode,
        keepX = keepX,
        keepY = keepY,
        max.iter = max.iter,
        tol = tol,
        all.outputs = all.outputs,
        DA = FALSE
    )
    
    # choose the desired output from 'result'
    out <- list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        study = result$study,
        mode = result$mode,
        keepX = result$keepX,
        keepY = result$keepY,
        variates = result$variates,
        loadings = result$loadings,
        variates.partial = result$variates.partial,
        loadings.partial = result$loadings.partial,
        names  =  result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        prop_expl_var = result$prop_expl_var
    )
    
    if (verbose.call) {
        c <- out$call
        out$call <- mget(names(formals()))
        out$call <- append(c, out$call)
        names(out$call)[1] <- "simple.call"
    }
    
    class(out) <- c("mint.spls","mixo_spls")
    return(invisible(out))
    
}

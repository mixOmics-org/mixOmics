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
# pls: perform a PLS
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

#' Partial Least Squares (PLS) Regression
#' 
#' Function to perform Partial Least Squares (PLS) regression.
#' 
#' \code{pls} function fit PLS models with \eqn{1, \ldots ,}\code{ncomp}
#' components. Multi-response models are fully supported. The \code{X} and
#' \code{Y} datasets can contain missing values.
#' 
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References). Different modes relate on how the Y matrix is deflated across
#' the iterations of the algorithms - i.e. the different components.
#' 
#' - Regression mode: the Y matrix is deflated with respect to the information
#' extracted/modelled from the local regression on X. Here the goal is to
#' predict Y from X (Y and X play an asymmetric role). Consequently the latent
#' variables computed to predict Y from X are different from those computed to
#' predict X from Y.
#' 
#' - Canonical mode: the Y matrix is deflated to the information
#' extracted/modelled from the local regression on Y. Here X and Y play a
#' symmetric role and the goal is similar to a Canonical Correlation type of
#' analysis.
#' 
#' - Invariant mode: the Y matrix is not deflated
#' 
#' - Classic mode: is similar to a regression mode. It gives identical results
#' for the variates and loadings associated to the X data set, but differences
#' for the loadings vectors associated to the Y data set (different
#' normalisations are used). Classic mode is the PLS2 model as defined by
#' Tenenhaus (1998), Chap 9.
#' 
#' Note that in all cases the results are the same on the first component as
#' deflation only starts after component 1.
#' 
#' @section missing values: 
#' The estimation of the missing values can be performed using the
#' \code{\link{impute.nipals}} function. Otherwise, missing values are handled
#' by element-wise deletion in the \code{pls} function without having to delete
#' the rows with missing data.
#' 
#' @section multilevel:
#' Multilevel (s)PLS enables the integration of data measured on two different
#' data sets on the same individuals. This approach differs from multilevel
#' sPLS-DA as the aim is to select subsets of variables from both data sets that
#' are highly positively or negatively correlated across samples. The approach
#' is unsupervised, i.e. no prior knowledge about the sample groups is included.
#' 
#' @section logratio and multilevel: 
#' logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#' 
#' @template arg/X.matrix
#' @template arg/Y.matrix
#' @template arg/ncomp
#' @templateVar scale.default TRUE
#' @template arg/scale
#' @templateVar modes \code{"regression"}, \code{"canonical"}, \code{"invariant"} or \code{"classic"}
#' @template arg/mode
#' @template arg/tol
#' @templateVar max.iter.default 100
#' @template arg/max.iter
#' @template arg/near.zero.var
#' @templateVar logratio.values.allowed ('none','CLR')
#' @template arg/logratio
#' @templateVar multilevel.example in \code{?spls}.
#' @template arg/multilevel
#' @template arg/all.outputs
#' @template arg/verbose.call
#' @return \code{pls} returns an object of class \code{"pls"}, a list that
#' contains the following components:
#' 
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{mode}{the algorithm used to fit the model.} \item{variates}{list
#' containing the variates.} \item{loadings}{list containing the estimated
#' loadings for the \eqn{X} and \eqn{Y} variates. The loading weights multiplied with their associated deflated (residual) matrix gives the variate.} \item{loadings.stars}{list containing the estimated
#' weighted loadings for the \eqn{X} and \eqn{Y} variates. The loading weights are projected so that when multiplied with their associated original matrix we obtain the variate.} \item{names}{list containing
#' the names to be used for individuals and variables.} \item{tol}{the
#' tolerance used in the iterative algorithm, used for subsequent S3 methods}
#' \item{iter}{Number of iterations of the algorithm for each component}
#' \item{max.iter}{the maximum number of iterations, used for subsequent S3
#' methods} \item{nzv}{list containing the zero- or near-zero predictors
#' information.} \item{scale}{whether scaling was applied per predictor.}
#' \item{logratio}{whether log ratio transformation for relative proportion
#' data was applied, and if so, which type of transformation.}
#' \item{prop_expl_var}{The proportion of the variance explained by each
#' variate / component divided by the total variance in the \code{data} (after
#' removing the possible missing values) using the definition of 'redundancy'.
#' Note that contrary to \code{PCA}, this amount may not decrease in the
#' following components as the aim of the method is not to maximise the
#' variance, but the covariance between data sets (including the dummy matrix
#' representation of the outcome variable in case of the supervised
#' approaches).}
#' \item{input.X}{numeric matrix of predictors in X that was input, before any
#' scaling / logratio / multilevel transformation.} \item{mat.c}{matrix of
#' coefficients from the regression of X / residual matrices X on the
#' X-variates, to be used internally by \code{predict}.}
#' \item{defl.matrix}{residual matrices X for each dimension.}
#' @author Sébastien Déjean, Ignacio González, Florian Rohart, Kim-Anh Lê Cao,
#' Al J Abadi
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}} and
#' http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' 
#' Abdi H (2010). Partial least squares regression and projection on latent
#' structure regression (PLS Regression). \emph{Wiley Interdisciplinary
#' Reviews: Computational Statistics}, 2(1), 97-106.
#' @keywords regression multivariate
#' @export
#' @examples
#' data(linnerud)
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.pls <- pls(X, Y, mode = "classic")
#' 
#' \dontrun{
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' toxicity.pls <- pls(X, Y, ncomp = 3)
#' }
#' 
pls <- function(X,
                Y,
                ncomp = 2,
                scale = TRUE,
                mode = c("regression", "canonical", "invariant", "classic"),
                tol = 1e-06,
                max.iter = 100,
                near.zero.var = FALSE,
                logratio = "none",
                # one of "none", "CLR"
                multilevel = NULL,
                all.outputs = TRUE,
                verbose.call = FALSE)
{
    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(
        X = X,
        Y = Y,
        ncomp = ncomp,
        scale = scale,
        near.zero.var = near.zero.var,
        mode = mode,
        max.iter = max.iter,
        tol = tol,
        logratio = logratio,
        multilevel = multilevel,
        DA = FALSE,
        all.outputs = all.outputs
    )
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        variates = result$variates,
        loadings = result$loadings,
        loadings.star = result$loadings.star,
        names = result$names,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = scale,
        logratio = logratio,
        prop_expl_var = result$prop_expl_var,
        input.X = result$input.X,
        mat.c = result$mat.c#,
        #defl.matrix = result$defl.matrix
    )
    
    if (verbose.call) {
        c <- out$call
        out$call <- mget(names(formals()))
        out$call <- append(c, out$call)
        names(out$call)[1] <- "simple.call"
    }
    
    class(out) = c("mixo_pls")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mixo_mlpls",class(out))
    }
    
    return(invisible(out))
    
}

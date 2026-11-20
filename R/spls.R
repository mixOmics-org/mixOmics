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
# spls: perform a sPLS
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

#' Sparse Partial Least Squares (sPLS)
#' 
#' Function to perform sparse Partial Least Squares (sPLS). The sPLS approach
#' combines both integration and variable selection simultaneously on two data
#' sets in a one-step strategy.
#' 
#' \code{spls} function fit sPLS models with \eqn{1, \ldots ,}\code{ncomp}
#' components. Multi-response models are fully supported. The \code{X} and
#' \code{Y} datasets can contain missing values.
#' 
#' @inherit pls
#' @templateVar multilevel.example .
#' @template arg/multilevel
#' @param keepX numeric vector of length \code{ncomp}, the number of variables
#' to keep in \eqn{X}-loadings. By default all variables are kept in the model.
#' @param keepY numeric vector of length \code{ncomp}, the number of variables
#' @return \code{spls} returns an object of class \code{"spls"}, a list that
#' contains the following components:
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized original response vector or matrix.}
#' \item{ncomp}{the number of components included in the model.}
#' \item{mode}{the algorithm used to fit the model.} \item{keepX}{number of
#' \eqn{X} variables kept in the model on each component.} \item{keepY}{number
#' of \eqn{Y} variables kept in the model on each component.}
#' \item{variates}{list containing the variates.} \item{loadings}{list
#' containing the estimated loadings for the \eqn{X} and \eqn{Y} variates.}
#' \item{names}{list containing the names to be used for individuals and
#' variables.} \item{tol}{the tolerance used in the iterative algorithm, used
#' for subsequent S3 methods} \item{iter}{Number of iterations of the algorithm
#' for each component} \item{max.iter}{the maximum number of iterations, used
#' for subsequent S3 methods} \item{nzv}{list containing the zero- or near-zero
#' predictors information.} \item{scale}{whether scaling was applied per
#' predictor.} \item{logratio}{whether log ratio transformation for relative
#' proportion data was applied, and if so, which type of transformation.}
#' \item{prop_expl_var}{Proportion of variance explained per component (note
#' that contrary to PCA, this amount may not decrease as the aim of the method
#' is not to maximise the variance, but the covariance between data sets).}
#' \item{input.X}{numeric matrix of predictors in X that was input, before any
#' saling / logratio / multilevel transformation.} \item{mat.c}{matrix of
#' coefficients from the regression of X / residual matrices X on the
#' X-variates, to be used internally by \code{predict}.}
#' \item{defl.matrix}{residual matrices X for each dimension.}
#' @author Sébastien Déjean, Ignacio González, Florian Rohart, Kim-Anh Lê Cao,
#' Al J abadi
#' @seealso \code{\link{pls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{cim}}, \code{\link{network}},
#' \code{\link{predict}}, \code{\link{perf}} and http://www.mixOmics.org for
#' more details.
#' @references Sparse PLS: canonical and regression modes:
#' 
#' Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009).
#' Sparse canonical methods for biological data integration: application to a
#' cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
#' 
#' Lê Cao, K.-A., Rossouw, D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. \emph{Statistical
#' Applications in Genetics and Molecular Biology} \bold{7}, article 35.
#' 
#' Sparse SVD: Shen, H. and Huang, J. Z. (2008). Sparse principal component
#' analysis via regularized low rank matrix approximation. \emph{Journal of
#' Multivariate Analysis} \bold{99}, 1015-1034.
#' 
#' PLS methods: Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic. Chapters 9 and 11.
#' 
#' Abdi H (2010). Partial least squares regression and projection on latent
#' structure regression (PLS Regression). \emph{Wiley Interdisciplinary
#' Reviews: Computational Statistics}, 2(1), 97-106.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' 
#' On multilevel analysis:
#' 
#' Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel
#' approach for biomarker selection and the integration of repeated measures
#' experiments from two platforms. \emph{BMC Bioinformatics} \bold{13}:325.
#' 
#' Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K.
#' (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA.
#' \emph{Metabolomics}, \bold{6}(1), 119-128.
#' @keywords regression multivariate
#' @export
#' @example ./examples/spls-examples.R
spls <- function(X,
                 Y,
                 ncomp = 2,
                 mode = c("regression", "canonical", "invariant", "classic"),
                 keepX,
                 keepY,
                 scale = TRUE,
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
    result <- internal_wrapper.mint(
        X = X,
        Y = Y,
        ncomp = ncomp,
        scale = scale,
        near.zero.var = near.zero.var,
        mode = mode,
        keepX = keepX,
        keepY = keepY,
        max.iter = max.iter,
        tol = tol,
        logratio = logratio,
        multilevel = multilevel,
        DA = FALSE,
        all.outputs = all.outputs
    )
    
    # choose the desired output from 'result'
    out <- list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        keepX = result$keepX,
        keepY = result$keepY,
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
        mat.c = result$mat.c
    )
    
    if (verbose.call) {
        c <- out$call
        out$call <- mget(names(formals()))
        out$call <- append(c, out$call)
        names(out$call)[1] <- "simple.call"
    }
    
    
    class(out) = c("mixo_spls")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mixo_mlspls",class(out))
    }
    
    return(invisible(out))
}



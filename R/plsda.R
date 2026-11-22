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
# plsda: perform a PLS-DA
# this function is a particular setting of internal_mint.block, the formatting of the input is checked in internal_wrapper.mint
# ========================================================================================================

#' Partial Least Squares Discriminant Analysis (PLS-DA).
#' 
#' Function to perform standard Partial Least Squares regression to classify
#' samples.
#' 
#' \code{plsda} function fit PLS models with \eqn{1,...,}\code{ncomp} components
#' to the factor or class vector \code{Y}. The appropriate indicator matrix is created.
#' 
#' Logratio transformation and multilevel analysis are
#' performed sequentially as internal pre-processing step, through
#' \code{\link{logratio.transfo}} and \code{\link{withinVariation}}
#' respectively. Logratio can only be applied if the data do not contain any 0 value (for
#' count data, we thus advise the normalise raw data with a 1 offset).
#' 
#' The type of deflation used is \code{'regression'} for discriminant algorithms.
#' i.e. no deflation is performed on Y.
#' 
#' @inheritParams pls
#' @param Y a factor or a class vector for the discrete outcome.
#' @param multilevel sample information for multilevel decomposition for
#' repeated measurements. A numeric matrix or data frame indicating the
#' repeated measures on each individual, i.e. the individuals ID. See examples
#' in \code{?splsda}.
#' @return \code{plsda} returns an object of class \code{"plsda"}, a list that
#' contains the following components:
#' 
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized indicator response vector or matrix.}
#' \item{ind.mat}{the indicator matrix.} \item{ncomp}{the number of components
#' included in the model.} \item{variates}{list containing the \code{X} and
#' \code{Y} variates.} \item{loadings}{list containing the estimated loadings
#' associated to each component/variate. The loading weights multiplied with the
#' deflated (residual) matrix gives the variate.} \item{loadings.stars}{list
#' containing the estimated loadings associated to each component/variate. The
#' loading weights are projected so that when multiplied with the original
#' matrix we obtain the variate.} \item{names}{list containing the names to be
#' used for individuals and variables.} \item{nzv}{list containing the zero- or
#' near-zero predictors information.} \item{tol}{the tolerance used in the
#' iterative algorithm, used for subsequent S3 methods} \item{max.iter}{the
#' maximum number of iterations, used for subsequent S3 methods}
#' \item{iter}{Number of iterations of the algorithm for each component}
#' \item{prop_expl_var}{The proportion of the variance explained by each
#' variate / component divided by the total variance in the \code{data} (after
#' removing the possible missing values) using the definition of 'redundancy'.
#' Note that contrary to \code{PCA}, this amount may not decrease in the
#' following components as the aim of the method is not to maximise the
#' variance, but the covariance between data sets (including the dummy matrix
#' representation of the outcome variable in case of the supervised
#' approaches).}\item{mat.c}{matrix of coefficients from the regression of X /
#' residual matrices X on the X-variates, to be used internally by
#' \code{predict}.} \item{defl.matrix}{residual matrices X for each dimension.}
#' @author Ignacio González, Kim-Anh Lê Cao, Florian Rohart, Al J Abadi
#' @seealso \code{\link{splsda}}, \code{\link{summary}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{mint.block.plsda}},
#' \code{\link{block.plsda}} and https://mixomics.org/ for more details.
#' @references On PLSDA: Barker M and Rayens W (2003). Partial least squares
#' for discrimination. \emph{Journal of Chemometrics} \bold{17}(3), 166-173.
#' Perez-Enciso, M. and Tenenhaus, M. (2003). Prediction of clinical outcome
#' with microarray data: a partial least squares discriminant analysis (PLS-DA)
#' approach. \emph{Human Genetics} \bold{112}, 581-592. Nguyen, D. V. and
#' Rocke, D. M. (2002). Tumor classification by partial least squares using
#' microarray gene expression data. \emph{Bioinformatics} \bold{18}, 39-50. On
#' log ratio transformation: Filzmoser, P., Hron, K., Reimann, C.: Principal
#' component analysis for compositional data with outliers. Environmetrics
#' 20(6), 621-632 (2009) Lê Cao K.-A., Costello ME, Lakis VA, Bartolo, F,Chua
#' XY, Brazeilles R, Rondeau P. MixMC: Multivariate insights into Microbial
#' Communities. PLoS ONE, 11(8): e0160169 (2016). On multilevel decomposition:
#' Westerhuis, J.A., van Velzen, E.J., Hoefsloot, H.C., Smilde, A.K.:
#' Multivariate paired data analysis: multilevel plsda versus oplsda.
#' Metabolomics 6(1), 119-128 (2010) Liquet, B., Lê Cao K.-A., Hocini, H.,
#' Thiebaut, R.: A novel approach for biomarker selection and the integration
#' of repeated measures experiments from two assays. BMC bioinformatics 13(1),
#' 325 (2012)
#' @keywords regression multivariate
#' @export
#' @example ./examples/plsda-examples.R
plsda <- function(X,
                  Y,
                  ncomp = 2,
                  scale = TRUE,
                  tol = 1e-06,
                  max.iter = 100,
                  near.zero.var = FALSE,
                  logratio = c("none", "CLR"),
                  multilevel = NULL,
                  all.outputs = TRUE)
{
    
    #-- validation des arguments --#
    # most of the checks are done in the wrapper.mint.spls.hybrid function
    if (is.null(multilevel))
    {
        if (is.null(Y))
            stop("'Y' has to be something else than NULL.")
        
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
            stop("'Y' should be a factor with more than one level")
        
        Y.mat = unmap(Y)
        colnames(Y.mat) = levels(Y)
        
    } else {
        # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
        multilevel = data.frame(multilevel)
        
        if ((nrow(X) != nrow(multilevel)))
            stop("unequal number of rows in 'X' and 'multilevel'.")
        
        if (ncol(multilevel) != 1)
            stop("'multilevel' should have a single column for the repeated measurements, other factors should be included in 'Y'.")
        
        if (!is.null(ncol(Y)) && !ncol(Y) %in% c(0,1,2))# multilevel 1 or 2 factors
            stop("'Y' should either be a factor, a single column data.frame containing a factor, or a 2-columns data.frame containing 2 factors.")
        multilevel = data.frame(multilevel, Y)
        multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
        
        Y.mat = NULL
    }
    
    logratio <- match.arg(logratio)
    
    # call to 'internal_wrapper.mint'
    result = internal_wrapper.mint(
        X = X,
        Y = Y.mat,
        ncomp = ncomp,
        scale = scale,
        near.zero.var = near.zero.var,
        mode = 'regression',
        max.iter = max.iter,
        tol = tol,
        logratio = logratio,
        multilevel = multilevel,
        DA = TRUE,
        all.outputs = all.outputs
    )
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY][[1]],
        Y = if (is.null(multilevel))
        {
            Y
        } else {
            result$Y.factor
        }, 
        ind.mat = result$A[result$indY][[1]],
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
        mat.c = result$mat.c
    )
    
    class(out) = c("mixo_plsda","mixo_pls","DA")
    # output if multilevel analysis
    if (!is.null(multilevel))
    {
        out$multilevel = multilevel
        class(out) = c("mixo_mlplsda",class(out))
    }
    
    return(invisible(out))
}


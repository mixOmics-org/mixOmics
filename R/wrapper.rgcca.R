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

#' mixOmics wrapper for Regularised Generalised Canonical Correlation Analysis
#' (rgcca)
#' 
#' Wrapper function to perform Regularized Generalised Canonical Correlation
#' Analysis (rGCCA), a generalised approach for the integration of multiple
#' datasets. For more details, see the \code{help(rgcca)} from the \pkg{RGCCA}
#' package.
#' 
#' This wrapper function performs rGCCA (see \pkg{RGCCA}) with \eqn{1, \ldots
#' ,}\code{ncomp} components on each block data set. A supervised or
#' unsupervised model can be run. For a supervised model, the
#' \code{\link{unmap}} function should be used as an input data set. More
#' details can be found on the package \pkg{RGCCA}.
#' 
#' @param X a list of data sets (called 'blocks') matching on the same samples.
#' Data in the list should be arranged in samples x variables. \code{NA}s are
#' not allowed.
#' @param design numeric matrix of size (number of blocks in X) x (number of
#' blocks in X) with values between 0 and 1. Each value indicates the strenght
#' of the relationship to be modelled between two blocks using sGCCA; a value
#' of 0 indicates no relationship, 1 is the maximum value. If \code{Y} is
#' provided instead of \code{indY}, the \code{design} matrix is changed to
#' include relationships to \code{Y}.
#' @param tau numeric vector of length the number of blocks in \code{X}. Each
#' regularization parameter will be applied on each block and takes the value
#' between 0 (no regularisation) and 1. If tau = "optimal" the shrinkage
#' paramaters are estimated for each block and each dimension using the Schafer
#' and Strimmer (2005) analytical formula.
#' @param ncomp the number of components to include in the model. Default to 1.
#' @param keepX A vector of same length as X.  Each entry keepX[i] is the
#' number of X[[i]]-variables kept in the model.
#' @param scale Logical. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE)
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var Logical, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Setting this argument to FALSE (when appropriate) will speed up the
#' computations. Default value is FALSE
#' @param all.outputs Logical. Computation can be faster when some specific
#' (and non-essential) outputs are not calculated. Default = \code{TRUE}.
#' @return \code{wrapper.rgcca} returns an object of class \code{"rgcca"}, a
#' list that contains the following components:
#' 
#' \item{data}{the input data set (as a list).} \item{design}{the input
#' design.} \item{variates}{the sgcca components.} \item{loadings}{the loadings
#' for each block data set (outer wieght vector).} \item{loadings.star}{the
#' laodings, standardised.} \item{tau}{the input tau parameter.}
#' \item{ncomp}{the number of components
#' included in the model for each block.} \item{crit}{the convergence
#' criterion.} \item{AVE}{Indicators of model quality based on the Average
#' Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner
#' model)..} \item{names}{list containing the names to be used for individuals
#' and variables.} More details can be found in the references.
#' Note that the argument 'scheme' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' @author Arthur Tenenhaus, Vincent Guillemot, Kim-Anh Lê Cao, 
#' Florian Rohart, Benoit Gautier
#' @seealso \code{\link{wrapper.rgcca}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{wrapper.sgcca}} and
#' \url{http://www.mixOmics.org} for more details.
#' @references
#' Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical
#' Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' 
#' Schafer J. and Strimmer K., (2005), A shrinkage approach to large-scale
#' covariance matrix estimation and implications for functional genomics.
#' Statist. Appl. Genet. Mol. Biol. 4:32.
#' @keywords multivariate
#' @export
#' @examples
#' data(nutrimouse)
#' # need to unmap the Y factor diet
#' Y = unmap(nutrimouse$diet)
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid, Y = Y)
#' # with this design, gene expression and lipids are connected to the diet factor
#' # design = matrix(c(0,0,1,
#' #                   0,0,1,
#' #                   1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#' 
#' # with this design, gene expression and lipids are connected to the diet factor
#' # and gene expression and lipids are also connected
#' design = matrix(c(0,1,1,
#' 1,0,1,
#' 1,1,0), ncol = 3, nrow = 3, byrow = TRUE)
#' #note: the tau parameter is the regularization parameter
#' wrap.result.rgcca = wrapper.rgcca(X = data, design = design, tau = c(1, 1, 0),
#' ncomp = 2)
#' #wrap.result.rgcca
#' 
wrapper.rgcca <-
    function(
        X,
        design = 1 - diag(length(X)),
        tau = rep(1, length(X)),
        ncomp = 1,
        keepX,
        scale = TRUE,
        tol = .Machine$double.eps,
        max.iter = 1000,
        near.zero.var = FALSE,
        all.outputs = TRUE)
    {

        # hardcode scheme and init
        scheme = "horst"
        init = "svd.single"
        
        check = Check.entry.rgcca(X = X, design = design, tau = tau, ncomp = ncomp, scheme = scheme, scale = scale,
                                  init = init, tol = tol, max.iter = max.iter, near.zero.var = near.zero.var,keepX = keepX)
        X = check$A
        ncomp = check$ncomp
        design = check$design
        init = check$init
        scheme = check$scheme
        nzv.A = check$nzv.A
        keepA = check$keepA
        
        keepA.save=keepA
        
        keepAA = vector("list", length = max(ncomp)) # one keepA per comp
        names(keepAA) = paste0("comp",1:max(ncomp))
        for(comp in 1:max(ncomp)) # keepA[[block]] [1:ncomp]
            keepAA[[comp]] = lapply(keepA, function(x) x[comp])
        
        keepA = lapply(keepAA, expand.grid)
        
        result.rgcca = internal_mint.block(A = X, design = design, tau = tau,
                                           ncomp = ncomp,
                                           scheme = scheme, scale = scale,
                                           init = init, tol = tol, keepA = keepA,
                                           max.iter = max.iter,
                                           study = factor(rep(1,nrow(X[[1]]))),#mint.rgcca not coded yet
                                           mode = "canonical",
                                           all.outputs = all.outputs
        )
        
        
        out = list(
            call = match.call(),
            X = result.rgcca$A,
            variates = result.rgcca$variates,
            loadings = result.rgcca$loadings,
            loadings.star = result.rgcca$loadings.star,
            keepX=keepA.save,
            design = result.rgcca$design,
            tau = result.rgcca$tau,
            scheme = result.rgcca$scheme,
            ncomp = result.rgcca$ncomp,
            crit = result.rgcca$crit,
            AVE = result.rgcca$AVE,
            names = result.rgcca$names,#names = list(indiv = rownames(X[[1]]), var = sapply(X, colnames)),
            init = result.rgcca$init,
            tol = result.rgcca$tol,
            iter = result.rgcca$iter,
            max.iter = result.rgcca$max.iter,
            nzv = result.rgcca$nzv,
            scale = result.rgcca$scale,
            design = result.rgcca$design,
            scheme = result.rgcca$scheme,
            prop_expl_var = result.rgcca$prop_expl_var
        )
        
        class(out) = c("sparse.rgcca","rgcca")
        return(invisible(out))
    }


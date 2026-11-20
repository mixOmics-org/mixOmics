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

#' mixOmics wrapper for Sparse Generalised Canonical Correlation Analysis
#' (sgcca)
#' 
#' Wrapper function to perform Sparse Generalised Canonical Correlation
#' Analysis (sGCCA), a generalised approach for the integration of multiple
#' datasets. For more details, see the \code{help(sgcca)} from the \pkg{RGCCA}
#' package.
#' 
#' This wrapper function performs sGCCA (see \pkg{RGCCA}) with \eqn{1, \ldots
#' ,}\code{ncomp} components on each block data set. A supervised or
#' unsupervised model can be run. For a supervised model, the
#' \code{\link{unmap}} function should be used as an input data set. More
#' details can be found on the package \pkg{RGCCA}.
#' 
#' Note that this function is the same as \code{\link{block.spls}} with
#' different default arguments.
#' 
#' More details about the PLS modes in \code{?pls}.
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
#' @param penalty numeric vector of length the number of blocks in \code{X}.
#' Each penalty parameter will be applied on each block and takes the value
#' between 0 (no variable selected) and 1 (all variables included).
#' @param ncomp the number of components to include in the model. Default to 1.
#' @param keepX A vector of same length as X.  Each entry keepX[i] is the
#' number of X[[i]]-variables kept in the model.
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
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
#' @return \code{wrapper.sgcca} returns an object of class \code{"sgcca"}, a
#' list that contains the following components:
#' 
#' \item{data}{the input data set (as a list).} \item{design}{the input
#' design.} \item{variates}{the sgcca components.} \item{loadings}{the loadings
#' for each block data set (outer wieght vector).} \item{loadings.star}{the
#' laodings, standardised.} \item{penalty}{the input penalty parameter.}
#' \item{scheme}{the input schme.} \item{ncomp}{the number of components
#' included in the model for each block.} \item{crit}{the convergence
#' criterion.} \item{AVE}{Indicators of model quality based on the Average
#' Variance Explained (AVE): AVE(for one block), AVE(outer model), AVE(inner
#' model)..} \item{names}{list containing the names to be used for individuals
#' and variables.} More details can be found in the references.
#' Note that the argument 'scheme' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' @author Arthur Tenenhaus, Vincent Guillemot, Kim-Anh Lê Cao, 
#' Florian Rohart, Benoit Gautier, Al J Abadi
#' @seealso \code{\link{wrapper.sgcca}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{wrapper.rgcca}} and
#' \url{http://www.mixOmics.org} for more details.
#' @references
#' Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical
#' Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' 
#' Tenenhaus A., Phillipe C., Guillemot, V., Lê Cao K-A., Grill J., Frouin, V.
#' Variable Selection For Generalized Canonical Correlation Analysis. 2013. (in
#' revision)
#' @keywords multivariate
#' @export
#' @examples
#' data(nutrimouse)
#' # need to unmap the Y factor diet if you pretend this is not a classification pb.
#' # see also the function block.splsda for discriminant analysis  where you dont
#' # need to unmap Y.
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
#' 
#' #note: the penalty parameters will need to be tuned
#' wrap.result.sgcca = wrapper.sgcca(X = data, design = design, penalty = c(.3,.5, 1),
#' ncomp = 2)
#' wrap.result.sgcca
#' #did the algo converge?
#' wrap.result.sgcca$crit  # yes
#' 
wrapper.sgcca <- 
    function(
        X,
        design = 1 - diag(length(X)),
        penalty = NULL,
        ncomp = 1,
        keepX,
        mode = "canonical",
        scale = TRUE,
        tol = .Machine$double.eps,
        max.iter = 1000,
        near.zero.var = FALSE,
        all.outputs = TRUE
    ){

        # hardcode scheme and init
        scheme = "horst"
        init = "svd.single"
        
        check=Check.entry.sgcca(X = X, design = design ,ncomp = ncomp , scheme = scheme , scale = scale,
                                init = init , tol = tol, mode = mode, max.iter = max.iter,near.zero.var = near.zero.var,keepX = keepX)
        
        
        A = check$A
        design = check$design
        ncomp = check$ncomp
        init = check$init
        scheme = check$scheme
        near.zero.var = check$near.zero.var
        keepA = check$keepA
        nzv.A = check$nzv.A
        
        keepAA = vector("list", length = max(ncomp)) # one keepA per comp
        names(keepAA) = paste0("comp",1:max(ncomp))
        for(comp in 1:max(ncomp)) # keepA[[block]] [1:ncomp]
            keepAA[[comp]] = lapply(keepA, function(x) x[comp])
        
        keepA = lapply(keepAA, expand.grid)
        
        result.sgcca = internal_mint.block(A = A, design = design, tau = NULL,
                                           ncomp = ncomp,
                                           scheme = scheme, scale = scale,
                                           init = init, tol = tol,
                                           keepA = keepA,
                                           max.iter = max.iter,
                                           study = factor(rep(1,nrow(A[[1]]))),
                                           mode = mode,penalty = penalty,
                                           all.outputs = all.outputs
        )
        
        
        out = list(
            call = match.call(),
            X = result.sgcca$A,
            variates = result.sgcca$variates,
            loadings = result.sgcca$loadings,
            loadings.star = result.sgcca$loadings.star,
            design = result.sgcca$design,
            penalty = penalty,
            scheme = result.sgcca$scheme,
            ncomp = result.sgcca$ncomp,
            crit = result.sgcca$crit,
            AVE = result.sgcca$AVE,
            names = result.sgcca$names,
            init = result.sgcca$init,
            tol = result.sgcca$tol,
            iter = result.sgcca$iter,
            max.iter = result.sgcca$max.iter,
            nzv = result.sgcca$nzv,
            scale = result.sgcca$scale,
            design = result.sgcca$design,
            scheme = result.sgcca$scheme,
            prop_expl_var = result.sgcca$prop_expl_var
        )
        
        class(out) = 'sgcca'
        return(invisible(out))
    }


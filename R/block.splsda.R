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

#' N-integration and feature selection with Projection to Latent Structures
#' models (PLS) with sparse Discriminant Analysis
#' 
#' Integration of multiple data sets measured on the same samples or
#' observations to classify a discrete outcome to classify a discrete outcome
#' and select features from each data set, ie. N-integration with sparse
#' Discriminant Analysis. The method is partly based on Generalised Canonical
#' Correlation Analysis.
#' 
#' 
#' \code{block.splsda} function fits a horizontal integration PLS-DA model with
#' a specified number of components per block). A factor indicating the
#' discrete outcome needs to be provided, either by \code{Y} or by its position
#' \code{indY} in the list of blocks \code{X}.
#' 
#' \code{X} can contain missing values. Missing values are handled by being
#' disregarded during the cross product computations in the algorithm
#' \code{block.pls} without having to delete rows with missing data.
#' Alternatively, missing data can be imputed prior using the
#' \code{\link{impute.nipals}} function.
#' 
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and \code{?pls} for more details).
#' 
#' Note that our method is partly based on sparse Generalised Canonical
#' Correlation Analysis and differs from the MB-PLS approaches proposed by
#' Kowalski et al., 1989, J Chemom 3(1), Westerhuis et al., 1998, J Chemom,
#' 12(5) and sparse variants Li et al., 2012, Bioinformatics 28(19); Karaman et
#' al (2014), Metabolomics, 11(2); Kawaguchi et al., 2017, Biostatistics.
#' 
#' Variable selection is performed on each component for each block of \code{X}
#' if specified, via input parameter \code{keepX}.
#' 
#' @inheritParams block.plsda
#' @inheritParams block.spls
#' @template arg/verbose.call
#' @return \code{block.splsda} returns an object of class \code{"block.splsda",
#' "block.spls"}, a list that contains the following components:
#' 
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{keepX}{Number of
#' variables used to build each component of each block} \item{variates}{list
#' containing the variates of each block of X.} \item{loadings}{list containing
#' the estimated loadings for the variates.} \item{names}{list containing the
#' names to be used for individuals and variables.} \item{nzv}{list containing
#' the zero- or near-zero predictors information.} \item{iter}{Number of
#' iterations of the algorithm for each component} \item{weights}{Correlation
#' between the variate of each block and the variate of the outcome. Used to
#' weight predictions.} \item{prop_expl_var}{Percentage of explained
#' variance for each component and each block}
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' Note that the argument 'scheme' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{plotIndiv}}, \code{\link{plotArrow}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.plsda}},
#' \code{\link{block.spls}} and https://mixomics.org/mixDIABLO/ for more
#' details and examples.
#' @references On multiple integration with sPLS-DA and 4 data blocks:
#' 
#' Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê
#' Cao K.A. (2016). DIABLO: multi omics integration for biomarker discovery.
#' BioRxiv available here:
#' \url{http://biorxiv.org/content/early/2016/08/03/067611}
#' 
#' On data integration:
#' 
#' Tenenhaus A., Philippe C., Guillemot V, Lê Cao K.A., Grill J, Frouin V.
#' Variable selection for generalized canonical correlation analysis.
#' \emph{Biostatistics}. kxu001
#' 
#' Gunther O., Shin H., Ng R. T. , McMaster W. R., McManus B. M. , Keown P. A.
#' , Tebbutt S.J. , Lê Cao K-A. , (2014) Novel multivariate methods for
#' integration of genomics and proteomics data: Applications in a kidney
#' transplant rejection study, OMICS: A journal of integrative biology, 18(11),
#' 682-95.
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate
#' @example ./examples/block.splsda-examples.R
#' @export
block.splsda <- function(X,
                         Y,
                         indY,
                         ncomp = 2,
                         keepX,
                         design,
                         scale = TRUE,
                         tol = 1e-06,
                         max.iter = 100,
                         near.zero.var = FALSE,
                         all.outputs = TRUE,
                         verbose.call = FALSE)
{
    # check inpuy 'Y' and transformation in a dummy matrix
    if (!missing(Y))
    {
        if (is.null(dim(Y)))
        {
            Y = factor(Y)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
            stop("'Y' should be a factor with more than one level")
        
        Y.input = Y
        Y = unmap(Y)
        colnames(Y) = levels(Y.input)
        rownames(Y) = rownames(X[[1]])
        
    } else if (!missing(indY)) {
        temp = X[[indY]]
        #not called Y to not be an input of the wrapper.sparse.mint.block
        if (is.null(dim(temp)))
        {
            temp = factor(temp)
        } else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(temp) == 1)
            stop("'X[[indY]]' should be a factor with more than one level")
        
        Y.input = temp
        X[[indY]] = unmap(temp)
        colnames(X[[indY]]) = levels(Y.input)
        rownames(X[[indY]]) = rownames(X[[ifelse(indY == 1, 2, 1)]])
        
    } else if (missing(indY)) {
        stop("Either 'Y' or 'indY' is needed")
        
    }
    
    # call to 'internal_wrapper.mint.block'
    result = internal_wrapper.mint.block(
        X = X,
        Y = Y,
        indY = indY,
        ncomp = ncomp,
        keepX = keepX,
        DA = TRUE,
        design = design,
        scheme = 'horst',
        mode = 'regression',
        scale = scale,
        init = 'svd.single',
        tol = tol,
        max.iter = max.iter,
        near.zero.var = near.zero.var,
        all.outputs = all.outputs
    )
    
    # calculate weights for each dataset
    weights = get.weights(result$variates, indY = result$indY)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY],
        Y = Y.input,
        ind.mat = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        keepX = result$keepX[-result$indY],
        keepY = result$keepX[result$indY][[1]],
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
        indY = result$indY,
        weights = weights,
        prop_expl_var = result$prop_expl_var
    )#[-result$indY])
    
    if (verbose.call) {
        c <- out$call
        out$call <- mget(names(formals()))
        out$call <- append(c, out$call)
        names(out$call)[1] <- "simple.call"
    }
    
    # give a class
    class(out) = c("block.splsda", "block.spls", "sgccda", "sgcca", "DA")
    return(invisible(out))
    
}

#' @rdname block.splsda
#' @export
wrapper.sgccda <- block.splsda

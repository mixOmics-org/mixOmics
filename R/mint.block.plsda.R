#############################################################################################################
# Author :
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 22-04-2015
# last modified: 05-10-2017
#
# Copyright (C) 2015
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# ========================================================================================================
# mint.block.plsda: perform a horizontal and vertical PLS-DA on a combination of datasets, input as a list in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint.block, which then call 'internal_mint.block'
# ========================================================================================================

#' NP-integration with Discriminant Analysis
#' 
#' Function to integrate data sets measured on the same samples (N-integration)
#' and to combine multiple independent studies measured on the same variables
#' or predictors (P-integration) using variants of multi-group and generalised
#' PLS-DA for supervised classification.
#' 
#' The function fits multi-group generalised PLS models with a specified number
#' of \code{ncomp} components. A factor indicating the discrete outcome needs
#' to be provided, either by \code{Y} or by its position \code{indY} in the
#' list of blocks \code{X}.
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
#' References and more details in \code{?pls}).
#' 
#' @inheritParams mint.block.pls
#' @param Y A factor or a class vector indicating the discrete outcome of each
#' sample.
#' @return \code{mint.block.plsda} returns an object of class
#' \code{"mint.plsda", "block.plsda"}, a list that contains the following
#' components:
#' 
#' \item{X}{the centered and standardized original predictor matrices.}
#' \item{Y}{the original factor indicating the discrete outcome of each
#' sample.}
#' \item{ind.mat}{the centered and standardized indicator matrix corresponding
#' to \code{Y}.}
#' \item{indY}{the position of the outcome Y in \code{variates} and
#' \code{loadings}.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{mat.c}{matrix of
#' coefficients from the regression of X / residual matrices X on the
#' X-variates, to be used internally by \code{predict}.} \item{variates}{list
#' containing the \eqn{X} and \eqn{Y} variates.} \item{loadings}{list
#' containing the estimated loadings for the variates.} \item{names}{list
#' containing the names to be used for individuals and variables.}
#' \item{nzv}{list containing the zero- or near-zero predictors information.}
#' \item{tol}{the tolerance used in the iterative algorithm, used for
#' subsequent S3 methods} \item{max.iter}{the maximum number of iterations,
#' used for subsequent S3 methods} \item{iter}{Number of iterations of the
#' algorithm for each component}
#' \item{weights}{Correlation between the variate of each block and the
#' variate of the outcome. Used to weight predictions.}
#' Note that the argument 'scheme' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{predict}}, \code{\link{perf}},
#' \code{\link{mint.block.spls}}, \code{\link{mint.block.plsda}},
#' \code{\link{mint.block.splsda}} and http://www.mixOmics.org/mixMINT for more
#' details.
#' @references On multi-group PLS:
#' 
#' Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017). MINT: A
#' multivariate integrative approach to identify a reproducible biomarker
#' signature across multiple experiments and platforms. BMC Bioinformatics
#' 18:128.
#' 
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms
#' for multi-group PLS. J. Chemometrics, 28(3), 192-201.
#' 
#' On multiple integration with PLSDA:
#' 
#' Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt S. and Lê
#' Cao K.A. (2016). DIABLO: multi omics integration for biomarker discovery.
#' BioRxiv available here:
#' \url{http://biorxiv.org/content/early/2016/08/03/067611} Tenenhaus A.,
#' Philippe C., Guillemot V, Lê Cao K.A., Grill J, Frouin V. Variable selection
#' for generalized canonical correlation analysis. \emph{Biostatistics}. kxu001
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
#' @examples
#' 
#' data(breast.TCGA)
#' 
#' # for the purpose of this example, we consider the training set as study1 and
#' # the test set as another independent study2.
#' study = c(rep("study1",150), rep("study2",70))
#' 
#' mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
#' mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
#' data = list(mrna = mrna, mirna = mirna)
#' 
#' Y = c(breast.TCGA$data.train$subtype, breast.TCGA$data.test$subtype)
#' 
#' res = mint.block.plsda(data,Y,study=study, ncomp=2)
#' 
#' res
#' @export
mint.block.plsda <- function(X,
                             Y,
                             indY,
                             study,
                             ncomp = 2,
                             design,
                             scale = TRUE,
                             tol = 1e-06,
                             max.iter = 100,
                             near.zero.var = FALSE,
                             all.outputs = TRUE)
{
    # checking that the outcome, either in Y or X[indY] is a factor with more than 1 level
    if(!missing(Y))
    {
        if (is.null(dim(Y)))
        {
            Y = as.factor(Y)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        
        if (nlevels(Y) == 1)
            stop("'Y' should be a factor with more than one level")
        
        Y.input=Y
        Y=unmap(Y)
        colnames(Y) = paste0("Y", 1:ncol(Y))
        rownames(Y) = rownames(X[[1]])
        
    }else if(!missing(indY)) {
        temp=X[[indY]] #not called Y to not be an input of the wrapper.sparse.mint.block
        if (is.null(dim(temp)))
        {
            temp = as.factor(temp)
        }  else {
            stop("'Y' should be a factor or a class vector.")
        }
        if (nlevels(temp) == 1)
            stop("'X[[indY]]' should be a factor with more than one level")
        
        Y.input=temp
        X[[indY]]=unmap(temp)
        rownames(X[[indY]]) = rownames(X[[ifelse(indY==1,2,1)]])
    } else if(missing(indY)) {
        stop("Either 'Y' or 'indY' is needed")
    }
    
    # call to 'internal_wrapper.mint.block'
    result = internal_wrapper.mint.block(
        X = X,
        Y = Y,
        indY = indY,
        study = study,
        ncomp = ncomp,
        design = design,
        scheme = 'horst',
        mode = 'regression',
        scale = scale,
        init = 'svd.single',
        tol = tol,
        max.iter = max.iter,
        near.zero.var = near.zero.var,
        all.outputs = all.outputs,
        DA = TRUE
    )
    
    # calculate weights for each dataset
    weights = get.weights(result$variates, indY = result$indY)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A[-result$indY],
        Y = Y.input,
        indY = result$indY,
        weights = weights,
        ind.mat = result$A[result$indY][[1]],
        ncomp = result$ncomp,
        mode = result$mode,
        study = result$study,
        variates = result$variates,
        loadings = result$loadings,
        variates.partial = result$variates.partial,
        loadings.partial = result$loadings.partial,
        names = result$names,
        init = result$init,
        tol = result$tol,
        iter = result$iter,
        max.iter = result$max.iter,
        nzv = result$nzv,
        scale = result$scale)
    
    class(out) = c("mint.block.plsda","mint.block.pls","block.plsda","block.pls","sgccda","sgcca","DA")
    return(invisible(out))
    
}




# ========================================================================================================
# mint.block.spls: perform a horizontal and vertical sPLS on a combination of datasets, input as a list in X
# this function is a particular setting of internal_mint.block,
# the formatting of the input is checked in internal_wrapper.mint.block, which then call 'internal_mint.block'
# ========================================================================================================

#' NP-integration for integration with variable selection
#' 
#' Function to integrate data sets measured on the same samples (N-integration)
#' and to combine multiple independent studies (P-integration) using variants
#' of sparse multi-group and generalised PLS with variable selection
#' (unsupervised analysis).
#' 
#' The function fits sparse multi-group generalised PLS models with a specified
#' number of \code{ncomp} components. An outcome needs to be provided, either
#' by \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
#' 
#' Multi (continuous)response are supported. \code{X} and \code{Y} can contain
#' missing values. Missing values are handled by being disregarded during the
#' cross product computations in the algorithm \code{block.pls} without having
#' to delete rows with missing data. Alternatively, missing data can be imputed
#' prior using the \code{nipals} function.
#' 
#' The type of algorithm to use is specified with the \code{mode} argument.
#' Four PLS algorithms are available: PLS regression \code{("regression")}, PLS
#' canonical analysis \code{("canonical")}, redundancy analysis
#' \code{("invariant")} and the classical PLS algorithm \code{("classic")} (see
#' References and more details in \code{?pls}).
#' 
#' @inheritParams mint.block.pls
#' @inheritParams block.spls
#' @return \code{mint.block.spls} returns an object of class \code{"mint.spls",
#' "block.spls"}, a list that contains the following components:
#' 
#' \item{X}{the centered and standardized original predictor and response matrices.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{design}{the design matrix indicating the relationships modelled
#' between the blocks.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{keepX}{Number of
#' variables used to build each component of each block} \item{keepY}{Number of
#' variables used to build each component of Y} \item{mat.c}{matrix of
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
#' \code{\link{mint.block.pls}}, \code{\link{mint.block.plsda}},
#' \code{\link{mint.block.splsda}} and http://www.mixOmics.org/mixMINT for more
#' details.
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#' 
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2014). Algorithms
#' for multi-group PLS. J. Chemometrics, 28(3), 192-201.
#' @keywords regression multivariate
#' @examples
#' data(breast.TCGA)
#' 
#' # for the purpose of this example, we create data that fit in the context of
#' # this function.
#' # We consider the training set as study1 and the test set as another
#' # independent study2.
#' 
#' study = c(rep("study1",150), rep("study2",70))
#' 
#' # to put the data in the MINT format, we rbind the two studies
#' mrna = rbind(breast.TCGA$data.train$mrna, breast.TCGA$data.test$mrna)
#' mirna = rbind(breast.TCGA$data.train$mirna, breast.TCGA$data.test$mirna)
#' 
#' # For the purpose of this example, we create a continuous response by
#' # taking the first mrna variable, and removing it from the data
#' Y = mrna[,1,drop=FALSE]
#' mrna = mrna[,-1]
#'
#' data = list(mrna = mrna, mirna = mirna)
#'
#' # we can now apply the function
#' res = mint.block.spls(data, Y, study=study, ncomp=2,
#' keepX = list(mrna=c(10,10), mirna=c(20,20)))
#' 
#' res
#' @export
mint.block.spls <- function(X,
                            Y,
                            indY,
                            study,
                            ncomp = 2,
                            keepX,
                            keepY,
                            design,
                            mode,
                            scale = TRUE,
                            tol = 1e-06,
                            max.iter = 100,
                            near.zero.var = FALSE,
                            all.outputs = TRUE)
{
    # call to 'internal_wrapper.mint.block'
    result = internal_wrapper.mint.block(
        X = X,
        Y = Y,
        indY = indY,
        study = study,
        ncomp = ncomp,
        keepX = keepX,
        keepY = keepY,
        design = design,
        scheme = 'horst',
        mode = mode,
        scale = scale,
        init = 'svd.single',
        tol = tol,
        max.iter = max.iter,
        near.zero.var = near.zero.var,
        all.outputs = all.outputs,
        DA = FALSE
    )
    
    # calculate weights for each dataset
    weights = get.weights(result$variates, indY = result$indY)
    
    # choose the desired output from 'result'
    out = list(
        call = match.call(),
        X = result$A,
        indY = result$indY,
        ncomp = result$ncomp,
        mode = result$mode,
        study = result$study,
        keepX = result$keepX[-result$indY],
        keepY = result$keepX[[result$indY]],
        design = result$design,
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
        scale = result$scale,
        weights = weights
    )
    
    class(out) = c("mint.block.spls","block.spls","sgcca")
    return(invisible(out))
    
}

# =============================================================================
# block.spls: perform a horizontal PLS-DA on a combination of datasets,
#   input as a list in X
#   this function is a particular setting of internal_mint.block,
#   the formatting of the input is checked in internal_wrapper.mint.block
# =============================================================================

#' N-integration and feature selection with sparse Projection to Latent
#' Structures models (sPLS)
#' 
#' Integration of multiple data sets measured on the same samples or
#' observations, with variable selection in each data set, ie. N-integration.
#' The method is partly based on Generalised Canonical Correlation Analysis.
#' 
#' \code{block.spls} function fits a horizontal sPLS model with a specified
#' number of components per block). An outcome needs to be provided, either by
#' \code{Y} or by its position \code{indY} in the list of blocks \code{X}.
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
#' References and \code{?pls} for more details).
#' 
#' Note that our method is partly based on sparse Generalised Canonical
#' Correlation Analysis and differs from the MB-PLS approaches proposed by
#' Kowalski et al., 1989, J Chemom 3(1), Westerhuis et al., 1998, J Chemom,
#' 12(5) and sparse variants Li et al., 2012, Bioinformatics 28(19); Karaman et
#' al (2014), Metabolomics, 11(2); Kawaguchi et al., 2017, Biostatistics.
#' 
#' Variable selection is performed on each component for each block of
#' \code{X}, and for \code{Y} if specified, via input parameter \code{keepX}
#' and \code{keepY}.
#' 
#' Note that if \code{Y} is missing and \code{indY} is provided, then variable
#' selection on \code{Y} is performed by specifying the input parameter
#' directly in \code{keepX} (no \code{keepY} is needed).
#' 
#' @inheritParams block.pls
#' @param Y Matrix response for a multivariate regression framework. Data
#' should be continuous variables (see \code{?block.splsda} for supervised
#' classification and factor response).
#' @param keepX A named list of same length as X. Each entry is the number of
#' variables to select in each of the blocks of X for each component. By
#' default all variables are kept in the model.
#' @param keepY Only if Y is provided (and not \code{indY}). Each entry is the number of variables to
#' select in each of the blocks of Y for each component.
#' @template arg/verbose.call
#' @return \code{block.spls} returns an object of class \code{"block.spls"}, a
#' list that contains the following components:
#' 
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{indY}{the position of the outcome Y in the output list X.}
#' \item{ncomp}{the number of components included in the model for each block.}
#' \item{mode}{the algorithm used to fit the model.} \item{keepX}{Number of
#' variables used to build each component of each block} \item{keepY}{Number of
#' variables used to build each component of Y} \item{variates}{list containing
#' the variates of each block of X.} \item{loadings}{list containing the
#' estimated loadings for the variates.} \item{names}{list containing the names
#' to be used for individuals and variables.} \item{nzv}{list containing the
#' zero- or near-zero predictors information.} \item{iter}{Number of iterations
#' of the algorithm for each component} \item{prop_expl_var}{Percentage of
#' explained variance for each component and each block after setting possible
#' missing values in the centered data to zero}
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' Note that the argument 'scheme' has now been hardcoded to 'horst' and 'init' to 'svd.single'. 
#' @author Florian Rohart, Benoit Gautier, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{plotIndiv}}, \code{\link{plotArrow}},
#' \code{\link{plotLoadings}}, \code{\link{plotVar}}, \code{\link{predict}},
#' \code{\link{perf}}, \code{\link{selectVar}}, \code{\link{block.pls}},
#' \code{\link{block.splsda}} and http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' 
#' Tenenhaus A. and Tenenhaus M., (2011), Regularized Generalized Canonical
#' Correlation Analysis, Psychometrika, Vol. 76, Nr 2, pp 257-284.
#' 
#' Tenenhaus A., Philippe C., Guillemot V, Lê Cao K.A., Grill J, Frouin V.
#' Variable selection for generalized canonical correlation analysis.
#' \emph{Biostatistics}. kxu001
#' @keywords regression multivariate
#' @example ./examples/block.spls-examples.R
#' @export
block.spls = function(X,
                      Y,
                      indY,
                      ncomp = 2,
                      keepX,
                      keepY,
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
    result = internal_wrapper.mint.block(
        X = X,
        Y = Y,
        indY = indY,
        ncomp = ncomp,
        keepX = keepX,
        keepY = keepY,
        design = design,
        scheme = "horst",
        mode = mode,
        scale = scale,
        init = "svd.single",
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
        weights = weights,
        prop_expl_var = result$prop_expl_var
    )
    
    if (verbose.call) {
        c <- out$call
        out$call <- mget(names(formals()))
        out$call <- append(c, out$call)
        names(out$call)[1] <- "simple.call"
    }
    
    # give a class
    class(out) = c("block.spls", "sgcca")
    return(invisible(out))
}




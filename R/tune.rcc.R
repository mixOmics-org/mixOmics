library(BiocParallel)

#' Estimate the parameters of regularization for Regularized CCA
#' 
#' Computes leave-one-out or M-fold cross-validation scores on a
#' two-dimensional grid to determine optimal values for the parameters of
#' regularization in \code{rcc}.
#' 
#' If \code{validation="Mfold"}, the samples are randomly split into the
#' number of folds specified by \code{folds}.
#'
#' If \code{validation="loo"}, leave-one-out cross-validation is performed by
#' leaving out each sample in turn. In this case \code{folds} is ignored.
#'
#' The \code{scale} argument precedes \code{validation}. Calls that previously
#' supplied \code{validation} or later arguments by position should now name
#' those arguments explicitly.
#'
#' Both validation methods centre held-out samples using the training-fold
#' means. With \code{scale = TRUE} (the default), standard deviations are also
#' estimated on the training fold only and applied to the held-out samples. Supply
#' unscaled data to avoid using held-out observations to estimate scaling
#' parameters, including when using \code{tune.rcc} inside nested
#' cross-validation. Use the same \code{scale} setting when fitting the final
#' \code{\link{rcc}} model.
#'
#' Ridge penalties depend on the variance scale of the input variables,
#' unlike classical CCA. The default grids assume variables of roughly unit
#' variance. By default, variables are standardised within each training fold.
#' If \code{scale = FALSE}, choose grids appropriate for the unscaled data.
#'
#' The estimation of the missing values can be performed by the reconstitution
#' of the data matrix using the \code{nipals} function. Otherwise, missing
#' values are handled by pairwise-complete covariance estimates in
#' \code{rcc}. Missing held-out values are replaced by zero after applying
#' the training-fold transformation, equivalent to training-mean imputation.
#' 
#' @param X numeric matrix or data frame \eqn{(n \times p)}, the observations
#' on the \eqn{X} variables. \code{NA}s are allowed.
#' @param Y numeric matrix or data frame \eqn{(n \times q)}, the observations
#' on the \eqn{Y} variables. \code{NA}s are allowed.
#' @param grid1,grid2 vector numeric defining the values of \code{lambda1} and
#' \code{lambda2} at which cross-validation score should be computed. Defaults
#' to \code{grid1=grid2=seq(0.001, 1, length=5)}.
#' @param scale logical. If \code{TRUE}, standardise both data matrices
#' within each training fold and apply the training standard deviations to
#' held-out samples. Defaults to \code{TRUE}; samples are always centred
#' using the training-fold means.
#' @param validation character string. What kind of (internal) cross-validation
#' method to use, (partially) matching one of \code{"loo"} (leave-one-out) or
#' \code{"Mfold"} (M-folds). See Details.
#' @param folds positive integer. Number of folds to use if
#' \code{validation="Mfold"}. Defaults to \code{folds=10}.
#' @param BPPARAM a BiocParallel parameter object; see \code{BiocParallel::bpparam} 
#' for details. Default is \code{SerialParam()} for serial processing.
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'.
#' @return The returned value is a list with components: \item{opt.lambda1,}{}
#' \item{opt.lambda2}{value of the parameters of regularization on which the
#' cross-validation method reached its optimal.} \item{opt.score}{the optimal
#' cross-validation score reached on the grid.} \item{grid1, grid2}{original
#' vectors \code{grid1} and \code{grid2}.} \item{mat}{matrix containing the
#' cross-validation score computed on the grid.}
#' @author Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{image.tune.rcc}} and http://www.mixOmics.org for more
#' details.
#' @keywords multivariate dplot
#' @export
#' @example ./examples/tune.rcc-examples.R

tune.rcc <-
    function(X, 
             Y, 
             grid1 = seq(0.001, 1, length = 5), 
             grid2 = seq(0.001, 1, length = 5), 
             scale = TRUE,
             validation = c("loo", "Mfold"), 
             folds = 10,
             BPPARAM = SerialParam(),
             seed = NULL)
    {
      if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
        stop("'scale' must be either TRUE or FALSE.", call. = FALSE)

      BPPARAM$RNGseed <- seed
      set.seed(seed)
      
        # validation des arguments #
        #--------------------------#
        if (length(dim(X)) != 2 || length(dim(Y)) != 2) 
            stop("'X' and/or 'Y' must be a numeric matrix.")
        
        X = as.matrix(X)
        Y = as.matrix(Y)
        
        if (!is.numeric(X) || !is.numeric(Y)) 
            stop("'X' and/or 'Y' must be a numeric matrix.")
        
        if (nrow(X) != nrow(Y)) 
            stop("unequal number of rows in 'X' and 'Y'.")
        
        validation = match.arg(validation)
        grid = expand.grid(grid1, grid2)
        
        if (validation == "loo")
        {
            M = nrow(X)
            folds = split(1:M, 1:M)
            cv.score = bplapply(1:nrow(grid), function(i) {
                lambda = as.numeric(grid[i, ])  # Ensure lambda is numeric
                Mfold(X, Y, lambda[1], lambda[2], folds, scale = scale)
            }, BPPARAM = BPPARAM)
            
        } else {
            nr = nrow(X)
            M = length(folds)
            if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > nr)
            {
                stop("Invalid number of folds.")
            } else {
                M = round(folds)
                folds = split(sample(1:nr), rep(1:M, length = nr))
            }
            cv.score = bplapply(1:nrow(grid), function(i) {
                lambda = as.numeric(grid[i, ])  # Ensure lambda is numeric
                Mfold(X, Y, lambda[1], lambda[2], folds, scale = scale)
            }, BPPARAM = BPPARAM)
        }
        
        cv.score = unlist(cv.score)
        cv.score.grid = cbind(grid, cv.score)
        mat = matrix(cv.score, nrow = length(grid1), ncol = length(grid2))
        
        opt = cv.score.grid[cv.score.grid[, 3] == max(cv.score.grid[, 3]), ]
        
        out = list(opt.lambda1 = opt[[1]], opt.lambda2 = opt[[2]], 
                   opt.score = opt[[3]], grid1 = grid1, grid2 = grid2, mat = mat)
        
        out$call = match.call()
        
        class(out) = "tune.rcc"
        return(invisible(out))
    }

Mfold = function(X, Y, lambda1, lambda2, folds, scale = TRUE)
{
    xscore = NULL
    yscore = NULL
    M = length(folds)
    
    for (m in 1:M)
    {
        omit = folds[[m]]
        result = rcc(X[-omit, , drop = FALSE], Y[-omit, , drop = FALSE],
                     ncomp = 1, lambda1 = lambda1, lambda2 = lambda2,
                     method = "ridge", scale = scale)
        X.test = base::scale(X[omit, , drop = FALSE],
                             center = result$center$X, scale = result$scale$X)
        Y.test = base::scale(Y[omit, , drop = FALSE],
                             center = result$center$Y, scale = result$scale$Y)
        X.test[is.na(X.test)] = 0
        Y.test[is.na(Y.test)] = 0
        xscore = c(xscore, X.test %*% result$loadings$X[, 1])
        yscore = c(yscore, Y.test %*% result$loadings$Y[, 1])
    }
    
    cv.score = cor(xscore, yscore, use = "pairwise")
    return(invisible(cv.score))
}

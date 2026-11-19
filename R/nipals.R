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

#' Non-linear Iterative Partial Least Squares (NIPALS) algorithm
#' 
#' This function performs NIPALS algorithm, i.e. the singular-value
#' decomposition (SVD) of a data table that can contain missing values.
#' 
#' The NIPALS algorithm (Non-linear Iterative Partial Least Squares) has been
#' developed by H. Wold at first for PCA and later-on for PLS. It is the most
#' commonly used method for calculating the principal components of a data set.
#' It gives more numerically accurate results when compared with the SVD of the
#' covariance matrix, but is slower to calculate.
#' 
#' This algorithm allows to realize SVD with missing data, without having to
#' delete the rows with missing data or to estimate the missing data.
#' 
#' @inheritParams pca
#' @return An object of class 'mixo_nipals' containing slots: 
#' \item{eig}{Vector containing the pseudo-singular values of \code{X}, of length
#' \code{ncomp}.}
#' \item{t}{Matrix whose columns contain the left singular vectors of \code{X}.
#' Note that for a complete data matrix X, the return values \code{eig},
#' \code{t} and \code{p} such that \code{X = t * diag(eig) * t(p)}.}
#' @author Sébastien Déjean, Ignacio González, Kim-Anh Le Cao, Al J Abadi
#' @seealso \code{\link{impute.nipals}}, \code{\link{svd}},
#'   \code{\link{princomp}}, \code{\link{prcomp}}, \code{\link{eigen}} and
#'   http://www.mixOmics.org for more details.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' 
#' Wold H. (1975). Path models with latent variables: The NIPALS approach. In:
#' Blalock H. M. et al. (editors). \emph{Quantitative Sociology: International
#' perspectives on mathematical and statistical model building}. Academic
#' Press, N.Y., 307-357.
#' @keywords algebra multivariate
#' @export
nipals <- function (X,
                    ncomp = 2,
                    max.iter = 500,
                    tol = 1e-06)
{
    #-- X matrix
    X <- .check_numeric_matrix(X)
    ncomp <- .check_ncomp(ncomp = ncomp, X = X, default = 2)
    if (any(colSums(!is.na(X)) == 0) | any(rowSums(!is.na(X)) == 0 ))
        stop("some rows or columns are entirely missing. ", 
             "Remove those before running pca.", call. = FALSE)
    nc = ncol(X)
    nr = nrow(X)

    #-- pca approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#

    #-- initialisation --#
    comp_names <- paste0('PC', seq_len(ncomp))
    p <- matrix(nrow = nc, ncol = ncomp, dimnames = list(colnames(X), comp_names))
    t.mat <- matrix(nrow = nr, ncol = ncomp, dimnames = list(rownames(X), comp_names))
    eig <- vector("numeric", length = ncomp)
    names(eig) <- comp_names
    nc.ones <- rep(1, nc)
    nr.ones <- rep(1, nr)
    is.na.X <- is.na(X)
    na.X <- any(is.na.X)
    X.iter <- X
    #-- loop on h --#
    for (h in 1:ncomp)
    {
        if (na.X)
        {
            ## initialise with minimum missing value column
            init.col <- which.min(colSums(is.na(X.iter)))
        } else {
            ## initialise with maximum variance column
            init.col <- which.max(apply(X.iter, 2, var, na.rm = TRUE))
        }
        
        th <- X.iter[, init.col]
        th[is.na(th)] <- 0
        ph.old <- 1/rep(sqrt(nc), nc)
        ph.new <- vector("numeric", length = nc)
        iter <- 1
        diff <- 1
        
        if (na.X)
        {
            # TODO is copying X.iter necessary here? Maybe we can just store the NA
            # indices and restore later (is that even necessary to restore them
            # now that X.iter is separate from X?). Simply set all NA to 0 (if any)
            # and see if results are different
            X.aux <- X.iter
            X.aux[is.na.X] <- 0
        }
        
        while (diff > tol & iter <= max.iter)
        {
            if (na.X)
            {
                ph.new <- crossprod(X.aux, th)
                Th <- drop(th) %o% nc.ones
                Th[is.na.X] <- 0
                th.cross <- crossprod(Th)
                ph.new <- ph.new / diag(th.cross)
            } else {
                ph.new <- crossprod(X.iter, th) / drop(crossprod(th))
            }
            
            ph.new <- ph.new / drop(sqrt(crossprod(ph.new)))
            
            if (na.X)
            {
                th <- X.aux %*% ph.new
                P <- drop(ph.new) %o% nr.ones
                P[t(is.na.X)] <- 0
                ph.cross <- crossprod(P)
                th <- th / diag(ph.cross)
            } else {
                th <- X.iter %*% ph.new / drop(crossprod(ph.new))
            }
            
            diff <- drop(sum((ph.new - ph.old)^2, na.rm = TRUE))
            ph.old <- ph.new
            iter <- iter + 1
        }
        
        if (iter > max.iter)
            message(paste("Maximum number of iterations reached for comp: ", h))
        
        X.iter <- X.iter - th %*% t(ph.new)
        p[, h] <- ph.new
        if (h > 1)
        {
            max.comp.cor <-  max(abs(cor(th, t.mat[,seq_len(h-1),drop=FALSE])))
            if (max.comp.cor >= 0.05)
            {
                message(sprintf("Component %s is not orthogonal to previous ones. (cor = %s).\nConsider filtering features with high rate of missing values or imputing the missing values.\n",
                                h, 
                                round(max.comp.cor, 3)))
            }
           
        }
        t.mat[, h] <- th
        eig[h] <- sum(th * th, na.rm = TRUE)
    }
    
    eig <- sqrt(eig)
    t.mat <- scale(t.mat, center = FALSE, scale = eig)
    attr(t.mat, "scaled:scale") <- NULL
    
    res <- list(eig = eig, p = p, t = t.mat)

    class(res) <- c('mixo_nipals')
    res
    return(invisible(res))
}

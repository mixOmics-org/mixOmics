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

#' Impute missing values using NIPALS algorithm
#'
#' This function uses \code{\link{nipals}} function to decompose
#' \code{X} into a set of components (\code{t}), (pseudo-) singular-values
#' (\code{eig}), and feature loadings (\code{p}). The original matrix is then
#' approximated/reconstituted using the following equation:
#' \deqn{\hat{X} = t * diag(eig) * t(p)} 
#' The missing values from \code{X} are then approximated from this matrix. It
#' is best to ensure enough number of components are used in order to best
#' impute the missing values.
#' 
#' @param X A numeric matrix containing missing values
#' @param ncomp Positive integer, the number of components to derive from
#'   \code{X} using the \code{\link{nipals}} function and reconstitute
#'   the original matrix
#' @param ... Optional arguments passed to \code{\link{nipals}}
#'
#' @return A numeric matrix with missing values imputed.
#' @seealso \code{\link{impute.nipals}}, \code{\link{pca}}
#' @author Al J Abadi
#' @examples
#' data("nutrimouse")
#' X <- data.matrix(nutrimouse$lipid)
#' ## add missing values to X to impute and compare to actual values
#' set.seed(42)
#' na.ind <- sample(seq_along(X), size = 10)
#' true.values <- X[na.ind]
#' X[na.ind] <- NA
#' X.impute <- impute.nipals(X = X, ncomp = 5)
#' ## compare
#' round(X.impute[na.ind], 2)
#' true.values
#' @export
impute.nipals <- function(X, ncomp, ...)
{
    if (!any(is.na(X)))
    {
        cat("no missing values in 'X' to impute \n")
        return(X)
    }
    
    ## check ncomp is high enough for reliable imputation
    if (ncomp < min(5, min(dim(X))))
        message("consider high 'ncomp' for more accurate ",
                "imputation of the missing values.")
    
    nipals.res <- nipals(X = X, ncomp = ncomp, ...)
    X.impute <- .impute.nipals(X = X, 
                               t = nipals.res$t,
                               eig = nipals.res$eig,
                               p = nipals.res$p
                               )
    return(X.impute)
}

#' Reconstitute a matrix
#' 
#' Using scores, singular values, and loadings from nipals
#' @noRd
#' @keywords Internal
.reconstitute.matrix <- function(t, eig, p)
{
    t %*% diag(eig) %*% t(p)
}

#' Impute missing values
#' 
#' Given scores, singular values, and loadings from nipals
#' @noRd
#' @keywords Internal
.impute.nipals <- function(X, t, eig, p)
{
    X.hat <- .reconstitute.matrix(t = t, eig = eig, p = p)
    X[is.na(X)] <- X.hat[is.na(X)]
    return(X)
}

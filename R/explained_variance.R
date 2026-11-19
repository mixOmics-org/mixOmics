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

# =============================================================================
# Calculate variance explained in one dataset based on its components
# =============================================================================

#' Calculates the proportion of explained variance of multivariate components
#' 
#' \code{explained_variance} calculates the proportion of variance explained by
#' a set of *orthogonal* variates / components and divides by the total variance
#' in \code{data} using the definition of 'redundancy'. This applies to any
#' component-based approaches where components are orthogonal. It is worth
#' noting that any missing values are set to zero (which is the column mean for
#' the centered input data) prior to calculation of total variance in the data.
#' Therefore, this function would underestimate the total variance in presence
#' of abundant missing values. One can use \code{\link{impute.nipals}} function
#' to impute the missing values to avoid such behaviour.
#' 
#' @param data numeric matrix of predictors
#' @param variates variates as obtained from a \code{pls} object for instance
#' @param ncomp number of components. Should be lower than the number of
#' columns of \code{variates}
#' @return \code{explained_variance} returns a named numeric vector containing
#'   the proportion of explained variance for each variate after setting all
#'   missing values in the data to zero (see details).
#' @details Variance explained by component \eqn{t_h} in \eqn{X} for dimension
#'   \eqn{h}: \deqn{Rd(X, t_h) = \frac{1}{p} \sum_{j = 1}^p \mbox{cor}^2(X^j,
#'   t_h)} where \eqn{X^j} is the variable centered and scaled, \eqn{p} is the
#'   total number of variables.
#' @references
#' Tenenhaus, M.,  La Régression PLS théorie et pratique (1998). Technip, Paris, chap2.
#' @author Florian Rohart, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{spls}}, \code{\link{splsda}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{cim}}, \code{\link{network}}.
#' @keywords regression multivariate
#' @examples
#' 
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' 
#' toxicity.spls <- spls(X, Y, ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))
#' 
#' ex = explained_variance(toxicity.spls$X, toxicity.spls$variates$X, ncomp =2)
#' 
#' # ex should be the same as
#' toxicity.spls$prop_expl_var$X
#' 
#' @export
explained_variance <- function(data, variates, ncomp)
{
  #check input data
  check = Check.entry.single(data, ncomp)
  data = check$X
  ncomp = check$ncomp
  ## pre-allocate output
  expl_var <- vector(mode = 'numeric', length = ncomp)
  names(expl_var) <- paste0('comp', seq_len(ncomp))

  data[is.na(data)] <- 0 ## if there is any -- no warning as explained in docs
  norm2.X <- norm(data, type='F')^2 # total variance in the data
  
  for (h in 1:ncomp)
  {
    a <- crossprod(variates[, h, drop=FALSE], data)
    # this is equivalent to calculate the redundancy as detailed in the help file
    expl_var[h] <- tcrossprod(a) / c(crossprod(variates[, h])) / norm2.X
  }
  
  return(expl_var)
}

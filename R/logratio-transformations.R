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
# SPDX-License-Identifier: AGPL-3.0-only AND GPL-2.0-or-later
#
# Code from other packages keeps its own licence. It is listed here and
# marked in the code with SPDX-SnippetBegin and SPDX-SnippetEnd:
#   - ilr.transfo(), clr.backtransfo() and clr.transfo(): adapted from
#     isomLR(), pcaCoDa() and cenLR() of the robCompositions package by
#     Matthias Templ, Karel Hron and Peter Filzmoser, GPL (>= 2) since
#     robCompositions 2.0.6.

#' Log-ratio transformation
#' 
#' This function applies a log transformation to the data, either CLR or ILR
#' 
#' \code{logratio.transfo} applies a log transformation to the data, either CLR
#' (centered log ratio transformation) or ILR (Isometric Log Ratio
#' transformation). In the case of CLR log-transformation, X needs to be a
#' matrix of non-negative values and \code{offset} is used to shift the values
#' away from 0, as commonly done with counts data.
#' 
#' @param X numeric matrix of predictors
#' @param logratio log-ratio transform to apply, one of "none", "CLR" or "ILR"
#' @param offset Value that is added to X for CLR and ILR log transformation.
#' Default to 0.
#' @return \code{logratio.transfo} simply returns the log-ratio transformed
#' data.
#' @author Florian Rohart, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{pca}}, \code{\link{pls}}, \code{\link{spls}},
#' \code{\link{plsda}}, \code{\link{splsda}}.
#' @references Kim-Anh Lê Cao, Mary-Ellen Costello, Vanessa Anne Lakis,
#' Francois Bartolo, Xin-Yi Chua, Remi Brazeilles, Pascale Rondeau mixMC: a
#' multivariate statistical framework to gain insight into Microbial
#' Communities bioRxiv 044206; doi: http://dx.doi.org/10.1101/044206
#' 
#' John Aitchison. The statistical analysis of compositional data. Journal of
#' the Royal Statistical Society. Series B (Methodological), pages 139-177,
#' 1982.
#' 
#' Peter Filzmoser, Karel Hron, and Clemens Reimann. Principal component
#' analysis for compositional data with outliers. Environmetrics,
#' 20(6):621-632, 2009.
#' @examples
#' 
#' data(diverse.16S)
#' CLR = logratio.transfo(X = diverse.16S$data.TSS, logratio = 'CLR')
#' # no offset needed here as we have put it prior to the TSS, see https://mixomics.org/mixmc/
#' @name logratio-transformations
NULL
#' @rdname logratio-transformations
#' @export
logratio.transfo <- function(X,
                             logratio = c('none','CLR','ILR'),
                             offset = 0)
{
    X <- as.matrix(X)
    if (!is.numeric(X))
        stop("X must be a numeric matrix", call. = FALSE)
    
    logratio <- match.arg(logratio)
    if (logratio == 'ILR')
    {
        if (!is(X, 'ilr'))
        {   # data are ilr transformed, then the data lose 1 variable, but we'll use V to reconstruct the matrix
            X = ilr.transfo(X, offset = offset)
        }
    } else if (logratio == 'CLR') {
        X = clr.transfo(X, offset = offset)
    }
    return(X)
}


# SPDX-SnippetBegin
# SPDX-SnippetCopyrightText: Matthias Templ, Karel Hron, Peter Filzmoser (robCompositions)
# SPDX-License-Identifier: GPL-2.0-or-later
# 1 - ilr transform of the data
# -----------------

# KA changed the function to add a min value when many zeroes in data (prob with log and division by 0 otherwise)
#' @param fast if TRUE, it is approx. 10 times faster but numerical problems may occur for high dimensional data
#' @noRd
ilr.transfo = function(x, fast = TRUE, offset = 0)
{
    if(any(x==0) & offset ==0)
        stop("make sure you use pseudo counts before normalisation to avoid 0 values with log ratio transformation")
    # ilr transformation
    x.ilr = matrix(NA, nrow = nrow(x), ncol = ncol(x)-1)
    D = ncol(x)
    # KA added: a little something to avoid 0 values
    if (fast)
    {
        for (i in 1 : ncol(x.ilr))
        {
            x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(((apply(as.matrix(x[, (i+1) : D, drop = FALSE]), 1, prod) + offset)^(1 / (D-i))) / (x[,i]+ offset))
            #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(((apply(as.matrix(x[,(i+1):D,drop = FALSE]),1,prod))^(1/(D-i)))/(x[,i]))
        }
    } else {
        for (i in 1 : ncol(x.ilr))
        {
            x.ilr[,i] = sqrt((D-i) / (D-i+1)) * log(apply(as.matrix(x[, (i+1):D]), 1, function(x){exp(log(x))})/(x[, i]+ offset) + offset)
            #x.ilr[,i] = sqrt((D-i)/(D-i+1))*log(apply(as.matrix(x[,(i+1):D]), 1, function(x){exp(log(x))})/(x[,i]))
        }
    }
    class(x.ilr) = 'ilr'
    return(as.matrix(x.ilr))
}

# 2 - back transformation from ilr to clr space
clr.backtransfo = function(x)
{
    # construct orthonormal basis
    V = matrix(0, nrow = ncol(x), ncol = ncol(x)-1)
    for( i in 1:ncol(V) )
    {
        V[1:i, i] = 1/i
        V[i+1, i] = (-1)
        V[, i] = V[, i] * sqrt(i/(i+1))
    }
    rownames(V) = colnames(x)
    return(V)
    
}


# CLR transformation
clr.transfo = function(x, offset = 0)
{
    if (any(is.na(x) | x < 0)) {
        stop('\nFor CLR transformation, data must be non-negative with no missing values\n', call. = FALSE)
    }
    if(any(x==0) & offset ==0)
        stop("make sure you use pseudo counts before normalisation to avoid 0 values with log ratio transformation")
    
    # KA added
    #offset = min(x[which(x != 0)])*0.01
    
    
    #if (dim(x)[2] < 2) stop("data must be of dimension greater equal 2")
    if (dim(x)[2] == 1)
    {
        res = list(x.clr = x, gm = rep(1, dim(x)[1]))
    } else{
        geometricmean = function (x) {
            #       if (any(na.omit(x == 0)))
            #         0
            #       else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
            #     }
            # KA changed to
            exp(mean(log(x + offset)))
        }
        gm = apply(x, 1, geometricmean)
        # KA changed
        x.clr = log((x + offset) / (gm))
        res = x.clr #list(x.clr = x.clr, gm = gm)
    }
    class(res) = "clr"
    return(res)
}
# SPDX-SnippetEnd


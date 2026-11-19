#############################################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#   Kim-Anh Le Cao, French National Institute for Agricultural Research and
# ARC Centre of Excellence ins Bioinformatics, Institute for Molecular Bioscience, University of Queensland, Australia
#  Leigh Coonan, Queensland Faculty for Advanced Bioinformatics, Australia
#  Fangzhou Yao, Queensland Faculty for Advanced Bioinformatics, Australia
#  Sebastien Dejean, Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
#
# created: 2009
# last modified:
#
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
#   - The centring, scaling, SVD and standard deviation steps of
#     pcasvd(): adapted from prcomp() of the R stats package by the
#     R Core Team, GPL (>= 2).
#############################################################################################################


pcasvd <- function(X,
                   ncomp = 3,
                   retx = TRUE,
                   center = TRUE,
                   scale = FALSE)
{
    X = as.matrix(X)
    
    # SPDX-SnippetBegin
    # SPDX-SnippetCopyrightText: The R Core Team (R stats package)
    # SPDX-License-Identifier: GPL-2.0-or-later
    result = svd(get("scale", envir = .GlobalEnv)(X, center = center, scale = scale), nu = 0)
    cen = attr(X, "scaled:center")
    sc = attr(X, "scaled:scale")
    if (any(sc == 0))
        stop("cannot rescale a constant/zero column to unit variance.")
    # SPDX-SnippetEnd
    
    # If ncomp is not defined by user at the beginning then this algorithm will calculate an ncomp to be used
    if (is.null(ncomp))
        ncomp = min(nrow(X), ncol(X))
    
    # SPDX-SnippetBegin
    # SPDX-SnippetCopyrightText: The R Core Team (R stats package)
    # SPDX-License-Identifier: GPL-2.0-or-later
    if (ncomp < ncol(X))
        result$v = result$v[, 1:ncomp, drop = FALSE]
    
    result$d = result$d/sqrt(max(1, nrow(X) - 1))
    # SPDX-SnippetEnd
    
    if(retx)
    {
        r = list(sdev = result$d, rotation = result$v,
                 X = as.matrix(X) %*% result$v,
                 center = if (is.null(cen)) FALSE else cen,
                 scale = if (is.null(sc)) FALSE else sc)
    }
    
    if(!retx)
        r = list(sdev = result$d, rotation = result$v)
    
    class(r) = c("pca", "prcomp")
    return(invisible(r))
}

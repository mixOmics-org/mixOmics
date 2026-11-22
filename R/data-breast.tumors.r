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

#' Human Breast Tumors Data
#' 
#' This data set contains the expression of 1,000 genes in 47 surgical
#' specimens of human breast tumours from 17 different individuals before and
#' after chemotherapy treatment.
#' 
#' This data consists of 47 breast cancer samples and 1753 cDNA clones
#' pre-selected by Perez-Enciso \emph{et al.} (2003) to draw their Fig. 1. The
#' authors selected 47 samples for which there was information at least before
#' or before and after chemotherapy treatment. There were 20 tumours that were
#' microarrayed both before and after treatment. For illustrative purposes we
#' then randomly selected 1000 cDNA clones for this data set.
#' 
#' @name breast.tumors
#' @docType data
#' @usage data(breast.tumors)
#' @format A list containing the following components: \describe{
#' \item{list("gene.exp")}{data matrix with 47 rows and 1000 columns. Each row
#' represents an experimental sample, and each column a single gene.}
#' \item{list("sample")}{a list containing two character vector components:
#' \code{name} the name of the samples, and \code{treatment} the treatment
#' status.} \item{list("genes")}{a list containing two character vector
#' components: \code{name} the name of the genes, and \code{description} the
#' description of each gene.} }
#' @return none
#' @references Perez-Enciso, M. and Tenenhaus, M. (2003). Prediction of
#' clinical outcome with microarray data: a partial least squares discriminant
#' analysis (PLS-DA) approach. \emph{Human Genetics} \bold{112}, 581-592.
#' 
#' Perou, C. M., Sorlie, T., Eisen, M. B., van de Rijn, M., Jeffrey, S. S.,
#' Rees, C. A., Pollack, J. R., Ross, D. T., Johnsen, H., Akslen, L. A., Fluge,
#' O., Pergamenschikov, A., Williams, C., Zhu, S. X., Lonning, P. E.,
#' Borresen-Dale, A. L., Brown, P. O. and Botstein, D. (2000). Molecular
#' portraits of human breast tumours. \emph{Nature} \bold{406}, 747-752.
#' @source The Human Breast Tumors dataset is a companion resource for the
#' paper of Perou \emph{et al.} (2000), \url{https://doi.org/10.1038/35021093}.
#' It was downloaded from the Stanford Genomics Breast Cancer Consortium
#' Portal, which is no longer online.
#' @keywords datasets
NULL

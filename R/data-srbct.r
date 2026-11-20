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

#' Small version of the small round blue cell tumors of childhood data
#' 
#' This data set from Khan \emph{et al.}, (2001) gives the expression measure
#' of 2308 genes measured on 63 samples.
#' 
#' 
#' @name srbct
#' @docType data
#' @usage data(srbct)
#' @format A list containing the following components: \describe{
#' \item{list("gene")}{data frame with 63 rows and 2308 columns. The expression
#' measure of 2308 genes for the 63 subjects.} \item{list("class")}{A class
#' vector containing the class tumour of each case (4 classes in total).}
#' \item{list("gene.name")}{data frame with 2308 rows and 2 columns containing
#' further information on the genes.} }
#' @return none
#' @references Khan et al. (2001). Classification and diagnostic prediction of
#' cancers using gene expression profiling and artificial neural networks.
#' \emph{Nature Medicine} \bold{7}, Number 6, June.
#' @source \url{https://www.research.nhgri.nih.gov/projects/Microarray/Supplement/index.html}
#' @keywords datasets
NULL

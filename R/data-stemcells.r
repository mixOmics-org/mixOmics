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

#' Human Stem Cells Data
#' 
#' This data set contains the expression of a random subset of 400 genes in 125
#' samples from 4 independent studies and 3 cell types.
#' 
#' This data set contains the expression of a random subset of 400 genes in 125
#' samples from 4 independent studies and 3 cell types. Those studies can be
#' combined and analysed using the MINT procedure.
#' 
#' @name stemcells
#' @docType data
#' @usage data(stemcells)
#' @format A list containing the following components: \describe{
#' \item{list("gene")}{data matrix with 125 rows and 400 columns. Each row
#' represents an experimental sample, and each column a single gene.}
#' \item{list("celltype")}{a factor indicating the cell type of each sample.}
#' \item{list("study")}{a factor indicating the study from which the sample was
#' extracted.} }
#' @return none
#' @references Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017).
#' MINT: A multivariate integrative approach to identify a reproducible
#' biomarker signature across multiple experiments and platforms. BMC
#' Bioinformatics 18:128.
#' @keywords datasets
NULL

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

#' Yeast metabolomic study
#' 
#' Two Saccharomyces Cerevisiae strains were compared under two different
#' environmental conditions, 37 metabolites expression are measured.
#' 
#' In this study, two Saccharomyces cerevisiae strains were used - wild-type
#' (WT) and mutant (MT), and were carried out in batch cultures under two
#' different environmental conditions, aerobic (AER) and anaerobic (ANA) in
#' standard mineral media with glucose as the sole carbon source. After
#' normalization and pre processing, the metabolomic data results in 37
#' metabolites and 55 samples which include 13 MT-AER, 14 MT-ANA, 15 WT-AER and
#' 13 WT-ANA samples
#' 
#' @name yeast
#' @docType data
#' @usage data(yeast)
#' @format A list containing the following components: \describe{
#' \item{list("data")}{data matrix with 55 rows and 37 columns. Each row
#' represents an experimental sample, and each column a single metabolite.}
#' \item{list("strain")}{a factor containing the type of strain (MT or WT).}
#' \item{list("condition")}{a factor containing the type of environmental
#' condition (AER or ANA).} \item{list("strain.condition")}{a crossed factor
#' between \code{strain} and \code{condition}.} }
#' @return none
#' @references Villas-Boas S, Moxley J, Akesson M, Stephanopoulos G, Nielsen J:
#' High-throughput metabolic state analysis (2005). The missing link in
#' integrated functional genomics. \emph{Biochemical Journal},
#' \bold{388}:669–677.
#' @keywords datasets
NULL

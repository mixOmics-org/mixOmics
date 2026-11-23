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

#' Liver Toxicity Data
#' 
#' This data set contains the expression measure of 3116 genes and 10 clinical
#' measurements for 64 subjects (rats) that were exposed to non-toxic,
#' moderately toxic or severely toxic doses of acetaminophen in a controlled
#' experiment.
#' 
#' The data come from a liver toxicity study (Bushel \emph{et al}., 2007) in
#' which 64 male rats of the inbred strain Fisher 344 were exposed to non-toxic
#' (50 or 150 mg/kg), moderately toxic (1500 mg/kg) or severely toxic (2000
#' mg/kg) doses of acetaminophen (paracetamol) in a controlled experiment.
#' Necropsies were performed at 6, 18, 24 and 48 hours after exposure and the
#' mRNA from the liver was extracted. Ten clinical chemistry measurements of
#' variables containing markers for liver injury are available for each subject
#' and the serum enzymes levels are measured numerically. The data were further
#' normalized and pre-processed by Bushel \emph{et al}. (2007).
#' 
#' @name liver.toxicity
#' @docType data
#' @usage data(liver.toxicity)
#' @format A list containing the following components: \describe{
#' \item{list("gene")}{data frame with 64 rows and 3116 columns. The expression
#' measure of 3116 genes for the 64 subjects (rats).}
#' \item{list("clinic")}{data frame with 64 rows and 10 columns, containing 10
#' clinical variables for the same 64 subjects.} \item{list("treatment")}{data
#' frame with 64 rows and 4 columns, containing the treatment information on
#' the 64 subjects, such as doses of acetaminophen and times of necropsies.}
#' \item{list("gene.ID")}{data frame with 3116 rows and 2 columns, containing
#' geneBank IDs and gene titles of the annotated genes}
#' 
#' }
#' @return none
#' @references Bushel, P., Wolfinger, R. D. and Gibson, G. (2007). Simultaneous
#' clustering of gene expression data with clinical chemistry and pathological
#' evaluations reveals phenotypic prototypes. \emph{BMC Systems Biology}
#' \bold{1}, Number 15.
#' 
#' Lê Cao, K.-A., Rossouw, D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. \emph{Statistical
#' Applications in Genetics and Molecular Biology} \bold{7}, article 35.
#' @source The two liver toxicity data sets are a companion resource for the
#' paper of Bushel \emph{et al.} (2007), and were downloaded from the
#' additional files of that paper:
#' 
#' \url{https://doi.org/10.1186/1752-0509-1-15}
#' @keywords datasets
NULL

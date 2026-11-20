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

#' Simulated data based on the vac18 study for multilevel analysis
#' 
#' Simulated data based on the vac18 study to illustrate the use of the
#' multilevel analysis for one and two-factor analysis with sPLS-DA. This data
#' set contains the expression simulated of 500 genes.
#' 
#' In this cross-over design, repeated measurements are performed 12
#' experiments units (or unique subjects) for each of the 4 stimulations.
#' 
#' The simulation study was based on a mixed effects model (see reference for
#' details). Ten clusters of 100 genes were generated. Amongt those, 4 clusters
#' of genes discriminate the 4 stimulations (denoted LIPO5, GAG+, GAG- and NS)
#' as follows: \ -2 gene clusters discriminate (LIPO5, GAG+) versus (GAG-, NS)
#' \ -2 gene clusters discriminate LIPO5 versus GAG+, while GAG+ and NS have
#' the same effect \ -2 gene clusters discriminate GAG- versus NS, while LIPO5
#' and GAG+ have the same effect \ -the 4 remaining clusters represent noisy
#' signal (no stimulation effect) \
#' 
#' Only a subset of those genes are presented here (to save memory space).
#' 
#' @name vac18.simulated
#' @docType data
#' @usage data(vac18.simulated)
#' @format A list containing the following components: \describe{
#' \item{list("genes")}{data frame with 48 rows and 500 columns. The simulated
#' expression of 500 genes for 48 subjects.} \item{list("sample")}{a vector
#' indicating the repeated measurements on each unique subject. See Details.}
#' \item{list("stimulation")}{a factor indicating the stimulation condition on
#' each sample.} \item{list("time")}{a factor indicating the time condition on
#' each sample.} }
#' @return none
#' @references Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. (2012). A
#' novel approach for biomarker selection and the integration of repeated
#' measures experiments from two platforms. \emph{BMC Bioinformatics}
#' \bold{13}:325.
#' @keywords datasets
NULL

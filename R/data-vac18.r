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

#' Vaccine study Data
#' 
#' The data come from a trial evaluating a vaccine based on HIV-1 lipopeptides
#' in HIV-negative volunteers. The vaccine (HIV-1 LIPO-5 ANRS vaccine) contains
#' five HIV-1 amino acid sequences coding for Gag, Pol and Nef proteins. This
#' data set contains the expression measure of a subset of 1000 genes from
#' purified in vitro stimulated Peripheral Blood Mononuclear Cells from 42
#' repeated samples (12 unique vaccinated participants) 14 weeks after
#' vaccination, , 6 hours after in vitro stimulation by either (1) all the
#' peptides included in the vaccine (LIPO-5), or (2) the Gag peptides included
#' in the vaccine (GAG+) or (3) the Gag peptides not included in the vaccine
#' (GAG-) or (4) without any stimulation (NS).
#' 
#' This is a subset of the original study for illustrative purposes.
#' 
#' @name vac18
#' @docType data
#' @usage data(vac18)
#' @format A list containing the following components: \describe{
#' \item{list("gene")}{data frame with 42 rows and 1000 columns. The expression
#' measure of 1000 genes for the 42 samples (PBMC cells from 12 unique
#' subjects).} \item{list("stimulation")}{is a fctor of 42 elements indicating
#' the type of in vitro simulation for each sample.} \item{list("sample")}{is a
#' vector of 42 elements indicating the unique subjects (for example the value
#' '1' correspond to the first patient PBMC cells). Note that the design of
#' this study is unbalanced.} \item{list("tab.prob.gene")}{is a data frame with
#' 1000 rows and 2 columns, indicating the Illumina probe ID and the gene name
#' of the annotated genes.} }
#' @return none
#' @references Salmon-Ceron D, Durier C, Desaint C, Cuzin L, Surenaud M,
#' Hamouda N, Lelievre J, Bonnet B, Pialoux G, Poizot-Martin I, Aboulker J,
#' Levy Y, Launay O, trial group AV: Immunogenicity and safety of an HIV-1
#' lipopeptide vaccine in healthy adults: a phase 2 placebo-controlled ANRS
#' trial. AIDS 2010, 24(14):2211-2223.
#' @keywords datasets
NULL

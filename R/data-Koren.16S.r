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

#' 16S microbiome atherosclerosis study
#' 
#' The 16S data come from Koren et al. (2011) and compared the bodysites oral,
#' gut and plaque microbial communities in patients with atherosclerosis. The
#' data can be analysed with our mixMC module. The data include 43 samples
#' measured on 980 OTUs.
#' 
#' The data are from Koren et al. (2011) who examined the link between oral,
#' gut and plaque microbial communities in patients with atherosclerosis and
#' controls. Only healthy individuals were retained in the analysis. This study
#' contained partially repeated measures from multiple sites including 15
#' unique patients samples from saliva and stool, and 13 unique patients only
#' sampled from arterial plaque samples and we therefore considered a non
#' multilevel analysis for that experimental design. After prefiltering, the
#' data included 973 OTU for 43 samples. We strongly recommend to use log ratio
#' transformations on the \code{data.TSS} normalisd data, as implemented in the
#' PLS and PCA methods, see details on \url{https://mixomics.org/mixmc/}.
#' 
#' The \code{data.raw} include a 1 offset in order to be log ratios transformed
#' after TSS normalisation. Consequently, the \code{data.TSS} are TSS
#' normalisation of \code{data.raw}. The CSS normalisation was performed on the
#' orignal data (including zero values)
#' 
#' @name Koren.16S
#' @docType data
#' @usage data(Koren.16S)
#' @format A list containing two data sets, \code{data.TSS} and \code{data.raw}
#' and some meta data information: \describe{ \item{list("data.TSS")}{data
#' frame with 43 rows (samples) and 980 columns (OTUs). The prefiltered
#' normalised data using Total Sum Scaling normalisation.}
#' \item{list("data.raw")}{data frame with 43 rows (samples) and 980 columns
#' (OTUs). The prefiltered raw count OTU data which include a 1 offset (i.e. no
#' 0 values).} \item{list("taxonomy")}{data frame with 980 rows (OTUs) and 7
#' columns indicating the taxonomy of each OTU.} \item{list("indiv")}{data
#' frame with 43 rows indicating sample meta data.}
#' \item{list("bodysite")}{factor of length 43 indicating the bodysite with
#' levels arterial plaque, saliva and stool.} }
#' @return none
#' @references Lê Cao K.-A., Costello ME, Lakis VA, Bartolo, F,Chua XY,
#' Brazeilles R, Rondeau P. MixMC: Multivariate insights into Microbial
#' Communities. PLoS ONE, 11(8): e0160169 (2016).
#' 
#' Koren, O., Spor, A., Felin, J., Fak, F., Stombaugh, J., Tremaroli, V., et
#' al.: Human oral, gut, and plaque microbiota in patients with
#' atherosclerosis. Proceedings of the National Academy of Sciences
#' 108(Supplement 1), 4592-4598 (2011)
#' @source The raw data were downloaded from the QIITA database. Filtering and
#' normalisation described in our website \url{https://mixomics.org/mixmc/}
#' @keywords datasets
NULL

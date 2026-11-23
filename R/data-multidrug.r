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

#' Multidrug Resistence Data
#' 
#' This data set contains the expression of 48 known human ABC transporters
#' with patterns of drug activity in 60 diverse cancer cell lines (the NCI-60)
#' used by the National Cancer Institute to screen for anticancer activity.
#' 
#' The data come from a pharmacogenomic study (Szakacs \emph{et al.}, 2004) in
#' which two kinds of measurements acquired on the NCI-60 cancer cell lines are
#' considered: \itemize{ \item the expression of the 48 human ABC transporters
#' measured by real-time quantitative RT-PCR for each cell line; \item the
#' activity of 1429 drugs expressed as \eqn{GI_{50}} which corresponds to the
#' concentration at which the drug induces \eqn{50\%} inhibition of cellular
#' growth for the cell line tested. }
#' 
#' The NCI- 60 panel includes cell lines derived from cancers of colorectal (7
#' cell lines), renal(8), ovarian(6), breast(8), prostate(2), lung(9) and
#' central nervous system origin(6), as well as leukemias(6) and melanomas(8).
#' It was set up by the Developmental Therapeutics Program of the National
#' Cancer Institute (NCI, one of the U.S. National Institutes of Health) to
#' screen the toxicity of chemical compound repositories. The expressions of
#' the 48 human ABC transporters is available as a supplement to the paper of
#' Szak?cs \emph{et al.} (2004).
#' 
#' The drug dataset consiste of 118 compounds whose mechanisms of action are
#' putatively classifiable (Weinstein \emph{et al.}, 1992) and a larger set of
#' 1400 compounds that have been tested multiple times and whose screening data
#' met quality control criteria described elsewhere (Scherf \emph{et al.},
#' 2000). The two were combined to form a joint dataset that included 1429
#' compounds.
#' 
#' @name multidrug
#' @docType data
#' @usage data(multidrug)
#' @format A list containing the following components: \describe{
#' \item{list("ABC.trans")}{data matrix with 60 rows and 48 columns. The
#' expression of the 48 human ABC transporters.} \item{list("compound")}{data
#' matrix with 60 rows and 1429 columns. The activity of 1429 drugs for the 60
#' cell lines.} \item{list("comp.name")}{character vector. The names or the NSC
#' No. of the 1429 compounds.} \item{list("cell.line")}{a list containing two
#' character vector components: \code{Sample} the names of the 60 cell line
#' which were analysed, and \code{Class} the phenotypes of the 60 cell lines.}
#' }
#' @return none
#' @references Scherf, U., Ross, D. T., Waltham, M., Smith, L. H., Lee, J. K.,
#' Tanabe, L., Kohn, K. W., Reinhold, W. C., Myers, T. G., Andrews, D. T.,
#' Scudiero, D. A., Eisen, M. B., Sausville, E. A., Pommier, Y., Botstein, D.,
#' Brown, P. O. and Weinstein, J. N. (2000). A Gene Expression Database for the
#' Molecular Pharmacology of Cancer. \emph{Nature Genetics}, \bold{24},
#' 236-244.
#' 
#' Szakacs, G., Annereau, J.-P., Lababidi, S., Shankavaram, U., Arciello, A.,
#' Bussey, K. J., Reinhold, W., Guo, Y., Kruh, G. D., Reimers, M., Weinstein,
#' J. N. and Gottesman, M. M. (2004). Predicting drug sensivity and resistance:
#' Profiling ABC transporter genes in cancer cells. \emph{Cancer Cell}
#' \bold{4}, 147-166.
#' 
#' Weinstein, J.N., Kohn, K.W., Grever, M.R., Viswanadhan, V.N., Rubinstein,
#' L.V., Monks, A.P., Scudiero, D.A., Welch, L., Koutsoukos, A.D., Chiausa,
#' A.J. et al. 1992. Neural computing in cancer drug development: Predicting
#' mechanism of action. \emph{Science} \bold{258}, 447-451.
#' @source The NCI dataset was downloaded from The Genomics and Bioinformatics
#' Group Supplemental Table S1 to the paper of Szakacs \emph{et al.} (2004),
#' \url{https://doi.org/10.1016/j.ccr.2004.06.026}
#' 
#' The two drug data sets are a companion resource for the paper of Scherf
#' \emph{et al.} (2000), \url{https://doi.org/10.1038/73439}. They were
#' downloaded from the Genomics and Bioinformatics Group's web site; that page
#' is no longer available.
#' @keywords datasets
NULL

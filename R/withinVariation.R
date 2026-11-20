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

#' Within matrix decomposition for repeated measurements (cross-over design)
#' 
#' This function is internally called by \code{pca}, \code{pls}, \code{spls},
#' \code{plsda} and \code{splsda} functions for cross-over design data, but can
#' be called independently prior to any kind of multivariate analyses.
#' 
#' 
#' \code{withinVariation} function decomposes the Within variation in the
#' \eqn{X} data set. The resulting \eqn{Xw} matrix is then input in the
#' \code{multilevel} function.
#' 
#' One or two-factor analyses are available.
#' 
#' @param X numeric matrix of predictors. \code{NA}s are allowed.
#' @param design a numeric matrix or data frame. The first column indicates the
#' repeated measures on each individual, i.e. the individuals ID. The 2nd and
#' 3rd columns are to split the variation for a 2 level factor.
#' @return \code{withinVariation} simply returns the \eqn{Xw} within matrix,
#' which can be input in the other multivariate approaches already implemented
#' in mixOmics (i.e. spls or splsda, see \code{multilevel}, but also pca or
#' ipca).
#' @author Benoit Liquet, Kim-Anh Lê Cao, Benoit Gautier, Ignacio González,
#' Florian Rohart, AL J Abadi
#' @seealso \code{\link{spls}}, \code{\link{splsda}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{cim}}, \code{\link{network}}.
#' @references On multilevel analysis:
#' 
#' Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel
#' approach for biomarker selection and the integration of repeated measures
#' experiments from two platforms. \emph{BMC Bioinformatics} \bold{13}:325.
#' 
#' Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K.
#' (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA.
#' \emph{Metabolomics}, \bold{6}(1), 119-128.
#' @keywords regression multivariate
#' @export
#' @examples
#' ## Example: one-factor analysis matrix decomposition
#' #--------------------------------------------------------------
#' data(vac18)
#' X <- vac18$genes
#' # in design we only need to mention the repeated measurements to split the one level variation
#' design <- data.frame(sample = vac18$sample)
#' 
#' Xw <- withinVariation(X = X, design = design)
#' # multilevel PCA
#' res.pca.1level <- pca(Xw, ncomp = 3)
#' 
#' # compare a normal PCA with a multilevel PCA for repeated measurements.
#' # note: PCA makes the assumptions that all samples are independent,
#' # so this analysis is flawed and you should use a multilevel PCA instead
#' res.pca <- pca(X, ncomp = 3)
#' 
#' # set up colors for plotIndiv
#' col.stim <- c("darkblue", "purple", "green4","red3")
#' col.stim <- col.stim[as.numeric(vac18$stimulation)]
#' 
#' # plotIndiv comparing both PCA and PCA multilevel
#' plotIndiv(res.pca, ind.names = vac18$stimulation, group = col.stim)
#' title(main = 'PCA ')
#' plotIndiv(res.pca.1level, ind.names = vac18$stimulation, group = col.stim)
#' title(main = 'PCA multilevel')
#' 
withinVariation <-
  function(X, design){ 
    
    # need a matrix for matrix calculations
    X = as.matrix(X)
    rep.measures = factor(design[, 1])
    factors = design[, -1, drop = FALSE] 
    
    if(any(summary(as.factor(rep.measures)) == 1))
      stop("Multilevel analysis can only be performed when all samples are repeated at least once.", call. = FALSE)
    
    # calculate the variation
    # ---------------------------
    # added condition for the spls case where the condition is not needed
    # all we need is the rep.measures
    if ((ncol(factors) == 0) | (ncol(factors) == 1))
    {
      message("Splitting the variation for 1 level factor.")
      
      # save sample names for the output
      indiv.names = rownames(X)
      rownames(X) = as.character(rep.measures)
      
      # compute the mean for each unique individual
      # dealing with specific case with only one subject (leave one out case during prediction)
      X.mean.indiv = matrix(apply(X, 2, tapply, rep.measures, mean, na.rm = TRUE), # to deal with only one subject
                            nrow = length(unique(rep.measures)), ncol = dim(X)[2], dimnames = list(levels(as.factor(rep.measures)), colnames(X)))
      
      # fill the between matrix with those means 
      Xb = X.mean.indiv[as.character(rep.measures), ]
      
      # compute the within matrix as a difference between the original data 
      # and the between matrix
      Xw = X - Xb
      
      # put dimnames
      dimnames(Xw) = list(indiv.names, colnames(X))    
    } else {  # for 2 levels split
      message("Splitting the variation for 2 level factors.")
      
      ###### off set term
      Xm = colMeans(X)
      
      ###### compute the mean within each subject
      Xs = apply(X, 2, tapply, rep.measures, mean, na.rm = TRUE)
      Xs = Xs[rep.measures, ]
      
      # for the first factor
      xbfact1 = apply(X, 2, tapply, paste0(rep.measures, factors[, 1]), mean, na.rm = TRUE)
      xbfact1 = xbfact1[paste0(rep.measures, factors[, 1]), ]
      
      # for the second factor
      xbfact2 = apply(X, 2, tapply, paste0(rep.measures, factors[, 2]), mean, na.rm = TRUE)
      xbfact2 = xbfact2[paste0(rep.measures, factors[, 2]), ]
      
      #### fixed effect
      # for the first factor
      Xmfact1 = apply(X, 2, tapply, factors[, 1], mean, na.rm = TRUE)
      Xmfact1 = Xmfact1[factors[, 1], ]
      
      # for the second factor
      Xmfact2 = apply(X, 2, tapply, factors[, 2], mean, na.rm = TRUE)
      Xmfact2 = Xmfact2[factors[, 2], ]
      
      # compute the within matrix 
      Xw = X + Xs - xbfact1 - xbfact2 + Xmfact1 + Xmfact2
      Xw = sweep(Xw, 2, 2*Xm)  # see formula in refernece Liquet et al.
      
      # put dimnames
      dimnames(Xw) = dimnames(X)  
    }
    
    return(invisible(Xw))
  }


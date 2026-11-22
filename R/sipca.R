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
# SPDX-License-Identifier: AGPL-3.0-only AND (GPL-2.0-only OR GPL-3.0-only)
#
# Code from other packages keeps its own licence. It is listed here and
# marked in the code with SPDX-SnippetBegin and SPDX-SnippetEnd:
#   - The checks of X, ncomp and w.init in sipca(): adapted from fastICA()
#     of the fastICA package (version 1.1-13) by J L Marchini, C Heaton
#     and B D Ripley, GPL-2 | GPL-3.

#' Independent Principal Component Analysis
#' 
#' Performs sparse independent principal component analysis on the given data
#' matrix to enable variable selection.
#' 
#' See Details of ipca.
#' 
#' Soft thresholding is implemented on the independent loading vectors to
#' obtain sparse loading vectors and enable variable selection.
#' 
#' @inheritParams ipca
#' @param keepX the number of variable to keep on each dimensions.
#' @return \code{pca} returns a list with class \code{"ipca"} containing the
#' following components: \item{ncomp}{the number of principal components used.}
#' \item{unmixing}{the unmixing matrix of size (ncomp x ncomp)}
#' \item{mixing}{the mixing matrix of size (ncomp x ncomp} \item{X}{the
#' centered data matrix} \item{x}{the principal components (with sparse
#' independent loadings)} \item{loadings}{the sparse independent loading
#' vectors} 
#' \item{kurtosis}{the kurtosis measure of the independent loading
#' vectors}
#' \item{prop_expl_var}{Proportion of the explained variance of derived
#' components, after setting possible missing values to zero.}
#' @author Fangzhou Yao, Jeff Coquery, Francois Bartolo, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{ipca}}, \code{\link{pca}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}} and https://mixomics.org/ for more details.
#' @references Yao, F., Coquery, J. and Lê Cao, K.-A. (2011) Principal
#' component analysis with independent loadings: a combination of PCA and ICA.
#' (in preparation)
#' 
#' A. Hyvarinen and E. Oja (2000) Independent Component Analysis: Algorithms
#' and Applications, \emph{Neural Networks}, \bold{13(4-5)}:411-430
#' 
#' J L Marchini, C Heaton and B D Ripley (2010). fastICA: FastICA Algorithms to
#' perform ICA and Projection Pursuit. R package version 1.1-13.
#' @keywords algebra
#' @export
#' @example ./examples/sipca-examples.R
sipca <-
    function (X,
              ncomp  = 3,
              mode = c("deflation", "parallel"),
              fun = c("logcosh", "exp"),
              scale = FALSE,
              max.iter = 200,
              tol = 1e-04,
              keepX = rep(50, ncomp),
              w.init = NULL)
    {
        
        # SPDX-SnippetBegin
        # SPDX-SnippetCopyrightText: J L Marchini, C Heaton, B D Ripley (fastICA)
        # SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
        dim_x <- dim(X)
        d <- dim_x[dim_x != 1]
        if (length(d) != 2)
            stop("data must be in a matrix form")
        X <- if (length(d) != length(dim_x))
        {matrix(X, d[1], d[2])}
        else {as.matrix(X)}
        # SPDX-SnippetEnd
        
        alpha <- 1
        
        mode <- match.arg(mode)
        fun <- match.arg(fun)
        
        X.names = dimnames(X)[[2]]
        if (is.null(X.names)) X.names = paste("X", 1:ncol(X), sep = "")
        
        ind.names = dimnames(X)[[1]]
        if (is.null(ind.names)) ind.names = 1:nrow(X)
        
        X <- scale(X, scale = FALSE)
        if (scale) {X=scale(X, scale=scale)}
        svd_mat <- svd(X)
        right_sing_vect <- svd_mat$v
        right_sing_vect <- scale(right_sing_vect, center=TRUE, scale=TRUE)
        n <- nrow(t(X))
        p <- ncol(t(X))
        
        # SPDX-SnippetBegin
        # SPDX-SnippetCopyrightText: J L Marchini, C Heaton, B D Ripley (fastICA)
        # SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
        if (ncomp > min(n, p)) {
            message("'ncomp' is too large: reset to ", min(n, p))
            ncomp <- min(n, p)
        }
        if(is.null(w.init))
            w.init <- matrix(1/sqrt(ncomp),ncomp,ncomp)
        else {
            if(!is.matrix(w.init) || length(w.init) != (ncomp^2))
                stop("w.init is not a matrix or is the wrong size")
        }
        # SPDX-SnippetEnd
        
        X1 <- t(right_sing_vect)[1:ncomp,]
        
        if (mode == "deflation") {
            unmix_mat <- ica.def(X1, ncomp, tol = tol, fun = fun,
                                 alpha = alpha, max.iter = max.iter, verbose = FALSE, w.init = w.init)
        }
        else if (mode == "parallel") {
            unmix_mat <- ica.par(X1, ncomp, tol = tol, fun = fun,
                                 alpha = alpha, max.iter = max.iter, verbose = FALSE, w.init = w.init)
        }
        w <- unmix_mat 
        independent_mat <- w %*% X1
        #==order independent_mat by kurtosis==#
        kurt <- vector(length=ncomp)
        independent_mat.new <- matrix(nrow = ncomp, ncol = n)
        for(h in 1:ncomp){
            kurt[h] <- (mean(independent_mat[h,]^4)-3*(mean(independent_mat[h,]^2))^2)
        }
        for(i in 1:ncomp){
            independent_mat.new[i,] <- independent_mat[order(kurt,decreasing=TRUE)[i],]
            independent_mat.new[i,] <- independent_mat.new[i,]/as.vector(crossprod(independent_mat.new[i,]))
        }
        
        #== variable selection==#  
        v.sparse=matrix(nrow = ncomp, ncol = n)
        for(i in 1:ncomp){
            nx <- n - keepX[i]
            v.sparse[i,] = ifelse(abs(independent_mat.new[i,]) > abs(independent_mat.new[i,][order(abs(independent_mat.new[i,]))][nx]), 
                                  (abs(independent_mat.new[i,]) - abs(independent_mat.new[i,][order(abs(independent_mat.new[i,]))][nx])) * sign(independent_mat.new[i,]), 0)
        }
        independent_mat.new = v.sparse
        
        
        mix_mat <- t(w) %*% solve(w %*% t(w))
        
        ipc_mat = matrix(nrow=p, ncol=ncomp)
        ipc_mat = X %*% t(independent_mat.new)        
        ##== force orthogonality ==##
        for(h in 1:ncomp){
            if(h==1){ipc_mat[,h]=X %*% (t(independent_mat.new)[,h])}
            if(h>1){ipc_mat[,h]=(lsfit(y=X%*%(t(independent_mat.new)[,h]), ipc_mat[,1:(h-1)],intercept=FALSE)$res)}
            ipc_mat[,h]=ipc_mat[,h]/as.vector(sqrt(crossprod(ipc_mat[,h])))
        }
        ##== force over ==##          
        # put rownames of loading vectors
        colnames(independent_mat.new) = colnames(X)
        
        cl = match.call()
        cl[[1]] = as.name('sipca')
        
        result = (list(call=cl, X = X, ncomp=ncomp, keepX=keepX, unmixing = t(unmix_mat), mixing = t(mix_mat), loadings = list(X=t(independent_mat.new)), rotation = t(independent_mat.new),
                       kurtosis = kurt[order(kurt,decreasing=TRUE)],names = list(X = X.names, sample = ind.names)))
        
        result$x = ipc_mat
        result$variates=list(X=ipc_mat)
        dimnames(result$x) = list(ind.names, paste("IPC", 1:ncol(result$rotation), sep = " "))
        
        class(result) = c("sipca","ipca","pca")
        
        #calcul explained variance
        explX=explained_variance(X,result$variates$X,ncomp)
        result$prop_expl_var = list(X = explX)
        
        
        
        return(invisible(result))
    } 

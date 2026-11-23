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

#' Sparse Principal Components Analysis
#' 
#' Performs a sparse principal component analysis for variable selection using
#' singular value decomposition and lasso penalisation on the loading vectors.
#' 
#' \code{scale= TRUE} is highly recommended as it will help obtaining orthogonal
#' sparse loading vectors.
#' 
#' \code{keepX} is the number of variables to select in each loading vector,
#' i.e. the number of variables with non zero coefficient in each loading
#' vector.
#' 
#' Note that data can contain missing values only when \code{logratio = 'none'}
#' is used. In this case, \code{center=TRUE} should be used to center the data
#' in order to effectively ignore the missing values. This is the default
#' behaviour in \code{spca}.
#' 
#' According to Filzmoser et al., a ILR log ratio transformation is more
#' appropriate for PCA with compositional data. Both CLR and ILR are valid.
#' 
#' Logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#' 
#' Logratio can only be applied if the data do not contain any 0 value (for
#' count data, we thus advise the normalise raw data with a 1 offset). For ILR
#' transformation and additional offset might be needed.
#' 
#' The principal components are not guaranteed to be orthogonal in sPCA. 
#' We adopt the approach of Shen and Huang 2008 (Section 2.3) to estimate 
#' the explained variance  in the case where the sparse loading vectors 
#' (and principal components) are not orthogonal. The data are projected 
#' onto the space spanned by the first loading vectors and the variance 
#' explained is then adjusted for potential correlation between PCs. 
#' Note that in practice, the loading vectors tend to be orthogonal if the
#' data are centered and scaled in sPCA.
#' 
#' @inheritParams pca
#' @param X a numeric matrix (or data frame) which provides the data for the sparse
#' principal components analysis. It should not contain missing values.
#' @param scale (Default=TRUE) Logical indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. 
#' @param keepX numeric vector of length \code{ncomp}, the number of variables to keep
#' in loading vectors. By default all variables are kept in the model. See
#' details.
#' @param logratio one of ('none','CLR'). Specifies the log ratio
#' transformation to deal with compositional values that may arise from
#' specific normalisation in sequencing data. Default to 'none'
#' @template arg/verbose.call
#' @return \code{spca} returns a list with class \code{"spca"} containing the
#' following components:
#' \describe{
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' \item{ncomp}{the number of components to keep in the
#' calculation.} 
#' \item{prop_expl_var}{the adjusted percentage of variance
#' explained for each component.} 
#' \item{cum.var}{the adjusted cumulative percentage of variances
#' explained.}
#' \item{keepX}{the number of variables kept in each loading
#' vector.} 
#' \item{iter}{the number of iterations needed to reach convergence
#' for each component.} 
#' \item{rotation}{the matrix containing the sparse
#' loading vectors.} 
#' \item{x}{the matrix containing the principal components.}
#' }
#' @author Kim-Anh Lê Cao, Fangzhou Yao, Leigh Coonan, Ignacio Gonzalez, Al J Abadi
#' @seealso \code{\link{pca}} and https://mixomics.org/ for more details.
#' @references Shen, H. and Huang, J. Z. (2008). Sparse principal component
#' analysis via regularized low rank matrix approximation. \emph{Journal of
#' Multivariate Analysis} \bold{99}, 1015-1034.
#' @keywords algebra
#' @export
#' @example ./examples/spca-examples.R
spca <- 
    function(X, 
             ncomp = 2,
             center = TRUE, 
             scale = TRUE,
             keepX = rep(ncol(X), ncomp),
             max.iter = 500, 
             tol = 1e-06,
             logratio = c('none','CLR'),
             multilevel = NULL,
             verbose.call = FALSE)
    {
        
        #-- checking general input parameters --------------------------------------#
        #---------------------------------------------------------------------------#
        
        #-- check that the user did not enter extra arguments
        arg.call = match.call()
        user.arg = names(arg.call)[-1]
        
        err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                       error = function(e) e)
        
        if ("simpleError" %in% class(err))
            stop(err[[1]], ".", call. = FALSE)
        
        #-- X matrix
        if (is.data.frame(X))
            X = as.matrix(X)
        
        if (!is.matrix(X) || is.character(X))
            stop("'X' must be a numeric matrix.", call. = FALSE)
        
        if (any(apply(X, 1, is.infinite)))
            stop("infinite values in 'X'.", call. = FALSE)
        
        #-- ncomp
        if (is.null(ncomp))
            ncomp = min(nrow(X),ncol(X))
        
        ncomp = round(ncomp)
        
        if ( !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
            stop("invalid value for 'ncomp'.", call. = FALSE)
        
        if (ncomp > min(ncol(X), nrow(X)))
            stop("use smaller 'ncomp'", call. = FALSE)
        
        #-- keepX
        if (length(keepX) != ncomp)
            stop("length of 'keepX' must be equal to ", ncomp, ".")
        if (any(keepX > ncol(X)))
            stop("each component of 'keepX' must be no greater than ", ncol(X), ".")
        
        #-- log.ratio
        choices = c('CLR','none')
        logratio = choices[pmatch(logratio, choices)]
        logratio <- match.arg(logratio)
        
        if (logratio != "none")
        {
            ## ensure data do not contain NA or negative values
            if (any(is.na(X)))
            {
                stop("data must not contain any missing or negative values when using logratio ",
                     "transformation. Consider imputing the missing values using ",
                     "mixOmics::nipals(..., reconst = TRUE)", call. = FALSE)
            } else if (any(X < 0))
            {
                stop("'X' must contain nonnegative values for log-transformation. ",
                     "See details.", call. = FALSE)
            }
        }

                
        #-- cheking center and scale
        if (!is.logical(center))
        {
            if (!is.numeric(center) || (length(center) != ncol(X)))
                stop("'center' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
                     call. = FALSE)
        }
        
        if (!is.logical(scale))
        {
            if (!is.numeric(scale) || (length(scale) != ncol(X)))
                stop("'scale' should be either a logical value or a numeric vector of length equal to the number of columns of 'X'.",
                     call. = FALSE)
        }
        
        #-- max.iter
        if (is.null(max.iter) || !is.numeric(max.iter) || max.iter < 1 || !is.finite(max.iter))
            stop("invalid value for 'max.iter'.", call. = FALSE)
        
        max.iter = round(max.iter)
        
        #-- tol
        if (is.null(tol) || !is.numeric(tol) || tol < 0 || !is.finite(tol))
            stop("invalid value for 'tol'.", call. = FALSE)
        
        #-- end checking --#
        #------------------#
        
        #-----------------------------#
        #-- logratio transformation --#
        X = logratio.transfo(X = X, logratio = logratio, offset = 0)#if(logratio == "ILR") {ilr.offset} else {0})
        
        #-- logratio transformation --#
        #-----------------------------#
        
        #---------------------------------------------------------------------------#
        #-- multilevel approach ----------------------------------------------------#
        
        if (!is.null(multilevel))
        {
            # we expect a vector or a 2-columns matrix in 'Y' and the repeated measurements in 'multilevel'
            multilevel = data.frame(multilevel)
            
            if ((nrow(X) != nrow(multilevel)))
                stop("unequal number of rows in 'X' and 'multilevel'.")
            
            if (ncol(multilevel) != 1)
                stop("'multilevel' should have a single column for the repeated measurements.")
            
            multilevel[, 1] = as.numeric(factor(multilevel[, 1])) # we want numbers for the repeated measurements
            
            Xw = withinVariation(X, design = multilevel)
            X = Xw
        }
        #-- multilevel approach ----------------------------------------------------#
        #---------------------------------------------------------------------------#
        
        
        #--scaling the data--#
        X=scale(X,center=center,scale=scale)
        cen = attr(X, "scaled:center")
        sc = attr(X, "scaled:scale")
        if (any(sc == 0))
            stop("cannot rescale a constant/zero column to unit variance.")
        
        
        #--initialization--#
        X <- as.matrix(X)

        if (isFALSE(center) && any(colMeans(X, na.rm = TRUE) != 0)) {
            warning("data contain missing values which will be ",
                    "set to zero for calculations. Consider using center = TRUE ",
                    "for better performance, or impute missing values using ",
                    "'impute.nipals' function.")
        }
        is.na.X <- is.na(X)
        X[is.na.X] <- 0
        X.temp <- X
        n=nrow(X)
        p=ncol(X)
        
        # put a names on the rows and columns
        X.names = dimnames(X)[[2]]
        if (is.null(X.names)) X.names = paste("X", 1:p, sep = "")
        colnames(X) = X.names
        
        ind.names = dimnames(X)[[1]]
        if (is.null(ind.names)) ind.names = 1:nrow(X)
        rownames(X) = ind.names
        
        vect.varX=vector(length=ncomp)
        names(vect.varX) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
        
        vect.iter=vector(length=ncomp)
        names(vect.iter) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
        
        vect.keepX = keepX
        names(vect.keepX) = paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
        
        # KA: to add if biplot function (but to be fixed!)
        #sdev = vector(length = ncomp)
        
        mat.u=matrix(nrow=n, ncol=ncomp)
        mat.v=matrix(nrow=p, ncol=ncomp)
        colnames(mat.u)=paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
        colnames(mat.v)=paste("PC", 1:ncomp, sep = "")#c(1:ncomp)
        rownames(mat.v)=colnames(X)
        
        #--loop on h--#
        for(h in 1:ncomp){
            
            #--computing the SVD--#
            svd.X = svd(X.temp, nu = 1, nv = 1)
            u = svd.X$u[,1]
            loadings = svd.X$v[,1]#svd.X$d[1]*svd.X$v[,1]
            v.stab = FALSE
            u.stab = FALSE
            iter = 0
            
            #--computing nx(degree of sparsity)--#
            nx = p-keepX[h]
            #vect.keepX[h] = keepX[h]
            
            u.old = u
            loadings.old = loadings
            #--iterations on v and u--#
            repeat{
                
                iter = iter +1
                
                loadings = t(X.temp) %*% u
                
                #--penalisation on loading vectors--#
                if (nx != 0) {
                    absa = abs(loadings)
                    if(any(rank(absa, ties.method = "max") <= nx)) {
                        loadings = ifelse(rank(absa, ties.method = "max") <= nx, 0, sign(loadings) * (absa - max(absa[rank(absa, ties.method = "max") <= nx])))
                    }
                }
                loadings = loadings / drop(sqrt(crossprod(loadings)))
                
                u = as.vector(X.temp %*% loadings)
                #u = u/sqrt(drop(crossprod(u))) # no normalisation on purpose: to get the same $x as in pca when no keepX.
                
                #--checking convergence--#
                if(crossprod(u-u.old)<tol){break}
                if(crossprod(loadings-loadings.old)<tol){break}
                
                if (iter >= max.iter)
                {
                    warning(paste("Maximum number of iterations reached for the component", h),call. = FALSE)
                    break
                }
                
                u.old = u
                loadings.old = loadings
                
            }
            
            
            
            #v.final = v.new/sqrt(drop(crossprod(v.new)))
            
            #--deflation of data--#
            c = crossprod(X.temp, u) / drop(crossprod(u))
            X.temp= X.temp - u %*% t(c)  #svd.X$d[1] * svd.X$u[,1] %*% t(svd.X$v[,1])
            
            
            vect.iter[h] = iter
            mat.v[,h] = loadings
            mat.u[,h] = u
            
            #--calculating the variance explained by projecting the data onto the space spanned by the h loading vectors--#
            # this is in case the loading vectors are not orthogonal (see Shen and Huang 2008 Section 2.3)
            X.var = X %*% mat.v[,1:h]%*%solve(t(mat.v[,1:h])%*%mat.v[,1:h])%*%t(mat.v[,1:h])
            vect.varX[h] = sum(X.var^2)
            
            # KA: to add if biplot function (but to be fixed!)
            #sdev[h] = sqrt(svd.X$d[1])
            
            
        }#fin h
        
        rownames(mat.u) = ind.names
        
        cl = match.call()
        cl[[1]] = as.name('spca')
        
        var.tot <- sum(X^2)
        cum.var <- vect.varX/var.tot
        ## calculate per-component explained variance from cum.var,
        # the variance is adjusted to account for potential correlation between PCs:                
        # prop_expl_var <- c(cum.var[1], diff(cum.var))
        prop_expl_var <- explained_variance(X, mat.u, ncomp)
        names(prop_expl_var) <- paste0('PC', seq_len(ncomp))
        
        ## return missing values for output
        X[is.na.X] <- NA_real_
        result = (list(call = cl, X = X,
                       ncomp = ncomp,	
                       #sdev = sdev,  # KA: to add if biplot function (but to be fixed!)
                       #center = center, # KA: to add if biplot function (but to be fixed!)
                       #scale = scale,   # KA: to add if biplot function (but to be fixed!)
                       prop_expl_var = list(X=prop_expl_var),
                       cum.var = cum.var,
                       keepX = vect.keepX,
                       iter = vect.iter,
                       rotation = mat.v,
                       x = mat.u,
                       names = list(X = X.names, sample = ind.names),
                       loadings = list(X = mat.v),
                       variates = list(X = mat.u)
        ))
        ## warning if components are not orthogonal
        .eval_non.orthogonality(variates = result$variates$X, scale=scale)
        
        if (verbose.call) {
            c <- result$call
            result$call <- mget(names(formals()))
            result$call <- append(c, result$call)
            names(result$call)[1] <- "simple.call"
        }
        
        class(result) = c("spca","pca", "prcomp")
        return(invisible(result))
    }

#' Evaluate the non-orthogonality of spca components
#'
#' @param variates variates from spca model
#' @param scale Logical
#' @param warn_threshold Threshold for soft warning
#' @param code_red_threshold Threshold for fail warning
#'
#' @return NULL, or a condition
#' @noRd
#'
#' @examples
#' \dontrun{
#' set.seed(2718)
#' spca.res <- spca(matrix(rnorm(1000), nrow = 20), ncomp = 5, keepX = rep(40, 5), scale = FALSE)
#' mixOmics:::.eval_non.orthogonality(spca.res$variates$X, warn_threshold = 1e-23, scale = FALSE)
#' }

.eval_non.orthogonality <- function(variates, 
                                    scale, 
                                    warn_threshold = 0.01, 
                                    code_red_threshold = 0.1)
{
    
    covmat <- crossprod(variates)
    ## get total variance
    trace <- sum(diag(covmat))
    ## normalise varcov matrix by trace so
    ## the trace sums to 1
    covmat.norm <- covmat/trace
    ## calculate the sum of off-diag (covariances) 
    non_orthog <- sum(abs(covmat.norm - diag(diag(covmat.norm))))
    ## if it's more than 'threshold', throw warning. If it's way too high
    ## warn more strongly (failure). 
    ## make scale=TRUE recommendations if that is scale=FALSE used
    
    recom_scale_msg <- ifelse(
        non_orthog >= warn_threshold & isFALSE(scale),
        "Consider using 'scale = TRUE' ",
        " "
    )
    
    if (non_orthog >= code_red_threshold)
    {
        warning("The sPCA algorithm failed to extract orthogonal components. ", 
                recom_scale_msg)
    }
    else if (non_orthog >= warn_threshold & non_orthog < code_red_threshold)
    {
        warning("The extracted components are not completely orthogonal. ", 
                "Regard the results with care. ",
                recom_scale_msg, "\n")
    }
    NULL
}

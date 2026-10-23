
#' Tune the number of principal components in PCA
#' 
#' \code{tune.pca} can be used to quickly visualise the proportion of explained
#' variance for a large number of principal components in PCA.
#' 
#' The calculation is done either by a singular value decomposition of the
#' (possibly centered and scaled) data matrix, if the data is complete or by
#' using the NIPALS algorithm if there is data missing. Unlike
#' \code{\link{princomp}}, the print method for these objects prints the
#' results in a nice format and the \code{plot} method produces a bar plot of
#' the percentage of variance explaned by the principal components (PCs).
#' 
#' When using NIPALS (missing values), we make the assumption that the first
#' (\code{min(ncol(X),} \code{nrow(X)}) principal components will account for
#' 100 \% of the explained variance.
#' 
#' Note that \code{scale= TRUE} cannot be used if there are zero or constant
#' (for \code{center = TRUE}) variables.
#' 
#' Components are omitted if their standard deviations are less than or equal
#' to \code{comp.tol} times the standard deviation of the first component. With
#' the default null setting, no components are omitted. Other settings for
#' \code{comp.tol} could be \code{comp.tol = sqrt(.Machine$double.eps)}, which
#' would omit essentially constant components, or \code{comp.tol = 0}.
#' 
#' logratio transform and multilevel analysis are performed sequentially as
#' internal pre-processing step, through \code{\link{logratio.transfo}} and
#' \code{\link{withinVariation}} respectively.
#' 
#' @param ncomp integer, the number of components to initially analyse in
#' \code{tune.pca} to choose a final \code{ncomp} for \code{pca}. If
#' \code{NULL}, function sets \code{ncomp = min(nrow(X), ncol(X))}
#' @inheritParams tune
#' @param scale (Default=TRUE) Logical indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. Alternatively,
#' a vector of length equal the number of columns of \code{X} can be supplied.
#' The value is passed to \code{\link{scale}}.
#' @param logratio one of ('none','CLR','ILR'). Default to 'none'
#' @param V Matrix used in the logratio transformation id provided.
#' @param multilevel Design matrix for multilevel analysis (for repeated
#' measurements).
#' @return \code{tune.pca} returns a list with class \code{"tune.pca"}
#' containing the following components: \item{sdev}{the square root of the
#' eigenvalues of the covariance/correlation matrix, though the calculation is
#' actually done with the singular values of the data matrix).}
#' \item{prop_expl_var}{The proportion of explained variance accounted for
#' by each principal component.}
#' \item{cum.var}{the cumulative proportion of explained variance accounted for
#' by the sequential accumulation of principal components is calculated using
#' the sum of the proportion of explained variance}
#' @author Ignacio González, Leigh Coonan, Kim-Anh Le Cao, Fangzhou Yao, 
#' Florian Rohart, Al J Abadi
#' @seealso \code{\link{nipals}}, \code{\link{biplot}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}} and http://www.mixOmics.org
#' for more details.
#' @keywords algebra
#' @export
#' @example ./examples/tune.pca-examples.R

tune.pca <-
    function(X,
             ncomp = NULL,
             center = TRUE, 	# sets the mean of the data to zero, ensures that the first PC describes the direction of the maximum variance
             scale = TRUE, 	# variance is unit across different units
             max.iter = 100,
             tol = 1e-09,
             logratio = c('none','CLR','ILR'),
             V = NULL,
             multilevel = NULL)
    {
        
        logratio <- match.arg(logratio)
        result = pca(X = X, ncomp = ncomp,
                     center = center, scale = scale,
                     max.iter = max.iter, tol = tol,
                     logratio = logratio, V = V,
                     multilevel = multilevel)
        
        is.na.X = is.na(X)
        na.X = FALSE
        if (any(is.na.X)) na.X = TRUE
        
        #  list eigenvalues, prop. of explained variance and cumulative proportion of explained variance
        prop.var = result$prop_expl_var$X
        cum.var = result$cum.var
        
        ind.show = min(10, ncomp)
        
        
        # Plot the principal components and explained variance
        # note: if NA values, we have an estimation of the variance using NIPALS
        if(!na.X)
        {
            ylab = "Proportion of Explained Variance"
        } else{
            ylab = "Estimated Proportion of Explained Variance"
        }
 
        result$call = match.call()
        
        class(result) = c("tune.pca", "pca")
        return(invisible(result))
    }

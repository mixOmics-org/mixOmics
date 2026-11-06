#' Principal Components Analysis
#' 
#' Performs a principal components analysis on the given data matrix that can
#' contain missing values. If data are complete 'pca' uses Singular Value
#' Decomposition, if there are some missing values, it uses the NIPALS
#' algorithm.
#' 
#' The calculation is done either by a singular value decomposition of the
#' (possibly centered and scaled) data matrix, if the data is complete or by
#' using the NIPALS algorithm if there is data missing. Unlike
#' \code{\link{princomp}}, the print method for these objects prints the
#' results in a nice format and the \code{plot} method produces a bar plot of
#' the percentage of variance explained by the principal components (PCs).
#' 
#' When using NIPALS (missing values), we make the assumption that the first
#' (\code{min(ncol(X),} \code{nrow(X)}) principal components will account for
#' 100 \% of the explained variance.
#' 
#' Note that \code{scale = TRUE} will throw an error if there are constant
#' variables in the data, in which case it's best to filter these variables
#' in advance.
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
#' @param X a numeric matrix (or data frame) which provides the data for the
#'   principal components analysis. It can contain missing values in which case
#'   \code{center = TRUE} is used as required by the
#'   \code{\link{nipals}} function.
#' @param ncomp Integer, if data is complete \code{ncomp} decides the number of
#' components and associated eigenvalues to display from the \code{pcasvd}
#' algorithm and if the data has missing values, \code{ncomp} gives the number
#' of components to keep to perform the reconstitution of the data using the
#' NIPALS algorithm. If \code{NULL}, function sets \code{ncomp = min(nrow(X),
#' ncol(X))}
#' @param center (Default=TRUE) Logical, whether the variables should be shifted
#'   to be zero centered. Only set to FALSE if data have already been centered.
#'   Alternatively, a vector of length equal the number of columns of \code{X}
#'   can be supplied. The value is passed to \code{\link{scale}}. If the data
#'   contain missing values, columns should be centered for reliable results.
#' @param scale (Default=FALSE) Logical indicating whether the variables should be
#' scaled to have unit variance before the analysis takes place. The default is
#' \code{FALSE} for consistency with \code{prcomp} function, but in general
#' scaling is advisable. Alternatively, a vector of length equal the number of
#' columns of \code{X} can be supplied. The value is passed to
#' \code{\link{scale}}.
#' @param max.iter Integer, the maximum number of iterations in the NIPALS
#' algorithm.
#' @param tol Positive real, the tolerance used in the NIPALS algorithm.
#' @param logratio (Default='none') one of ('none','CLR','ILR'). Specifies the log ratio
#' transformation to deal with compositional values that may arise from
#' specific normalisation in sequencing data. Default to 'none'
#' @param ilr.offset (Default=0.001) When logratio is set to 'ILR', an offset must be input to
#' avoid infinite value after the logratio transform.
#' @param V Matrix used in the logratio transformation if provided.
#' @param multilevel sample information for multilevel decomposition for
#' repeated measurements.
#' @template arg/verbose.call
#' @return \code{pca} returns a list with class \code{"pca"} and \code{"prcomp"}
#' containing the following components: 
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' \item{X}{The input data matrix, possibly scaled and centered.}
#' \item{ncomp}{The number of principal components used.}
#' \item{center}{The centering used.}
#' \item{scale}{The scaling used.}
#' \item{names}{List of row and column names of data.}
#' \item{sdev}{The eigenvalues of the covariance/correlation matrix, though
#' the calculation is actually done with the singular values of the data
#' matrix or by using NIPALS.}
#' \item{loadings}{A length one list of matrix of variable loadings for X (i.e.,
#' a matrix whose columns contain the eigenvectors).} 
#' \item{variates}{Matrix containing the coordinate values corresponding to the
#' projection of the samples in the space spanned by the principal components.
#' These are the dimension-reduced representation of observations/samples.}
#' \item{var.tot}{Total variance in the data.}
#' \item{prop_expl_var}{Proportion of variance explained per
#' component after setting possible missing values in the data to zero (note
#' that contrary to PCA, this amount may not decrease as the aim of the method
#' is not to maximise the variance, but the covariance between X and the
#' dummy matrix Y).}
#' \item{cum.var}{The cumulative explained variance for components.}
#' \item{Xw}{If multilevel, the data matrix with within-group-variation removed.}
#' \item{design}{If multilevel, the provided design.}
#' @author Florian Rohart, Kim-Anh Lê Cao, Ignacio González, Al J Abadi
#' @seealso \code{\link{nipals}}, \code{\link{prcomp}}, \code{\link{biplot}},
#' \code{\link{plotIndiv}}, \code{\link{plotVar}} and http://www.mixOmics.org
#' for more details.
#' @references On log ratio transformations: Filzmoser, P., Hron, K., Reimann,
#' C.: Principal component analysis for compositional data with outliers.
#' Environmetrics 20(6), 621-632 (2009) Lê Cao K.-A., Costello ME, Lakis VA,
#' Bartolo, F,Chua XY, Brazeilles R, Rondeau P. MixMC: Multivariate insights
#' into Microbial Communities. PLoS ONE, 11(8): e0160169 (2016). On multilevel
#' decomposition: Westerhuis, J.A., van Velzen, E.J., Hoefsloot, H.C., Smilde,
#' A.K.: Multivariate paired data analysis: multilevel plsda versus oplsda.
#' Metabolomics 6(1), 119-128 (2010) Liquet, B., Lê Cao, K.-A., Hocini, H.,
#' Thiebaut, R.: A novel approach for biomarker selection and the integration
#' of repeated measures experiments from two assays. BMC bioinformatics 13(1),
#' 325 (2012)
#' @keywords algebra
#' @export
#' @example ./examples/pca-examples.R
pca <- function(X,
                ncomp = 2,
                center = TRUE,
                scale = FALSE,
                max.iter = 500,
                tol = 1e-09,
                logratio = c('none','CLR','ILR'),
                ilr.offset = 0.001,
                V = NULL,
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
    X <- .check_numeric_matrix(X)
    
    #-- put a names on the rows and columns of X --#
    X.names = colnames(X)
    if (is.null(X.names))
        X.names = paste("V", 1:ncol(X), sep = "")
    
    ind.names = rownames(X)
    if (is.null(ind.names))
        ind.names = 1:nrow(X)
    
    #-- ncomp
    if (is.null(ncomp))
        ncomp = min(nrow(X), ncol(X))
    
    ncomp = round(ncomp)
    
    if (!is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
        stop("invalid value for 'ncomp'.", call. = FALSE)
    
    if (ncomp > min(ncol(X), nrow(X)))
        stop("use smaller 'ncomp'", call. = FALSE)
    
    #-- log.ratio
    logratio_choice = c('CLR', 'ILR', 'none')
    logratio = logratio_choice[pmatch(logratio, logratio_choice)]
    logratio <- match.arg(logratio)

    if (logratio != "none" && any(X < 0))
        stop("'X' contains negative values, you can not log-transform your data"
        )
    
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
    
    if (is.null(V) & logratio == "ILR") # back-transformation to clr-space, will be used later to recalculate loadings etc
        V = clr.backtransfo(X)
    
    X = logratio.transfo(X = X, logratio = logratio, offset = if(logratio == "ILR") {ilr.offset} else {0})
    
    #as X may have changed
    if (ncomp > min(ncol(X), nrow(X)))
        stop("use smaller 'ncomp'", call. = FALSE)
    
    #-- logratio transformation --#
    #-----------------------------#
    
    #---------------------------------------------------------------------------#
    #-- Start: multilevel approach ----------------------------------------------------#
    
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
    #-- End: multilevel approach ----------------------------------------------------#
    
    .check_zero_var_columns(X, scale = scale)
    X <- scale(X, center = center, scale = scale)
    sc <- attr(X, 'scaled:scale')
    cen <- attr(X, 'scaled:center')
    result <- list(call = match.call(), 
                   X = X, 
                   ncomp = ncomp,
                   scale = if (is.null(sc)) FALSE else sc,
                   center = if (is.null(cen)) FALSE else cen,
                   names = list(X = X.names, sample = ind.names))
    
    #-- pca approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    variates <- NULL ## only changed in ILR case, otherwise calculated using X and loadings
    if (any(is.na(X))) {
        res <- 
            nipals(X, ncomp = ncomp,
                   max.iter = max.iter, tol = tol)
        sdev <- res$eig
        loadings <- res$p
    } else {
        if (logratio %in% c('CLR', 'none')) {
            #-- if data is complete use singular value decomposition
            #-- borrowed from 'prcomp' function (R stats package, GPL-2 | GPL-3)
            res = svd(X, nu = 0)
            sdev <- res$d
            loadings <- res$v
        } else {
            # if 'ILR', transform data and then back transform in clr space (from RobCompositions package)
            # data have been transformed above
            res = svd(X, nu = max(1, nrow(X) - 1))
            sdev = res$d[1:ncomp] / sqrt(max(1, nrow(X) - 1))  # Note: what differs with RobCompo is that they use: cumsum(eigen(cov(X))$values)/sum(eigen(cov(X))$values)
            # calculate loadings using back transformation to clr-space
            loadings = V %*% res$v[, 1:ncomp, drop = FALSE]
            # extract component score from the svd, multiply matrix by vector using diag, NB: this differ from our mixOmics PCA calculations
            # NB: this differ also from Filmoser paper, but ok from their code: scores are unchanged
            variates = res$u[, 1:ncomp, drop = FALSE] %*% diag(res$d[1:ncomp, drop = FALSE])
        }
        
    }
    rm(res)
    
    result <-
        .add_sdev_loadings_and_variates(
            result = result,
            sdev = sdev,
            loadings = loadings,
            variates = variates
        )
    ## add explained/cum/total variance
    result <- c(result, .get_var_stats(X = result$X, sdev = result$sdev))
    expected_output_names <- c("call", "X", "ncomp", "center", "scale", "names", 
                         "sdev", "loadings", "variates", "prop_expl_var", "var.tot",
                         "cum.var", "rotation", "x")
    
    if (names(result) %!=% expected_output_names)
    {
        stop("Unexpected error. Please submit an issue at\n",
             "https://github.com/mixOmicsTeam/mixOmics/issues/new/choose", call. = FALSE)
    }
    result <- result[expected_output_names]
    # output multilevel if needed
    if(!is.null(multilevel)) # TODO include in docs returns
        result=c(result, list(Xw = Xw, design = multilevel))
    
    
    
    if (verbose.call) {
        c <- result$call
        result$call <- mget(names(formals()))
        result$call <- append(c, result$call)
        names(result$call)[1] <- "simple.call"
    }
    
    class(result) = c("pca") 
    if(!is.null(multilevel))
        class(result)=c("mlpca",class(result))
    
    return(invisible(result))
}

#' Add sdev, loadings and components to result
#'
#' @param result Result list. $sdev, $loadings and $variates will be added to it.
#' @param X A numeric matrix.
#' @param ncomp Integer, number of components.
#' @param sdev numeric vector - sqrt(eig) from different analyses
#' @param loadings Numeric matrix of loadings from different analyses
#' @param variates Numeric matrix of variates from ILR case, or NULL.
#' @param row_names Vector of sample names
#' @param col_names Vector of feature names
#' @noRd
#' @return A modified list
.add_sdev_loadings_and_variates <-
    function(result,
             sdev,
             loadings,
             variates = NULL) {
        
        ncomp <- result$ncomp
        pc_names <- paste("PC", seq_len(ncomp), sep = "")
        
        X <- result$X
        X[is.na(X)] <- 0 ## if any
        
        sdev = sdev[seq_len(ncomp)] / sqrt(max(1, nrow(X) - 1))
        loadings = loadings[, seq_len(ncomp), drop = FALSE]
        if (is.null(variates))
            variates = X %*% loadings
        
        names(sdev) = pc_names
        dimnames(loadings) = list(colnames(X), pc_names)
        dimnames(variates) = list(rownames(X), pc_names)
        
        result[c('sdev', 'loadings', 'variates', 'x', 'rotation')] <- list(sdev, 
                                                                      list(X=loadings), 
                                                                      list(X=variates),
                                                                      variates,
                                                                      loadings)
        
        result
    }


#' Add cumulative and per-component explained variance 
#'
#' @param X Numeric matrix
#' @param sdev Numeric sd vector
#' 
#'
#' @return A list including total variance, loadings and
#' variates, as well as per-component and
#' cumulative explained variance.
#' @noRd
.get_var_stats <- function(X, sdev) {
    var.tot=sum(X^2, na.rm = TRUE) / max(1, nrow(X) - 1)
    
    # calculate explained variance
    prop_expl_var <- sdev^2 / var.tot
    cum.var = cumsum(prop_expl_var)
    prop_expl_var = list(X = prop_expl_var)
    list(prop_expl_var = prop_expl_var,
         var.tot = var.tot, 
         cum.var = cum.var)
}


#' Regularized Canonical Correlation Analysis
#' 
#' The function performs the regularized extension of the Canonical Correlation
#' Analysis to seek correlations between two data matrices.
#' 
#' The main purpose of Canonical Correlations Analysis (CCA) is the exploration
#' of sample correlations between two sets of variables \eqn{X} and \eqn{Y}
#' observed on the same individuals (experimental units) whose roles in the
#' analysis are strictly symmetric.
#' 
#' The \code{cancor} function performs the core of computations but additional
#' tools are required to deal with data sets highly correlated (nearly
#' collinear), data sets with more variables than units by example.
#' 
#' The \code{rcc} function, the regularized version of CCA, is one way to deal
#' with this problem by including a regularization step in the computations of
#' CCA. Such a regularization in this context was first proposed by Vinod
#' (1976), then developped by Leurgans \emph{et al.} (1993). It consists in the
#' regularization of the empirical covariances matrices of \eqn{X} and \eqn{Y}
#' by adding a multiple of the matrix identity, that is, Cov\eqn{(X)+ \lambda_1
#' I} and Cov\eqn{(Y)+ \lambda_2 I}.
#' 
#' When \code{lambda1=0} and \code{lambda2=0}, \code{rcc} performs a classical
#' CCA, if possible (i.e. when \eqn{n > p+q}.
#'
#' Unlike classical CCA, ridge regularisation is not invariant to rescaling
#' individual variables. The ridge parameters are on the variance scale of
#' the input data. With \code{scale = TRUE}, both matrices are centred and
#' scaled to unit variance before the covariance matrices are calculated.
#' This makes the penalties comparable across variables with different units.
#' The default grid in \code{\link{tune.rcc}} assumes variables of roughly
#' unit variance. Ridge fits are standardised by default with
#' \code{scale = TRUE}. Set \code{scale = FALSE} to retain historical unscaled
#' ridge fits. Constant columns and columns with fewer than two non-missing
#' values must be removed before using \code{scale = TRUE}.
#'
#' The \code{scale} argument precedes \code{verbose.call}. Calls that previously
#' passed \code{verbose.call} as the seventh positional argument must now name
#' that argument explicitly.
#'
#' The shrinkage estimates \code{method = "shrinkage"} can be used to bypass
#' \code{\link{tune.rcc}} to choose the shrinkage parameters - which can be
#' long and costly to compute with very large data sets. Note that both
#' functions \code{\link{tune.rcc}} (which uses cross-validation) and the
#' shrinkage parameters (which uses the formula from Schafer and Strimmer, see the corpcor package \code{\link{estimate.lambda}} ) may
#' output different results.
#' 
#' Note: when \code{method = "shrinkage"} the parameters are estimated using \code{\link{estimate.lambda}} 
#' from the corpcor package. Data are then centered to calculate
#' the regularised variance-covariance matrices in \code{\link{rcc}}.
#' 
#' Missing values are handled in the function, except when using \code{method = "shrinkage"}.
#' In that case the estimation of the missing values can be performed by the reconstitution
#' of the data matrix using the \code{nipals} function. 
#' 
#' @aliases rcc rcc.default
#' @param X numeric matrix or data frame \eqn{(n \times p)}, the observations
#' on the \eqn{X} variables. \code{NA}s are allowed.
#' @param Y numeric matrix or data frame \eqn{(n \times q)}, the observations
#' on the \eqn{Y} variables. \code{NA}s are allowed.
#' @param method One of "ridge" or "shrinkage". If "ridge", \code{lambda1} and
#' \code{lambda2} need to be supplied (see also our function tune.rcc); if
#' "shrinkage", parameters are directly estimated with Strimmer's formula, see
#' below and reference.
#' @param ncomp the number of components to include in the model. Default to 2.
#' @param lambda1,lambda2 a non-negative real. The regularization parameter for
#' the \emph{X} and \emph{Y} data. Defaults to \code{lambda1=lambda2=0}. Only
#' used if \code{method="ridge"}
#' @param scale logical. If \code{TRUE}, standardise both \code{X} and
#' \code{Y} to unit variance before fitting with \code{method = "ridge"}.
#' Defaults to \code{TRUE}. This argument is ignored for
#' \code{method = "shrinkage"}, which retains its own covariance estimation.
#' @template arg/verbose.call
#' @return \code{rcc} returns a object of class \code{"rcc"}, a list that
#' contains the following components: \item{X}{the \eqn{X} data, centred and
#' scaled when \code{method = "ridge"} and \code{scale = TRUE}; otherwise
#' the original data.}
#' \item{Y}{the \eqn{Y} data, transformed in the same way as \code{X}.}
#' \item{center}{list of training column means for \code{X} and \code{Y}.}
#' \item{scale}{list of training standard deviations for \code{X} and
#' \code{Y}, or \code{FALSE} for each matrix if unit-variance scaling was not
#' applied. Use these with \code{center} to transform new samples before
#' multiplying by \code{loadings}.}
#' \item{cor}{a vector containing the
#' canonical correlations.} \item{lambda}{a vector containing the
#' regularization parameters whether those were input if ridge method or
#' directly estimated with the shrinkage method.} \item{loadings}{list
#' containing the estimated coefficients used to calculate the canonical
#' variates in \eqn{X} and \eqn{Y}.} \item{variates}{list containing the
#' canonical variates.} \item{names}{list containing the names to be used for
#' individuals and variables.}
#' \item{prop_expl_var}{Proportion of the explained variance of derived
#' components, after setting possible missing values to zero.}
#' \item{call}{if \code{verbose.call = FALSE}, then just the function call is returned.
#' If \code{verbose.call = TRUE} then all the inputted values are accessable via
#' this component}
#' @author Sébastien Déjean, Ignacio González, Francois Bartolo, Kim-Anh Lê Cao, 
#' Florian Rohart, Al J Abadi
#' @seealso \code{\link{summary}}, \code{\link{tune.rcc}},
#' \code{\link{plot.rcc}}, \code{\link{plotIndiv}}, \code{\link{plotVar}},
#' \code{\link{cim}}, \code{\link{network}} and http://www.mixOmics.org for
#' more details.
#' @references González, I., Déjean, S., Martin, P. G., and Baccini, A. (2008).
#' CCA: An R package to extend canonical correlation analysis. Journal of
#' Statistical Software, 23(12), 1-14.
#' 
#' González, I., Déjean, S., Martin, P., Goncalves, O., Besse, P., and Baccini,
#' A. (2009). Highlighting relationships between heterogeneous biological data
#' through graphical displays based on regularized canonical correlation
#' analysis. Journal of Biological Systems, 17(02), 173-199.
#' 
#' Leurgans, S. E., Moyeed, R. A. and Silverman, B. W. (1993). Canonical
#' correlation analysis when the data are curves. \emph{Journal of the Royal
#' Statistical Society. Series B} \bold{55}, 725-740.
#' 
#' Vinod, H. D. (1976). Canonical ridge and econometrics of joint production.
#' \emph{Journal of Econometrics} \bold{6}, 129-137.
#' 
#' Opgen-Rhein, R., and K. Strimmer. 2007. Accurate ranking of differentially
#' expressed genes by a distribution-free shrinkage approach. Statist.
#' emphAppl. Genet. Mol. Biol. \bold{6}:9.
#' (http://www.bepress.com/sagmb/vol6/iss1/art9/)
#' 
#' Sch"afer, J., and K. Strimmer. 2005. A shrinkage approach to large-scale
#' covariance estimation and implications for functional genomics. Statist.
#' emphAppl. Genet. Mol. Biol.  \bold{4}:32.
#' (http://www.bepress.com/sagmb/vol4/iss1/art32/)
#' @keywords multivariate
#' @export
#' @examples
#' 
#' ## Classic CCA
#' data(linnerud)
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.res <- rcc(X, Y)
#' 
#' \dontrun{
#' ## Regularized CCA
#' data(nutrimouse)
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' nutri.res1 <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
#' 
#' ## using shrinkage parameters
#' nutri.res2 <- rcc(X, Y, ncomp = 3, method = 'shrinkage')
#' nutri.res2$lambda # the shrinkage parameters
#' }
rcc <-
  function(X,
           Y,
           ncomp = 2,
           method = c("ridge", "shrinkage"),
           lambda1 = 0,
           lambda2 = 0,
           scale = TRUE,
           verbose.call = FALSE
  )
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
    
    #-- data set names --#
    data.names = c(deparse(substitute(X)), deparse(substitute(Y)))
    
    #-- method
    method = match.arg(method)

    if (!is.logical(scale) || length(scale) != 1L || is.na(scale))
      stop("'scale' must be either TRUE or FALSE.", call. = FALSE)
    
    #-- X matrix
    if (is.data.frame(X)) X = as.matrix(X)
    
    if (!is.matrix(X) || is.character(X))
      stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite)))
      stop("infinite values in 'X'.", call. = FALSE)
    
    if (method == "shrinkage"){
      if (any(is.na(X)))
        stop("missing values in 'X' matrix. NAs not are allowed if method = 'shrinkage'.", call. = FALSE)
    }
    #-- Y matrix
    if (is.data.frame(Y)) Y = as.matrix(Y)
    
    if (!is.matrix(Y) || is.character(Y))
      stop("'Y' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(Y, 1, is.infinite)))
      stop("infinite values in 'Y'.", call. = FALSE)
    
    if (method == "shrinkage")
      if (any(is.na(Y)))
        stop("missing values in 'Y' matrix. NAs not are allowed if method = 'shrinkage'.", call. = FALSE)
    
    #-- equal number of rows in X and Y
    if ((n = nrow(X)) != nrow(Y))
      stop("unequal number of rows in 'X' and 'Y'.", call. = FALSE)
    
    p = ncol(X)
    q = ncol(Y)
    
    #-- put a names on the columns of X and Y --#
    X.names = colnames(X)
    if (is.null(X.names)) X.names = paste("X", 1:ncol(X), sep = "")
    
    Y.names = colnames(Y)
    if (is.null(Y.names)) Y.names = paste("Y", 1:ncol(Y), sep = "")
    
    #-- put a names on the samples --#
    ind.names = dimnames(X)[[1]]
    if (is.null(ind.names)) ind.names = dimnames(Y)[[1]]
    if (is.null(ind.names)) ind.names = 1:n
    
    #-- ncomp
    if (is.null(ncomp) || ncomp < 1 || !is.finite(ncomp))
      stop("invalid value for 'ncomp'.", call. = FALSE)
    
    ncomp = round(ncomp)
    
    if (ncomp > min(p, q))
      stop("'comp' must be <= ", min(p, q), ".",
           call. = FALSE)
    
    #-- lambda1
    if (!is.finite(lambda1) || is.null(lambda1))
      stop("invalid value for 'lambda1'.", call. = FALSE)
    
    if(lambda1 < 0)
      stop("'lambda1' must be a non-negative value.", call. = FALSE)
    
    #-- lambda2
    if (!is.finite(lambda2) || is.null(lambda2))
      stop("invalid value for 'lambda2'.", call. = FALSE)
    
    if(lambda2 < 0)
      stop("'lambda2' must be a non-negative value.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    
    #-- rcc approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- covariance matrices regularization --#
    if (method == "ridge") {
      if (scale) {
        if (any(colSums(!is.na(X)) < 2L))
          stop("each column of 'X' must have at least two non-missing values to scale.", call. = FALSE)
        if (any(colSums(!is.na(Y)) < 2L))
          stop("each column of 'Y' must have at least two non-missing values to scale.", call. = FALSE)
        .check_zero_var_columns(X, block_name = "X")
        .check_zero_var_columns(Y, block_name = "Y")
        X = base::scale(X, center = TRUE, scale = TRUE)
        Y = base::scale(Y, center = TRUE, scale = TRUE)
      }
      Cxx = var(X, na.rm = TRUE, use = "pairwise") + diag(lambda1, ncol(X))
      Cyy = var(Y, na.rm = TRUE, use = "pairwise") + diag(lambda2, ncol(Y))
      Cxy = cov(X, Y, use = "pairwise")
    }
    else { # if method == 'shrinkage'
      Cxx = cov.shrink(X, verbose = FALSE)
      Cyy = cov.shrink(Y, verbose = FALSE)
      
      lambda.x = attr(Cxx, "lambda")
      lambda.y = attr(Cyy, "lambda")
      
      sc.x = sqrt(var.shrink(X, verbose = FALSE))
      sc.y = sqrt(var.shrink(Y, verbose = FALSE))
      
      w = rep(1/n, n)
      xs = wt.scale(X, w, center = TRUE, scale = TRUE)
      ys = wt.scale(Y, w, center = TRUE, scale = TRUE)
      
      #-- bias correction factor
      h1 = n / (n - 1)
      
      #-- unbiased empirical estimator
      Cxy = h1 * crossprod(sweep(sweep(xs, 1, sqrt((1 - lambda.x) * w), "*"), 2, sc.x, "*"),
                           sweep(sweep(ys, 1, sqrt((1 - lambda.y) * w), "*"), 2, sc.y, "*"))
    }
    
    #-- calculation of the canonical correlations and canonical variables --#
    Cxx.fac = chol(Cxx)
    Cyy.fac = chol(Cyy)
    Cxx.fac.inv = solve(Cxx.fac)
    Cyy.fac.inv = solve(Cyy.fac)
    mat = t(Cxx.fac.inv) %*% Cxy %*% Cyy.fac.inv
    
    if (p >= q) {
      result = svd(mat, nu = ncomp, nv = ncomp)
      cor = result$d
      xcoef = Cxx.fac.inv %*% result$u
      ycoef = Cyy.fac.inv %*% result$v
    }
    else {
      result = svd(t(mat), nu = ncomp, nv = ncomp)
      cor = result$d
      xcoef = Cxx.fac.inv %*% result$v
      ycoef = Cyy.fac.inv %*% result$u
    }
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    names(cor) = 1:length(cor)
    if (method == "ridge" && scale) {
      X.aux = X
      Y.aux = Y
      scaling = list(X = attr(X, "scaled:scale"), Y = attr(Y, "scaled:scale"))
    } else {
      X.aux = base::scale(X, center = TRUE, scale = FALSE)
      Y.aux = base::scale(Y, center = TRUE, scale = FALSE)
      scaling = list(X = FALSE, Y = FALSE)
    }
    centering = list(X = attr(X.aux, "scaled:center"), Y = attr(Y.aux, "scaled:center"))
    X.aux[is.na(X.aux)] = 0
    Y.aux[is.na(Y.aux)] = 0
    
    U = X.aux %*% xcoef
    V = Y.aux %*% ycoef
    
    cl = match.call()
    cl[[1]] = as.name('rcc')
    
    if (method == "ridge"){
      lambda = c("lambda1" = lambda1, "lambda2" = lambda2)
    } else { # if method == 'shrinkage')
      lambda = c("lambda1" = lambda.x, "lambda2" = lambda.y)
    }
    
    result = list(call = cl,
                  X = X,
                  Y = Y,
                  ncomp = ncomp,
                  method = method,
                  cor = cor,
                  loadings = list(X = xcoef, Y = ycoef),
                  variates = list(X = U, Y = V),
                  names = list(sample = ind.names, colnames = list(X=colnames(X),Y=colnames(Y)), blocks = c("X","Y"),#list(X = X.names, Y = Y.names, indiv = ind.names,
                               data = data.names),
                  lambda = lambda)
    
    #calcul explained variance
    explX=explained_variance(result$X,result$variates$X,ncomp)
    explY=explained_variance(result$Y,result$variates$Y,ncomp)
    result$prop_expl_var=list(X=explX,Y=explY)
    result$center = centering
    result$scale = scaling
    
    if (verbose.call) {
      c <- result$call
      result$call <- mget(names(formals()))
      result$call <- append(c, result$call)
      names(result$call)[1] <- "simple.call"
    }
    
    
    class(result) = "rcc"
    return(invisible(result))
  }

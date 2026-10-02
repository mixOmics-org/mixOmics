## TODO use these utils throughout the package
## ------------------------ .check_numeric_matrix  ------------------------ ##
#' check if x is a valid numeric matrix -- possibly including NAs
#' 
#' Coerces to numeric matrix if necessary and ensures only numeric (including
#' NA) values are present.
#'@noRd
#'@examples
#' .check_numeric_matrix(mtcars)
.check_numeric_matrix <- function(X, block_name = 'X')
{
    err_msg <- sprintf("'%s' must be a numeric matrix (possibly including NA's) with finite values", block_name)
    X <- tryCatch(data.matrix(X, rownames.force = TRUE), 
                  error = function(e) stop(err_msg, call. = FALSE))
    if (!all(is.finite(X) | is.na(X)))
        stop(err_msg, call. = FALSE)
    X
}

## ----------------------------- .check_cv_args ----------------------------- ##
#' Check validation and nrepeat arguments
#'
#' Check that validation arg is one of c('Mfold', 'loo') and nrepeat is
#' valid. i.e., it is (set to) 1 if validation == 'loo'. In this case folds is also set to N The inputs are not
#' case-sensitive but the output will be one of c('Mfold','loo').
#' @param validation validation arg
#' @param nrepeat nrepeat arg
#' @param folds folds arg. All function must have a NULL default.
#' @param N Number of samples/rows
#'
#' @return A list of valid validation and nrepeat args (possibly with a warning), 
#' or a condition
#' @keywords Internal
#' @noRd
#' @examples
#'  dput(.check_cv_args(validation = 'mfold', nrepeat = 2, folds = 3, N = 10))
#'  dput(.check_cv_args(validation = 'LOO', nrepeat = 1, folds = NULL, N = 10))
#' \dontrun{
#'  ## invalid folds
#'  dput(.check_cv_args(validation = 'LOO', nrepeat = 1, folds = 3, N = 10))
#'  ## invalid nrepeat
#'  dput(.check_cv_args(validation = 'LOO', nrepeat = 2, folds = NULL, N = 10))
#' }
.check_cv_args <- function(validation, nrepeat, folds, N)
{
    validation <- match.arg(tolower(validation), choices = tolower(c("Mfold", "loo")))
    
    if (validation == 'loo')
    {
        if (nrepeat != 1)
            warning("Leave-One-Out cross-validation does not need to be repeated. ",
                    "'nrepeat' is set to '1'.", call. = FALSE)
        nrepeat <- 1
        
        if (!missing(folds) && folds != N)
            warning("'folds' is set to ", N, " for Leave-One-Out cross-validation.", call. = FALSE)
        folds <- N
    } else {
        if (missing(folds) || !is.numeric(folds) || (folds < 2 | folds > N))
            stop("'folds' should be an integer > 2 and < ", N, call. = FALSE)
        
        validation = 'Mfold'
    }
    list(validation = validation, nrepeat = as.integer(nrepeat), folds = as.integer(folds))
}

## ------------------------------ .check_folds ----------------------------- ##
#' Check folds supplied as a list
#'
#' Validates a list of test-set indices for M-fold cross-validation, as
#' supplied through the \code{folds} argument. Each element must contain the
#' row indices (integers between 1 and \code{n}) of the samples in one test
#' fold. Together the elements must cover every sample exactly once.
#'
#' @param folds a list of integer vectors, one per fold
#' @param n number of samples/rows
#' @param Y optional factor of class labels. When supplied, an error is raised
#'   if a class is absent from the training set of any fold: the model fitted
#'   on that training set would not know the class, and the downstream
#'   aggregation of the predictions assumes all levels of \code{Y}.
#'
#' @return The validated folds as an unnamed list of integer vectors, or a
#'   condition
#' @keywords Internal
#' @noRd
#' @examples
#'  dput(.check_folds(split(1:10, rep(1:2, 5)), n = 10))
#' \dontrun{
#'  ## repeated sample
#'  .check_folds(list(1:6, 5:10), n = 10)
#'  ## missing sample
#'  .check_folds(list(1:4, 5:9), n = 10)
#'  ## a list of lists
#'  .check_folds(list(list(1:5, 6:10), list(1:5, 6:10)), n = 10)
#' }
.check_folds <- function(folds, n, Y = NULL)
{
    if (!is.list(folds) || any(vapply(folds, is.list, logical(1))))
        stop("'folds' must be a list of integer vectors, one vector per fold.",
             call. = FALSE)

    if (length(folds) < 2 || length(folds) > n)
        stop("Invalid number of folds.", call. = FALSE)

    if (any(lengths(folds) == 0))
        stop("Invalid folds. Each fold must contain at least one sample.", call. = FALSE)

    ind <- unlist(folds, use.names = FALSE)

    if (!is.numeric(ind) || any(!is.finite(ind)) || any(ind != round(ind)) ||
        any(ind < 1) || any(ind > n))
        stop("Invalid folds. Each fold must contain integer indices between 1 and ",
             n, ".", call. = FALSE)

    if (length(ind) != n)
        stop("Invalid folds. The total number of samples in folds must be equal to ",
             n, ".", call. = FALSE)

    if (length(unique(ind)) != n)
        stop("Invalid folds. Repeated samples in folds.", call. = FALSE)

    folds <- lapply(unname(folds), as.integer)

    if (!is.null(Y))
    {
        Y <- as.factor(Y)
        for (i in seq_along(folds))
        {
            absent <- names(which(table(Y[-folds[[i]]]) == 0))
            if (length(absent) > 0)
                stop("Invalid folds. Class '", paste(absent, collapse = "', '"),
                     "' is not represented in the training set of fold ", i,
                     ": each class must be present in the training set of every fold.",
                     call. = FALSE)
        }
    }

    folds
}

## -------------------------- .check_test.keepX --------------------------- ##
#' Check test.keepX
#'
#' Check test.keepX is valid when X is either matrix or list
#' @param test.keepX test.keepX
#' @param X X input from mixOmics tune models
#' @param indY indY
#' @param already.tested.X already.tested.X
#'
#' @return test.keepX, possibly re-ordered by names for list X
#' @noRd
#' @keywords Internal
#' 
.check_test.keepX <- function(test.keepX, 
                              X,
                              indY = NULL, # TODO
                              already.tested.X = NULL # TODO
)
{
    ## -- checker for a pair of test.keepX and X
    .check_test.keepX_helper <- function(test.keepX_, X_)
    {
        if (is.data.frame(X_))
        {
            X_ <- as.matrix(X_)
        }
        
        if (mode(test.keepX_) != 'numeric' || !all(test.keepX_ %in% seq_len(ncol(X_))))
        {
            stop( "'test.keepX' values should be positive whole numbers < ncol(X) = ", 
                  ncol(X_), call. = FALSE)
        }
        invisible(NULL)
    }
    
    ## ------- X matrix ------- ##
    if (!is.null(dim(X)))
    {
        .check_test.keepX_helper(test.keepX, X)
    } 
    ## -------- X list -------- ##
    else
    {
        ## names
        if (! (is.list(test.keepX) && setequal(names(test.keepX), names(X)) ) )
        {
            stop("'test.keepX' must be a named list with names: names(X) = ", names(X))
        }
        else
        {
            test.keepX <- test.keepX[names(X)]
            mapply(z = test.keepX, w = X, FUN = function(z, w) .check_test.keepX_helper(z, w)) 
        }
        
    }
    invisible(test.keepX)
}

## ----------------------------- .check_ncomp ----------------------------- ##

#' Check ncomp
#'
#' Check that ncomp is positive integer <= smallest dim in the data
#' @param ncomp ncomp arg
#' @param X A matrix or a list of matrices with \code{dim} method
#' @param default Integer, default value if \code{ncomp} is \code{NULL}
#'
#' @return Integer, or a condition
#' @keywords Internal
#' @noRd
#' @examples
#' \dontrun{
#' .check_ncomp(300, X = mtcars)
#' #> Error in .check_ncomp(300, X = mtcars) : 
#' #>     'ncomp' must be smaller than or equal to the smallest dimension in X: 11
#' .check_ncomp(NULL, X = mtcars, default = 3)
#' #> 3
#' }
#' 
.check_ncomp <- function(ncomp, X, default = 2)
{
    formals(stop)$call. <- FALSE
    if (!is.null(dim(X))) X <- list(X)
    min.dim <- min(sapply(X, function(z) min(dim(z))))
    err_msg <- paste0("'ncomp' must be a positive integer smaller than ",
                      "or equal to the smallest dimenion in X: ", min.dim, "\n")
    if ( is.null(ncomp) )
        ncomp <- default
    if (mode(ncomp) != 'numeric' || ncomp%%1 != 0 | ncomp < 1 | ncomp > min.dim)
        stop(err_msg)

    return(ncomp)
}

## ---------------------------- .check_logical ---------------------------- ##
#' Check logical arguments
#'
#' @param arg A function argument provided as name
#'
#' @return arg or condition, invisibly
#'
#' @noRd
#' @keywords Internal
#' @examples
#' \dontrun{
#' foo <- function(a=TRUE) .check_logical(a)
#' foo(a = 1)
#' #> Error: 'a' must be a class logical (TRUE or FALSE).
#' }
.check_logical <- function(arg)
{
    if (!is.logical(arg))
    {
        stop(sQuote(deparse(substitute(arg))), 
             " must be a class logical (TRUE or FALSE).", 
             call. = FALSE)
    }
    invisible(arg)
}

## ----------------------------- .check_cpus ------------------------------ ##
#' Check cpus argument
#'
#' @param cpus Integer, the number of cpus for parallel processing.
#'
#' @return The number of cpus available/applicable, or condition.
#' @author Al J Abadi
#'
#' @keywords Internal
#' @noRd
.check_cpus <- function(cpus) {
    if (!is.numeric(cpus) || cpus <= 0)
        stop("Number of CPUs should be a positive integer.", call. = FALSE)
    
    if (cpus < 0 || is.null(cpus))
        cpus <- 1
    
    if (cpus > parallel::detectCores())
        message(sprintf("\nOnly %s CPUs available for parallel processing.\n", parallel::detectCores()))
    return(cpus)
    
}

## ----------------------------- .check_alpha ----------------------------- ##
#' Check significance threshold sanity
#'
#' @param alpha numeric, significance threshold for t.test
#'
#' @return NULL if acceptable value, otherwise condition
#' @author Al J Abadi
#'
#' @keywords Internal
#' @noRd
.check_alpha <- function(alpha=NULL) {
    if (is.null(alpha))
        alpha <- 0.01
    
    if (!is.numeric(alpha) || alpha < 0 || alpha > 1)
        stop("invalid 'signif.threshold'. Use 0.01 or 0.05.", call. = FALSE)
    alpha
}

## ----------------------------- .check_comp ----------------------------- ##
#' Check comp 
#'
#' Check comp is an integer of length (or 3)
#' @param comp comp arg
#' @param ncomp maximum comp allowed
#' @param style character, one of c('2d', '3d')
#'
#' @return Integer vector, or a condition
#' @keywords Internal
#' @noRd
#' @examples
#' \dontrun{
#' .check_comp(c(1,2), ncomp = 4)
#' .check_comp(c(1,5), ncomp = 4)
#' .check_comp(c(1,3,4), ncomp = 4)
#' .check_comp(c(1,3,4), ncomp = 4, style = '3d')
#' }
.check_comp <- function(comp, ncomp, style = c('2d', '3d'))
{
    formals(stop)$call. <- FALSE
    style <- match.arg(style)
    
    if (style == '2d')
        len <- 2
    else
        len <- 3
    
    err_msg <- paste0("'comp' must be a positive integer of length ", len,
                      " with values > 1 and <= ", ncomp, "\n")
    
    if (mode(comp) != 'numeric' || length(comp) != len | any(comp%%1 != 0) | any(comp < 1) | any(comp > ncomp))
        stop(err_msg)

    return(comp)
}

## --------------------------- .check_character --------------------------- ##
#' Check character input
#'
#' Check character and its length
#' @param arg character arg
#' @param len expected length
#'
#' @return character, or condition
#' @keywords Internal
#' @noRd
#' @examples
#' \dontrun{
#' .check_character(arg = 'foo', len = 1)
#' bar <- 'foo'
#' .check_character(arg = bar, len = 2)
#' .check_character(arg = NULL, default = 'foo')
#' .check_character(arg = NULL, default = c('foo', 'bar'))
#' .check_character(arg = 11)
#' }
.check_character <- function(arg, len = 1, default = NULL)
{ # TODO use this throughout for title etc
    formals(stop)$call. <- FALSE
    err_msg <- sprintf("'%s' must be character of length %s\n", deparse(substitute(arg)), len)
    if (is.null(arg))
    {
        arg <- default
        len <- length(arg)
    }
    
    if (is.character(arg))
    {
        if (length(arg) != len)
            stop(err_msg)
    } else {
        stop(err_msg)
    }
    
    return(arg)
}

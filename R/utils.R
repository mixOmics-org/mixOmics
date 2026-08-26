## ------------------------ .plotLoadings_barplot ------------------------- ##
#' plotLoadings helper
#'
#' @noRd
# .plotLoadings_barplot <- function(height, col, names.arg, cex.name, border, xlim) {
#     tryCatch({barplot(height, horiz = TRUE, las = 1, col = col, axisnames = TRUE, names.arg = names.arg, #names.arg = row.names(df),
#                       cex.names = cex.name, cex.axis = 0.7, beside = TRUE, border = border, xlim = xlim)},
#              error = function(e){
#                  if ( grepl(pattern = "figure margins too large", e) ){
#                      stop("\nplotLoadings encountered margin errors. Ensure feature names are not too long (see 'name.var' argument) and the 'Plots' pane is cleared and enlargened.\n", call. = FALSE)
#                  } else {
#                      stop(e$message, call. = FALSE)
#                  }
#              })
# }

.plotLoadings_barplot <- function(height, col, names.arg, cex.name, border, xlim, xlab, ylab, cex.lab = 1, cex.axis = size.axis) {
  tryCatch({
    barplot(height, horiz = TRUE, las = 1, col = col, names.arg = names.arg, cex.names = cex.name, cex.axis = cex.axis, beside = TRUE, 
            border = border, xlim = xlim, xlab = xlab, ylab = ylab, cex.lab = cex.lab)
  }, error = function(e) {
    if (grepl('figure margins too large', e$message)) {
      stop('plotLoadings encountered margin errors. Ensure feature names are not too long and the "Plots" pane is enlarged.', call. = FALSE)
    } else {
      stop(e$message, call. = FALSE)
    }
  })
}

## --------------------------- .unexpected_err ---------------------------- ##
#' Unexpected error handler for the package
#'
#' To be used in unexperimented situations where handlers fail
#' @param trying_to character, the context in which unexpected error occurred.
#'
#' @return Error message
#' @author Al J Abadi
#'
#' @noRd
.unexpected_err <- function(trying_to = NULL) {
    trying_to <- ifelse(is.null(trying_to), "", sprintf(" while trying to %s", trying_to))
    msg <- sprintf("Unexpected error%s. Please check the inputs and if problem persists submit an issue to https://github.com/mixOmicsTeam/mixOmics/issues", trying_to)
    stop(msg, call. = FALSE)
}

## ------------------------------- .on_unix ------------------------------- ##
#' Check OS type for parallel processing
#'
#' @return Logical, FALSE if windows OS, TRUE if unix OS.
#' @author Al J Abadi
#'
#' @noRd
.onUnix <- function() {
    return(ifelse(.Platform$OS.type == "unix", TRUE, FALSE))
}

## ----------------------- .unlist_repeat_cv_output ----------------------- ##
#' repeat_cv_perf.diablo helper to unlist internal outputs to
#' previous list format for downstream work
#'
#' @return list of outputs needed for downstream work
#' @author Al J Abadi
#'
#' @noRd
.unlist_repeat_cv_output <- function(list_rep=NULL) {
    
    ## so lapply will return a named list
    names_list <- as.list(names(list_rep[[1]]))
    names(names_list) <- names(list_rep[[1]])
    
    list_nrep <- lapply(names_list, function(x){
        lapply(list_rep, function(y) {y[[x]]})
    })
    
    
    return(list_nrep)
}

## ------------------------ stratified_subsampling ------------------------ ##
#' Perform stratified subsampling for cross-validation
#'
#' @param Y A factor or a class vector for the discrete outcome
#' @param folds  Integer, if validation = Mfold, how many folds?
#' Mainly to be used for reproducibility and unit tests.
#'
#' @return A list containing:
#' \item{SAMPLE} List of samples for each fold
#' \item{stop} A warning signal for when there is not enough samples in a 
#' class
#' @author Florian Rohart, Al J Abadi
#' @noRd
stratified.subsampling <- function(Y, folds = 10)
{
    
    stop <- 0
    for (i in seq_len(nlevels(Y)))
    {
        ai <- sample(which(Y == levels(Y)[i]), replace = FALSE) # random sampling of the samples from level i
        aai <- suppressWarnings(split(ai, factor(seq_len(
            min(folds, length(ai))
        ))))                       # split of the samples in k-folds
        if (length(ai) < folds)
            # if one level doesn't have at least k samples, the list is completed with "integer(0)"
        {
            for (j in (length(ai) + 1):folds)
                aai[[j]] <- integer(0)
            stop <- stop + 1
        }
        assign(paste("aa", i, sep = "_"), sample(aai, replace = FALSE))         # the `sample(aai)' is to avoid the first group to have a lot more data than the rest
    }
    
    # combination of the different split aa_i into SAMPLE
    SAMPLE <- list()
    for (j in seq_len(folds))
    {
        SAMPLE[[j]] <- integer(0)
        for (i in seq_len(nlevels(Y)))
        {
            SAMPLE[[j]] <- c(SAMPLE[[j]], get(paste("aa", i, sep = "_"))[[j]])
        }
    }# SAMPLE is a list of k splits
    
    ind0 <- sapply(SAMPLE, length)
    if (any(ind0 == 0))
    {
        SAMPLE <- SAMPLE [-which(ind0 == 0)]
        message(
            "Because of a too high number of 'folds' required, ",
            length(which(ind0 == 0)),
            " folds were randomly assigned no data: the number of 'folds' is reduced to ",
            length(SAMPLE)
        )
    }
    
    return(list(SAMPLE = SAMPLE, stop = stop))
}

## ------------------------------ .name_list ------------------------------ ##
#' Create a named list
#'
#' Creates a named list of a character vector where names and values
#' are the same, useful for lapply family.
#' @param char Character vector.
#'
#' @return A named list.
#' @author Al J Abadi
#'
#' @noRd
.name_list <- function(char, names=NULL) {
    out <- as.list(char)
    names <- if (is.null(names)) char else names
    names(out) <- names
    return(out)
}

## ------------------------------ .mixo_rng ------------------------------- ##
#' RNG used for package tests
#'
#' Creating a function so that it can easily be changed.
#' It's mainly used in form of the given example in unit tests.
#'
#' @return An R version
#' @author Al J Abadi
#' @example RNGversion(.mixo_rng())
#' @noRd
.mixo_rng <- function() {"3.6.0"}

## --------------------------- .change_if_null ---------------------------- ##
#' Set default value if NULL
#'
#' For arguments with default NULL value, sets the desired default value
#' 
#' @param arg The function argument with NULL/no default value
#' @param default The desired default value
#'
#' @return The desired default value
#' @author Al J Abadi
#'
#' @examples
#' .change_if_null(NULL, 10)
#' #> 10
#' .change_if_null(3, 10)
#' #> 3
#' @noRd
.change_if_null <- function(arg, default) {
    return(
        if ( missing(arg) || is.null(arg) )
            default
        else
            arg
    )
}

## ------------------------ .add_average_blocks ------------------------- ##
#' Add average blocks to DIABLO object
#'
#' For plotIndiv(..., blocks = "average")
#' 
#' @param block_object A diablo object
#' @param average_blocks One or both of c('average', 'weighted.average')
#'
#' @return A diablo object with average blocks added for smooth input into plotIndiv.sgccda
#' @author Al J Abadi
#' @examples 
#' \dontrun{
#' data(nutrimouse)
#' Y = nutrimouse$diet
#' data = list(gene = nutrimouse$gene, lipid = nutrimouse$lipid)
#' design1 = matrix(c(0,1,0,1), ncol = 2, nrow = 2, byrow = TRUE)
#' 
#' nutrimouse.sgccda1 <- wrapper.sgccda(X = data,
#'                                      Y = Y,
#'                                      design = design1,
#'                                      ncomp = 2,
#'                                      keepX = list(gene = c(10,10), lipid = c(15,15)),
#'                                      scheme = "centroid")
#' nutrimouse.sgccda.average <- mixOmics:::.add_average_blocks(nutrimouse.sgccda1, average_blocks = "average")
#' names(nutrimouse.sgccda.average$X)
#> "gene"      "lipid"     "average"
#' }
#' @noRd
.add_average_blocks <- function(block_object, average_blocks = c('average', 'weighted.average')) {
    X_blocks <- with(block_object, names$blocks[-which(names$block == 'Y')])
    average_blocks <- match.arg(average_blocks, several.ok = TRUE)
    
    .get_average_variates <- function(object, X_blocks, weighted = FALSE) {
        
        arrays <- object$variates
        arrays <- arrays[X_blocks]
        if (isTRUE(weighted)) {
            if (!is(block_object, "sgccda")) {
                if ('weighted.average' %in% average_blocks ) {
                    stop("'weighted.average' plots are only available for block.splsda objects ")
                }
            }
            weights <- object$weights
        } else {
            ## if unweighted, create a weight data.frame of 1 to avoid code complications
            weights <- matrix(rep(1, length(X_blocks)*object$ncomp[1]), nrow = length(X_blocks))
            dimnames(weights) <- list(X_blocks, paste0("comp", seq_len(object$ncomp[1])))
            weights <- as.data.frame(weights)
        }
        block_names <- .name_list(names(arrays))
        
        weighted_arrays <- lapply(block_names, function(x){
            variates <- arrays[[x]]
            weights <- diag(weights[x, ])
            weighted_variates <- variates %*% weights
            dimnames(weighted_variates) <- dimnames(variates)
            weighted_variates
        })
        wtd_sum <- Reduce(f = '+', weighted_arrays)
        ## weighted mean = weighted sum / sum(weight)
        sweep(wtd_sum, MARGIN = 2, colSums(weights), FUN = "/")
    }
    
    for (average_block in average_blocks) {
        block_object$X[[average_block]] <-  0
        if (average_block == "weighted.average") {
            block_object$variates[[average_block]] <-  .get_average_variates(object = block_object, X_blocks = X_blocks, weighted = TRUE)
        }
        if (average_block == "average") {
            block_object$variates[[average_block]] <-  .get_average_variates(object = block_object, X_blocks = X_blocks, weighted = FALSE)
        }
        block_object$names$blocks <- c(block_object$names$blocks, average_block)
        block_object$ncomp[average_block] <- block_object$ncomp[1]
        block_object$prop_expl_var[average_block] <- 0
        
    }
    
    block_object
}

## ------------------------------- mat.rank ------------------------------- ##
#' Matrix Rank
#' 
#' This function estimate the rank of a matrix.
#' 
#' \code{mat.rank} estimate the rank of a matrix by computing its singular
#' values \eqn{d[i]} (using \code{nipals}). The rank of the matrix can be
#' defined as the number of singular values \eqn{d[i] > 0}.
#' 
#' If \code{tol} is missing, it is given by
#' \code{tol=max(dim(mat))*max(d)*.Machine$double.eps}.
#' 
#' @param mat a numeric matrix or data frame that can contain missing values.
#' @param tol positive real, the tolerance for singular values, only those with
#' values larger than \code{tol} are considered non-zero.
#' @return The returned value is a list with components: \item{rank}{a integer
#' value, the matrix rank.} \item{tol}{the tolerance used for singular values.}
#' @author Sébastien Déjean, Ignacio González, Al J Abadi
#' @seealso \code{\link{nipals}}
#' @keywords algebra
#' @export
#' @examples
#' 
#' ## Hilbert matrix
#' hilbert <- function(n) { i <- 1:n; 1 / outer(i - 1, i, "+") }
#' mat <- hilbert(16)
#' mat.rank(mat)
#' 
#' \dontrun{
#' ## Hilbert matrix with missing data
#' idx.na <- matrix(sample(c(0, 1, 1, 1, 1), 36, replace = TRUE), ncol = 6)
#' m.na <- m <- hilbert(9)[, 1:6]
#' m.na[idx.na == 0] <- NA
#' mat.rank(m)
#' mat.rank(m.na)
#' }

mat.rank = function (mat, tol)
{
    
    if (length(dim(mat)) != 2) 
        stop("'mat' must be a numeric matrix.")
    
    mat = as.matrix(mat)
    
    if (!is.numeric(mat)) 
        stop("'mat' must be a numeric matrix.")
    
    d = nipals(mat)$eig
    max.d = d[1]
    min.d = d[length(d)]
    if (missing(tol)) 
        tol = max(dim(mat)) * max.d * .Machine$double.eps
    
    r = sum(d > tol)
    
    return(list(rank = r, tol = tol))
}

## ----------------------------- nearZeroVar ------------------------------ ##
#' Identification of zero- or near-zero variance predictors
#' 
#' Borrowed from the \pkg{caret} package. It is used as an internal function in
#' the PLS methods, but can also be used as an external function, in
#' particular when the data contain a lot of zeroes values and need to be
#' pre-filtered beforehand.
#' 
#' This function diagnoses predictors that have one unique value (i.e. are zero
#' variance predictors) or predictors that are have both of the following
#' characteristics: they have very few unique values relative to the number of
#' samples and the ratio of the frequency of the most common value to the
#' frequency of the second most common value is large.
#' 
#' For example, an example of near zero variance predictor is one that, for
#' 1000 samples, has two distinct values and 999 of them are a single value.
#' 
#' To be flagged, first the frequency of the most prevalent value over the
#' second most frequent value (called the ``frequency ratio'') must be above
#' \code{freqCut}. Secondly, the ``percent of unique values,'' the number of
#' unique values divided by the total number of samples (times 100), must also
#' be below \code{uniqueCut}.
#' 
#' In the above example, the frequency ratio is 999 and the unique value
#' percentage is 0.0001.
#' 
#' @param x a numeric vector or matrix, or a data frame with all numeric data.
#' @param freqCut the cutoff for the ratio of the most common value to the
#' second most common value.
#' @param uniqueCut the cutoff for the percentage of distinct values out of the
#' number of total samples.
#' @return \code{nearZeroVar} returns a list that contains the following
#' components:
#' 
#' \item{Position}{a vector of integers corresponding to the column positions
#' of the problematic predictors that will need to be removed.}
#' \item{Metrics}{a data frame containing the zero- or near-zero predictors
#' information with columns: \code{freqRatio}, the ratio of frequencies for the
#' most common value over the second most common value and,
#' \code{percentUnique}, the percentage of unique data points out of the total
#' number of data points.}
#' @author Max Kuhn, Allan Engelhardt, Florian Rohart, Benoit Gautier, AL J Abadi
#' for mixOmics
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}
#' @keywords utilities
#' @export
#' @examples
#' 
#' data(diverse.16S)
#' nzv = nearZeroVar(diverse.16S$data.raw)
#' length(nzv$Position) # those would be removed for the default frequency cut
nearZeroVar = function (x, freqCut = 95/5, uniqueCut = 10)
{
    
    if (is.vector(x))
        x = matrix(x, ncol = 1)
    
    freqRatio = apply(x, 2, function(data)
    {
        data = na.omit(data)
        
        if (length(unique(data)) == length(data))
        { # No duplicate
            return(1)
        } else if (length(unique(data)) == 1) { # Same value
            return(0)
        } else {
            t = table(data)
            return(max(t, na.rm = TRUE)/max(t[-which.max(t)], na.rm = TRUE))
        }
    })
    
    lunique = apply(x, 2, function(data) length(unique(data[!is.na(data)])))
    
    percentUnique = 100 * lunique/nrow(x)
    zeroVar = (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
    
    out = list()
    out$Position = which((freqRatio > freqCut & percentUnique <= uniqueCut) | zeroVar)
    names(out$Position) = NULL
    out$Metrics = data.frame(freqRatio = freqRatio, percentUnique = percentUnique)
    out$Metrics = out$Metrics[out$Position, ]
    return(out)
}


## ----------------------------- %=% ----------------------------- ##
## check if LHS and RHS contain exactly the same elements without caring about order
## TODO use this throughout for col.per.group etc.
#'@noRd
#'@examples
#' letters %=% rev(letters)
#' #> TRUE
'%=%' <- function(LHS, RHS)
{
    diffs <- c(setdiff(LHS, RHS), setdiff(RHS, LHS))
    if (length(diffs) | (length(LHS) != length(RHS)))
    {
        return(FALSE)
    }
    return(TRUE)
}
## NOT '%=%'
#' @noRd
#' @examples
#' letters %!=% rev(letters)
#' #> FALSE
#' letters %!=% LETTERS
#' #> TRUE
'%!=%' <- function(LHS, RHS)
{
    ! (LHS %=% RHS)
}

## ---------------------------- .create_design ---------------------------- ##
#' create design matrix for block.(s)pls(da)
#'
#' @param X Named list of datasets. When
#' type = 'null', for non-DA analyses the first one is 
#' taken to be the response matrix which is fully connected to others.
#' @param design Character, one of c('full', 'null'). If 'full' all blocks will
#'   be connected, otherwise only the first block is connected. Alternatively, a
#'   numeric between 0 and 1 which specifies all off-diagonal elements of a
#'   fully connected matrix. Default is 'full'.
#' @param indY Integer, index of Y dataset. It could also be NUL/missing.
#' @return A matrix whose both dimension names match the names of \code{X} while
#' after adding a placeholder 'Y' to it.
#' @noRd
#' @examples
#' \dontrun{
#' .create_design(list(rna = cars, met = mtcars, acc = faithful), design = 'full')
#' .create_design(list(rna = cars, met = mtcars, acc = faithful), design = 'null')
#' .create_design(list(rna = cars, met = mtcars, acc = faithful), design = 'null', indY = 2)
#' .create_design(list(rna = cars, met = mtcars, acc = faithful), design = 0.3)
#' .create_design(list(rna = cars, met = mtcars, acc = faithful), design = 0.3, indY = 2)
#' 
#' data("breast.TCGA")
#' data = list(mrna = breast.TCGA$data.train$mrna, mirna = breast.TCGA$data.train$mirna,
#'             protein = breast.TCGA$data.train$protein)
#' diag(design) =  0
#' ncomp = c(2)
#' list.keepX = list(mrna = rep(20, 2), mirna = rep(10,2), protein = rep(10, 2))
#' 
#' block.splsda(X = data, Y = breast.TCGA$data.train$subtype, ncomp = ncomp, keepX = list.keepX, design = 'null')$design
#' block.splsda(X = data, Y = breast.TCGA$data.train$subtype, ncomp = ncomp, keepX = list.keepX, design = 'full')$design
#' block.splsda(X = data, Y = breast.TCGA$data.train$subtype, ncomp = ncomp, keepX = list.keepX, design = 0.1)$design
#' block.spls(X = data, Y = data$mrna, ncomp = ncomp, keepX = list.keepX, design = 'null')$design
#' block.spls(X = data, Y = data$protein, ncomp = ncomp, keepX = list.keepX, design = 'full')$design
#' block.spls(X = data, Y = data$protein, ncomp = ncomp, keepX = list.keepX, design = 0.1)$design
#'}
.create_design <- function(X, design = 'full', indY = NULL) {
    if (!all(is.list(X) && length(unique(names(X))) == length(X)))
        stop("'X' must be a named list. See documentation.", call. = FALSE)
    
    if ( missing(indY) || is.null(indY) || isTRUE(tryCatch(is.integer(as.integer(indY)) && indY > length(X))) ){
        indY <- length(X) + 1
        X <- c(X, list(Y = matrix())) ## just so we have Y in X for design, ind
    }
    
    blocks <- names(X)
    ## diag will be 0 at the end, specify off-diag elements
    off_diag <- 1
    if (is.character(design))
    {
        design <- match.arg(tolower(design), choices = c('full', 'null'))
        if (design == 'null')
            off_diag <- 0
    } else if (isTRUE(tryCatch(design<= 1 & design >= 0)))
    {
        off_diag <- design
    } else {
        stop("'design' must be a matrix, or one of c('full', 'null'), or a numeric ",
             "between 0 and 1. See documentation for details.")
    }
    
    design <-  matrix(off_diag, nrow = length(blocks), ncol = length(blocks) ,dimnames = rep(list(blocks), 2))
    
    if (!(missing(indY) || is.null(indY)) ) {
        if (isTRUE(tryCatch(is.integer(as.integer(indY)) && indY <= length(X))))
        {
            indY <- as.integer(indY)
        } else 
        {
            stop("'indY' must be an integer from 1:length(X):", seq_along(X))
        }
        design[,indY] <-  design[indY,] <- 1
    }
    
    diag(design) <- 0
    return(design)
}

## ----------------------- .check_zero_var_columns ------------------------ ##
#' Check if scaling can be performed (no constant vectors)
#' @noRd
#' @examples
#' \dontrun
#' {
#' .check_zero_var_columns(matrix(rep(c(1, 2, 3), 2), nrow=2, byrow=TRUE))
#' }
# TODO use this through package
.check_zero_var_columns <- function(X, scale = TRUE, block_name = 'X')
{
    if (!isFALSE(scale))
    {
        zero_var_cols <- which(colVars(X, na.rm = TRUE) == 0)
        #' @importFrom matrixStats colVars
        if (length(zero_var_cols) > 0)
            stop("columns with zero variance in '", 
                 block_name, 
                 "': ",
                 paste(zero_var_cols, collapse = ','),
                 ". Remove these columns before scaling.\n",
                 call. = FALSE)
    }
    NULL
}

## ----------- .stop ----------- 
## custom stop to define specific error classes for testthat
.stop <- function(message, .subclass = NULL,  call = NULL, ...) {
    formals(stop)$call. <- FALSE
    msg <- structure(
        list(
            message = message,
            call = call,
            ...
        ),
        class = c(.subclass, "error", "condition")
    )
    stop(msg)
}

## ----------- .warning ----------- 
## custom warning to define specific warning classes for testthat and custom suppression
.warning <- function(message, .subclass = NULL) {
    formals(warning)$call. <- FALSE
    msg <- structure(
        list(
            message = message
        ),
        class = c(.subclass, "character")
    )
    warning(msg)
}

## ----------- .message ----------- 
## custom message to define specific message classes for testthat and custom suppression
.message <- function(message, .subclass = NULL) {
    msg <- structure(
        list(
            message = message
        ),
        class = c(.subclass)
    )
    message(msg)
}

## ----------------------- .progressBar ------------------------ ##
#' Nice progress bar
#' @noRd
#' @examples
#' \dontrun
#' {
#' for (i in 1:30)
#' {
#'     Sys.sleep(0.02)
#'     progress(i, max = 30)
#' }
#' }
# TODO use this through package
.progressBar <- function (x, max = 1, title = NULL, breakline = TRUE) {
    percent <- x / max * 100
    if (!is.null(title))
        cat(sprintf("\033[33m  %s\033[39m", title))
        
    cat(sprintf('\033[33m\r[%-50s] %d%%\033[39m',
                paste(rep('=', percent / 2), collapse = ''),
                floor(percent)))
    if (x == max & breakline)
        cat('\n')
    invisible(NULL)
}

## ----------------------- .relist ------------------------ ##
#' Swap list hierarchy
#'
#' For a nested list of say indices within categories, swap the hierarchy
#' to categories within indices. Used to processs CV results. See examples.
#' @param lst 
#'
#' @return A re-ordered list
#' @noRd
#' @keywords Internal
#' @examples
#' .relist(list(cat1 = list('1' = 11, '2' =  12),
#'                           cat2 = list('1' = 21, '2' = 22)))
#' #> list('1' = list(cat1 = 11, cat2 = 21), '2' = list(cat1 = 12, 
#' #>                                                   cat2 = 22))
.relist <- function(lst) {
    lapply(.name_list(names(lst[[1]])), function(entry) {
        lapply(.name_list(names(lst)), function(x) lst[[x]][[entry]])
    })
}
## ---------------------------- mixo gg theme ----------------------------- ##

#' ggplot2 theme
#'
#' @param cex text size
#' @param x.angle x text angle
#' @param background.fill background fill for plot
#'
#' @noRd
mixo_gg.theme <- function(cex, x.angle = 90, background.fill = 'grey97', subtitle.cex = NULL) {
    subtitle.cex <- .change_if_null(subtitle.cex, cex)
    text.size <- as.integer(cex*10)
    subtitle.text.size <- as.integer(subtitle.cex*10)
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(linewidth = 0.5, linetype = "solid",
                                   colour = "black"),
          panel.background = element_rect(fill = background.fill),
          
          axis.text = element_text( size = text.size ),
          axis.text.x = element_text( size = text.size, angle = x.angle, hjust = x.angle/90),
          axis.title = element_text( size = text.size),
          legend.text = element_text( size = text.size ),
          legend.title =  element_text( size = text.size),
          plot.title = element_text(hjust = 0.5),
          # subtitles
          strip.text = element_text(size = subtitle.text.size, 
                                    face = 'bold')
    )
}

#' Format text with colour to input into cat
#' 
#' @example cat(.colour.txt('foo'))
#' @noRd
#' @keywords Internal
.colour.txt <- function(char, sQuote = TRUE) {
    if (sQuote) ## add single quotes around char
        sprintf("\033[33m'%s'\033[39m", char)
    else
        sprintf("\033[33m%s\033[39m", char)
}

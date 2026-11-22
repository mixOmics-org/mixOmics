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

.are.colors <- function(cols) {
    #TODO use this to deduplicate code in internals
    ## check cols, if all valid colors, return cols, else, throw error
    all.cols <- sapply(cols, function(X) {
        tryCatch(is.matrix(col2rgb(X)), 
                 error = function(e) FALSE)
    })
    if (any(!all.cols))
    {
        stop('invalid colour character/vector.\n', call. = FALSE)
    }
    cols
}

.get.colors <- function(col, n_col=0) {
    #TODO use this to deduplicate code in internals
    ## ensure 'col' is valid color(s) and of length 'col.length', repeat it if
    ## necessary
    ## return a vector of lenfth col.length, or throw condition
    if (is.null(col))
    {
        col <- 'black'
    }
    col <- .are.colors(col)
    if (length(col) == 1) 
    {
        col <- rep(col, n_col)
    }
    else if (length(col == n_col))
    {
        col <- col
    }
    else
    {
        stop("'col' must be character of length 1 or ", paste0(n_col))
    }
    col
}

.get.character.vector <- function(arg, vec) {
    #TODO use this to deduplicate code in internals
    ## if 'arg' is TRUE, return vec, if FALSE, return NULL, if character, ensure
    ## it's of the same length as vec
    ## return a character of the same length as vec
    stopifnot(is.character(vec))
    fail.msg <- paste0(sQuote(deparse(substitute(arg))), 
                       " must be either logical or a character vector of length ", length(vec))
    
    if (isTRUE(arg)) 
    {
        res <- vec
    } 
    else if (isFALSE(arg))
    {
        res <- NULL
    }
    else if (is.character(arg))
    {
        if (length(arg) != length(vec))
        {
            stop(fail.msg, call. = FALSE)
        }
        else 
        {
            res <- arg
        }
    } else
    {
        res <- tryCatch(as.character(arg))
        if (length(res) != length(vec)) 
        {
            stop(fail.msg, call. = FALSE)
        }
        
    }
    res
}

#' @noRd
#' @examples 
#' .get.group(group = NULL, object = NULL, n_ind = 10)
#' .get.group(group = letters[1:10], object = NULL, n_ind = 10)
#' .get.group(group = NULL, object = structure(list(Y = 10:19), class='DA'), n_ind = 10)
#' \dontrun{
#' .get.group(group = letters[1:11], object = NULL, n_ind = 10)
#' .get.group(group = NULL, object = structure(list(Y = 10:19), class='DA'), n_ind = 10)
#' }
.get.group <- function(group, object=NULL, n_ind=0) {
    #TODO use this to deduplicate code in internals
    ## check group, if and if it's not NULL return a factor
    ## if object is of class 'DA', return object$Y which should be factor
    out <- group
    ## if no group given
    if (is.null(out))
    {
        if (is (object, 'DA'))
        {
            out <- object$Y
            if (n_ind != length(out))
            {
                stop("Unexpected error. Please submit an issue with reproducible example to\n",
                     "https://github.com/mixOmics-org/mixOmics/issues")
            }
        }
        else
        {
            ## for consistency of plot so group is always an aes
            out <- factor(rep(1, n_ind))
        }
        
    } 
    ## group is factor and does not drop NAs
    if (!is.factor(out))
    {
        out[is.na(out)] <- "missing group (NA)"
        out <- as.factor(out)
    }
    ## group length is valid
    if (length(out) != n_ind)
    {
        stop("'group' must be a factor of length ", n_ind)
    }
    
    return(out)
}

#'@noRd
#'@examples
#' .get.cols.and.group(n_ind = 3)
#' .get.cols.and.group(col = 'blue', n_ind = 3)
#' .get.cols.and.group(group = c(2,2,1))
#' .get.cols.and.group(col.per.group = c('1'='red', '2'='blue'), group = c(2,2,1))
.get.cols.and.group <- function(
    col.per.group=NULL, 
    group=NULL, 
    col=NULL, 
    object=NULL, 
    n_ind=ifelse(is.null(group), 0, length(group))
) {
    #TODO use this to deduplicate code in internals
    ## check group, if and if it's not NULL return a factor
    ## if object is of class 'DA', return object$Y which should be factor
    
    ## group factor
    group <- .get.group(group, object, n_ind)
   ## create col.per.group
    if (is.null(col.per.group))
    {
        if (nlevels(group) == 1)
        {
            col.per.group <- ifelse(is.null(col), color.mixo(1), .are.colors(col))
        } else if (nlevels(group) < 10)
        { 
            col.per.group <- color.mixo(seq_len(nlevels(group)))
        } else if (nlevels(group) >= 10)
        { 
            col.per.group <- color.jet(nlevels(group))
        }
       
        names(col.per.group) <- levels(group)
    }
    ## valid col.per.group
    if (!all(names(col.per.group) %in% levels(group)))
    {
        stop("'col.per.group' should be a named vector whose names match group (",
             levels(group), ")")
    }
    ## return group and col.per.group
    return(list(
        group = group,
        col.per.group = col.per.group
    ))
    
}
.get.ind.colors <- function(group=NULL, col=NULL, col.per.group=NULL, n_ind=0) {
    #TODO use this to deduplicate code in internals
    ## consider all valid combinations and return a vector of colors
    ## if group != NULL, return a vector whose names are colors and values the
    ## group levels for ggplot
    ## ------------- group == NULL
    if (is.null(group))
    {
        col <- .get.colors(col = col, n_col = n_ind)
    }
    else
    {
        col.per.group <- 
            .change_if_null(col.per.group, 
                         color.mixo(seq_len(1+nlevels(group))[-3] ## bc 3rd one is grey
                         )
            )
        names(col.per.group) <- levels(group)
        
        if (!all(names(col.per.group) %in% levels(group)))
        {
            stop("'col.per.group' should be a named vector whose names match group (",
                 levels(group), ")")
        }
        col <- group
        names(col) <- col.per.group[group]
    }
    return(col)
}

#' @noRd
#' @examples
#' .get.pch(pch = NULL, pch.levels = NULL)
#' .get.pch(pch = 6, pch.levels = NULL)
#' .get.pch(pch = mtcars$cyl, pch.levels = NULL)
.get.pch <- function(pch = NULL,
                     pch.levels = NULL,
                     pch.legend = TRUE,
                     n_ind = length(pch))
{
    ## NULL -> 16
    if (is.null(pch))
    {
        pch <- 16
    }
    
    if (length(pch) == 1)
    {
        names(pch.levels) <- pch.levels <-  pch
        pch <- rep(pch, n_ind)
    }
    
    pch <- factor(pch)
    
    if (is.null(pch.levels))
    {
        all.pch.levels <- 0:24
        main.pch.levels <- c(16, 15, 17, 3, 4)
        ref.pch.levels <- c(main.pch.levels, all.pch.levels[-which(all.pch.levels %in% main.pch.levels)])
        
        pch.levels <- ref.pch.levels[seq_len(nlevels(pch))]
        names(pch.levels) <- levels(pch)
        
    }
    else
    {
        # if ( length(setdiff(names(pch.levels),  levels(pch)) ) > 0 |
        #      length(setdiff(levels(pch),  names(pch.levels)) ) > 0 )
        # {
        #     stop("'pch.levels' must be a vector of pch values with names: ", 
        #          levels(pch), call. = FALSE)
        # }
    }
pch.legend <- pch.legend & ifelse(length(pch.levels) == 1, FALSE, TRUE)
    return(list(pch = pch, 
                pch.levels = pch.levels,
                pch.legend = pch.legend))
}

.check.cutoff <- function(cutoff = 0, 
                          style = 'ggplot2') {
    ncomp.vis <- ifelse(style != '3d', yes = 2, no = 3)
    cutoff <- .change_if_null(cutoff, 0)
    ## check that cutoff value is valid
    ## it could be of length 2 when different cutoff for different components
    ## are needed (useful when most of selected variables tend to contribute
    ## to one component only so we can tailor the cutoff)
    err <- FALSE
    if (length(cutoff) == 1)
    {
        if (! (is.numeric(cutoff) | (cutoff > 1) | (cutoff < 0)) )
        {
            err <- TRUE   
        } else
        {
            cutoff <- rep(cutoff, ncomp.vis)
        }
    }
    else if (length(cutoff) == ncomp.vis)
    {
        if (! (all(is.numeric(cutoff) | (cutoff > 1) | (cutoff < 0)) ))
            err <- TRUE
    }
    else
    {
        err <- TRUE
    }
    
    if (isTRUE(err))
    {
        stop("cutoff must be a numeric of length 1 (or 'ncomp') between 0 and 1", call. = FALSE)
    }
    return(cutoff)
}

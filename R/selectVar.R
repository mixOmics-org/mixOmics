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

# ========================================================================================================
# selectVar: output the variables that were selected on 'comp', for 'block'
# ========================================================================================================

#' Output of selected variables
#' 
#' This function outputs the selected variables on each component for the
#' sparse versions of the approaches (was also generalised to the non sparse
#' versions for our internal functions).
#' 
#' \code{selectVar} provides the variables selected on a given component. \
#' \describe{ \item{list("name")}{outputs the name of the selected variables
#' (provided that the input data have column names) ranked in decreasing order of
#' importance.} \item{list("value")}{outputs the loading value for each
#' selected variable, the loadings are ranked according to their absolute
#' value.} } These functions are only implemented for the sparse versions.
#' 
#' @aliases selectVar selectVar.mixo_pls selectVar.mixo_spls selectVar.pca
#' selectVar.sgcca selectVar.rgcca select.var
#' @param object object of class inherited from \code{"pls"}, \code{"spls"},
#'   \code{"plsda"},\code{"splsda"},\code{"sgcca"}, \code{"rgcca"},
#'   \code{"pca"}, \code{"spca"}, \code{"sipca"}.
#' @param comp integer value indicating the component of interest.
#' @param block for an object of class \code{"sgcca"}, the block data sets can
#' be specified as an input vector, for example \code{c(1,2)} for the first two
#' blocks. Default to NULL (all block data sets)
#' @param ... other arguments.
#' @return none
#' @author Kim-Anh Lê Cao, Florian Rohart, Al J Abadi
#' @export
#' @example ./examples/selectVar-examples.R
selectVar <-
    function(...)
        UseMethod("selectVar")

## ------------------------------- Methods -------------------------------- ##

#' @rdname selectVar
#' @method selectVar mixo_pls
#' @export
selectVar.mixo_pls  <- function(object, comp =1, block=NULL, ...)
{
    
    # check arguments
    # -----------------
    if (length(comp) > 1)
        stop("Expecting one single value for 'comp'")
    
    if (is.null(block))
    {
        if (any(comp > object$ncomp))
            stop("'comp' is greater than the number of components in the fitted model")
        null.block=TRUE
        block=1:length(object$loadings)
        
    }else{
        if (any(class(object)%in%c("pca")))
            object$names$blocks="X"
        
        if (is.numeric(block))
        {
            if (any(block>length(object$names$blocks)))
                stop("'block' needs to be lower than the number of blocks in the fitted model, which is length(object$names$blocks)")
            
        }else if (is.character(block) & sum(!is.na(match(block,object$names$blocks)))==0) {
            stop("No entry of 'block'  match object$names$blocks")
            
        }else if (is.character(block) & sum(is.na(match(block,object$names$blocks)))>0) {
            warning("At least one entry of 'block' does not match object$names$blocks")
        }
        
        if (length(object$ncomp)>1)
        {
            if (any(comp > object$ncomp[block]))
                stop("'comp' is greater than the number of components in the fitted model for the block you specified. See object$ncomp")
            
        }else{
            if (any(comp > object$ncomp))
                stop("'comp' is greater than the number of components in the fitted model")
        }
        
        null.block=FALSE
    }
    
    # main function: get the names and values of the non zero loadings
    # -----------------
    out = lapply(object$loadings[block],get.name.and.value,comp=comp)
    
    
    # outputs
    # ----------
    #if all blocks are considered by default (null.block=TRUE) and it's a DA analysis, then we don't show Y
    if (null.block)
    {
        if (any(class(object)%in%c("block.plsda","block.splsda")))# the position of Y is in indY
        {
            out=out[-object$indY] #remove Y
        }else if (any(class(object)%in%c("mint.plsda","mint.splsda","mixo_plsda","mixo_splsda"))) {
            # Y is always in second position
            out=out[[1]]
        }else if (any(class(object)%in%c("pca"))) { #keep the result as a list
            out=out[[1]]
        }
        
    } else {
        if (length(grep("pca",class(object)))>0)
            out=out
    }
    
    #we add comp as an output
    out$comp=comp
    
    return(out)
}

#' @rdname selectVar
#' @method selectVar mixo_spls
#' @export
selectVar.mixo_spls <- selectVar.mixo_pls

#' @rdname selectVar
#' @method selectVar pca
#' @export
selectVar.pca <- selectVar.mixo_pls

#' @rdname selectVar
#' @method selectVar sgcca
#' @export
selectVar.sgcca <- selectVar.mixo_pls

#' @rdname selectVar
#' @method selectVar rgcca
#' @export
selectVar.rgcca <- selectVar.mixo_pls

## -------------------------------- helper -------------------------------- ##
get.name.and.value <- function(x,comp)
{
    value <- data.frame(value.var = x[,comp])
    rownames(value) <- rownames(x)
    value <- value[abs(value$value.var) > .Machine$double.eps,,drop=FALSE]
    value <- value[order(-abs(value$value.var)),,drop=FALSE]
    
    name.var <- rownames(value)
    return(list(name = name.var, value = value))
}


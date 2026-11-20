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

#----------------------------------------------------------------------------------------------------------#
#-- Includes plotArrow for PLS, sPLS, rCC,  rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#
#' Arrow sample plot
#' 
#' Represents samples from multiple coordinates to assess the alignment in the
#' latent space.
#' 
#' Graphical of the samples (individuals) is displayed in a superimposed manner
#' where each sample will be indicated using an arrow. The start of the arrow
#' indicates the location of the sample in \eqn{X} in one plot, and the tip the
#' location of the sample in \eqn{Y} in the other plot. Short arrows indicate a
#' strong agreement between the matching data sets, long arrows a disagreement
#' between the matching data sets. The representation space is scaled using the
#' range of coordinates so minimum and maximum values are equal for all blocks.
#' Since the algorithm maximises the covariance of these components, the
#' absolute values do not affect the alignment.
#' 
#' For objects of class \code{"GCCA"} and if there are more than 2 blocks, the
#' start of the arrow indicates the centroid between all data sets for a given
#' individual and the tips of the arrows the location of that individual in
#' each block.
#' 
#' @inheritParams biplot
#' @param object object of class inheriting from \pkg{mixOmics}: \code{PLS,
#' sPLS, rCC, rGCCA, sGCCA, sGCCDA}
#' @param pch plot character. A character string or a named vector of single
#' characters or integers whose names match those of \code{object$variates}.
#' @param ind.names.position One of c('start', 'end') indicating where to show
#'   the ind.names . Not used in block analyses, where centroids are used.
#' @param arrow.alpha Numeric between 0 and 1 determining the opacity of arrows.
#' @param arrow.size Numeric, variable arrow head size.
#' @param arrow.length Numeric, length of the arrow head in 'cm'.
#' @param ... Not currently used.
#' sample size to display sample names.
#' @return A ggplot object
#' @author Al J Abadi
#' @seealso \code{\link{arrows}}, \code{\link{text}}, \code{\link{points}} and
#' http://mixOmics.org/graphics for more details.
#' @references Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P.
#' (2009). Sparse canonical methods for biological data integration:
#' application to a cross-platform study. \emph{BMC Bioinformatics}
#' \bold{10}:34.
#' @keywords multivariate hplot dplot
#' @export
#' @example ./examples/plotArrow-examples.R
plotArrow <- function(object,
                      comp = c(1,2),
                      ind.names = TRUE,
                      group = NULL,
                      col = NULL,
                      ind.names.position = c('start', 'end'), 
                      ind.names.size = 2,
                      pch = NULL,
                      pch.size = 2,
                      # arrow.col = 'grey50',
                      arrow.alpha = 0.6,
                      arrow.size = 0.5,
                      arrow.length = 0.2,
                      legend = if (is.null(group)) FALSE else TRUE,
                      legend.title = NULL,
                      ...
)
{
    # TODO allow for arbitrary blocks (e.g. Y and X[[i]] in block.spls)
    ## ------------- checks
    class.object = class(object)
    object.pls=c("mixo_pls", "mixo_plsda", "mixo_spls", "mixo_splsda", "rcc")
    object.blocks=c("sgcca", "sgccda", "rgcca")
    
    if (! any(class.object %in% c(object.pls,object.blocks)))
        stop( " 'plotArrow' is only implemented for the following objects: pls, spls, rcc, sgcca, sgccda, rgcca", call.=FALSE)
    
    if ("DA" %in% class.object && !any(object.blocks %in% class.object)) {
        stop("'plotArrow' not implemented for (s)PLSDA or MINT sPLSDA", call.=FALSE)
    }
    
    ind.names.position <- match.arg(ind.names.position)
    is.multiblock <- ifelse(is.list(object$X), TRUE, FALSE)
    
    col.per.group <- col
    
    # check for inappropriate args
    extra_args <- list(...)
    # check for deprecated args
    if ("col.per.group" %in% names(extra_args)) {
      warning("'col.per.group' is deprecated, please use 'col' to specify colours for each group")
    }
    if ("pch.levels" %in% names(extra_args)) {
        warning("'pch.levels' is deprecated, please use 'pch' to specify point types")
    }
    
    ## keep the plot components only and name columns as x and y
    variates <- mapply(arr = object$variates, block = names(object$variates), FUN = function(arr, block) {
        df <- data.frame(arr[,comp])
        colnames(df) <- c('x', 'y')
        df
    }, SIMPLIFY = FALSE)
    ## remove the Y block and DA analyses
    if (is(object, 'DA'))
    {
        variates$Y <- NULL
    }
    
    ## standardise x and y for all blocks
    variates <- lapply(variates, function(df){
        x_min <- min(df[,1])
        x_max <- max(df[,1])
        y_min <- min(df[,2])
        y_max <- max(df[,2])
        df[,1] <- (df[,1] - x_min)/ (x_max - x_min) - 0.5
        df[,2] <- (df[,2] - y_min)/ (y_max - y_min) - 0.5
        df
    })
    
    blocks <- names(variates)
    if ('centroid' %in%  names(blocks))
        stop("'centroid' is a reserved name. Please choose another name for blocks")

    if (is.multiblock) {
        ## calculate and add centroids
        centroid_variates <- lapply(c(x='x', y='y'), function(w){
            xORy <- lapply(variates, function(v) v[,w, drop=FALSE])
            xORy <- Reduce(x = xORy, f = cbind)
            xORy <- rowMeans(xORy)
        })
        stopifnot(! ('centroid' %in% names(variates)))
        variates$centroid <- data.frame(centroid_variates)
    }

    ## prefix colnames by block for cbind() and then ggplot
    variates <- mapply(df = variates, block = names(variates), FUN = function(df, block) {
        colnames(df) <- sprintf("%s_%s", colnames(df), block)
        df
    }, SIMPLIFY = FALSE)
    
    # cbind
    variates <- Reduce(x = variates, f = cbind)
    ind.names <- .get.character.vector(ind.names, rownames(variates))
    ## axes labels
    xylabels <- paste0('Dimension ', comp)
    
    ## colours
    if (is.null(group) && is(object = object, 'DA'))
    {
        group <- object$Y
    }
    ## ouputs col.per.group and group based on possible inputs
    col.group <-
        .get.cols.and.group(
            col.per.group = col.per.group,
            group = group,
            col = col,
            object = object,
            n_ind = nrow(variates)
        )
    group <- col.group$group
    if(!missing(col.per.group)){
      col.per.group <- col.group$col.per.group}
    pch.legend <- TRUE
    
    if (is(object, 'DA') && is.null(pch))
    {
        pch <- seq_len(length(blocks)+1)[-8][seq_along(blocks)] ## keep pch = 8 for centroids
        names(pch) <- blocks
    }
    
    if (length(col.per.group) == 1)
    {
        legend <- FALSE
    }
    
    if (is.null(pch))
    {
        pch <- 1
    } else if (length(pch) == length(blocks))
    {
        if (names(pch) %!=% blocks)
        {
            stop("'pch' must be either a single value, or a vector whose names are the block names ",
                 "and whose values are the plot characters") 
        }
    } else if ( length(pch) != 1 )
    {
        stop("'pch' must be either a single value, or a vector whose names are the block names ",
             "and whose values are the plot characters") 
    }
    if (length(pch) == 1)
    {
        pch <- rep(pch, length(blocks))
        names(pch) <- blocks
        pch.legend <- FALSE
    }
        
    if (is.multiblock)
    {
        if (any(pch == 8)) 
            stop("'pch=8' is a reserved value. Please use another values")
        pch <- c(pch, 'centroid' = 8)
    }
    
    ## ------------- outline
    variates$group <- group
    p <- ggplot(variates) + 
        theme_classic() +
        labs(x = xylabels[1], y = xylabels[2])# +

    if (is(object, 'DA'))
    {
        legend.title <- .change_if_null(legend.title, 'Y')
    }
    legend.title <- .change_if_null(legend.title, as.character(as.list(match.call())['group']))
    geom_point_blocks <- if(is.multiblock) c('centroid', blocks) else blocks
    for (block in geom_point_blocks)
    {
        x <- paste0('x_', block)
        y <- paste0('y_', block)
        df <- variates[,c(x, y, 'group')]
        df$pch <- block
        # Convert variable names from strings to symbols
        p <- p + geom_point(data = df, aes(x = !!sym(x), 
                                               y = !!sym(y),
                                               col = !!sym('group'),
                                               shape = !!sym('pch')),
                                size = pch.size)
    }
    p <- p + scale_color_manual(values = col.per.group, guide = if (isTRUE(legend)) guide_legend(title = legend.title) else NULL)
    p <- p +  scale_shape_manual(values = pch, guide = if (pch.legend) guide_legend(title = 'Block') else NULL)
    
    ## ind.names
    
    if (!is.null(ind.names))
    {
        if (is.multiblock) {
            x_label <- 'x_centroid'
            y_label <- 'y_centroid'
        } else {
            if (ind.names.position == 'start') {
                x_label <- 'x_Y'
                y_label <- 'y_Y'
            } else {
                x_label <- 'x_X'
                y_label <- 'y_X'
            }
        }
        
        
        p <- p + 
            geom_text(mapping = aes(x = !!sym(x_label),
                                    y = !!sym(y_label),
                                    label = !!sym('ind.names')
            ), 
            size = ind.names.size,
            show.legend = FALSE, vjust=-0.1, hjust=-0.1)
        
    }
    p
    ## arrows
    arrow.offset <- 0.04*pch.size/2
    if (is.multiblock) {
        for (block in blocks)
        {
            cols <- c('x_centroid', 'y_centroid', paste0(c('x_', 'y_'), block), 'group')
            df <- variates[,cols]
            colnames(df) <- c('xs', 'ys', 'xe', 'ye', 'group')
            xs <- ys <- xe <- ye <- NULL ## avoid check warnings
            p <- p + geom_segment(data = df,
                                  aes(
                                      x = xs + (xe - xs) * arrow.offset,
                                      y = ys + (ye - ys) * arrow.offset,
                                      xend = xe -  (xe - xs) * arrow.offset,
                                      yend = ye - (ye - ys) * arrow.offset,
                                      col = group
                                  ),
                                  alpha = arrow.alpha,
                                  arrow = arrow(type='closed', angle = 30, length = unit(arrow.length, "cm"), ends = 'last'),
                                  linewidth = arrow.size,
                                  show.legend = FALSE
            )
        }
    } else {
        x_X <- y_X <- x_Y <- y_Y <- NULL ## avoid check warnings
        p <- p + geom_segment(
            aes(
                x = x_X + (x_Y - x_X)*arrow.offset,
                y = y_X+ (y_Y - y_X)*arrow.offset,
                xend = x_Y - (x_Y - x_X)*arrow.offset,
                yend = y_Y - (y_Y - y_X)*arrow.offset,
                col = group
            ),
            alpha = arrow.alpha,
            arrow = arrow(type='closed', angle = 30, length = unit(arrow.length, "cm"), ends = 'last'),
            linewidth = arrow.size,
            show.legend = FALSE
        ) 
    } 
    
    ## no axes values
    p <- p + scale_y_continuous(breaks = NULL) +
        scale_x_continuous(breaks = NULL) 
    return(p)
}

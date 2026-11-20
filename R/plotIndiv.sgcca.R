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
#-- Includes plotIndiv rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#
#' @rdname plotIndiv
#' @method plotIndiv sgcca
#' @export
plotIndiv.sgcca <-
  function(object,
           comp = NULL,
           blocks = NULL,
           style = "ggplot2",
           ind.names = TRUE,
           group,
           col,
           ellipse = FALSE,
           ellipse.level = 0.95,
           centroid = FALSE,
           star = FALSE,
           pch,
           title = NULL,
           subtitle,
           legend = FALSE,
           legend.title = "Legend",
           legend.title.pch = "Legend",
           legend.position = "right",
           X.label = NULL,
           Y.label = NULL,
           Z.label = NULL,
           abline = FALSE,
           xlim = NULL,
           ylim = NULL,
           axes.box = "box",
           cex,
           alpha = 0.2,
           point.lwd = 1,
           size.title = rel(2),
           size.subtitle = rel(1.5),
           size.xlabel = rel(1),
           size.ylabel = rel(1),
           size.axis = rel(0.8),
           size.legend = rel(1),
           size.legend.title = rel(1.1),
           ...
           
  )
  {

    # check for inappropriate args
    extra_args <- list(...)
    if ("rep.space" %in% names(extra_args)) {
        warning("'rep.space' is not used for multiblock models, use 'blocks' to specify which data blocks to plot")
    }
    if ("study" %in% names(extra_args) || "layout" %in% names(extra_args)) {
        warning("'study' and 'layout' arguments are only used for MINT models")
    }
    # check for deprecated args
    if ("col.per.group" %in% names(extra_args)) {
        warning("'col.per.group' is deprecated, please use 'col' to specify colours for each group")
    }      
    if ("pch.levels" %in% names(extra_args)) {
         warning("'pch.levels' is deprecated, please use 'pch' to specify point types")
    }
    if ("cols" %in% names(extra_args)) {
         warning("'cols' is not a valid argument, did you mean 'col' ?")
    }
    

    plot_parameters = list(
      size.title = size.title,
      size.subtitle = size.subtitle,
      size.xlabel = size.xlabel,
      size.ylabel = size.ylabel,
      size.axis = size.axis,
      size.legend = size.legend,
      size.legend.title = size.legend.title,
      legend.title = legend.title,
      legend.title.pch = legend.title.pch,
      legend.position = legend.position,
      point.lwd = point.lwd,
      alpha = alpha
    )
    
    if (any(class(object)%in%c("mint.block.pls", "mint.block.spls", "mint.block.plsda", "mint.block.splsda")))
      stop("No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.")
    
    #-- rep.space
    rep.space = "multi" # rep.space is not used afterwards, put to "multi" to plot all blocks
    layout = NULL
    input_average_blocks <- NULL
    valid_average_blocks <- c('average', 'weighted.average')
    
    if (is.null(blocks))
    {
      blocks = names(object$X)#names$blocks
    } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <=  length(object$names$blocks)) {
      blocks = object$names$blocks[blocks]
    } else if (is.character(blocks)) {
      if (any(grepl(blocks, pattern = "average", ignore.case = TRUE))) {
        input_average_blocks <- match.arg(blocks, choices = valid_average_blocks, several.ok = TRUE)
        supported.classes <- c("sgcca", "rgcca", "sgccda")
        if (! any(supported.classes %in% class(object)) ){
          stop(sprintf("average plots are supported for objects of classes: %s", paste(supported.classes, collapse = ", ")), call. = FALSE)
        }
        object <- .add_average_blocks(block_object = object, average_blocks = input_average_blocks)
      }
      invalid_blocks <- setdiff(blocks, object$names$blocks)
      if (length(invalid_blocks) > 0  ) {
        valid_blocks <- object$names$blocks
        stop(sprintf("Block(s) not found: %s. Blocks must be a selection of: %s, average, or weighted.average", paste(invalid_blocks , collapse = ', '), paste(valid_blocks, collapse = ', ')))
      }
    } else {
      stop("Incorrect value for 'blocks'", call. = FALSE)
    }
    #object$variates = object$variates[names(object$variates) %in% blocks] # reduce the variate to the 'blocks' we are looking at
    object$variates = object$variates[match(blocks, names(object$variates))] # reduce the variate to the 'blocks' we are looking at
    if (any(object$ncomp[blocks] ==  1))
      stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "), "' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
    ncomp = object$ncomp[blocks]
    
    
    if(length(blocks)!= length(unique(blocks)))
      stop("Duplicate in 'blocks' not allowed")
    
    if(!missing(subtitle))
    {
      if(length(subtitle)!= length(blocks) | length(subtitle)!= length(unique(subtitle)))
        stop("'subtitle' indicates the subtitle of the plot for each 'blocks'; it needs to be the same length as 'blocks' and duplicate are not allowed.")
    }
    
    
    
    #-- check inputs
    check = check.input.plotIndiv(
      object = object,
      comp = comp,
      blocks = blocks,
      ind.names = ind.names,
      style = style,
      ellipse = ellipse,
      ellipse.level = ellipse.level,
      centroid = centroid,
      star = star,
      legend = legend,
      X.label = X.label,
      Y.label = Y.label,
      Z.label = Z.label,
      abline = abline,
      xlim = xlim,
      ylim = ylim,
      alpha = alpha,
      axes.box = axes.box,
      plot_parameters = plot_parameters
    )
    #-- retrieve outputs from the checks
    axes.box = check$axes.box
    comp = check$comp
    xlim = check$xlim
    ylim = check$ylim
    ind.names = check$ind.names
    display.names = check$display.names
    
    
    #-- get the variates
    variate = internal_getVariatesAndLabels(
      object,
      comp,
      blocks = blocks,
      style = style,
      X.label = X.label,
      Y.label = Y.label,
      Z.label = Z.label,
      rep.space = rep.space
    )
    #-- retrieve outputs
    x = variate$x
    y = variate$y
    z = variate$z
    X.label = variate$X.label
    Y.label = variate$Y.label
    Z.label = variate$Z.label
    
    n = nrow(object$X[[1]])
    
    # create data frame df that contains (almost) all the ploting information
    out = shape.input.plotIndiv(
      object = object,
      n = n,
      blocks = blocks,
      x = x,
      y = y,
      z = z,
      ind.names = ind.names,
      group,
      col.per.group = col.per.group,
      style = style,
      study = "global",
      ellipse = ellipse,
      ellipse.level = ellipse.level,
      centroid = centroid,
      star = star,
      title = title,
      xlim = xlim,
      ylim = ylim,
      col = col,
      cex = cex,
      pch = pch,
      display.names = display.names,
      plot_parameters = plot_parameters
    )
    #-- retrieve outputs
    df = out$df
    df.ellipse = out$df.ellipse
    col.per.group = out$col.per.group
    title = out$title
    display.names = out$display.names
    xlim = out$xlim
    ylim = out$ylim
    #missing.col = out$missing.col
    ellipse = out$ellipse
    centroid = out$centroid
    star = out$star
    plot_parameters = out$plot_parameters
    
    if (!is.null(input_average_blocks)) {
      ## ------ drop 'Block: ' from average labels
      df$Block <- as.character(df$Block)
      
      df$Block[df$Block == 'Block: weighted.average'] <- 'average (weighted)'
      df.ellipse$Block[df.ellipse$Block == 'Block: weighted.average'] <- 'average (weighted)'
      
      df$Block[df$Block == 'Block: average'] <- 'average'
      df.ellipse$Block[df.ellipse$Block == 'Block: average'] <- 'average'
      
      df$Block <- factor(df$Block)
      
      ## ------ drop explained variance on axes if only average asked
      if (length(unique(df$Block)) == 1) {
        if (grepl("% expl. var", X.label)) { X.label <- 'variate 1 - average' }
        if (grepl("% expl. var", Y.label)) { Y.label <- 'variate 2 - average' }
        if (style=="3d") { if (grepl("% expl. var", Z.label)) { Z.label <- 'variate 3 - average' } }
      }
    }
    
    # change the levels of df.final$Block to "subtitle"
    if(!missing(subtitle))
    {
      df$Block = factor(df$Block, labels = subtitle)
      if(ellipse)
        df.ellipse$Block = factor(df.ellipse$Block, labels = subtitle)
    }
    
    #call plot module (ggplot2, lattice, graphics, 3d)
    res = internal_graphicModule(
      df = df,
      centroid = centroid,
      col.per.group = col.per.group,
      title = title,
      X.label = X.label,
      Y.label = Y.label,
      Z.label = Z.label,
      xlim = xlim,
      ylim = ylim,
      class.object = class(object),
      display.names = display.names,
      legend = legend,
      abline = abline,
      star = star,
      ellipse = ellipse,
      df.ellipse = df.ellipse,
      style = style,
      layout = layout,
      #missing.col = missing.col,
      axes.box = axes.box,
      plot_parameters = plot_parameters,
      alpha = alpha
    )
    
    return(invisible(list(df = df, df.ellipse = df.ellipse, graph = res)))
    
    
  }

#' @rdname plotIndiv
#' @method plotIndiv rgcca
#' @export
plotIndiv.rgcca <- plotIndiv.sgcca

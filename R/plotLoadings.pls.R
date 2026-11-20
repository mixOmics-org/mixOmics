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

## --------------------- PLS, sPLS, rCC, rGCCA, sGCCA --------------------- ##

#' @rdname plotLoadings
#' @method plotLoadings mixo_pls
#' @export
plotLoadings.mixo_pls <- 
  function(object,
           comp = 1,
           style = 'graphics',
           ndisplay = NULL,
           xlim = NULL,
           layout = NULL,
           col = NULL,
           border = NA,
           name.var = NULL,
           size.name = 0.7,
           title = NULL,
           subtitle,
           size.title = 2,
           size.subtitle = rel(1.5),
           size.axis = 0.7,
           X.label = NULL,
           Y.label = NULL,
           size.labs = 1,
           block,
           ...) {

    # Input checks
    if (!is.numeric(comp) || length(comp) != 1 || comp <= 0)
      stop("'comp' must be a positive integer.")
    if (!style %in% c('graphics', 'ggplot2'))
      stop("'style' must be either 'graphics' or 'ggplot2'.")
    
    if (!is.null(col) && (length(col) != 1 || !col %in% colors()))
      stop("'col' must be a single valid color.")
    if (!is.null(ndisplay) && (!is.numeric(ndisplay) || length(ndisplay) != 1 || ndisplay <= 0))
      stop("'ndisplay' must be a positive integer.")
    
    if (!is.null(title) && !is.character(title))
      stop("'title' must be NULL or a character string.")
    if (!is.null(xlim) && (!is.numeric(xlim) || length(xlim) != 2))
      stop("'xlim' must be a numeric vector of length 2.")
    if (!is.null(X.label) && !is.character(X.label))
      stop("'X.label' must be NULL or a character string.")
    if (!is.null(Y.label) && !is.character(Y.label))
      stop("'Y.label' must be NULL or a character string.")
    
    if (!is.numeric(size.name) || size.name <= 0)
      stop("'size.name' must be a positive numeric value.")
    if (!is.numeric(size.title) || size.title <= 0)
      stop("'size.title' must be a positive numeric value.")
    if (!is.numeric(size.subtitle) || size.subtitle <= 0)
      stop("'size.subtitle' must be a positive numeric value.")
    if (!is.numeric(size.labs) || size.labs <= 0)
      stop("'size.labs' must be a positive numeric value.")
    if (!is.numeric(size.axis) || size.axis <= 0)
      stop("'size.axis' must be a positive numeric value.")
    
    # check for inappropriate args
    extra_args <- list(...)
    if ("name.var.complete" %in% names(extra_args)) {
      warning("'name.var.complete' argument is deprecated")
    }
    name.var.complete <- FALSE
    
    # -- input checks
    check = check.input.plotLoadings(object = object, block = block, subtitle = subtitle, size.name = size.name, title = title, col = col, name.var = name.var, xlim = xlim)
    
    col = check$col
    size.name = check$size.name
    block = check$block
    xlim = check$xlim
    
    # -- layout
    res = layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = block)
    reset.mfrow = res$reset.mfrow
    opar = res$opar
    omar = par("mar") #reset mar at the end
    
    if (length(block) == 1 & !is.null(name.var))
      name.var = list(name.var = name.var)
    
    contrib.df <- list()
    
    # create a list to hold ggplot objects if needed
    plot_list <- list()
    
    for (i in 1 : length(block))
    {
      res = get.loadings.ndisplay(object = object, comp = comp, block = block[i], name.var = name.var[[i]], name.var.complete = name.var.complete, ndisplay = ndisplay)
      X = res$X
      names.block = res$names.block
      colnames.X = res$colnames.X
      value.selected.var = res$value.selected.var
      
      df <- data.frame(importance = value.selected.var) # contribution of the loading

      if (style == 'graphics') {
        # barplot with contributions using base graphics
        if (!is.null(title) & length(block) > 1)
        {
          par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 6, 2))
        } else {
          par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 4, 2))
        }
        
        .plotLoadings_barplot(height = df$importance, col = col, names.arg = colnames.X, cex.name = size.name, border = border, xlim = xlim[i, ],
                              xlab = X.label, ylab = Y.label, cex.lab = size.labs, cex.axis = size.axis)
        
        if ( (length(block) == 1 & is.null(title)) | (length(block) > 1 & missing(subtitle)))
        {
          title(paste0('Loadings on comp ', comp, "\nBlock '", names.block,"'"), line=0, cex.main = size.title)
        } else if (length(block) == 1) {
          title(paste(title), line=0, cex.main = size.title)
        } else if (length(block) > 1 & !missing(subtitle)) {
          title(paste(subtitle[i]), line=0, cex.main = size.subtitle)
        }
        
      } else if (style == 'ggplot2') {
        
          if (missing(subtitle)) {
            # Combine title and block info in one title line when subtitle is missing
            plot.title <- paste0("Loadings on comp ", comp, "\nBlock '", names.block, "'")
          } else {
            # If subtitle is provided for multiple blocks, use it per block.
            plot.title <- subtitle[i]
          }

        df$names <- colnames.X
        
        p <- ggplot(df, aes(x = reorder(names, -abs(as.numeric(importance))), y = importance)) +
          geom_bar(stat = "identity", aes(fill = col), color = border) +
          scale_fill_identity() +  # This ensures the colors are used as-is
          theme_minimal() +
          theme(axis.text.y = element_text(size = size.name * 8),
                axis.text.x = element_text(size = size.axis * 8),
                axis.title.x = element_text(size = size.labs * 8),
                axis.title.y = element_text(size = size.labs * 8),
                plot.title = element_text(face = "bold", hjust = 0.5, size = size.title * 8)) +
          labs(title = plot.title, y = X.label, x = Y.label)
        
        # optionally control x axis
        if (!is.null(xlim)) {p <- p + scale_y_continuous(limits = xlim[1,], expand = c(0,0))}
        
        # flip and store it
        p <- p + coord_flip()
        plot_list[[i]] <- p
      }
      
      contrib.df <- c(contrib.df, list(df))
    }
    
    names(contrib.df) <- block
    
    ## Final additions to plots
    if (style == "graphics" & length(block) > 1 & !is.null(title)){
      title(title, outer=TRUE, line = -2, cex.main = size.title)
      if (reset.mfrow)
        par(opar)
      par(mar = omar) #reset mar
    }
    
    
    if (style == 'ggplot2') {
      # Add overall plot title if set
      # If there is more than one plot, arrange them side by side using gridExtra
      if (length(plot_list) > 1) {
        grid::grid.newpage() # clear previous grids
        if (is.null(layout)) {layout <- c(1, length(plot_list))} # default layout if one isn't specified
        if(is.null(title)){
          gridExtra::grid.arrange(
            grobs = plot_list, 
            layout_matrix = matrix(seq_along(plot_list), nrow = layout[1], ncol = layout[2], byrow = TRUE)
          )
        } else {
          title_grob <- grid::textGrob(title, gp = grid::gpar(fontsize = size.title * 8, fontface = "bold"))
          plot_grobs <- gridExtra::arrangeGrob(
            grobs = plot_list,
            layout_matrix = matrix(seq_along(plot_list), nrow = layout[1], ncol = layout[2], byrow = TRUE)
          )
          combined <- gridExtra::arrangeGrob(title_grob, plot_grobs, ncol = 1, heights = c(0.1, 1))
          grid::grid.draw(combined)
        }
        
        # if just one plot
      } else {
        print(plot_list[[1]])
      }
      
    }
    
    
    # return the contribution matrix
    return(invisible(contrib.df))
  }


#' @rdname plotLoadings
#' @method plotLoadings mixo_spls
#' @export
plotLoadings.mixo_spls <- plotLoadings.mixo_pls

#' @rdname plotLoadings
#' @method plotLoadings rcc
#' @export
plotLoadings.rcc <- plotLoadings.mixo_pls

#' @rdname plotLoadings
#' @method plotLoadings sgcca
#' @export
plotLoadings.sgcca <- plotLoadings.mixo_pls

#' @rdname plotLoadings
#' @method plotLoadings rgcca
#' @export
plotLoadings.rgcca <- plotLoadings.mixo_pls
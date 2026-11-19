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

## ---------------------------- MINT.(s)PLSDA ----------------------------- ##

#' @rdname plotLoadings
#' @method plotLoadings mint.plsda
#' @export
plotLoadings.mint.plsda <- 
    function(object,
             comp = 1,
             style = "graphics",
             ndisplay = NULL,
             xlim = NULL,
             layout = NULL,
             border = NA,
             name.var = NULL,
             size.name = 0.7,
             title = NULL,
             subtitle,
             size.title = 2,
             size.subtitle = 1.7,
             size.axis = 0.7,
             X.label = NULL,
             Y.label = NULL,
             size.labs = 1,
             contrib = NULL,  # choose between 'max" or "min", NULL does not color the barplot
             method = "mean", # choose between 'mean" or "median"
             show.ties = TRUE,
             col.ties = "white",
             legend = TRUE,
             legend.color = NULL,
             legend.title = 'Outcome',
             size.legend = 0.8,
             study = "global",
             ...
    ) {

        ## Check input args
        # Input checks
        if (!is.numeric(comp) || length(comp) != 1 || comp <= 0)
            stop("'comp' must be a positive integer.")
        if (!style %in% c('graphics', 'ggplot2'))
            stop("'style' must be either 'graphics' or 'ggplot2'.")
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
        if ("block" %in% names(extra_args)) {
            warning("'block' argument is not used for mint.plsda or mint.splsda objects")
        }
        block <- "X"
        
        if(any(study == "global"))
        {
            # if study == "global" then we plot the results on the concatenated data
            plotLoadings.mixo_plsda(object = object, 
                                  style = style,
                                  contrib = contrib, 
                                  method = method, 
                                  block = "X", 
                                  comp = comp, 
                                  ndisplay = ndisplay,
                                  size.name = size.name,
                                  size.legend = size.legend,
                                  name.var = name.var,
                                  legend = legend,
                                  legend.color = legend.color,
                                  title = if(!is.null(title)){title}else{paste0('Contribution on comp ', comp, "\n All studies")},
                                  subtitle = subtitle,
                                  legend.title = legend.title,
                                  xlim = xlim,
                                  layout = layout,
                                  size.title = size.title,
                                  size.subtitle = size.subtitle,
                                  border = border,
                                  col.ties = col.ties,
                                  show.ties = show.ties,
                                  size.axis = size.axis,
                                  size.labs = size.labs,
                                  X.label = X.label,
                                  Y.label = Y.label)
            
        } else {
            # if study != "global" then we plot the results on each study
            
            # -- input checks
            check = check.input.plotLoadings(object = object, block = "X", size.name = size.name, size.legend = size.legend,
                                             title = title, col = NULL, name.var = name.var, contrib = contrib)
            
            size.name = check$size.name
            size.legend = check$size.legend
            block = "X"  # Always use block "X"
            
            #study needs to be either: from levels(object$study), numbers from 1:nlevels(study) or "global"
            if (any(!study%in%c(levels(object$study), "global" , "all.partial")))
                stop("'study' must from one of 'object$study', 'global' or 'all.partial', see help file.")
            
            study.init = unique(study)
            # replace "all.partial" by all levels of object$study
            ind.all.partial = which(study.init == "all.partial")
            if (length(ind.all.partial) > 0)
            {
                if (ind.all.partial > 1 & ind.all.partial < length(study.init))
                {
                    # there are things before and after "all.partial"
                    study.init = c(study.init[1:(ind.all.partial-1)], levels(object$study), study.init[(ind.all.partial+1) : length(study.init)])
                } else if (ind.all.partial == 1 & ind.all.partial < length(study.init)) {
                    # there are only things after "all.partial"
                    study.init = c(levels(object$study), study.init[(ind.all.partial+1) : length(study.init)])
                } else if (ind.all.partial > 1 & ind.all.partial == length(study.init)) {
                    # there are things only before "all.partial"
                    study.init = c(study.init[1:(ind.all.partial-1)], levels(object$study))
                } else if (ind.all.partial == 1 & ind.all.partial == length(study.init)) {
                    # there's only "all.partial"
                    study.init = levels(object$study)
                }
            }
            study.init = unique(study.init) #once again cause we added studies if "all.partial"
            study = study.init
            
            if (!missing(subtitle))
            {
                if (length(subtitle)!=length(study))
                    stop("'subtitle' indicates the subtitle of the plot for each study and it needs to be the same length as 'study' (", length(study),"), which includes: ", paste(study, collapse = ", "))
            }
            
            # check xlim, has to be a matrix with number of rows=number of studies, or a vector of two values
            if(length(study) == 1 & !is.null(xlim))
            {
                if(length(xlim) !=2)
                    stop("'xlim' must be a vector of length 2")
                
                xlim = matrix(xlim, nrow = 1)
            }
            
            if(length(study)>1 & !is.null(xlim))
            {
                if(is.matrix(xlim) && ( !nrow(xlim) %in%c(1, length(study))  | ncol(xlim) != 2 ))
                    stop("'xlim' must be a matrix with ",length(study)," rows (length(study)) and 2 columns")
                
                if(is.vector(xlim))
                {
                    if(length(xlim) !=2)
                        stop("'xlim' must be a matrix with ",length(study)," rows (length(study)) and 2 columns")
                    
                    xlim = matrix(xlim, nrow = 1)
                }
                
                if(nrow(xlim) != length(study)) # we complete xlim to have one xlim per block
                    xlim = matrix(rep(xlim, length(study)), nrow = length(study), byrow=TRUE)
            }
            
            # method
            # ----
            if (length(method) !=1 || !method %in% c("mean","median"))
            {
                method = "median"
                warning("'method' should be either 'mean' or 'median', set to 'median' by default")
            }
            
            # get the selected variables on the concatenated data
            res = get.loadings.ndisplay(object = object, comp = comp, block = block, name.var = name.var, name.var.complete = name.var.complete, ndisplay = ndisplay)
            X = res$X
            colnames.X = res$colnames.X
            name.selected.var = res$name.selected.var
            value.selected.var = res$value.selected.var
            
            # -- layout
            res = layout.plotLoadings(layout = layout, plot = TRUE, legend = legend, block = study.init)
            reset.mfrow = res$reset.mfrow
            opar = res$opar
            omar = par("mar") #reset mar at the end
            
            # Set up layout for multiple plots
            if (length(study.init) > 1) {
                if (is.null(layout)) {
                    # Default layout: arrange in a grid that's as square as possible
                    n = length(study.init)
                    layout = c(ceiling(sqrt(n)), ceiling(n/ceiling(sqrt(n))))
                }
                
                if (style == "graphics") {
                    # Set up the plotting area
                    par(mfrow = layout)
                    # Adjust margins for better spacing
                    par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 4, 2))
                }
            }
            
            # swap loadings partial for loadings
            object$loadings.global = object$loadings
            if (!block %in% names(object$loadings.partial)) {
                stop("Block '", block, "' not found in object$loadings.partial. Available blocks: ", 
                     paste(names(object$loadings.partial), collapse = ", "))
            }
            object$loadings = object$loadings.partial[[block]]
            object$names$block = levels(object$study)
            
            X.study = study_split(X, study = object$study)
            Y = object$Y #v6: all $Y are factors for DA methods
            Y.study = study_split(Y, study = object$study)
            
            df.final = list()
            plot_list = list() # to store ggplot objects if style is ggplot2
            
            for (i in 1 : length(study.init))
            {
                value.selected.var = object$loadings.partial[[block]][[study.init[i]]][, comp] [name.selected.var]
                
                #legend.color
                #-----
                if (!is.null(legend.color) & (length(legend.color) != nlevels(Y)))
                {
                    warning('legend.color must be the same length than the number of group, by default set to default colors')
                    legend.color = color.mixo(1:10)  # by default set to the colors in color.mixo (10 colors)
                }
                if (is.null(legend.color))
                    legend.color = color.mixo(1:10)[1:nlevels(Y)] # by default set to the colors in color.mixo (10 colors)
                
                if (col.ties%in%legend.color[1:nlevels(Y)])
                    stop("'col.ties' should not be in 'legend.color'")
                
                if(!is.null(contrib))
                {
                    df = get.contrib.df(Y = factor(Y.study[[study.init[i]]]), X = X.study[[study.init[i]]], method = method, contrib = contrib, value.selected.var = value.selected.var, colnames.X = colnames.X, name.selected.var = name.selected.var, legend.color = legend.color, col.ties = col.ties)
                    # when working with sparse counts in particular and using the median to measure contribution
                    # ties to determine the contribution of a variable may happen, in that case remove them, otherwise they are showns as blank
                    if (show.ties == FALSE)
                    {
                        df = df[!df$color %in% col.ties, ]
                        colnames.X = rownames(df)
                    }
                    
                } else {
                    # if contrib is NULL, then we plot the loadings without colors
                    df = data.frame(importance = value.selected.var, color = "white", stringsAsFactors = FALSE) # contribution of the loading
                    border = TRUE
                }
                
                if (style == "graphics") {
                    #display barplot with names of variables
                    #added condition if all we need is the contribution stats
                    if (!is.null(title) & length(study.init) > 1)
                    {
                        par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 6, 2))
                    } else {
                        par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/2), 4, 2))
                    }
                    
                    .plotLoadings_barplot(height = df$importance, 
                                        col = df$color, 
                                        names.arg = colnames.X, 
                                        cex.name = size.name, 
                                        border = border, 
                                        xlim = xlim[i, ],
                                        xlab = X.label, ylab = Y.label, cex.lab = size.labs, cex.axis = size.axis)
                    
                    if ( length(study.init) == 1 & is.null(title) )
                    {
                        title(paste0('Contribution on comp ', comp, "\nStudy '", study.init[i],"'"), line=0, cex.main = size.title)
                    } else if (length(study.init) == 1) {
                        title(paste(title), line=0, cex.main= size.title)
                    } else if ((length(study.init) > 1 & missing(subtitle))) {
                        title(paste0('Contribution on comp ', comp, "\nStudy '", study.init[i],"'"), line=0, cex.main = size.subtitle)
                    } else if (length(study.init) > 1 & !missing(subtitle)) {
                        title(paste(subtitle[i]), line=0, cex.main = size.subtitle)
                    }
                    
                    if (legend && !is.null(contrib))
                    {
                        par(mar = c(5, 0, 4, 3) + 0.1)
                        plot(1,1, type = "n", axes = FALSE, ann = FALSE)
                        legend(0.8, 1, col = legend.color[1:nlevels(Y)], legend = levels(Y), pch = 19,
                               title = paste(legend.title),
                               cex = size.legend)
                    }
                    
                } else if (style == "ggplot2") {
                    # Create ggplot version
                    df$names <- colnames.X
                    
                    # Create the base plot
                    p <- ggplot(df, aes(x = reorder(names, -abs(importance)), y = importance)) +
                        geom_bar(stat = "identity", aes(fill = color), color = border) +
                        theme_minimal() +
                        theme(axis.text.y = element_text(size = size.name * 8),
                              axis.text.x = element_text(size = size.axis * 8),
                              axis.title.x = element_text(size = size.labs * 8),
                              axis.title.y = element_text(size = size.labs * 8),
                              plot.title = element_text(face = "bold", hjust = 0.5, size = size.title * 8)) +
                        labs(title = if(length(study.init) == 1 & is.null(title)) {
                            paste0('Contribution on comp ', comp, "\nStudy '", study.init[i],"'")
                        } else if(length(study.init) == 1) {
                            title
                        } else if(length(study.init) > 1 & missing(subtitle)) {
                            paste0('Contribution on comp ', comp, "\nStudy '", study.init[i],"'")
                        } else {
                            subtitle[i]
                        },
                        y = X.label, x = Y.label)
                    
                    # Control x axis limits if specified
                    if (!is.null(xlim)) {
                        p <- p + scale_y_continuous(limits = xlim[i,], expand = c(0,0))
                    }
                    
                    # Flip coordinates for horizontal bar plot
                    p <- p + coord_flip()
                    
                    # Add legend if needed
                    if (legend && !is.null(contrib)) {
                        # Create a named vector for the legend
                        legend_values <- c(legend.color[1:nlevels(Y)], setNames(col.ties, "ties"))
                        legend_labels <- c(levels(Y), "ties")
                        
                        p <- p + 
                            scale_fill_identity(guide = "legend",
                                              breaks = legend_values,
                                              labels = legend_labels,
                                              name = legend.title)
                    } else {
                        p <- p + scale_fill_identity()
                    }
                    
                    plot_list[[i]] <- p
                }
                
                df.final[[i]] = df
            }
            names(df.final) = study.init
            
            if (style == "graphics") {
                if (length(study.init) > 1 & !is.null(title))
                    title(title, outer=TRUE, line = -2, cex.main = size.title)
                
                if (reset.mfrow)
                    par(opar)
                
                par(mar = omar) #reset mar
            } else if (style == "ggplot2") {
                # Add overall plot title if set
                # If there is more than one plot, arrange them side by side using gridExtra
                if (length(plot_list) > 1) {
                    grid::grid.newpage() # clear previous grids
                    if (is.null(layout)) {
                        # Default layout: arrange in a grid that's as square as possible
                        n = length(plot_list)
                        layout = c(ceiling(sqrt(n)), ceiling(n/ceiling(sqrt(n))))
                    }
                    
                    if(is.null(title)){
                        gridExtra::grid.arrange(
                            grobs = plot_list, 
                            layout_matrix = matrix(seq_along(plot_list), nrow = layout[1], ncol = layout[2], byrow = TRUE),
                            widths = rep(1, layout[2]),  # Equal widths
                            heights = rep(1, layout[1])  # Equal heights
                        )
                    } else {
                        title_grob <- grid::textGrob(title, gp = grid::gpar(fontsize = size.title * 8, fontface = "bold"))
                        plot_grobs <- gridExtra::arrangeGrob(
                            grobs = plot_list,
                            layout_matrix = matrix(seq_along(plot_list), nrow = layout[1], ncol = layout[2], byrow = TRUE),
                            widths = rep(1, layout[2]),  # Equal widths
                            heights = rep(1, layout[1])  # Equal heights
                        )
                        combined <- gridExtra::arrangeGrob(title_grob, plot_grobs, ncol = 1, heights = c(0.1, 1))
                        grid::grid.draw(combined)
                    }
                } else {
                    print(plot_list[[1]])
                }
            }
            
            return(invisible(df.final))
        }
    }

#' @rdname plotLoadings
#' @method plotLoadings mint.splsda
#' @export
plotLoadings.mint.splsda <- plotLoadings.mint.plsda 
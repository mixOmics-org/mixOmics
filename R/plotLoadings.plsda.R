## ----------------------- PLS-DA, SPLS-DA, sGCCDA ------------------------ ##

#' @rdname plotLoadings
#' @method plotLoadings mixo_plsda
#' @export
plotLoadings.mixo_plsda <- 
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
             size.subtitle = 1.6,
             size.axis = 0.7,
             X.label = NULL,
             Y.label = NULL,
             size.labs = 1,
             block, # only for sgccda object
             contrib = NULL,  # choose between 'max" or "min", NULL does not color the barplot
             method = "mean", # choose between 'mean" or "median"
             show.ties = TRUE,
             col.ties = "white",
             legend = TRUE,
             legend.color = NULL,
             legend.title = 'Outcome',
             size.legend = 0.8,
             ...
    ) {
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
      
        # -- input checks
        check = check.input.plotLoadings(object = object, block = block, subtitle = subtitle, size.name = size.name, size.legend = size.legend,
                                         title = title, col = NULL, contrib = contrib, name.var = name.var, xlim = xlim)
        
        size.name = check$size.name
        size.legend = check$size.legend
        block = check$block
        xlim = check$xlim
        
        # contrib
        # --
        
        # if contrib is NULL, then we switch to the classical plotLoadings (without contribution/colors)
        if(is.null(contrib))
        {
            xlim <- xlim[1,]
            plotLoadings.mixo_pls(object = object, block = block, comp = comp, ndisplay = ndisplay,
                                      size.name = size.name, name.var = name.var,
                                      title = title, subtitle = subtitle, size.title = size.title, size.subtitle = size.subtitle,
                                      xlim = xlim, layout = layout, size.axis = size.axis,
                                      X.label = X.label, Y.label = Y.label, size.labs = size.labs,
                                      border = TRUE, col = "white", style = style)

        } else {
            
            # -- layout
            res = layout.plotLoadings(layout = layout, plot = plot, legend = legend, block = block)
            reset.mfrow = res$reset.mfrow
            opar = res$opar
            omar = par("mar") #reset mar at the end
            
            # method
            # ----
            if (length(method) !=1 || !method %in% c("mean","median"))
            {
                method = "median"
                warning("'method' should be either 'mean' or 'median', set to 'median' by default")
            }
            
            if (length(block) == 1 & !is.null(name.var))
                name.var = list(name.var = name.var)
            
            
            contrib.df <- list()
            plot_list <- list() # to store ggplot objects if style is ggplot2
            
            for (i in 1 : length(block))
            {
                res = get.loadings.ndisplay(object = object, comp = comp, block = block[i], name.var = name.var[[i]], name.var.complete = name.var.complete, ndisplay = ndisplay)
                X = res$X
                names.block = res$names.block
                colnames.X = res$colnames.X
                name.selected.var = res$name.selected.var
                value.selected.var = res$value.selected.var
                
                Y = object$Y #v6: all $Y are factors for DA methods
                
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
                
                #  determine the colors/groups matching max contribution
                df = get.contrib.df(Y = Y, X = X, method = method, contrib = contrib, value.selected.var = value.selected.var, colnames.X = colnames.X, name.selected.var = name.selected.var, legend.color = legend.color, col.ties = col.ties)
                
                # when working with sparse counts in particular and using the median to measure contribution
                # ties to determine the contribution of a variable may happen, in that case remove them, otherwise they are showns as blank
                if (show.ties == FALSE)
                {
                    df = df[!df$color %in% col.ties, ]
                    colnames.X = rownames(df)
                }
                
                if (style == "graphics") {
                    # display barplot with names of variables
                    if (!is.null(title) & length(block) > 1)
                    {
                        par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 6, 2))
                    } else {
                        par(mar = c(4, max(7, max(sapply(colnames.X, nchar),na.rm = TRUE)/3), 4, 2))
                    }
                    
                    .plotLoadings_barplot(height = df$importance, col = df$color, names.arg = colnames.X, 
                    cex.name = size.name, border = border, xlim = xlim[i, ], 
                    xlab = X.label, ylab = Y.label, cex.lab = size.labs, cex.axis = size.axis)
                    
                    if ( length(block) == 1 & is.null(title) )
                    {
                        title(paste0('Contribution on comp ', comp), line=0, cex.main = size.title)
                    } else if (length(block) == 1) {
                        title(paste(title), line=1, cex.main = size.title)
                    } else if ((length(block) > 1 & missing(subtitle))) {
                        title(paste0('Contribution on comp ', comp, "\nBlock '", names.block,"'"), line=0, cex.main = size.subtitle)
                    } else if (length(block) > 1 & !missing(subtitle)) {
                        title(paste(subtitle[i]), line=1, cex.main = size.subtitle)
                    }
                    
                    if (legend)
                    {
                        par(mar = c(5, 0, 4, 3) + 0.1)
                        plot(1,1, type = "n", axes = FALSE, ann = FALSE)
                        legend(0.8, 1.1, col = legend.color[1:nlevels(Y)], legend = levels(Y), pch = 19,
                               title = paste(legend.title),
                               cex = size.legend)
                    }
                } else if (style == "ggplot2") {
                    # Create ggplot version
                    df$names <- colnames.X
                    
                    # Create the base plot
                    p <- ggplot(df, aes(x = reorder(names, -abs(importance)), y = importance)) +
                        geom_bar(stat = "identity", aes(fill = color), color = border) +
                        scale_fill_identity() +  # This ensures the colors are used as-is
                        theme_minimal() +
                        theme(axis.text.y = element_text(size = size.name * 8),
                              axis.text.x = element_text(size = size.axis * 8),
                              axis.title.x = element_text(size = size.labs * 8),
                              axis.title.y = element_text(size = size.labs * 8),
                              plot.title = element_text(face = "bold", hjust = 0.5, size = size.title * 8)) +
                        labs(title = if(length(block) == 1 & is.null(title)) {
                            paste0('Contribution on comp ', comp)
                        } else if(length(block) == 1) {
                            title
                        } else if(length(block) > 1 & missing(subtitle)) {
                            paste0('Contribution on comp ', comp, "\nBlock '", names.block,"'")
                        } else {
                            subtitle[i]
                        },
                        y = X.label, x = Y.label)
                    
                    # Add legend if requested
                    if (legend) {
                        # Create invisible points for the legend
                        legend_data <- data.frame(
                            group = levels(Y),
                            color = legend.color[1:nlevels(Y)]
                        )
                        p <- p + 
                            geom_point(data = legend_data, aes(x = 0, y = 0, color = group)) +
                            scale_color_manual(name = legend.title,
                                             values = legend.color[1:nlevels(Y)],
                                             labels = levels(Y))
                    }
                    
                    # Control x axis limits if specified
                    if (!is.null(xlim)) {
                        p <- p + scale_y_continuous(limits = xlim[i,], expand = c(0,0))
                    }
                    
                    # Flip coordinates for horizontal bar plot
                    p <- p + coord_flip()
                    
                    plot_list[[i]] <- p
                }

                contrib.df <- c(contrib.df, list(df))
            }
            
            names(contrib.df) <- block
            
            if (style == "graphics") {
                # legend
                if (length(block) > 1 & !is.null(title))
                    title(title, outer=TRUE, line = -2, cex.main = size.title)
                
                if (reset.mfrow)
                    par(opar)#par(mfrow = omfrow)
                
                par(mar = omar) #reset mar
            } else if (style == "ggplot2") {
                # Arrange multiple plots if needed
                if (length(plot_list) > 1) {
                    grid::grid.newpage()
                    if (is.null(layout)) {
                        layout <- c(1, length(plot_list))
                    }
                    if(is.null(title)) {
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
                } else {
                    print(plot_list[[1]])
                }
            }
            
            # return the contribution matrix
            return(invisible(contrib.df)) # df
        }# end contrib missing
    }


#' @rdname plotLoadings
#' @method plotLoadings mixo_splsda
#' @export
plotLoadings.mixo_splsda <- plotLoadings.mixo_plsda

#' @rdname plotLoadings
#' @method plotLoadings sgccda
#' @export
plotLoadings.sgccda <- plotLoadings.mixo_plsda
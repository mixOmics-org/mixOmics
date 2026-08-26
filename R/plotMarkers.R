#' Plot the values for multivariate markers in block analyses
#'
#' Plots the standardised values (after centring and/or scaling) for the
#' selected variables for a given block on a given component. Only applies to
#' \code{block.splsda} or \code{block.spls}.
#' 
#' @param object An object of class \code{block.splsda} or \code{block.spls}
#' @param block Name or index of the block to use
#' @param comp Integer, the component to use
#' @param markers Character or integer, only include these markers. If integer,
#'   the top 'markers' features are shown
#' @param group Factor, the grouping variable (only required for
#'   \code{block.spls} objects)
#' @template arg/col.per.group
#' @param global Logical indicating whether to show the global plots (TRUE) or
#'   segregate by feature (FALSE). Only available when \code{object$scale=TRUE}
#' @param title The plot title
#' @param violin (if global = FALSE) Logical indicating whether violin plots should also be shown
#' @param boxplot.width Numeric, adjusts the width of the box plots
#' @param violin.width Numeric, adjusts the width of the violin plots
#'
#' @return A ggplot object
#' @seealso \code{\link{plotLoadings}}, \code{\link{block.splsda}}, \code{\link{block.spls}}
#' @export
#'
#' @examples
#' # see ?block.splsda and ?block.spls
plotMarkers <-
    function(object,
             block,
             markers = NULL,
             comp = 1,
             group = NULL,
             col.per.group = NULL,
             global = FALSE,
             title = NULL,
             violin = TRUE,
             boxplot.width = NULL,
             violin.width = 0.9
             )
    {
        boxplot.width <- .change_if_null(boxplot.width, ifelse(violin, 0.1, 0.5))
        if (object$scale == FALSE && isTRUE(global))
        {
            cat("\n'global' plots are only available with 'scale=TRUE'. Using 'global=FALSE'\n")
            global <- FALSE
        }
        
        blocks <- names(object$X)
        if (is.numeric(block))
        {
            blocks <- blocks[block]
        }
        if (!block %in% blocks)
            stop(message = sprintf("block must be one of: %s", paste0(blocks, collapse = ', ')))
        df <- data.frame(object$X[[block]], check.names = FALSE)
        
        ## group factor
        group <- .get.group(group, object, n_ind = nrow(df))
        
        col.group <- .get.cols.and.group(col.per.group = col.per.group, 
                                         group = group)
        group <- col.group$group
        col.per.group <- col.group$col.per.group
        vars <- selectVar(object, block=block, comp=comp)[[1]]$value
        vars <- data.frame(feature = rownames(vars), loading = vars$value.var)
        var.names <- vars$feature
        
        df <- df[,var.names]
        df$group <- group
        df <-
            melt(df,
                 id.vars = 'group',
                 variable.name = 'feature',
                 value.name = 'value')
        df <- merge(df, vars)
        df$feature <- factor(df$feature, levels = var.names, ordered = TRUE)
        df$sign <- ifelse(df$loading > 0, 'Positive Loading', 'Negative Loading')
        # to show +ves on top
        df$sign <- factor(df$sign, levels = c('Positive Loading', 'Negative Loading'), ordered = TRUE)
        if (global == TRUE)
        {
            feature <- value <- NULL
            #' @importFrom dplyr group_by summarise
            df <- group_by(.data = df, feature, group, sign)
            df <- summarise(.data = df, value = median(value, na.rm = TRUE))
        }
        if (!is.null(markers))
        {
            if (is.numeric(markers))
            {
                markers <- selectVar(object = object, comp = comp)[[block]]$name
            }
            if ( is.character(markers))
            {
                invalid.markers <- setdiff(markers, df$feature)
                if (length(invalid.markers) > 0)
                    stop("invalid feature names: ", paste0(invalid.markers, collapse = ", "), call. = FALSE)
                #' @importFrom dplyr filter
                feature <- NULL
                df <- filter(df, feature %in% markers)
            }
        }
        if (global)
        {
            df$feature <- 'plotMarkers'
        }
        
        df$facet <- paste0(df$feature, "\n", df$value)
        if (is.null(title))
            title <- sprintf("Block: %s | Component: %s", block, comp)
        y.title <- sprintf("%svalue%s",
                        ifelse(isTRUE(global), "median ", ""),
                        ifelse(isTRUE(object$scale), " (standardised)", "")
                        )
        p <- ggplot(df, aes(.data[["group"]], .data[["value"]], fill = .data[["group"]]))
        
        if (violin)
        {
            p <- p + geom_violin(adjust = violin.width, scale = "count")
        }
        p <- p +
            geom_boxplot(width=boxplot.width) + 
            scale_fill_manual(values = col.per.group) + 
            theme_classic() + 
            labs(x='', 
                 y=y.title, 
                 title = title) + 
            theme(legend.position = 'none', 
                  axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), 
                  plot.title = element_text(hjust = 0.5, colour = "grey40"), 
                  strip.background = element_rect(colour="black", fill="grey80"))
        if (global)
        {
            p <- p + facet_grid(sign~feature, scales = 'free') 
        } else {
            p <- p + facet_grid(.~feature, scales = 'free') 
        }
        p
    }

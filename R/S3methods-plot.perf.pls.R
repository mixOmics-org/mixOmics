## --------------------------- plot.perf.(s)pls --------------------------- ##
#' Plot for model performance for PLS analyses
#' 
#' Function to plot performance criteria, such as MSEP, RMSEP, \eqn{R^2},
#' \eqn{Q^2} for s/PLS methods as a function of the number of components.
#' 
#' \code{plot.perf} creates one plot for each response variable in the model,
#' laid out in a multi-panel display. See ?perf for examples.
#' 
#' @param x an \code{perf.pls} object.
#' @param criterion character string. What type of validation criterion to plot
#'   for \code{pls} or \code{spls}. One of \code{"MSEP"}, \code{"RMSEP"},
#'   \code{"R2"} or \code{"Q2"}. More measures available for pls2 methods. See
#'   \code{\link{perf}}.
#' @param xlab,ylab titles for \eqn{x} and \eqn{y} axes.  Typically character
#' strings, but can be expressions (e.g., \code{expression(R^2)}).
#' @param LimQ2 numeric value. Signification limit for the components in the
#' model. Default is \code{LimQ2 = 0.0975}.
#' @param LimQ2.col character string specifying the color for the \code{LimQ2}
#' line to be plotted. If \code{"none"} the line will not be plotted.
#' @param sd If 'nrepeat' was used in the call to 'perf', error bar shows the
#' standard deviation if sd=TRUE. For mint objects sd is set to FALSE as the
#' number of repeats is 1.
#' @param pch Plot character to use.
#' @param pch.size Plot character size to use.
#' @param cex A numeric which adjusts the font size in the plot.
#' @param col Character. Colour to be used for data points.
#' @param title Character, Plot title. Not used by PLS2 feature-wise measure
#'   plots.
#' @param ... Not used.
#' @return none
#' @author Al J Abadi
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{perf}}.
#' @references
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate hplot
#' @name plot.perf.pls
NULL
## --------------------------- plot.perf.(s)pls --------------------------- ##
#' @method plot perf.pls.mthd
#' @rdname plot.perf.pls
#' @export
plot.perf.pls.mthd <-
    function (x,
              criterion = "MSEP", # TODO homogenise the use of criterion/measure in methods and functions
              xlab = "Number of components",
              ylab = NULL,
              LimQ2 = 0.0975,
              LimQ2.col = 'grey30',
              sd = NULL,
              pch = 1,
              pch.size = 3,
              cex = 1.2,
              col = color.mixo(1),
              title = NULL, ## not used by feature-wise measures
              ...
    )
    {
        ## fix check issues
        comp <- lower <- upper <- NULL
        if (length(criterion) > 1 || !(criterion %in% names(x$measures) ))
            stop("'criterion' must be one of names(", 
                 deparse(substitute(x)),"$measures): ", "\n", 
                 paste(names(x$measures), collapse = ', '), call. = FALSE)
        df = x$measures[[criterion]]$summary
        df <- df[,c('feature', 'comp', 'mean', 'sd')]
        repeated <- all(!is.na(df$sd))
        sd <- .change_if_null(sd, default = ifelse(repeated, TRUE, FALSE))
        if (is.null(ylab)) {
            ylab <- criterion
            if (ylab == 'R2')
                ylab <-  expression(R^~2)
            else if (ylab == 'Q2')
                ylab <- expression(Q^~2)
            else if (ylab == 'Q2.total')
                ylab <- bquote(Q[total]^~2)
            else if(grepl('pred', ylab)) ## subscript
            { ## {RSS/cor}{.}{t/u}{pred}
                measure <- substr(ylab, start = 0, stop = 3)
                v <- substr(ylab, start = 5, stop = 5)
                ylab <- bquote(.(measure)[pred]^.(v))
                
            }
        }
        features <- as.character(unique(df$feature))
        if (length(features) == 1)
        {
            ## will be used as subtitle by facet_wrap
            title <- .check_character(arg = title, len = 1, default = '')
            df$feature <- factor(title)
        }
        
        # TODO Q2 with pls2
        df$upper <- df$mean + df$sd
        df$lower <- df$mean - df$sd
        
        p <- ggplot(df, aes(comp, mean)) + theme_bw(base_size = as.integer(cex*10)) +
            labs(x = xlab, y = ylab) #+ mixo_gg.theme(cex = cex, x.angle = 0, subtitle.cex = 1.2*cex)
        ## discrete x
        p <- p + theme(panel.grid.minor = element_blank())
        
        if (length(features) == 1 && title == '')
            p <- p +  theme(strip.background = element_blank())
        if (grepl('Q2', criterion))
        {
            # min.y <- ifelse(repeated, min(df$lower), min(df$mean))
            # max.y <- ifelse(repeated, max(df$upper), max(df$mean))
            # y.breaks <- sort(c(round(seq(min.y, max.y, length.out=5), 2), LimQ2))
            p <- p + geom_hline(yintercept = LimQ2, col = LimQ2.col, linewidth  = 1.3) #+
            # geom_text(x = max(df$comp), y = 0.0975, label = 'Q2 = 0.0975', hjust=1, vjust=-0.5)
            # scale_y_continuous(breaks = y.breaks)
        }
        
        
        p <- p + geom_point(shape = pch, col = col, size = pch.size)  +
            scale_x_continuous(breaks = as.integer(sort(unique(df$comp))))
        
        if (grepl('^cor', criterion)) ## correlation
            p <- p + ylim(c(min(0, df$lower), max(1, df$upper)))
        if (sd)
            if (any (is.na(df$sd)))
                message("error bars cannot be calculated as nrepeat < 3")
        else
            p <- p + geom_errorbar(aes(ymin = lower, ymax = upper),
                                   width = 0.08,
                                   col = color.mixo(2))
        
        p <- p + facet_wrap(.~feature)
        
        p
    }


#' @method plot perf.spls.mthd
#' @rdname plot.perf.pls
#' @export
plot.perf.spls.mthd <- plot.perf.pls.mthd

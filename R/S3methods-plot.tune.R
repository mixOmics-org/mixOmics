## ------------------------------------------------------------------------ ##
###                               plot.tune                                ###
## ------------------------------------------------------------------------ ##
#' Plot model performance
#' 
#' Function to plot performance criteria, such as classification error rate or
#' correlation of cross-validated components for different models.
#' 
#' \code{plot.tune.splsda} plots the classification error rate or the balanced
#' error rate from x$error.rate, for each component of the model. A lozenge
#' highlights the optimal number of variables on each component.
#' 
#' \code{plot.tune.block.splsda} plots the classification error rate or the
#' balanced error rate from x$error.rate, for each component of the model. The
#' error rate is ordered by increasing value, the yaxis shows the optimal
#' combination of keepX at the top (e.g. `keepX on block 1'_`keepX on block
#' 2'_`keepX on block 3')
#' 
#' \code{plot.tune.spls} plots either the correlation of cross-validated
#' components or the Residual Sum of Square (RSS) values for these components
#' against those from the full model for both \code{t} (X components) and
#' \code{u} (Y components). The optimal number of features chosen are indicated
#' by squares.
#' 
#' If neither of the \code{object$test.keepX} or \code{object$test.keepY} are
#' fixed, a dot plot is produced where a larger size indicates the strength of
#' the measure (higher correlation or lower RSS). Otherwise, the measures are
#' plotted against the number of features selected. In both cases, the colour
#' shows the dispersion of the values across repeated cross validations.
#' 
#' \code{plot.tune.spca} plots the correlation of cross-validated components from
#' the \code{tune.spca} function with respect to the full model.
#' @inheritParams plotIndiv
#' @param x a \code{tune} object. See details for supported objects.
#' @param comp Integer of length 2 denoting the components to plot.
#' @param measure Character. Measure used for plotting a \code{tune.spls} object.
#' One of c('cor', 'RSS').
#' @param optimal If TRUE, highlights the optimal keepX per component
#' @param sd If \eqn{nrepeat >= 3} was used in the call, error bar
#' shows the standard deviation if sd=TRUE. Note that the values might exceeed
#' the valid performance measures (such as [0, 1] for accuracy)
#' @param col character (or symbol) color to be used, possibly vector. One
#' colour per component.
#' @param title Plot title.
#' @param size.range Numeric vector of length 2. Range of sizes used in plot.
#' @param ... Not currently used.
#' @return none
#' @author Kim-Anh Lê Cao, Florian Rohart, Francois Bartolo, Al J Abadi
#' @seealso \code{\link{tune.mint.splsda}}, \code{\link{tune.splsda}},
#'   \code{\link{tune.block.splsda}}, \code{\link{tune.spca}} and
#'   http://www.mixOmics.org for more details.
#' @keywords regression multivariate hplot
#' @name plot.tune
#' @example ./examples/plot.tune-examples.R
NULL
## --------------------------- plot.tune.(s)pls --------------------------- ##
#' @method plot tune.spls
#' @rdname plot.tune
#' @section plot arguments for pls2 tuning:
#' For tune.spls objects where tuning is performed on both X and Y, arguments
#' 'col.low.sd' and 'col.high.sd' can be used to indicate a low and high sd, 
#' respectively. Default to 'blue' & 'red'.
#' @export
plot.tune.spls <-
    function(x, measure = NULL, comp = c(1,2), pch = 16, cex = 1.2, title = NULL, size.range = c(3,10), sd = NULL,...)
    {
        
        ## if measure not given, use object's 'measure.tune' for spls
        if (is.null(measure) & is(x, 'tune.spls') )
        {
            measure <- x$call$measure
        } else {
            measure <- match.arg(measure, c('cor', 'RSS'))    
        }
        df <- x$measure.pred[x$measure.pred$measure == measure,]
        values <- grepl('value', colnames(df))
        df <- df[,!values]
        df$comp <- paste0('comp ', df$comp)
        na.opt <- is.na(df$optimum.keepA)
        if (any(na.opt)) ## for plot
            df$optimum.keepA[na.opt] <- FALSE
        
        col.low.sd <- .change_if_null(arg = list(...)$col.low.sd, default = 'blue')
        col.high.sd <- .change_if_null(arg = list(...)$col.high.sd, default = 'red')
        ## R CMD check
        keepX <- keepY <- NULL
        ggplot_pls2 <- function(df, col.low.sd, col.high.sd, title = NULL) {
            
            p <- ggplot(df, aes(factor(keepX), factor(keepY))) + 
                geom_point(aes(size = .data[["mean"]], col = .data[["sd"]]), shape = pch) +
                scale_color_gradient(low = col.low.sd, 
                                     high = col.high.sd, 
                                     na.value = color.mixo(1))
            
            ## optimal keepX/keepY
            # opt.size.coef <- ifelse(measure == 'cor', 2, 0.00001)
            # df$mean <- df$mean * opt.size.coef #> we'll have cor > 1 in legend
            if (any(!is.na(df$optimum.keepA))) ## to make it possible to plot the unused measure too
                p <- p + geom_point(data = df[df$optimum.keepA,], 
                                    aes(factor(keepX), factor(keepY), size = mean), 
                                    shape = 0, 
                                    col = 'green', 
                                    stroke = 1.3,
                                    show.legend = FALSE)
            
            p <- p + labs(x = 'keepX', y = 'keepY', size = 'mean', col = 'SD') +
                facet_grid(V~comp)
            
            p <- p + scale_size_continuous(range = if (measure == 'RSS') rev(size.range) else size.range)
            
            p <- p + guides(colour = guide_legend(order=2, override.aes = list(size=2)),
                            size = guide_legend(order=1))
            
            list(gg.plot = p, df= df)
        }
        
        ## this should not bee needed at all as tune.pls1 is different
        ## and uses plot.tune.pls1
        ggplot_pls1 <- function(df, ## from .get_ut_df
                                title = NULL, ## title
                                keepA = 'keepX',
                                sd,
                                cex) ## which keepA is not fixed?
            {
            # TODO do we need this one now that we have plot.tune.spls1? Is it used?
            ## fix check issues
            comp <- lower <- upper <- NULL
            
            # if sd is NULL & sd values are present, set it to TRUE
            sd <- .change_if_null(sd, default = !any(is.na(df$sd)))
            sd <- .check_logical(sd)
            if (isTRUE(sd) & any(is.na(df$sd)))
            {
                cat("the model is not repeated > 2 times. setting 'sd' to FALSE")
                sd <- FALSE
            }
                
            if (sd)
            {
                df$lower <- df$mean - df$sd
                df$upper <- df$mean + df$sd
            }
            ## keepX or keepY must be removed before running this
            p <- ggplot(df, aes(x = .data[[keepA]], y = .data[["mean"]], col = .data[["comp"]])) +
                geom_point(shape = pch, size = cex) +
                geom_line(show.legend = FALSE)
            if (sd)
            {
                p <- p + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.04)
            }
            p <- p + labs(x = keepA, y = measure, col = title)
            
            ## optimal
            if (any(!is.na(df$optimum.keepA))) ## to make it possible to plot the unused measure too
                p <- p + geom_point(data = df[df$optimum.keepA,], 
                                    aes(x = .data[[keepA]], y = .data[["mean"]]),
                                    size = cex*1.2,
                                    shape = 0,
                                    col = 'green',
                                    show.legend = FALSE)
                
            p <- p + guides(fill = guide_legend( override.aes = list(size=2)))
            
            list(gg.plot = p, df= df)
        }
        
        if ( length(unique(df$keepY)) > 1 & length(unique(df$keepX)) > 1)
            res <- ggplot_pls2(df, title = title, 
                               col.low.sd = col.low.sd, 
                               col.high.sd = col.high.sd)
        else if ( length(unique(df$keepY)) == 1 & length(unique(df$keepX)) > 1)
        {
            df <- df[df$V == 't',]
            res <- ggplot_pls1(df, title = title, sd = sd, keepA = 'keepX', cex = 2*cex)
        }
           
        else if ( length(unique(df$keepY)) > 1 & length(unique(df$keepX)) == 1)
        {
            df <- df[df$V == 'u',]
            res <- ggplot_pls1(df, title = title, sd = sd, keepA = 'keepY', cex = 2*cex)
        }
        else
            .stop('Unexpected error. Inavlid keepX and/or keepY.')
        
        text.size = as.integer(cex*10)
        
        if (is.null(title))
        {
            title <- sprintf("measure = '%s'", measure)
            if (measure != x$call$measure)
                title <- sprintf("%s (tune.measure = '%s')", title, x$call$measure)
                
        }
        
        res$gg.plot <- res$gg.plot + mixo_gg.theme(cex = cex) +
            labs(title = title)
        
        res$gg.plot
    }

## ------------------------ plot.tune.block.(s)plsda ---------------------- ##
#' @importFrom gridExtra grid.arrange
#' @rdname plot.tune
#' @method plot tune.block.splsda
#' @export
plot.tune.block.splsda =
    function(x, sd = NULL, col, ...)
    {
        
        # R check
        error.sd=NULL
        
        sd = .change_if_null(sd, !is.null(x$error.rate.sd))
        error <- x$error.rate
        if(sd & !is.null(x$error.rate.sd))
        {
            error.rate.sd = x$error.rate.sd
            ylim = range(c(error + error.rate.sd), c(error - error.rate.sd))
        } else {
            error.rate.sd = NULL
            ylim = range(error)
        }
        select.keepX <- x$choice.keepX
        comp.tuned = length(select.keepX[[1]])
        
        if (length(select.keepX) < 10)
        {
            #only 10 colors in color.mixo
            if(missing(col))
                col = color.mixo(seq_len(comp.tuned))
        } else {
            #use color.jet
            if(missing(col))
                col = color.jet(comp.tuned)
        }
        if(length(col) != comp.tuned)
            stop("'col' should be a vector of length ", comp.tuned,".")
        
        legend=NULL
        measure = x$measure
        
        
        if(measure == "overall")
        {
            ylab = "Classification error rate"
        } else if (measure == "BER")
        {
            ylab = "Balanced error rate"
        }
        
        # if(FALSE)
        # {
        #     # not ordered graph
        #     
        #     # creating one dataframe with all the comp
        #     error.plot = data.frame(comp = rep(colnames(error), each = nrow(error)), names = do.call("rbind", as.list(rownames(error))), error = do.call("rbind", as.list(error)), error.sd = do.call("rbind", as.list(error.rate.sd)), color = rep(col, each = nrow(error)))
        #     
        #     #    p = ggplot(error.plot, aes(x=reorder(names, -error), y=error)) +
        #     p = ggplot(error.plot, aes(x=names, y=error)) + 
        #         theme_minimal() +
        #         geom_bar(stat="identity", fill = error.plot$color)
        #     if(sd) p = p + geom_errorbar(aes(ymin=error-error.sd, ymax = error+error.sd), width=0.4)
        #     
        #     p= p +
        #         ylab(ylab)+
        #         xlab("Number of selected features for each block")+
        #         coord_flip()+
        #         facet_grid(~comp,scales='free')
        #     p
        # }
        
        pp=list()
        for(comp in seq_len(comp.tuned))
        {
            # order error per comp
            so = sort(error[,comp], index.return=TRUE, decreasing = TRUE)
            
            error.ordered = so$x
            error.sd.ordered = error.rate.sd[so$ix,comp]
            
            error.plot = data.frame (names = names(error.ordered), error = error.ordered, error.sd = error.sd.ordered, color = col[comp])
            
            ## ggplot
            p = ggplot(error.plot, aes(x=reorder(names, -error), y=error)) +
                theme_classic() +
                geom_bar(stat="identity", fill = error.plot$color)
            if(sd) p = p + geom_errorbar(aes(ymin=error-error.sd, ymax = error+error.sd), width=0.4)
            
            p= p +
                ylab(ylab)+
                xlab("Number of selected features for each block")+
                ggtitle(colnames(error)[comp])+
                coord_flip()
            
            
            if(comp==1)
                p1=p
            if(comp==2)
                p2=p
            #+theme(axis.text.x = element_text(angle = 90, hjust = 1))
            
            pp[[comp]] = p#assign(paste0("p", colnames(error)[comp]), p)
            
        }
        
        do.call("grid.arrange", c(pp, nrow=ceiling(comp.tuned/3)))
        
        
    }

## --------------------------- plot.tune.spca --------------------------- ##
#' @rdname plot.tune
#' @method plot tune.spca
#' @export
plot.tune.spca <-
    function(x, optimal = TRUE, sd = NULL, col=NULL, ...)
    {
        ncomp <- length(x$cor.comp)
        nrepeat <- x$call$nrepeat
        
        if (nrepeat <= 2 & isTRUE(sd))
        {
            cat("nrepeat < 2 so no SD can be calculated.",
                "setting sd to FALSE")
            sd <- FALSE
        } else if (nrepeat > 2 & is.null(sd))
        {
            sd <- TRUE
        }
        
        cors <- mapply(z=x$cor.comp, w=seq_len(ncomp), FUN = function(z, w){
            z$comp = w
            rownames(z) <- NULL
            if (nrepeat > 2)
            {
                z$corQ1 = z$cor.mean - z$cor.sd
                z$corQ3 = z$cor.mean + z$cor.sd
            }
            z
        }, SIMPLIFY = FALSE)
        cors <- Reduce(rbind, cors)
        cors$comp <- factor(cors$comp)
        
        if (is.null(col))
        {
            col <- color.mixo(seq_len(ncomp))
        }
        names(col) <- seq_len(ncomp)
        p <- ggplot(cors, aes(.data[["keepX"]], .data[["cor.mean"]], col = .data[["comp"]])) +
            theme_minimal() +
            geom_line() +
            geom_point() +
            scale_x_continuous(trans='log10', breaks = cors$keepX) +
            ylim(c(min(cors$corQ1, 0), max(cors$corQ3, 1))) +
            labs(x= 'Number of features selected',
                 y = 'Correlation of components',
                 col = 'Comp') +
            scale_color_manual(values = col)
        
        if (nrepeat > 2) {
            p <- p +  geom_point(data=cors[!is.na(cors$opt.keepX),], 
                                 aes(.data[["keepX"]], .data[["cor.mean"]], col = .data[["comp"]]),
                                 size=6, shape = 18, show.legend = FALSE)
        }
        
        
        ## ----- error bars
        if (isTRUE(sd))
        {
            p <- p + geom_errorbar(aes(ymin = .data[["corQ1"]], ymax = .data[["corQ3"]]),
                                   # position = position_dodge(0.02),
                                   width = 0.04,
                                   ...)
            ## suppress "position_dodge requires non-overlapping x intervals"
            suppressWarnings(print(p))
            return(invisible(p))
        } else
        {
            return(p)
        }
        
    }

## --------------------------- plot.tune.(s)pls --------------------------- ##
#' Plot for model performance
#' 
#' Function to plot performance criteria, such as classification error rate or
#' balanced error rate on a tune.splsda result.
#' 
#' \code{plot.tune.splsda} plots the classification error rate or the balanced
#' error rate from x$error.rate, for each component of the model. A lozenge
#' highlights the optimal number of variables on each component.
#' 
#' \code{plot.tune.block.splsda} plots the classification error rate or the
#' balanced error rate from x$error.rate, for each component of the model. The
#' error rate is ordered by increasing value, the yaxis shows the optimal
#' combination of keepX at the top (e.g. `keepX on block 1'_`keepX on block
#' 2'_`keepX on block 3')
#' 
#' @param x an \code{tune.splsda} object.
#' @param optimal If TRUE, highlights the optimal keepX per component
#' @param sd If 'nrepeat' was used in the call to 'tune.splsda', error bar
#' shows the standard deviation if sd=TRUE
#' @param col character (or symbol) color to be used, possibly vector. One
#' colour per component.
#' @return none
#' @author Kim-Anh Lê Cao, Florian Rohart, Francois Bartolo, AL J Abadi
#' @seealso \code{\link{tune.mint.splsda}}, \code{\link{tune.splsda}}
#' \code{\link{tune.block.splsda}} and http://www.mixOmics.org for more
#' details.
#' @keywords regression multivariate hplot
#' @name plot.tune
#' @method plot tune.spls1
#' @importFrom reshape2 melt
#' @export
plot.tune.spls1 <-
    function(x, optimal = TRUE, sd = NULL, col, ...)
    {
        # TODO add examples
        # to satisfy R CMD check that doesn't recognise x, y and group (in aes)
        y = Comp = lwr = upr = NULL
        
        if (!is.logical(optimal))
            stop("'optimal' must be logical.", call. = FALSE)
        sd = .change_if_null(sd, !is.null(x$error.rate.sd))
        error <- x$error.rate
        error.rate.sd = x$error.rate.sd # for LOGOCV and nrepeat=1, will be NULL
        
        #
        if (is.null(x$error.rate.sd)) {
            message("Note: sd bars cannot be calculated when nrepeat = 1.\n")
            ylim = range(error)
        } else {
            ylim = range(c(error + error.rate.sd), c(error - error.rate.sd))
        }
        
        optimal <- optimal && (any(grepl('mint', class(x))) || x$call$nrepeat > 2)
        select.keepX <- x$choice.keepX[colnames(error)]
        comp.tuned = length(select.keepX)
        
        legend=NULL
        measure = x$measure
        
        if (length(select.keepX) < 10)
        {
            #only 10 colors in color.mixo
            if(missing(col))
                col = color.mixo(seq_len(comp.tuned))
        } else {
            #use color.jet
            if(missing(col))
                col = color.jet(comp.tuned)
        }
        if(length(col) != comp.tuned)
            stop("'col' should be a vector of length ", comp.tuned,".")
        
        if(measure == "overall")
        {
            ylab = "Classification error rate"
        } else if (measure == "BER")
        {
            ylab = "Balanced error rate"
        } else if (measure == "MSE"){
            ylab = "MSE"
        }else if (measure == "MAE"){
            ylab = "MAE"
        }else if (measure == "Bias"){
            ylab = "Bias"
        }else if (measure == "R2"){
            ylab = "R2"
        }else if (measure == "AUC"){
            ylab = "AUC"
        }
        
        #legend
        names.comp = substr(colnames(error),5,10) # remove "comp" from the name
        if(length(x$choice.keepX) == 1){
            #only first comp tuned
            legend = "1"
        } else if(length(x$choice.keepX) == comp.tuned) {
            # all components have been tuned
            legend = c("1", paste("1 to", names.comp[-1]))
        } else {
            #first components were not tuned
            legend = paste("1 to", names.comp)
        }
        
        
        # creating data.frame with all the information
        df = melt(error)
        colnames(df) = c("x","Comp","y")
        df$Comp = factor(df$Comp, labels=legend)
        
        p = ggplot(df, aes(x = x, y = y, color = Comp)) +
            labs(x = "Number of selected features", y = ylab) +
            theme_bw() +
            geom_line()+ geom_point()
        p = p+ scale_x_continuous(trans='log10') +
            scale_color_manual(values = col)
        
        # error bar
        if(!is.null(error.rate.sd) && sd)
        {
            dferror = melt(error.rate.sd)
            df$lwr = df$y - dferror$value
            df$upr = df$y + dferror$value
            
            #adding the error bar to the plot
            p = p + geom_errorbar(data=df,aes(ymin=lwr, ymax=upr), width = 0.04)
        }
        
        if(optimal)
        {
            index = NULL
            for(i in seq_len(comp.tuned))
                index = c(index, which(df$x == select.keepX[i] & df$Comp == levels(df$Comp)[i]))
            
            # adding the choseen keepX to the graph
            p = p + geom_point(data=df[index,],size=7, shape = 18)
            p = p + guides(color = guide_legend(override.aes =
                                                    list(size=0.7,stroke=1)))
        }
        
        p
    }

## -------------------------- plot.tune.splsda -------------------------- ##
#' @rdname plot.tune
#' @method plot tune.splsda
#' @export
plot.tune.splsda <- plot.tune.spls1
# TODO add examples


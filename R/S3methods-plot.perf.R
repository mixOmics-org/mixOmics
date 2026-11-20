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

## --------------------------- plot.perf.(s)plsda --------------------------- ##
#' Plot for model performance for PSLDA analyses
#' 
#' Function to plot classification performance for supervised
#' methods, as a function of the number of components.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017).
#' See ?perf for examples.
#' 
#' @author Ignacio González, Florian Rohart, Francois Bartolo, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{perf}}.
#' @references
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords regression multivariate hplot
#' @name plot.perf
#' @return none
NULL

## -------------------------- plot.perf.(s)plsda -------------------------- ##
#' @param x an \code{perf.plsda} object.
#' @param dist prediction method applied in \code{perf} for \code{plsda} or
#' \code{splsda}. See \code{\link{perf}}.
#' @param measure Two misclassification measure are available: overall
#' misclassification error \code{overall} or the Balanced Error Rate \code{BER}
#' @param col character (or symbol) colour to be used, possibly vector. One
#' color per distance \code{dist}.
#' @param xlab,ylab titles for \eqn{x} and \eqn{y} axes.  Typically character
#' strings, but can be expressions (e.g., \code{expression(R^2)}).
#' @param overlay parameter to overlay graphs; if 'all', only one graph is
#' shown with all outputs; if 'measure', a graph is shown per distance; if
#' 'dist', a graph is shown per measure.
#' @param legend.position position of the legend, one of "vertical" (only one
#' column) or "horizontal" (two columns).
#' @param sd If 'nrepeat' was used in the call to 'perf', error bar shows the
#' standard deviation if sd=TRUE. For mint objects sd is set to FALSE as the
#' number of repeats is 1.
#' @param ... Not used.
#' @method plot perf.plsda.mthd
#' @rdname plot.perf
#' @importFrom methods hasArg
#' @export
plot.perf.plsda.mthd <- 
    function (x,
              dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
              measure = c("all", "overall", "BER"),
              col,
              xlab = NULL,
              ylab = NULL,
              overlay = c("all", "measure", "dist"),
              legend.position = c("vertical", "horizontal"),
              sd = TRUE,
              ...
    )
    {
        # maybe later, so far we set type = "l"
        type = "l"
        
        if (hasArg(pred.method))
            stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
        pred.method = NULL # to pass R CMD check
        
        
        if (any(measure == "all"))
            measure = names(x$error.rate)
        
        if (is.null(measure) || !any(measure %in% names(x$error.rate)))
            stop("'measure' should be among the ones used in your call to 'perf': ", paste(names(x$error.rate),collapse = ", "),".")
        
        if (any(dist == "all"))
            dist = colnames(x$error.rate[[1]])
        
        
        if (is.null(dist) || !any(dist %in% colnames(x$error.rate[[1]])))
            stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$error.rate[[1]]),collapse = ", "),".")
        
        if(missing(col)) #one col per distance
        {
            col = color.mixo(1:length(dist))
        } else {
            if(length(col) != length(dist))
                stop("'col' should be a vector of length ", length(dist),".")
        }
        
        if (is.null(ylab))
            ylab = 'Classification error rate'
        
        if (is.null(xlab))
            xlab = 'Component'
        
        if(length(overlay) >1 )
            overlay = overlay[1]
        
        if(length(legend.position) >1 )
            legend.position = legend.position[1]
        
        # error.rate is a list [[measure]]
        # error.rate[[measure]] is a matrix of dist columns and ncomp rows
        # same for error.rate.sd, if any
        error.rate = x$error.rate
        if(sd)
        {
            error.rate.sd = x$error.rate.sd
        } else {
            error.rate.sd = NULL
        }
        def.par = par(no.readonly = TRUE)
        
        internal_graphic.perf(error.rate = error.rate, error.rate.sd = error.rate.sd,
                              overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                              xlab = xlab, ylab = ylab, sd = sd, color = col, ...)
        
        par(def.par)
        # error.bar(out,as.vector(mat.error.plsda),as.vector(cbind(x$error.rate.sd$overall,x$error.rate.sd$BER)))
        
        return(invisible())
        
    }


#' @method plot perf.splsda.mthd
#' @rdname plot.perf
#' @export
plot.perf.splsda.mthd <- plot.perf.plsda.mthd

## ----------------------- plot.perf.mint.(s)plsda ------------------------ ##
#' @param study Indicates which study-specific outputs to plot. A character
#' vector containing some levels of \code{object$study}, "all.partial" to plot
#' all studies or "global" is expected. Default to "global".
#' @rdname plot.perf
#' @method plot perf.mint.plsda.mthd
#' @importFrom methods hasArg
#' @export
plot.perf.mint.plsda.mthd <-
    function (x,
              dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
              measure = c("all", "overall", "BER"),
              col,
              xlab = NULL,
              ylab = NULL,
              study = "global",
              overlay = c("all", "measure", "dist"),
              legend.position = c("vertical", "horizontal"),
              ...
    )
    {
        if (isTRUE(list(...)$sd))
            message("'sd' not applicable to perf.mint.plsda objects. See ?plot.perf.")
        # maybe later, so far we set type = "l"
        type = "l"
        
        if (hasArg(pred.method))
            stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
        pred.method = NULL # to pass R CMD check
        
        
        if (any(measure == "all"))
            measure = c("BER","overall")
        
        if (is.null(measure) || !any(measure %in% c("BER","overall")))
            stop("'measure' should be among the ones used in your call to 'perf': ", paste(c("BER","overall"),collapse = ", "),".")
        
        
        if (any(dist == "all"))
            dist = colnames(x$global.error[[1]])
        
        if(length(overlay) >1 )
            overlay = overlay[1]
        
        if(length(legend.position) >1 )
            legend.position = legend.position[1]
        
        if(missing(col)) #one col per distance
        {
            col = color.mixo(1:length(dist))
        } else {
            if(length(col) != length(dist))
                stop("'col' should be a vector of length ", length(dist),".")
        }
        
        
        if(any(study == "global"))
        {
            if (is.null(dist) || !any(dist %in% colnames(x$global.error[[1]])))
                stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$global.error[[1]]),collapse = ", "),".")
            
            
            if (is.null(ylab))
                ylab = 'Classification error rate'
            
            if (is.null(xlab))
                xlab = 'Component'
            
            # error.rate is a list [[measure]]
            # error.rate[[measure]] is a matrix of dist columns and ncomp rows
            # same for error.rate.sd, if any
            
            error.rate = x$global.error
            
            def.par = par(no.readonly = TRUE)
            
            internal_graphic.perf(error.rate = error.rate, error.rate.sd = NULL,
                                  overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                                  xlab = xlab, ylab = ylab, color = col, ...)
            
            par(def.par)
            
        } else {
            
            def.par = par(no.readonly = TRUE)
            
            
            if (any(study == "all.partial"))
                study = 1:length(x$study.specific.error)
            
            
            if (any(dist == "all"))
                dist = colnames(x$study.specific.error[[1]][[1]])
            
            if((length(study) >1) & (overlay != "all"))
                stop("When more than one study is plotted, overlay must be 'all'")
            
            
            if (is.null(dist) || !any(dist %in% colnames(x$study.specific.error[[1]][[1]])))
                stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$study.specific.error[[1]][[1]]),collapse = ", "),".")
            
            if (is.null(ylab))
                ylab = 'Classification error rate'
            
            if (is.null(xlab))
                xlab = 'Component'
            
            
            if(overlay=="all")
            {
                par(mfrow=c(1,length(study)))
                
            } else if(overlay=="measure") {
                par(mfrow=c(length(study),length(dist)))
            } else if(overlay=="dist") {
                par(mfrow=c(length(study),length(measure)))
            }
            
            
            for(stu in study)
            {
                error.rate = x$study.specific.error[[stu]]
                
                internal_graphic.perf(error.rate = error.rate, error.rate.sd = NULL,
                                      overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                                      xlab = xlab, ylab = ylab, color = col, ...)
                
                if (overlay == "all")
                    title(stu, line = 1)
            }
            
            if((length(study)==1) & (length(measure) > 1) & overlay != "all")
                title(stu, outer=TRUE, line = -1)#,...)
            
            
            par(def.par)
            
        }
        return(invisible())
        
    }

#' @method plot perf.mint.splsda.mthd
#' @rdname plot.perf
#' @export
plot.perf.mint.splsda.mthd <- plot.perf.mint.plsda.mthd

## --------------------------- plot.perf.sgccda --------------------------- ##
#' @method plot perf.sgccda.mthd 
#' @importFrom methods hasArg
#' @param weighted plot either the performance of the Majority vote or the
#' Weighted vote.
#' @rdname plot.perf
#' @export
plot.perf.sgccda.mthd <-
    function (x,
              dist = c("all","max.dist","centroids.dist","mahalanobis.dist"),
              measure = c("all","overall","BER"),
              col,
              weighted = TRUE,
              xlab = NULL,
              ylab = NULL,
              overlay= c("all", "measure", "dist"),
              legend.position=c("vertical","horizontal"),
              sd = TRUE,
              ...)
    {
        # maybe later, so far we set type = "l"
        type = "l"
        
        if (hasArg(pred.method))
            stop("'pred.method' argument has been replaced by 'dist' to match the 'tune' and 'perf' functions")
        pred.method = NULL # to pass R CMD check
        
        
        measure.input = measure
        measure = NULL
        if(any(measure.input == "all"))
            measure.input = c("BER", "overall")
        
        if(any(measure.input == "BER"))
            measure = c(measure, "Overall.BER")
        
        if (any(measure.input == "overall"))
            measure = c(measure, "Overall.ER")
        
        if(!all(measure.input %in% c("all", "overall", "BER")))
            stop("'measure' must be 'all', 'overall' or 'BER'")
        
        if (any(dist == "all"))
            dist = colnames(x$error.rate[[1]])
        
        if(length(overlay) >1 )
            overlay = overlay[1]
        
        if(length(legend.position) >1 )
            legend.position = legend.position[1]
        
        
        if (is.null(dist) || !any(dist %in% colnames(x$error.rate[[1]])))
            stop("'dist' should be among the ones used in your call to 'perf': ", paste(colnames(x$error.rate[[1]]),collapse = ", "),".")
        
        if(missing(col)) #one col per distance
        {
            col = color.mixo(1:length(dist))
        } else {
            if(length(col) != length(dist))
                stop("'col' should be a vector of length ", length(dist),".")
        }
        
        if (is.null(ylab))
            ylab = 'Classification error rate'
        
        if (is.null(xlab))
            xlab = 'Component'
        
        if(weighted == TRUE)
        {
            perfo = "WeightedVote.error.rate"
            perfo.sd = "WeightedVote.error.rate.sd"
        } else {
            perfo = "MajorityVote.error.rate"
            perfo.sd = "MajorityVote.error.rate.sd"
        }
        
        if(sd == TRUE)
        {
            if(is.null(x[[perfo.sd]]))
                sd = FALSE
        }
        
        # error.rate is a list [[measure]]
        # error.rate[[measure]] is a matrix of dist columns and ncomp rows
        # same for error.rate.sd, if any
        error.rate = error.rate.sd = list()
        for(mea in measure)
        {
            error.temp = error.temp.sd = NULL
            for(di in dist)
            {
                temp = t(x[[perfo]][[di]][mea, , drop=FALSE])
                colnames(temp) = di
                error.temp = cbind(error.temp, temp)
                if(sd)
                {
                    temp.sd = t(x[[perfo.sd]][[di]][mea, , drop=FALSE])
                    colnames(temp.sd) = di
                    error.temp.sd = cbind(error.temp.sd, temp.sd)
                }
                
            }
            error.rate[[mea]] = error.temp
            if(sd)
            {
                error.rate.sd[[mea]] = error.temp.sd
            } else {
                error.rate.sd = NULL
            }
        }
        
        
        
        def.par = par(no.readonly = TRUE)
        
        internal_graphic.perf(error.rate = error.rate, error.rate.sd = error.rate.sd,
                              overlay = overlay, type = type, measure = measure, dist = dist, legend.position = legend.position,
                              xlab = xlab, ylab = ylab, color = col, ...)
        
        par(def.par)
        return(invisible())
        
    }



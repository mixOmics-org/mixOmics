################################################################################
# Author:
#   Francois Bartolo,
#   Florian Rohart,
#
# created: 29-08-2016
# last modified: 30-08-2016
#
# Copyright (C) 2016
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################

# -----------------------------------------------------------------------------
# Internal helpers functions to run "plot.perf" functions
# -----------------------------------------------------------------------------


internal_graphic.perf <-
    function (error.rate,
              error.rate.sd,
              overlay,
              type,
              measure,
              dist,
              legend.position,
              xlab,
              ylab,
              color,
              ylim = NULL,
              ...
    )
{
    # error.rate is a list [[measure]]
    # error.rate[[measure]] is a matrix of dist columns and ncomp rows
    # color is a vector of length dist
    
    # concatenation of the error in a single matrix
    error.rate.concat = matrix(nrow = nrow(error.rate[[1]]), ncol = 0)
    for(mea in measure)
    {
        temp = error.rate[[mea]][, dist, drop = FALSE]
        colnames(temp) = paste(mea, colnames(temp),sep="_")
        error.rate.concat = cbind(error.rate.concat, temp)
    }
    
    if(!is.null(error.rate.sd))
    {
        error.rate.sd.concat = matrix(nrow = nrow(error.rate.sd[[1]]), ncol = 0)
        for(mea in measure)
        {
            temp = error.rate.sd[[mea]][, dist, drop = FALSE]
            colnames(temp) = paste(mea, colnames(temp),sep="_")
            error.rate.sd.concat = cbind(error.rate.sd.concat, temp)
        }
        if(is.null(ylim))
            ylim = range(c(error.rate.concat + error.rate.sd.concat),
                         c(error.rate.concat - error.rate.sd.concat))
        
    } else {
        error.rate.sd.concat = NULL
        if(is.null(ylim))
            ylim = range(error.rate.concat)
        
    }
    
    ## extract component numbers from rownames(error.rate.concat)
    
    rownames(error.rate.concat) = gsub("comp","",rownames(error.rate.concat))
    
    component = as.numeric(rownames(error.rate.concat))
    
    if(overlay == "all")
    {
        out<-matplot(component, error.rate.concat, type = type, lty =
                         rep(c(1:length(measure)), each = length(dist)),
                     col = rep(color, length(measure)),
                     lwd = 2, xlab = xlab, ylab = ylab, xaxt="n", ylim = ylim)
        
        axis(1, rownames(error.rate.concat))#, rownames(error.rate.concat))
        
        if(legend.position == "vertical")
        {
            legend('topright', legend = c(measure, dist),
                   lty = c(1:length(measure), rep(NA, length(dist))),
                   pch = c(rep(NA, length(measure)), rep(16, length(dist))),
                   col = c(rep('black', length(measure)), color), ncol = 1, lwd = 2)
        } else if(legend.position == "horizontal") {
            legend('topright', legend = c(measure, "" , dist),
                   lty = c(1:length(measure), rep(NA, (length(dist)+1))),
                   pch = c(rep(NA, (length(measure)+1)), rep(16, length(dist))),
                   col = c(rep('black', length(measure)), NA, color), ncol = 2,
                   lwd = 2)
            
        }
        if(!is.null(error.rate.sd.concat))
        {
            col = matrix(rep(color, length(measure)),
                         nrow = nrow(error.rate.concat), ncol = length(measure)*length(dist),
                         byrow=TRUE)
            
            for(j in 1:ncol(error.rate.concat))
                plot_error_bar(x = component, y = error.rate.concat[, j],
                               uiw=error.rate.sd.concat[, j], add=TRUE, col = col[, j])
        }
    } else if(overlay == "measure") {
        
        # overlaying all measure on one graph, a graph per distance
        par(mfrow=c(1, length(dist)))
        for(di in dist)
        {
            new_mat.error =
                error.rate.concat[, grep(di, colnames(error.rate.concat)),
                                  drop = FALSE]
            
            out<-matplot(new_mat.error, type = type, lty = c(1:length(measure)),
                         col = rep(color[which(di == dist)]),
                         lwd = 2, xlab = xlab, ylab = ylab, xaxt="n", ylim = ylim)
            
            axis(1, rownames(error.rate.concat))#, rownames(error.rate.concat))
            #axis(2)
            
            if(legend.position == "vertical")
            {
                legend('topright', legend = measure, lty = 1:length(measure),
                       col = rep(color[which(di == dist)]), ncol = 1, lwd = 2)
            } else if(legend.position == "horizontal") {
                legend('topright', legend = c(measure), lty = 1:length(measure),
                       col = rep(color[which(di == dist)]), ncol = 2, lwd = 2)
            }
            if(!is.null(error.rate.sd.concat))
            {
                new_sd.error = error.rate.sd.concat[, grep(di,
                                                           colnames(error.rate.sd.concat)), drop = FALSE]
                for(j in 1:ncol(new_mat.error))
                    plot_error_bar(x = component, y = new_mat.error[, j],
                                   uiw=new_sd.error[, j], add=TRUE, col = color[which(di == dist)])
            }
            #add title for each subgraph
            title(di, line = 1)
        }
        
    } else if(overlay == "dist") {
        
        # overlaying all distance on one graph, a graph per measure
        
        par(mfrow=c(1, length(measure)))
        for(mea in measure)
        {
            new_mat.error=error.rate.concat[, grep(mea,
                                                   colnames(error.rate.concat)), drop = FALSE]
            
            out<-matplot(new_mat.error, type = type, lty = c(1:length(dist)),
                         col = color,
                         lwd = 2, xlab = xlab, ylab = ylab, xaxt="n", ylim = ylim)
            
            axis(1, rownames(error.rate.concat))#, rownames(error.rate.concat))
            #axis(2)
            
            if(legend.position == "vertical")
            {
                legend('topright', legend = dist, lty = 1:length(dist),
                       col = color, ncol = 1, lwd = 2)
            } else if(legend.position == "horizontal") {
                legend('topright', dist, lty = 1:length(dist),
                       col = color, ncol = 2, lwd = 2)
            }
            if(!is.null(error.rate.sd.concat))
            {
                col = matrix(color, nrow = nrow(error.rate.concat),
                             ncol = length(dist), byrow=TRUE)
                
                new_sd.error=error.rate.sd.concat[,
                                                  grep(mea, colnames(error.rate.sd.concat)), drop = FALSE]
                for(j in 1:ncol(new_mat.error))
                    plot_error_bar(x = component, y = new_mat.error[, j],
                                   uiw=new_sd.error[, j], add=TRUE, col = col[, j])
            }
            title(mea, line = 1)
        }
        
    }
    
}


plot_error_bar <- function (x, y = as.numeric(x), uiw, add=FALSE, col = "black")
{
    # x are the xaxis values
    # y are the y axis values
    arglist <- list()#...)
    liw = uiw
    #y <- as.numeric(x)
    #x <- seq(along = x)
    z <- y
    ui <- z + uiw
    li <- z - liw
    arglist$xlab <- deparse(substitute(x))
    arglist$ylab <- deparse(substitute(y))
    arglist$ylim <- range(c(y, ui, li), na.rm = TRUE)
    arglist$xlim <- range(c(x, ui, li), na.rm = TRUE)
    plotpoints <- TRUE
    ul <- c(li, ui)
    pin <- par("pin")
    usr <- par("usr")
    x.to.in <- pin[1]/diff(usr[1:2])
    y.to.in <- pin[2]/diff(usr[3:4])
    gap <- rep(0, length(x)) * diff(par("usr")[3:4])
    smidge <- par("fin")[1] * 0.01
    nz <- abs(li - pmax(y - gap, li)) * y.to.in > 0.001
    scols <- col
    
    arrow.args <- c(list(lty = par("lty"), angle = 90, length = smidge,
                         code = 1, col = scols))
    
    do.call("arrows", c(list(x[nz], li[nz], x[nz], pmax(y -
                                                            gap, li)[nz]), arrow.args))
    
    nz <- abs(ui - pmin(y + gap, ui)) * y.to.in > 0.001
    #scols <- rep(scol, length.out = length(x))[nz]
    arrow.args$col <- scols
    
    do.call("arrows", c(list(x[nz], ui[nz], x[nz], pmin(y +
                                                            gap, ui)[nz]), arrow.args))
    
    do.call("points", c(list(x, y, bg = par("bg"), col =col)))
    
    invisible(list(x = x, y = y))
}

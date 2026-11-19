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
#-- Includes plotIndiv for PCA, sPCA, IPCA, sIPCA --#
#----------------------------------------------------------------------------------------------------------#
#' @rdname plotIndiv
#' @method plotIndiv pca
#' @export
plotIndiv.pca <- 
    function(object,
             comp = NULL,
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
             legend = FALSE,
             legend.title.pch = "Legend",
             legend.title = "Legend",
             legend.position = "right",
             X.label = NULL,
             Y.label = NULL,
             Z.label = NULL,
             xlim = NULL,
             ylim = NULL,
             axes.box = "box",
             abline = FALSE,
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
            warning("'rep.space' is not used for PCA, sPCA, IPCA, sIPCA")
        }
        if ("study" %in% names(extra_args) || "layout" %in% names(extra_args)) {
            warning("'study' and 'layout' arguments are only used for MINT models")
        }
        if ("blocks" %in% names(extra_args)) {
            warning("'blocks' argument is only used for multiblock models")
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
            point.lwd = point.lwd
        )
        
        blocks = "X"
        rep.space = "X-variate"
        layout = NULL
        
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
        
        # retrieve outputs from the checks
        axes.box = check$axes.box
        comp = check$comp
        xlim = check$xlim
        ylim = check$ylim
        ind.names = check$ind.names
        display.names = check$display.names
        
        
        #-- Get variates
        x = y = z = list()
        x[[1]] = object$variates$X[, comp[1]]
        y[[1]] = object$variates$X[, comp[2]]
        if(style == "3d") z[[1]] = object$variates$X[, comp[3]]
        
        
        #-- Variance explained on X, Y and Z labels
        
        if (style ==  "3d")
        {
            inf = object$prop_expl_var$X[c(comp[1], comp[2], comp[3])]
            inf = round(inf, 2)
        } else {
            inf = object$prop_expl_var$X[c(comp[1], comp[2])]
            inf = round(inf, 2)}
        
        
        if (is.null(X.label))
        {
            X.label = paste("PC", comp[1], sep = '')
            percentage = paste0(inf[1]*100, "% expl. var")
            X.label = paste(X.label, percentage, sep = ": ")
        }
        if (is.null(Y.label))
        {
            Y.label = paste("PC", comp[2], sep = '')
            percentage = paste0(inf[2]*100, "% expl. var")
            Y.label = paste(Y.label, percentage, sep = ": ")
        }
        if (is.null(Z.label)&&style == "3d")
        {
            Z.label = paste("PC", comp[3], sep = '')
            percentage = paste0(inf[3]*100, "% expl. var")
            Z.label = paste(Z.label, percentage, sep = ": ")
        }
        
        
        n = nrow(object$X)
        
        # create data frame df that contains (almost) all the ploting information
        out = shape.input.plotIndiv(
            object = object,
            n = n,
            blocks = blocks,
            x = x,
            y = y,
            z = z,
            ind.names = ind.names,
            group = group,
            col.per.group = col,
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


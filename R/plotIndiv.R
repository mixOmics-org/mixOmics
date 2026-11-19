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

#' Plot of Individuals (Experimental Units)
#' 
#' This function provides scatter plots for individuals (observations)
#' representation in (s)(I)PCA, (r)CCA, (s)PLS(DA), multiblock (s)PLS(DA) and
#' MINT models. 
#' 
#' \code{plotIndiv} method makes scatter plot for individuals representation
#' depending on the subspace of projection. Each point corresponds to an
#' individual.
#' 
#' If \code{ind.names=TRUE} and row names is \code{NULL}, then
#' \code{ind.names=1:n}, where \code{n} is the number of individuals. Also, if
#' \code{pch} is an input, then \code{ind.names} is set to FALSE as we do not
#' show both names and shapes.
#' 
#' The outputs of \code{"ellipse", "ellipse.level", "centroid", "star"} are:
#' 1) based on the sample group membership (specified with the \code{group}) for unsupervised
#' models such as \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA,
#' sGCCA} or 2) based on the outcome \code{Y} for supervised models such as  \code{PLS-DA, SPLS-DA,sGCCDA}
#' 
#' #' In the specific case of a single `omics supervised model
#' (\code{\link{plsda}}, \code{\link{splsda}}), users can overlay prediction
#' results to sample plots in order to visualise the prediction areas of each
#' class, via the \code{background} input parameter. Note that this
#' functionality is only available for models with less than 2 components as
#' the surfaces obtained for higher order components cannot be projected onto a
#' 2D representation in a meaningful way. For more details, see
#' \code{\link{background.predict}}
#' 
#' The \code{group} parameter is used to colour the samples by group membership, and by 
#' default the shape of samples will also correspond to these groups. Sample shapes can
#' be further customised using the \code{pch} parameter. This can be done in 3 ways:
#' 1) If \code{pch} is one numeric value, shape of all points will be this. 
#' 2) If \code{pch} is a vector of numeric values length equal to the number of groups 
#' set in \code{group}, then each group will have the given different shape. 
#' 3) If \code{pch} is a character vector of length equal to the number of samples, then the
#' shape of each sample will correspond to this grouping, i.e. this allows users to visualise
#' second group membership. In case 3), the elements of the vector should be the names of the 
#' second group membership, as these will appear in the legend. The shapes of the second group
#' membership will be 1, 2, 3, etc by default, control of the shape is only possible by
#' reordering the levels of the second group membership.
#' See examples.
#' 
#' For customized plots (i.e. adding points, text), use the style = 'graphics'
#' (default is ggplot2).
#' 
#' Note: the ellipse options were borrowed from the \pkg{ellipse}.
#' 
#' @param object object of class inherited from any \pkg{mixOmics}: \code{PLS,
#' sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA, MINT}
#' @param comp integer vector of length two (or three to 3d). The components
#' that will be used on the horizontal and the vertical axis respectively to
#' project the individuals.
#' @param rep.space - only for objects of class \code{"pca", "plsda", "plsda", "rcc"}. Can be 
#' \code{"X-variate", "Y-variate", "XY-variate", "multi"}.
#' @param blocks - only for multiblock objects. Integer value or name(s) of block(s) to be plotted OR
#' "average" (averages the components from all blocks to produce a consensus plot) OR
#' "weighted.average" (weighted average of the components according to their correlation with the
#' outcome Y). See examples.
#' @param study - only for MINT models. Indicates whether to 1) plot all studies together
#' \code{"global"} (default), 2) facet plot by study to show all at once \code{"all.partial"}
#' or 3) plot specififice studies (see examples). 
#' @param layout layout parameter passed to mfrow. Only for MINT model and only 
#' used when \code{study} is "all.partial"
#' @param style argument to be set to either \code{'graphics'},
#' \code{'lattice'}, \code{'ggplot2'} or \code{'3d'} for a style of plotting.
#' Default set to 'ggplot2'. See details. \code{3d} is not available for MINT objects.
#' @param ind.names either a character vector of names for the individuals to
#' be plotted, or \code{FALSE} for no names. If \code{TRUE}, the row names of
#' the first (or second) data matrix is used as names (see Details). If `pch` is set this
#' will overwrite the names as shapes. Default is \code{TRUE}. 
#' Not avaliable for MINT objects. 
#' @param group factor indicating the group membership for each sample, useful
#' for colouring samples by groups, adding ellipses, centroids and stars.
#' Coded as default for the supervised methods \code{PLS-DA,
#' SPLS-DA,sGCCDA}, but needs to be input for the unsupervised methods
#' \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA}. 
#' Order of levels will reflect the order the groups appear in legends and 
#' correspond to the order of set colours in \code{col}. 
#' @param col character (or symbol) color to be used. If \code{group} provided, should be a 
#' vector of the same length as groups, order of colours will be respected to correspond to 
#' order of \code{group} levels. 
#' @param ellipse Logical indicating if ellipse plots should be plotted. See details. 
#' @param ellipse.level Numerical value indicating the confidence level of
#' ellipse being plotted when \code{ellipse =TRUE} (i.e. the size of the
#' ellipse). The default is set to 0.95, for a 95\% region. See details. 
#' @param centroid Logical indicating whether centroid points should be
#' plotted. See details. 
#' @param star Logical indicating whether a star plot should be plotted, with
#' arrows starting from the centroid (see argument \code{centroid}, and ending
#' for each sample belonging to each group or outcome. See details. 
#' @param background color the background by the predicted class, see
#' \code{\link{background.predict}}
#' @param pch value or vector, controls shape of points to be plotted. Can be used 
#' to control the shape of the points directly or to create a second group membership
#' and legend. See details.
#' @param title set of characters indicating the title plot.
#' @param subtitle subtitle for each plot, only used when several \code{block}
#' or \code{study} are plotted.
#' @param legend Logical. Whether the legend should be added. Default is FALSE.
#' @param legend.title title of the legend
#' @param legend.title.pch title of the second legend created by \code{pch}, if
#' any.
#' @param legend.position position of the legend, one of "bottom", "left",
#' "top" and "right".
#' @param X.label x axis titles.
#' @param Y.label y axis titles.
#' @param Z.label z axis titles (when style = '3d').
#' @param axes.box for style '3d', argument to be set to either \code{'axes'},
#' \code{'box'}, \code{'bbox'} or \code{'all'}, defining the shape of the box.
#' @param abline should the vertical and horizontal line through the center be
#' plotted? Default set to \code{FALSE}
#' @param xlim,ylim numeric list of vectors of length 2 and length
#' =length(blocks), giving the x and y coordinates ranges.
#' @param cex numeric character (or symbol) expansion, possibly vector.
#' @param alpha Semi-transparent colors (0 < \code{'alpha'} < 1)
#' @param point.lwd \code{lwd} of the points, used when \code{ind.names =
#' FALSE}
#' @param size.title size of the title
#' @param size.subtitle size of the subtitle
#' @param size.xlabel size of xlabel
#' @param size.ylabel size of ylabel
#' @param size.axis size of the axis
#' @param size.legend size of the legend
#' @param size.legend.title size of the legend title
#' @param ... Optional arguments or type par can be added with \code{style =
#' 'graphics'}
#' @return none
#' @author Ignacio González, Benoit Gautier, Francois Bartolo, Florian Rohart,
#' Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{text}}, \code{\link{background.predict}},
#' \code{\link{points}} and http://mixOmics.org/graphics for more details.
#' @keywords multivariate hplot dplot
#' @export
#' @example ./examples/plotIndiv-examples.R 
plotIndiv  <- function(object, ...) {
  UseMethod("plotIndiv")
}

# --------------------------------------------------------------------------------------
# Internal helpers functions to run "plotIndiv" functions
# --------------------------------------------------------------------------------------

check.input.plotIndiv <-
  function(object,
           comp = NULL,
           blocks = NULL, # to choose which block data to plot, when using GCCA module
           ind.names = TRUE,
           style = "ggplot2", # can choose between graphics, 3d, lattice or ggplot2
           #study = "global",
           ellipse = FALSE,
           ellipse.level = 0.95,
           centroid = FALSE,
           star = FALSE,
           legend = FALSE,
           X.label = NULL,
           Y.label = NULL,
           Z.label = NULL,
           abline = FALSE,
           xlim = NULL,
           ylim = NULL,
           alpha = 0.2,
           axes.box = "box",
           plot_parameters
  )
  {
    
    
    # --------------------------------------------------------------------------------------
    #           independent from class.object
    # --------------------------------------------------------------------------------------
    
    ### Start: Validation of arguments
    ncomp = object$ncomp
    
    #-- style
    if (!style %in% c("ggplot2", "lattice", "graphics", "3d"))
      stop("'style' must be one of 'ggplot2', 'lattice', 'graphics' or '3d' .", call. = FALSE)
    
    #-- axes.box
    choices = c("box", "bbox", "both")
    axes.box = choices[pmatch(axes.box, choices)]
    
    if (is.na(axes.box))
      stop("'axes.box' should be one of 'box', 'bbox' or 'both'.", call. = FALSE)
    
    #-- ellipse.level
    if ((ellipse.level > 1) | (ellipse.level < 0))
      stop("The value taken by 'ellipse.level' must be between 0 and 1")
    
    #-- legend
    if (length(legend) !=  1 || !is.logical(legend))
      stop("'legend' must be a logical value.", call. = FALSE)
    
    #-- alpha correlation
    if (!is.numeric(alpha) | (alpha > 1) | (alpha < 0))
      stop("The value taken by 'alpha' must be between 0 and 1", call. = FALSE)
    
    
    #-- comp
    if (is.null(comp))
    {
      if (style == "3d")
      {
        comp = c(1:3)
      } else {
        comp = c(1:2)
      }
    }
    if (length(comp) !=  2 && !(style == "3d"))
    {
      stop("'comp' must be a numeric vector of length 2.", call. = FALSE)
    } else if (length(comp) !=  3 && (style == "3d")) {
      stop("'comp' must be a numeric vector of length 3.", call. = FALSE)
    }
    
    if (!is.numeric(comp))
      stop("Invalid vector for 'comp'.")
    
    if (any(ncomp < max(comp)))
      stop(paste0("The number of components of the object to be plotted (ncomp = ", max(object$ncomp), ") is smaller than 'comp' (", paste(comp, collapse = ", "), "). Please increase ncomp  or decrease 'comp'"), call. = FALSE)
    
    comp1 = round(comp[1])
    comp2 = round(comp[2])
    if (style == "3d")
    {
      comp3 = round(comp[3])
    } else {
      comp3 = NULL
    }
    
    
    #ellipse
    if (!is.logical(ellipse))
      stop("'ellipse' must be either TRUE or FALSE", call. = FALSE)
    
    #centroid
    if (!is.logical(centroid))
      stop("'centroid' must be either TRUE or FALSE", call. = FALSE)
    
    #star
    if (!is.logical(star))
      stop("'star' must be either TRUE or FALSE", call. = FALSE)
    
    #legend
    if (!is.logical(legend))
      stop("'legend' must be either TRUE or FALSE", call. = FALSE)
    
    # abline
    if (!is.logical(abline))
      stop("'abline' must be either TRUE or FALSE", call. = FALSE)
    
    #X.label, Y.label, Z.label
    if (!is.null(X.label))
    {
      if (length(X.label)!= 1 | !is.vector(X.label))
        stop("'X.label' must be a vector of length 1", call. = FALSE)
    }
    
    if (!is.null(Y.label))
    {
      if (length(Y.label)!= 1 | !is.vector(Y.label))
        stop("'Y.label' must be a vector of length 1", call. = FALSE)
    }
    
    if (!is.null(Z.label))
    {
      if (style!= "3d")
        warning("'Z.label' is not used as style!= '3d'")
      if (length(Z.label)!= 1 | !is.vector(Z.label))
        stop("'Z.label' must be a vector of length 1", call. = FALSE)
    }
    
    # plot_parameters
    #plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel, size.axis = size.axis,
    #size.legend = size.legend, size.legend.title = size.legend.title, legend.position = legend.position)
    size.title = plot_parameters$size.title
    size.subtitle = plot_parameters$size.subtitle
    size.xlabel = plot_parameters$size.xlabel
    size.ylabel = plot_parameters$size.ylabel
    size.axis = plot_parameters$size.axis
    size.legend = plot_parameters$size.legend
    size.legend.title = plot_parameters$size.legend.title
    legend.title = plot_parameters$legend.title
    legend.title.pch = plot_parameters$legend.title.pch
    legend.position = plot_parameters$legend.position
    point.lwd = plot_parameters$point.lwd
    
    if (!is.numeric(size.title) || length(size.title)>1 || size.title<0)
      stop("'size.title' needs to be a non negative number")
    
    if (!is.numeric(size.subtitle) || length(size.subtitle)>1 || size.subtitle<0)
      stop("'size.subtitle' needs to be a non negative number")
    
    if (!is.numeric(size.xlabel) || length(size.xlabel)>1 || size.xlabel<0)
      stop("'size.xlabel' needs to be a non negative number")
    
    if (!is.numeric(size.ylabel) || length(size.ylabel)>1 || size.ylabel<0)
      stop("'size.ylabel' needs to be a non negative number")
    
    if (!is.numeric(size.axis) || length(size.axis)>1 || size.axis<0)
      stop("'size.axis' needs to be a non negative number")
    
    if (!is.numeric(size.legend) || length(size.legend)>1 || size.legend<0)
      stop("'size.legend' needs to be a non negative number")
    
    if (!is.numeric(size.legend.title) || length(size.legend.title)>1 || size.legend.title<0)
      stop("'size.legend.title' needs to be a non negative number")
    
    if (length(legend.position)>1 || !legend.position%in%c("bottom", "left", "right", "top"))
      stop('"legend.position" needs to be one of "bottom", "left", "right" or "top"')
    
    if (length(legend.title)>1)
      stop("'legend.title' needs to be a single value (length 1)")
    
    if (length(legend.title.pch)>1)
      stop("'legend.title.pch' needs to be a single value (length 1)")
    
    if (!is.numeric(point.lwd) || length(point.lwd)>1 || point.lwd<0)
      stop("'point.lwd' needs to be a non negative number")
    
    if (is.logical(ind.names) & isTRUE(ind.names))
      ind.names = object$names$sample
    if (length(ind.names) > 1)
    {
      if (length(ind.names) !=  length(object$names$sample))
        stop("'ind.names' must be a character vector of length ", length(object$names$sample), " or a Logical atomic vector.")
    }
    
    
    display.names = FALSE
    if (length(ind.names) == length(object$names$sample))
      display.names = TRUE
    
    # --------------------------------------------------------------------------------------
    #           need blocks
    # --------------------------------------------------------------------------------------
    #-- xlim, ylim
    if(style %in% c("lattice", "graphics"))
    {
      if (!is.null(xlim))
      {
        if (length(blocks) == 1) # a single graph is plotted, xlim needs to be a vector of length 2
        {
          if (!is.numeric(xlim) || length(xlim) !=  2)
            stop("'xlim' must be a vector of length 2.", call. = FALSE)
          
          xlim = list(xlim)
          
        } else { # multiple graphs are plotted, xlim needs to be a list of vectors of length 2
          
          if (!is.list(xlim) || length(xlim) !=  length(blocks) || length(unlist(xlim)) !=  2 * length(blocks))
            stop("'xlim' must be a list of ", length(blocks), " vectors of length 2.", call. = FALSE)
        }
      }
      
      if (!is.null(ylim))
      {
        if (length(blocks) == 1) # a single graph is plotted, ylim needs to be a vector of length 2
        {
          if (!is.numeric(ylim) || length(ylim) !=  2)
            stop("'ylim' must be a vector of length 2.", call. = FALSE)
          
          ylim = list(ylim)
          
        } else { # multiple graphs are plotted, ylim needs to be a list of vectors of length 2
          
          if (!is.list(ylim) || length(ylim) !=  length(blocks) || length(unlist(ylim)) !=  2 * length(blocks))
            stop("'ylim' must be a list of ", length(blocks), " vectors of length 2.", call. = FALSE)
        }
      }
    }else if(style =="ggplot2") { #xlim/ylim needs to be a vector: same limits for all graphs
      if (!is.null(xlim))
      {
        if (!is.numeric(xlim) || length(xlim) !=  2)
          stop("'xlim' must be a vector of length 2.", call. = FALSE)
      }
      if (!is.null(ylim))
      {
        if (!is.numeric(ylim) || length(ylim) !=  2)
          stop("'ylim' must be a vector of length 2.", call. = FALSE)
      }
      
    }# for style = 3d, no xlim, ylim used
    
    out = list(axes.box = axes.box, comp = c(comp1, comp2, comp3), xlim = xlim, ylim = ylim, ind.names = ind.names, display.names = display.names)
    
  }


internal_getVariatesAndLabels <- 
  function(object,
           comp,
           blocks.init,
           blocks,
           rep.space,
           style,
           X.label,
           Y.label,
           Z.label
  )
  {
    
    class.object = class(object)
    object.pls = c("mixo_pls", "mixo_spls", "mixo_mlspls", "rcc")
    object.blocks = c("sgcca", "rgcca")
    object.mint = c("mint.pls", "mint.spls", "mint.plsda", "mint.splsda")
    
    #-- Start: Retrieve variates from object
    x = y = z = list()
    if (any(class.object %in%  c(object.pls, object.blocks, object.mint)))
    {
      
      
      x = lapply(object$variates, function(x){x[, comp[1], drop = FALSE]})
      y = lapply(object$variates, function(x){x[, comp[2], drop = FALSE]})
      if (style == "3d")
      {
        z = lapply(object$variates, function(x){x[, comp[3], drop = FALSE]})
      } else {
        z = NULL
      }
      
      
      #-- Variance explained for pls and block.pls objects

    # we display explained variance when only 1 block is plotted
    if (length(blocks) == 1 && rep.space != "XY-variate") {
      if (!is.null(object$prop_expl_var)) {
        if (style == "3d") {
          inf = 100 * round(object$prop_expl_var[[blocks]][c(comp[1], comp[2], comp[3])], 2)
        } else {
          if (any(class.object %in% object.mint)) {
            if (blocks %in% c("X", "Y")) {  # means that study == "global"
              inf = 100 * round(object$prop_expl_var[[blocks]]$"all data"[c(comp[1], comp[2])], 2)
            } else {
              ## PATCH for MINT models: safely coerce blocks to character
              blocks_char <- as.character(blocks)
              if (!is.null(object$prop_expl_var[[blocks.init]])) {
                # Check if we can access the explained variance
                if (blocks_char %in% names(object$prop_expl_var[[blocks.init]])) {
                  var_data <- object$prop_expl_var[[blocks.init]][[blocks_char]]
                  if (!is.null(var_data) && is.numeric(var_data)) {
                    inf = 100 * round(var_data[c(comp[1], comp[2])], 2)
                  } else {
                    inf = c(NA, NA)
                  }
                } else {
                  inf = c(NA, NA)
                }
              } else {
                inf = c(NA, NA)
              }
            }
          } else {
            inf = 100 * round(object$prop_expl_var[[blocks]][c(comp[1], comp[2])], 2)
          }
        }
      } else {
        inf = rep("NC", 3)
      }

        
        # for future development: if a label with explained variance for each blocks :
        #   inf = lapply(object$prop_expl_var, function(x){round(x[c(comp[1], comp[2])], 2)})
        #   lapply(inf, function(x){paste0("variate ", comp[2], ": ", x[comp[2]], "% expl. var")})}
        
        
        if (is.null(X.label))
        {
          percentage = paste0(inf[1], "% expl. var")
          if (rep.space == "multi")
            X.label = paste0("variate ", comp[1], ": ", percentage)
          if (rep.space == "X-variate")
            X.label = paste0("X-variate ", comp[1], ": ", percentage)
          if (rep.space == "Y-variate")
            X.label = paste0("Y-variate ", comp[1], ": ", percentage)
          if (rep.space == "XY-variate")
            X.label = paste0("XY-variate ", comp[1], ": ", percentage)
        }
        if (is.null(Y.label))
        {
          percentage = paste0(inf[2], "% expl. var")
          if (rep.space == "multi")
            Y.label = paste0("variate ", comp[2], ": ", percentage)
          if (rep.space == "X-variate")
            Y.label = paste0("X-variate ", comp[2], ": ", percentage)
          if (rep.space == "Y-variate")
            Y.label = paste0("Y-variate ", comp[2], ": ", percentage)
          if (rep.space == "XY-variate")
            Y.label = paste0("XY-variate ", comp[2], ": ", percentage)
        }
        if (is.null(Z.label)&&style == "3d")
        {
          percentage = paste0(inf[3], "% expl. var")
          if (rep.space == "multi")
            Z.label = paste0("variate ", comp[3], ": ", percentage)
          if (rep.space == "X-variate")
            Z.label = paste0("X-variate ", comp[3], ": ", percentage)
          if (rep.space == "Y-variate")
            Z.label = paste0("Y-variate ", comp[3], ": ", percentage)
          if (rep.space == "XY-variate")
            Z.label = paste0("XY-variate ", comp[3], ": ", percentage)
        }
      }
      
      
      
      if (is.null(X.label))
      {
        if (rep.space == "multi")
          X.label = paste0("variate ", comp[1])
        if (rep.space == "X-variate")
          X.label = paste("X-variate", comp[1])
        if (rep.space == "Y-variate")
          X.label = paste("Y-variate", comp[1])
        if (rep.space == "XY-variate")
          X.label = paste("XY-variate", comp[1])
      }
      if (is.null(Y.label))
      {
        if (rep.space == "multi")
          Y.label = paste("variate", comp[2])
        if (rep.space == "X-variate")
          Y.label = paste("X-variate", comp[2])
        if (rep.space == "Y-variate")
          Y.label = paste("Y-variate", comp[2])
        if (rep.space == "XY-variate")
          Y.label = paste("XY-variate", comp[2])
      }
      if (is.null(Z.label)&&style == "3d")
      {
        if (rep.space == "multi")
          Z.label = paste("variate", comp[3])
        if (rep.space == "X-variate")
          Z.label = paste("X-variate", comp[3])
        if (rep.space == "Y-variate")
          Z.label = paste("Y-variate", comp[3])
        if (rep.space == "XY-variate")
          Z.label = paste("XY-variate", comp[3])
      }
      
    }
    #-- End: Retrieve variates from object
    out = list(x = x, y = y, z = z, X.label = X.label, Y.label = Y.label, Z.label = Z.label)
    
    
  }


#' @importFrom ellipse ellipse
shape.input.plotIndiv <- 
  function(object,
           n,
           #number of total samples
           blocks = NULL,
           # to choose which block data to plot, when using GCCA module
           x,
           y,
           z,
           ind.names = TRUE,
           group,
           # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
           col.per.group,
           style = "ggplot2",
           # can choose between graphics, 3d, lattice or ggplot2
           study = "global",
           ellipse = FALSE,
           ellipse.level = 0.95,
           centroid = FALSE,
           star = FALSE,
           title = NULL,
           xlim = NULL,
           ylim = NULL,
           col,
           cex,
           pch,
           display.names,
           plot_parameters
  )
  {
    
    class.object = class(object)
    object.mint = c("mint.pls", "mint.spls", "mint.plsda", "mint.splsda")
    
    size.title = plot_parameters$size.title
    size.subtitle = plot_parameters$size.subtitle
    size.xlabel = plot_parameters$size.xlabel
    size.ylabel = plot_parameters$size.ylabel
    size.axis = plot_parameters$size.axis
    size.legend = plot_parameters$size.legend
    size.legend.title = plot_parameters$size.legend.title
    legend.title = plot_parameters$legend.title
    legend.title.pch = plot_parameters$legend.title.pch
    legend.position = plot_parameters$legend.position
    point.lwd = plot_parameters$point.lwd
    
    # --------------------------------------------------------------------------------------
    #           need class.object whether it's DA
    # --------------------------------------------------------------------------------------
    
    #-- Define group
    missing.group = TRUE
    if (missing(group) & is(object, "DA"))
    {
      group = object$Y#factor(map(object$ind.mat), labels = object$names$colnames$Y)
      object$ind.mat = unmap(group) # added in v6 cause $ind.mat is the scaled (if scale = TRUE) version of ind.mat( = unmap(Y))
      missing.group = FALSE #not user defined
    } else if (!missing(group)) {
      missing.group = FALSE
      if (!is.factor(group)) {
        group[is.na(group)] <- "missing group (NA)" ## create a character so it remains in factor levels
        group = as.factor(group)
      }
      
      
      object$ind.mat = unmap(group)
      
      if (length(group) !=  n)
        stop("Length of 'group' should be of length ", n, ", the sample size of your data")
    } else {
      if (star || centroid || ellipse)
        warning('star , ellipse and centroid work only if !group == NULL')
      star = centroid = ellipse = FALSE
      group = factor(rep("No group", n))
      object$ind.mat = unmap(group)
    }
    
    
    
    # --------------------------------------------------------------------------------------
    #           independent from class.object
    # --------------------------------------------------------------------------------------
    
    #at this stage, we have a 'group' - user defined or DA, or by default 1 single group
    
    # if there are groups
    if(!missing.group){
      # if col set, we use col as a col.per.group
      if(!missing(col)) {
        if(length(col) !=  nlevels(group))
          stop("Length of 'col' should be of length ", nlevels(group), " (the number of groups).")
        col.per.group = col
        # if col not set, use default colours
        } else {
            if (nlevels(group) < 10) {
            # only 10 colors in color.mixo
            col.per.group = color.mixo(1:nlevels(group))
            } else {
            # use color.jet
            col.per.group = color.jet(nlevels(group))
            }
        }
      
      # make a vector of one color per sample from col.per.group
      levels.color = vector(, n)
      if (length(col.per.group) !=  n)
      {
        for (i in 1 : nlevels(group))
          levels.color[group == levels(group)[i]] = col.per.group[i]
      } else {
        levels.color = col.per.group
      }
      
    } else { #missing group, we require a 'col' of length n (or repeated to length n) and not a 'col.per.group'
      # col creates a group argument, which creates a col.per.group (levels from 'col')

      if(!missing(col))
      {
        if (length(col) > n)
          stop("Length of 'col' should be of length inferior or equal to ", n, ".")
        
        col = rep(col, ceiling(n/length(col)))[1 : n]
        group = factor(col)
        levels.color = col
        col.per.group = levels(group)
        object$ind.mat = unmap(group)
      } else { # no group, no color => a single color by default
        col.per.group = color.mixo(1)
        levels.color = rep(col.per.group, n)
        col = levels.color
      }
      
    }
    # 'group' and 'col' are always the same factor, but different values
    # here we have a group, a col.per.group (length nlevels(group)) and a levels.color (length n)
    
    
    #-- cex argument
    if (missing(cex))
    {
      if (style == "ggplot2")
      {
        cex = rep(3, n)
        cex = cex[as.factor(group)]
      } else {
        cex = rep(1, n)
        cex = cex[as.factor(group)]
      }
    } else {
      if (length(cex) == 1)
      {
        cex = rep(cex, n)
        cex = cex[as.factor(group)]
      } else if (length(cex) == length(unique(group))) {
        cex = cex[as.factor(group)]
      } else {
        if(length(unique(group))>1)
        {
          stop("Length of 'cex' should be either 1 or ",length(unique(group)), ".")
        }else{ # one group
          stop("Length of 'cex' should be 1.")
        }
        #cex = rep(cex, ceiling(n/length(cex)))[1 : n]
      }
    }
    
    if (ellipse)
    {
      # removing calculation for classes with only one sample
      
      nlevels.class = 1 : ncol(object$ind.mat)
      ind.unique = which(apply(object$ind.mat, 2, sum) == 1)
      if(length(ind.unique) > 0)
      {
        nlevels.class = nlevels.class[-ind.unique]
      }
      
      #-- Start: Computation ellipse
      min.ellipse = max.ellipse = xlim.min = xlim.max = ylim.min = ylim.max = list()
      ind.gp = matrice = cdg = variance = list()
      ind.gp = lapply(1 : ncol(object$ind.mat), function(x){which(object$ind.mat[, x] == 1)})
      matrice = lapply(1 : length(x), function(z1) {lapply(ind.gp, function(z2){matrix(c(x[[z1]][z2], y[[z1]][z2]), ncol = 2)})})
      cdg = lapply(1 : length(x), function(z){ lapply(matrice[[z]], colMeans)})
      variance = lapply(1 : length(x), function(z){lapply(matrice[[z]], var)})
      coord.ellipse = lapply(1 : length(x), function(z1){ lapply(nlevels.class, function(z2){ellipse(variance[[z1]][[z2]],
                                                                                                     centre = cdg[[z1]][[z2]],
                                                                                                     level = ellipse.level)})})
      max.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, max)})})
      min.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, min)})})
      #-- End: Computation ellipse
      
      if (is.null(xlim))
        xlim = lapply(1 : length(x), function(z) {c(min(x[[z]], min.ellipse[[z]][1, ]), max(x[[z]], max.ellipse[[z]][1, ]))})
      if (is.null(ylim))
        ylim = lapply(1 : length(x), function(z) {c(min(y[[z]], min.ellipse[[z]][2, ]), max(y[[z]], max.ellipse[[z]][2, ]))})
      if(style == "ggplot2") # no lists, a single vector of two values is expected
      {
        temp = matrix(unlist(xlim),ncol=2,byrow=TRUE)
        xlim = c(min(temp[,1]),max(temp[,2]))
        temp = matrix(unlist(ylim),ncol=2,byrow=TRUE)
        ylim = c(min(temp[,1]),max(temp[,2]))
      }
      
    }
    # no need for xlim and ylim as ggplot2, lattice and graphics are good without by default
    #else {
    #  if (is.null(xlim))
    #    xlim = lapply(1 : length(x), function(z) {c(min(x[[z]]), max(x[[z]]))})
    #    if (is.null(ylim))
    #    ylim = lapply(1 : length(x), function(z) {c(min(y[[z]]), max(y[[z]]))})
    #}
    
    
    # --------------------------------------------------------------------------------------
    #           not independent from class.object: for the title of the plot, either "PlotIndiv" or "block:.."
    # --------------------------------------------------------------------------------------

    #-- pch argument
    if (missing(pch) & !any(class.object%in%object.mint))
    {
      
      if (style == "3d")
      {
        pch = unlist(lapply(1 : length(length(levels(group))), function(x){rep(c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")[x], length(group == x))}))
        
      } else {
        pch = as.numeric(group)
      }
      
    }else if (any(class.object%in%object.mint)) {
      if (missing(pch))
      {
        # a pch per study, forced
        pch = as.numeric(object$study)
        pch.levels = pch
      } else {
        if (length(pch)!= length(object$study))
          stop("'pch' needs to be of length 'object$study' as each of 'pch' represents a specific study", call. = FALSE)
      }
      
    } else {
      if (style == "3d")
      {
        if (!all(unlist(pch) %in% c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")))
          stop("pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.",
               call. = FALSE)
      }
      
      if (length(pch) == 1)
      {
        pch = rep(pch, n)
        
      } else if (length(pch) > n) {
        stop("Length of 'pch' should be of length inferior or equal to ", n, ".")
      } else if (length(pch) == length(unique(group)) & length(pch)!=n ) { # prevent from reordering pch when 1 group per sample (length(pch)=length(group)=n)
        pch = pch[as.factor(group)]
        
      } else {
        pch = rep(pch, ceiling(n/length(pch)))[1 : n]
      }
      # if pch is given and ind.names is TRUE, pch takes over
      display.names = FALSE
    }
    
    
    
    # constructing data.frame df
    if (any(study == "global"))# | length(study) == 1)
    {
      #-- Start: data set
      df = list()
      if (style == "3d")
      {
        for (i in 1 : length(x))
          df[[i]] = data.frame(x = x[[i]], y = y[[i]], z = z[[i]], group = group)
      } else {
        for (i in 1 : length(x))
          df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group)
      }
      
      title.save = title # to use for ellipse
      #if (any(class.object %in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda")) &
      if(length(blocks) == 1 & !any(class(object)%in%c(object.mint, "sgcca", "rgcca"))) # add blocks == 1 to allow "multi" with plsda
      {
        if (is.null(title))
        {
          df = data.frame(do.call(rbind, df), "Block" = "PlotIndiv")
          if (style %in%c("graphics"))
            title = "PlotIndiv" # to add title to graphics
          
        } else {
          df = data.frame(do.call(rbind, df), "Block" = title)
          if (style %in%c("ggplot2", "lattice"))
            title = NULL # to avoid double title
          
        }
        
        # no subtitle with these objects
        if(size.title != rel(2)) # rel(2) is the default
          size.subtitle = size.title
        
        
        df$Block = as.factor(df$Block)
      } else {
        df = data.frame(do.call(rbind, df), "Block" = paste0("Block: ", unlist(lapply(1 : length(df), function(z){rep(blocks[z], nrow(df[[z]]))}))))
        df$Block = factor(df$Block, levels = unique(df$Block))
      }
      
      if (style == "3d")
      {
        names(df)[1:3] = c("x", "y", "z")
      } else {
        names(df)[1:2] = c("x", "y")
      }
      
      if (display.names)
        df$names = rep(ind.names, length(x))
      
      df$pch = pch; df$pch.levels = pch
      df$cex = cex
      #df$col.per.group = levels.color#[group] #FR: don't understand what is that changing as levels.color is already group?
      df$col = levels.color#as.character(col)
      
      
      if (centroid == TRUE || star == TRUE)
      {
        df = cbind(df, rep(0, nrow(df)))
        n = ncol(df)
        df = cbind(df, rep(0, nrow(df)))
        for (i in 1:nlevels(group))
        {
          if (length(x)>1)
          {
            for (k in 1 : length(x))
            {
              x0 = mean(df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "x"])
              y0 = mean(df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]) , "y"])
              df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), n] = x0
              df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), n+1] = y0
              names(df)[c(ncol(df)-1, ncol(df))] = c("x0", "y0")
            }
          } else {
            x0 = mean(df[df$group == levels(group)[i] , "x"])
            y0 = mean(df[df$group == levels(group)[i]  , "y"])
            df[df$group == levels(group)[i] , n] = x0
            df[df$group == levels(group)[i] , n+1] = y0
            names(df)[c(ncol(df)-1, ncol(df))] = c("x0", "y0")
          }
          
        }
      }
      
      if (ellipse == TRUE)
      {
        df.ellipse = data.frame(do.call("rbind", lapply(1 : length(x), function(k){do.call("cbind", coord.ellipse[[k]])})), "Block" = paste0("Block: ", rep(blocks, each = 100)))
        if(length(ind.unique) > 0)
        {
          names(df.ellipse)[1 : (2*length(nlevels.class))] = paste0("Col", c(1 : (2*nlevels(group)))[-c(2*(ind.unique-1)+1,2*ind.unique)])
        } else {
          names(df.ellipse)[1 : (2*nlevels(group))] = paste0("Col", c(1 : (2*nlevels(group))))
          
        }
        df.ellipse$ellipse.level = ellipse.level
      } else {
        df.ellipse = NULL
      }
      
      if(ellipse == TRUE && length(blocks) == 1 & !any(class(object)%in%c(object.mint, "sgcca", "rgcca"))) # add blocks == 1 to allow "multi" with plsda
        #if (ellipse == TRUE && any(class.object %in% c("ipca", "sipca", "pca", "spca", "prcomp", "splsda", "plsda"))& length(blocks) == 1& !any(class(object)%in%object.mint))
      {
        if (is.null(title.save))
        {
          df.ellipse$Block = "PlotIndiv"
        } else {
          df.ellipse$Block = title.save
        }
      }
      
      
      
      # pch.legend = NULL
      # for (i in 1:nlevels(group))
      # {
      #   pch.legend = c(pch.legend, df[df$group == levels(group)[i], ]$pch)
      # }
      df$pch.legend <- df$pch
      
    } else {
      
      #mint object
      #display.names = FALSE # so far ggplot and lattice require a unique vector of names. when the code changes, we can use ind.names (list)
      # Handle both numeric position and study names
      if (is.numeric(study)) {
        # Convert numeric position to actual study name
        if (any(study > length(levels(object$study)))) {
          stop("'study' position out of bounds. Maximum study position is ", length(levels(object$study)))
        }
        study <- levels(object$study)[study]
      } else if (!all(study %in% levels(object$study))) {
        stop("'study' must be either a numeric position or one of: ", paste(levels(object$study), collapse = ", "))
      }

      # For MINT models, we need to handle both blocks and study
      if (length(blocks) == 1) {
        # Single block case
        group.mint = split(group, object$study)[study]
        group = as.factor(unlist(group.mint))
        pch = as.vector(unlist(split(pch, object$study)[study]))
        cex = as.vector(unlist(split(cex, object$study)[study]))
        pch.levels = as.vector(unlist(split(pch.levels, object$study)[study]))
        col.per.group.mint = as.vector(unlist(split(levels.color, object$study)[study]))
      } else {
        # Multiple blocks case
        group.mint = lapply(blocks, function(block) {
          split(group, object$study)[study]
        })
        group = as.factor(unlist(lapply(group.mint, unlist)))
        pch = as.vector(unlist(lapply(blocks, function(block) {
          unlist(split(pch, object$study)[study])
        })))
        cex = as.vector(unlist(lapply(blocks, function(block) {
          unlist(split(cex, object$study)[study])
        })))
        pch.levels = as.vector(unlist(lapply(blocks, function(block) {
          unlist(split(pch.levels, object$study)[study])
        })))
        col.per.group.mint = as.vector(unlist(lapply(blocks, function(block) {
          unlist(split(levels.color, object$study)[study])
        })))
      }
      
      #-- Start: data set
      df = list()
      if (style == "3d") {
        for (i in 1 : length(x)) {
          if (length(blocks) == 1) {
            df[[i]] = data.frame(x = x[[i]], y = y[[i]], z = z[[i]], group = group)
          } else {
            df[[i]] = data.frame(x = x[[i]], y = y[[i]], z = z[[i]], group = group.mint[[i]])
          }
        }
      } else {
        for (i in 1 : length(x)) {
          if (length(blocks) == 1) {
            df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group)
          } else {
            if(length(study)>1){
              df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group.mint[[1]][[i]])
            } else {df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group.mint[[i]])}
            
          }
        }
      }
      
      df = data.frame(do.call(rbind, df), "Block" = paste0("Study: ", unlist(lapply(1 : length(df), function(z){rep(blocks[z], nrow(df[[z]]))}))))
      df$Block = factor(df$Block, levels = unique(df$Block))
      
      #print(df)
      
      if (style == "3d")
      {
        names(df)[1:3] = c("x", "y", "z")
      } else {
        names(df)[1:2] = c("x", "y")
      }
      
      # no names for MINT object
      #if (display.names)
      #df$names = rep(ind.names, length(x))
      
      # have added if statements here to deal with MINT cases when length of pch, pch.level, cex and col.per.group.mint dont match number of data points (rows in df)
      # in case when one study is plotted in 'multi' rep space, these params are half length of data points (which are across X and Y variate space), so need to be repeated twice
      # in case where multiple studies are plotted, these params are too long (nrow(df) repeated by number of studies) so need to be shorterned
      # any other deviation of length of these vectors as compared to number of data points is caught as an error
      if (length(pch) != nrow(df)){
        if (length(pch) == nrow(df)/2){pch <- c(pch, pch)}
        else if (length(pch) == length(study)*nrow(df) & length(pch) %% length(study) == 0 & all(pch[1:(length(pch)/length(study))] == pch[(length(pch)/length(study) + 1):length(pch)]))
        {pch <- pch[1:(length(pch) / length(study))]} 
        else {stop("unexpected error occured! 'pch' length does not match number of samples")}
      }
      df$pch = pch
      
      if (length(pch.levels) != nrow(df)){
        if (length(pch.levels) == nrow(df)/2){pch.levels <- c(pch.levels, pch.levels)}
        else if (length(pch.levels) == length(study)*nrow(df) & length(pch.levels) %% length(study) == 0 & all(pch.levels[1:(length(pch.levels)/length(study))] == pch.levels[(length(pch.levels)/length(study) + 1):length(pch.levels)])){
          pch.levels <- pch.levels[1:(length(pch.levels) / length(study))]
        } else {stop("unexpected error occured! 'pch.levels' length does not match number of samples")}
      }
      df$pch.levels = pch.levels
      
      if (length(cex) != nrow(df)){
        if (length(cex) == nrow(df)/2){cex <- c(cex, cex)}
        else if (length(cex) == length(study)*nrow(df) & length(cex) %% length(study) == 0 & all(cex[1:(length(cex)/length(study))] == cex[(length(cex)/length(study) + 1):length(cex)])){
          cex <- cex[1:(length(cex) / length(study))]
        } else {stop("unexpected error occured! 'cex' length does not match number of samples")}
      }
      df$cex = cex
      
      if (length(col.per.group.mint) != nrow(df)){
        if (length(col.per.group.mint) == nrow(df)/2){col.per.group.mint <- c(col.per.group.mint, col.per.group.mint)}
        else if (length(col.per.group.mint) == length(study)*nrow(df) & length(col.per.group.mint) %% length(study) == 0 & all(col.per.group.mint[1:(length(col.per.group.mint)/length(study))] == col.per.group.mint[(length(col.per.group.mint)/length(study) + 1):length(col.per.group.mint)])){
          col.per.group.mint <- col.per.group.mint[1:(length(col.per.group.mint) / length(study))]
        } else {stop("unexpected error occured! 'col.per.group.mint' length does not match number of samples")}
      }
      df$col = col.per.group.mint#as.character(col)
      
      
      pch.legend = NULL
      for (i in 1:nlevels(group))
        pch.legend = c(pch.legend, df[df$group == levels(group)[i], ]$pch)
      
      df$pch.legend = pch.legend
      df.ellipse = NULL # no ellipse so far
    }
    
    if (any(study == "global"))
      study = levels(object$study)
    
    # match study with names of the study
    study.ind = match(study, levels(object$study))
    
    #print(df)
    plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel,
                           size.axis = size.axis, size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
                           legend.title.pch = legend.title.pch, legend.position = legend.position, point.lwd = point.lwd)
    
    out = list(df = df, study.ind = study.ind, df.ellipse = df.ellipse, col.per.group = col.per.group, title = title, display.names = display.names, xlim = xlim, ylim = ylim, ellipse = ellipse, centroid = centroid, star = star, plot_parameters = plot_parameters)
  }


# --------------------------------------------------------------------------------------
# Internal helpers functions to run some plots functions
# --------------------------------------------------------------------------------------

#-- Function to display an error message (used for the parameters var.names, cex, col, pch and font)
stop.message <- function(argument, data){
  if (length(data) == 1) {
    count.data = sapply(data, length)
  } else {
    count.data = paste(paste(sapply(data[-length(data)], length), collapse =  ", "), length(data[[length(data)]]), sep = " and ")
  }
  stop(argument, " must be either a vector of length ", length(data),
       " or a list of ", length(data), " vector components of length ", count.data, " respectively.",call.= FALSE)
}

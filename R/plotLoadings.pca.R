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

## ------------------------ PCA, sPCA, IPCA, sIPCA ------------------------ ##

#' @rdname plotLoadings
#' @method plotLoadings pca
#' @export
plotLoadings.pca <- 
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
           size.title = 2, 
           size.axis = 0.7,
           X.label = NULL,
           Y.label = NULL,
           size.labs = 1,
           ...) {
    
    # Input checks
    if (!is.list(object) || is.null(object$X))
      stop("The object must be a valid PCA object containing the 'X' matrix.")
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
    if (!is.numeric(size.labs) || size.labs <= 0)
      stop("'size.labs' must be a positive numeric value.")
    if (!is.numeric(size.axis) || size.axis <= 0)
      stop("'size.axis' must be a positive numeric value.")
    
    if (style == 'ggplot2' && !is.null(layout))
      warning("layout is ignored for PCA objects when style is set to 'ggplot2'.")
    
    # check for inappropriate args
    extra_args <- list(...)
    if ("name.var.complete" %in% names(extra_args)) {
      warning("'name.var.complete' argument is deprecated")
    }
    name.var.complete <- FALSE
  
  # -- input checks
  object$names$blocks <- 'X'
  check <- check.input.plotLoadings(object = object, block = 'X', size.name = size.name, title = title, col = col, name.var = name.var, xlim = xlim)
  
  col <- check$col
  size.name <- check$size.name
  block <- 1
  object$X <- list(X = object$X)
  xlim <- check$xlim
  
  # -- layout
  res <- layout.plotLoadings(layout = layout, plot = TRUE, legend = FALSE, block = block)
  reset.mfrow <- res$reset.mfrow
  opar <- res$opar
  omar <- par('mar')
  
  # -- title
  default.title <- paste0("Loadings on comp ", comp)
  plot.title <- ifelse(is.null(title), default.title, title)
  
  res <- get.loadings.ndisplay(object = object, comp = comp, block = block, name.var = name.var, name.var.complete = name.var.complete, ndisplay = ndisplay)
  X <- res$X
  colnames.X <- res$colnames.X
  value.selected.var <- res$value.selected.var
  
  df <- data.frame(importance = value.selected.var, names = colnames.X)
  
  # -- Plotting based on style
  if (style == 'ggplot2') {
    p <- ggplot(df, aes(x = reorder(names, -abs(importance)), y = importance)) +
      geom_bar(stat = 'identity', fill = col, color = border) +
      theme_minimal() +
      theme(axis.text.y = element_text(size = size.name * 8),
            axis.text.x = element_text(size = size.axis * 8),
            axis.title.x = element_text(size = size.labs * 8),
            axis.title.y = element_text(size = size.labs * 8),
            plot.title = element_text(face = "bold", hjust = 0.5, size = size.title * 8)) +
      labs(title = plot.title, y = X.label, x = Y.label)
    
    # optionally control x axis
    if (!is.null(xlim)) {p <- p + scale_y_continuous(limits = xlim, expand = c(0,0))}
    
    # flip and print plot
    p <- p + coord_flip()
    print(p)
    
  } else {
    par(mar = c(4, max(7, max(sapply(colnames.X, nchar), na.rm = TRUE) / 3), 4, 2))
    .plotLoadings_barplot(height = df$importance, col = col, names.arg = colnames.X, cex.name = size.name, border = border, 
                          xlim = xlim, xlab = X.label, ylab = Y.label, cex.lab = size.labs, cex.axis = size.axis)
    title(title %||% paste0('Loadings on comp ', comp), cex.main = size.title)
  }
  
  if (reset.mfrow) par(opar)
  par(mar = omar)
  
  return(invisible(df))
}
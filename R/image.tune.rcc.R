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

#' Plot the cross-validation score.
#' 
#' This function provide a image map (checkerboard plot) of the
#' cross-validation score obtained by the \code{tune.rcc} function.
#' 
#' \code{plot.tune.rcc} creates an image map of the matrix \code{object$mat}
#' containing the cross-validation score obtained by the \code{tune.rcc}
#' function. Also a color scales strip is plotted.
#' 
#' @param x object returned by \code{tune.rcc}.
#' @param col a character string specifying the colors function to use:
#' \code{\link{terrain.colors}}, \code{\link{topo.colors}},
#' \code{\link{rainbow}} or similar functions. Defaults to
#' \code{\link{heat.colors}}.
#' @param ... not used currently.
#' @return none
#' @author Sébastien Déjean, Ignacio González, Kim-Anh Le Cao, Al J Abadi
#' @seealso \code{\link{tune.rcc}}, \code{\link{image}}.
#' @keywords dplot hplot
#' @examples
#' 
#' data(nutrimouse)
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' 
#' ## this can take some seconds
#' cv.score <- tune.rcc(X, Y, validation = "Mfold")
#' plot(cv.score)
#' 
#' # image(cv.score) # same result as plot()
#' @rdname image.tune.rcc
#' @method image tune.rcc
#' @export
image.tune.rcc <- function(x, col = heat.colors, ...) 
    
{
    
    opar = par(no.readonly = TRUE)
    
    grid1 = x$grid1
    grid2 = x$grid2
    mat = x$mat
    nlevel = min(255, length(c(unique(mat))))
    if (nlevel / 2 - trunc(nlevel / 2) == 0 & nlevel > 6) nlevel = nlevel + 1
    col = col(nlevel)
    
    def.par = par(no.readonly = TRUE)
    layout(matrix(c(2, 1), ncol = 2, nrow = 1, byrow = FALSE),
           widths = c(1, 0.18))
    
    #-- layout 1 --#
    min.mat = min(mat)
    max.mat = max(mat)
    par(pty = "m", mai = c(1.1, 0.1, 0.95, 0.5))
    z = seq(min.mat, max.mat, length = nlevel)
    breaks = seq(min.mat, max.mat, length = nlevel + 1)
    z = matrix(z, nrow = 1)
    
    image(z, col = col, 
          zlim = c(min.mat, max.mat), oldstyle = TRUE,
          xaxt = "n", yaxt = "n", breaks = breaks)
    box()
    
    par(usr = c(0, 1, min.mat, max.mat))
    binwidth = (max.mat - min.mat) / nlevel
    midpoints = seq(min.mat + binwidth/2, max.mat - binwidth/2, by = binwidth)
    
    if (nlevel <= 6)
    {
        axis(4, at = midpoints, labels = round(midpoints, 3))
    }else{
        binwidth = seq(1, nlevel, by = round(nlevel / 5))
        if (binwidth[length(binwidth)] != nlevel){
            binwidth[length(binwidth) + 1] = nlevel
        }
        
        axis(4, at = midpoints[binwidth], labels = round(midpoints[binwidth],
                                                         3))
        
    }
    
    
    
    #-- layout 2 --#
    par(mai = c(0.9, 0.85, 0.75, 0.2))
    image(grid1, grid2, mat, col = col,
          xlab = expression(lambda[1]), ylab=expression(lambda[2]),
          main = expression(CV(lambda[1], lambda[2])), axes = FALSE,
          zlim = c(min.mat, max.mat), oldstyle = TRUE)
    
    
    if (length(grid1) > 10) {
        grid1 = seq(min(grid1), max(grid1), length = 11)
    }
    if (length(grid2) > 10) {
        grid2 = seq(min(grid2), max(grid2), length = 11)
    }
    
    axis(1, at = grid1, labels = as.character(round(grid1, 4)))
    axis(2, at = grid2, labels = as.character(round(grid2, 4)))
    box()
    
    par(opar) #reset par
    
}

#' @rdname image.tune.rcc
#' @method plot tune.rcc
#' @export
plot.tune.rcc <- image.tune.rcc

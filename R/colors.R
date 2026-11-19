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

## --------------------------- mixOmics Colors ---------------------------- ##
#' Color Palette for mixOmics
#' 
#' The functions create a vector of \code{n} "contiguous" colors (except the
#' \code{color.mixo} which are colors used internally to fit our logo colors).
#' 
#' The function \code{color.jet(n)} create color scheme, beginning with dark
#' blue, ranging through shades of blue, cyan, green, yellow and red, and
#' ending with dark red. This colors palette is suitable for displaying ordered
#' (symmetric) data, with \code{n} giving the number of colors desired.
#' 
#' @param n an integer, the number of colors \eqn{(\geq 1)} to be in the
#' palette.
#' @param alpha a numeric value between 0 and 1 for alpha channel (opacity).
#' @param num.vector for \code{color.mixo} an integer vector specifying which
#' colors to use in the mixOmics palette (there are only 10 colors available.
#' @return For \code{color.jet(n)}, \code{color.spectral(n)},
#' \code{color.GreenRed(n)} a character vector, \code{cv}, of color names. This
#' can be used either to create a user-defined color palette for subsequent
#' graphics by \code{palette(cv)}, a \code{col=} specification in graphics
#' functions or in \code{par}.
#' 
#' For \code{color.mixo}, a vector of colors matching the mixOmics logo (10
#' colors max.)
#' @author Ignacio Gonzalez, Kim-Anh Lê Cao, Benoit Gautier, Al J Abadi
#' @seealso \code{\link{colorRamp}}, \code{\link{palette}},
#' \code{\link{colors}} for the vector of built-in "named" colors;
#' \code{\link{hsv}}, \code{\link{gray}}, \code{\link{rainbow}},
#' \code{\link{terrain.colors}}, ... to construct colors; and
#' \code{\link{heat.colors}}, \code{\link{topo.colors}} for images.
#' @keywords color
#' @examples
#' 
#' # -----------------------
#' # jet colors
#' # ----------------------
#' par(mfrow = c(3, 1))
#' z <- seq(-1, 1, length = 125)
#' for (n in c(11, 33, 125)) {
#' image(matrix(z, ncol = 1), col = color.jet(n),
#' xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
#' box()
#' par(usr = c(-1, 1, -1, 1))
#' axis(1, at = c(-1, 0, 1))
#' }
#' 
#' \dontrun{
#' # -----------------------
#' # spectral colors
#' # ----------------------
#' par(mfrow = c(3, 1))
#' z <- seq(-1, 1, length = 125)
#' for (n in c(11, 33, 125)) {
#' image(matrix(z, ncol = 1), col = color.spectral(n),
#' xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
#' box()
#' par(usr = c(-1, 1, -1, 1))
#' axis(1, at = c(-1, 0, 1))
#' }
#' 
#' # -----------------------
#' # GreenRed colors
#' # ----------------------
#' par(mfrow = c(3, 1))
#' z <- seq(-1, 1, length = 125)
#' for (n in c(11, 33, 125)) {
#' image(matrix(z, ncol = 1), col = color.GreenRed(n),
#' xaxt = 'n', yaxt = 'n', main = paste('n = ', n))
#' box()
#' par(usr = c(-1, 1, -1, 1))
#' axis(1, at = c(-1, 0, 1))
#' }
#' 
#' # # --------------------------------
#' # mixOmics colors
#' # # -------------------------------
#' data(nutrimouse)
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
#' 
#' my.colors = color.mixo(1:5)
#' my.pch = ifelse(nutrimouse$genotype == 'wt', 16, 17)
#' #plotIndiv(nutri.res, ind.names = FALSE, group = my.colors, pch = my.pch, cex = 1.5)
#' }
#' @name colors
NULL
#' @export
#' @rdname colors
color.mixo <- function(num.vector)
{
    
    if (is.factor(num.vector)) num.vector=as.numeric(num.vector)
    
    if (!is.numeric(num.vector)){
        stop(paste("num.vector has to be numeric", call. = FALSE))
    }
    
    # these are the colors in the logo (the first 3)
    mixo.gray = gray.colors(1, start = 0.76, gamma = 1)
    
    mixo.col = c('#388ECC', # mixOmics logo blue
                 '#F68B33', # mixOmics logo orange
                 mixo.gray, # mixOmics logo grey
                 '#009E73', # shiny dark green
                 '#CC79A7', # shiny purple/pink
                 '#F0E442', #shiny yellow
                 'black',
                 '#D55E00', #shiny dark orange
                 '#0072B2', #shiny dark blue
                 '#999999'  # shiny grey
    )
    
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    n = length(num.vector)
    #-- n: check that there are more colors available than requested
    if (isTRUE(num.vector) > length(mixo.col)){
        stop(paste("We only have a few mix.colors available, n <= ",
                   length(mixo.col)), call. = FALSE)
    }
    
    if (isTRUE(!is.finite((num.vector))) ||  (n < 1)){
        stop("'num.vector' must be an integer vector with positive values.",
             call. = FALSE)
    }
    #-- end checking --#
    #------------------#
    
    return(mixo.col[num.vector])
}

## ---------------------- Green Black red Gradients ----------------------- ##
#' @export
#' @rdname colors
color.GreenRed <-
    function (n, alpha = 1)
    {
        #-- checking general input parameters -------------------------------------#
        #--------------------------------------------------------------------------#
        
        #-- n
        if (length(n) > 1 || !is.finite(n))
            stop("'n' must be an integer positive value.", call. = FALSE)
        
        if (n < 1)
            stop("'n' must be an integer positive value.", call. = FALSE)
        
        #-- alpha
        if (length(alpha) > 1 || !is.finite(alpha))
            stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
        
        if (alpha < 0 || alpha > 1)
            stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
        
        alpha = round(255 * alpha)
        
        #-- end checking --#
        #------------------#
        
        ramp = colorRampPalette(c("green", "darkgreen", "black", "darkred", "red"))
        ramp = ramp(101)
        green = ramp[1:43]
        red = ramp[59:101]
        ramp = colorRamp(c(green, "black", red), space = "Lab")
        
        rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
    }
## ------------------------------ Jet Colors ------------------------------ ##
#' @export
#' @rdname colors
color.jet <-
    function (n, alpha = 1)
    {
        #-- checking general input parameters -------------------------------------#
        #--------------------------------------------------------------------------#
        
        #-- n
        if (length(n) > 1 || !is.finite(n))
            stop("'n' must be an integer positive value.", call. = FALSE)
        
        if (n < 1)
            stop("'n' must be an integer positive value.", call. = FALSE)
        
        #-- alpha
        if (length(alpha) > 1 || !is.finite(alpha))
            stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
        
        if (alpha < 0 || alpha > 1)
            stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
        
        alpha = round(255 * alpha)
        
        #-- end checking --#
        #------------------#
        
        ramp = colorRamp(c("#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
                           "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
                           "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
                           "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
                           "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
                           "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
                           "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
                           "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
                           "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
                           "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
                           "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
                           "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
                           "#AF0000", "#9F0000", "#8F0000", "#800000"),
                         space = "Lab")
        
        rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
    }
## --------------------------- Spectral Colors ---------------------------- ##
#' @export
#' @rdname colors
color.spectral <-
    function (n, alpha = 1) 
    {
        #-- checking general input parameters -------------------------------------#
        #--------------------------------------------------------------------------#
        
        #-- n
        if (length(n) > 1 || !is.finite(n))
            stop("'n' must be an integer positive value.", call. = FALSE)
        
        if (n < 1)
            stop("'n' must be an integer positive value.", call. = FALSE)
        
        #-- alpha
        if (length(alpha) > 1 || !is.finite(alpha))
            stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
        
        if (alpha < 0 || alpha > 1)
            stop("'alpha' must be an numeric value in the range [0, 1].", call. = FALSE)
        
        alpha = round(255 * alpha)
        
        #-- end checking --#
        #------------------#
        
        ramp = colorRamp(c("#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598", 
                           "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", 
                           "#9E0142"), space = "Lab")
        
        rgb(ramp(seq(0, 1, length = n)), alpha = alpha, maxColorValue = 255)
    }

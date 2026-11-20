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

# -----------------------------------------
# --- small example but see help file
# ----------------------------------------

#data(liver.toxicity)
#X <- liver.toxicity$gene
#Y <- liver.toxicity$clinic
#imgCor(X=X[,1:10],Y=Y[,1:4])

#' Image Maps of Correlation Matrices between two Data Sets
#' 
#' Display two-dimensional visualizations (image maps) of the correlation
#' matrices within and between two data sets.
#' 
#' If \code{type="combine"}, the correlation matrix is computed of the combined
#' matrices \code{cbind(X, Y)} and then plotted. If \code{type="separate"},
#' three correlation matrices are computed, \code{cor(X)}, \code{cor(Y)} and
#' \code{cor(X,Y)} and plotted separately on a device. In both cases, a color
#' correlation scales strip is plotted.
#' 
#' The correlation matrices are pre-processed before calling the \code{image}
#' function in order to get, as in the numerical representation, the diagonal
#' from upper-left corner to bottom-right one.
#' 
#' Missing values are handled by casewise deletion in the \code{imgCor}
#' function.
#' 
#' If \code{X.names = FALSE}, the name of each X-variable is hidden. Default
#' value is \code{TRUE}.
#' 
#' If \code{Y.names = FALSE}, the name of each Y-variable is hidden. Default
#' value is \code{TRUE}.
#' 
#' @param X numeric matrix or data frame \eqn{(n} x \eqn{p)}, the observations
#' on the \eqn{X} variables. \code{NA}s are allowed.
#' @param Y numeric matrix or data frame \eqn{(n} x \eqn{q)}, the observations
#' on the \eqn{Y} variables. \code{NA}s are allowed.
#' @param type character string, (partially) maching one of \code{"combine"} or
#' \code{"separated"}, determining the kind of plots to be produced. See
#' Details.
#' @param X.var.names,Y.var.names logical, should the name of \eqn{X}- and/or
#' \eqn{Y}-variables be shown? If \code{TRUE} (defaults) \code{object$names$X}
#' and/or \code{object$names$Y} are used. Possible character vector with
#' \eqn{X}- and/or \eqn{Y}-variable labels to use.
#' @param sideColors character vector of length two. The color name for
#' horizontal and vertical side bars that may be used to annotate the \eqn{X}
#' and \eqn{Y} correlation matrices.
#' @param interactive.dev Logical. The current graphics device that will be
#' opened is interactive?
#' @param title logical, should the main titles be shown?
#' @param color,xlab,ylab arguments passed to \code{cim}.
#' @param row.cex,col.cex positive numbers, used as \code{cex.axis} in for the
#' row or column axis labeling. The defaults currently only use number of rows
#' or columns, respectively.
#' @param symkey Logical indicating whether the color key should be made
#' symmetric about 0. Defaults to \code{TRUE}.
#' @param keysize positive numeric value indicating the size of the color key.
#' @param margins numeric vector of length two containing the margins (see
#' \code{\link{par}(mar)}) for column and row names respectively.
#' @param lhei,lwid arguments passed to \code{layout} to divide the device up
#' into two rows and two columns, with the row-heights \code{lhei} and the
#' column-widths \code{lwid}.
#' @return NULL (invisibly)
#' @author Ignacio González, Kim-Anh Lê Cao, Florian Rohart, Al J Abadi
#' @seealso \code{\link{cor}}, \code{\link{image}}, \code{\link{color.jet}}.
#' @keywords multivariate dplot
#' @examples
#' 
#' data(nutrimouse)
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' 
#' ## 'combine' type plot (default)
#' imgCor(X, Y)
#' 
#' \dontrun{
#' ## 'separate' type plot
#' 
#' imgCor(X, Y, type = "separate")
#' 
#' ## 'separate' type plot without the name of datas
#' imgCor(X, Y, X.var.names = FALSE, Y.var.names = FALSE, type = "separate")
#' }
#' @export
imgCor <-
    function(X,
             Y,  
             type = "combine", 
             X.var.names = TRUE, 
             Y.var.names = TRUE,  	 
             sideColors = TRUE,
             interactive.dev = TRUE,
             title = TRUE,
             color, row.cex, col.cex,symkey, keysize, 
             xlab, ylab, margins, lhei, lwid) 
    {
        #-- checking general input arguments ---------------------------------------#
        #---------------------------------------------------------------------------#
        
        #-- check that the user did not enter extra arguments
        arg.call = match.call()
        user.arg = names(arg.call)[-1]
        
        err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())), 
                       error = function(e) e)
        
        if ("simpleError" %in% class(err))
            stop(err[[1]], ".", call. = FALSE)
        
        default.arg = c("color", "row.cex", "col.cex","symkey", "keysize", 
                        "xlab", "ylab", "margins", "lhei", "lwid")
        function.arg = c(names(mget(names(formals()), sys.frame(sys.nframe()))),
                         default.arg)
        not.arg = !(user.arg %in% function.arg)
        
        if (any(not.arg)) {
            unused.arg = user.arg[not.arg]
            not.arg = which(not.arg) + 1
            output = rep("", length(not.arg))
            
            for (i in 1:length(not.arg)) {
                output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
            }
            
            output = paste0("(",paste(output, collapse = ", "), ").")
            msg = "unused argument "
            if (length(not.arg) > 1) msg ="unused arguments "  
            stop(msg, output, call. = FALSE)
        }
        
        #-- X and Y matrices
        if (is.data.frame(X)) X = as.matrix(X)
        if (is.data.frame(Y)) Y = as.matrix(Y)
        
        if (!is.matrix(X) || !is.matrix(Y)) 
            stop("'X' and/or 'Y' must be a numeric matrix.", call. = FALSE)
        
        if (!is.numeric(X) || !is.numeric(Y)) 
            stop("'X' and/or 'Y' must be a numeric matrix.", call. = FALSE)
        
        p = ncol(X)
        q = ncol(Y)  
        
        #-- X.var.names
        if (!is.logical(X.var.names)) 
            stop("'X.var.names' must be a logical constant (TRUE or FALSE).",
                 call. = FALSE)
        
        #-- Y.var.names
        if (!is.logical(Y.var.names)) 
            stop("'Y.var.names' must be a logical constant (TRUE or FALSE).",
                 call. = FALSE)
        
        #-- type
        choices = c("combine", "separate")
        type = choices[pmatch(type, choices)]
        
        if (is.na(type)) 
            stop("'type' should be one of 'combine' or 'separate'.", 
                 call. = FALSE)
        
        #-- sideColors
        if (is.logical(sideColors)) {
            if(isTRUE(sideColors)) sideColors = c("blue", "red") else sideColors = NULL
        }
        else {
            if (length(sideColors) != 2)
                stop("'sideColors' must be a character vector of length ", 2, ".", 
                     call. = FALSE)
        }
        
        isColor = sapply(sideColors, function(x) { 
            tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })
        
        if (any(!isColor)) 
            stop("'sideColors' must be a character vector of recognized colors.", 
                 call. = FALSE)
        
        #-- interactive.dev
        if (!is.logical(interactive.dev))
            stop("'interactive.dev' must be a logical constant (TRUE or FALSE).",
                 call. = FALSE)
        
        #-- title
        if (!is.logical(title))
            stop("'title' must be a logical constant (TRUE or FALSE).",
                 call. = FALSE)
        
        #-- end checking --#
        #------------------#
        
        
        #-- calling the cim function for mapping -----------------------------------#
        #---------------------------------------------------------------------------#
        
        #-- if type = "combine" --#
        if (type == "combine") {		
            
            if (isTRUE(X.var.names)) X.var.names = colnames(X)
            else X.var.names = rep("", p)
            
            if (isTRUE(Y.var.names)) Y.var.names = colnames(Y)
            else Y.var.names = rep("", q)
            
            mat.lab = c(X.var.names, Y.var.names)
            
            matcor = cor(cbind(X, Y), use = "pairwise")
            colnames(matcor) = rownames(matcor) = mat.lab
            matcor=matcor[(p + q):1,]
            
            row.sideColors = col.sideColors = NULL
            
            if (!is.null(sideColors)) {
                bg.col = par("bg")
                row.sideColors = c(rep(sideColors[2], q),rep(sideColors[1], p))
                col.sideColors = rev(row.sideColors)
            }
            
            cim(matcor, cluster = "none", 
                title = if(title) "[X,Y] correlation matrix" else NULL,
                col.sideColors = col.sideColors,
                row.sideColors = row.sideColors,
                row.names = TRUE, col.names = TRUE)
        }
        
        # representation des matrices de correlation de #
        # X, Y et entre X et Y                          #
        #-----------------------------------------------#
        if (type == "separate") {
            Xcor = cor(X, use = "pairwise")
            Ycor = cor(Y, use = "pairwise")
            XYcor = cor(X, Y, use = "pairwise")
            
            Xcor = Xcor[p:1,]#t(Xcor[p:1, ])
            Ycor = Ycor[q:1,]#t(Ycor[q:1, ])
            XYcor = XYcor[p:1,]
            
            if (interactive.dev) {   
                cim(Xcor, cluster = "none",
                    title = if(title) "X correlation matrix" else NULL,
                    row.names = X.var.names, col.names = X.var.names)
                
                devAskNewPage(TRUE)
                cim(Ycor, cluster = "none", 
                    title = if(title) "Y correlation matrix" else NULL,
                    row.names = Y.var.names, col.names = Y.var.names)
                
                cim(XYcor, cluster = "none",
                    title = if(title) "XY correlation matrix" else NULL,			
                    row.names = X.var.names, col.names = Y.var.names)
            } else {
                GD = getOption("device")
                if (is.character(GD)) {
                    if (GD == "RStudioGD") GD = FALSE
                }
                else 
                    GD = TRUE
                
                if (isTRUE(GD)) getOption("device")()
                cim(Xcor, cluster = "none",
                    title = if(title) "X correlation matrix" else NULL,
                    row.names = X.var.names, col.names = X.var.names)
                
                if (isTRUE(GD)) getOption("device")()
                cim(Ycor, cluster = "none", 
                    title = if(title) "Y correlation matrix" else NULL,
                    row.names = Y.var.names, col.names = Y.var.names)
                
                if (isTRUE(GD)) getOption("device")()
                cim(XYcor, cluster = "none",
                    title = if(title) "XY correlation matrix" else NULL,
                    row.names = X.var.names, col.names = Y.var.names)
            }
        }
        return(invisible(NULL))
    }

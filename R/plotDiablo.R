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

# ========================================================================================================
# functions
#     1) plotIndiv_diablo; 2 sub-functions; splotMatPlot() and panel.ellipses
#     2) circosPlot_diablo
#     3) heatmap_diablo
#     4) enrichPathwayNetwork_diablo
#
# ========================================================================================================

#' Graphical output for the DIABLO framework
#' 
#' Function to visualise correlation between components from different data
#' sets
#' 
#' The function uses a plot.data.frame to plot the component \code{ncomp}
#' calculated from each data set to visualise whether DIABLO (block.splsda) is
#' successful at maximising the correlation between each data sets' component.
#' The lower triangular panel indicated the Pearson's correlation coefficient,
#' the upper triangular panel the scatter plot.
#'
#' @param object,x object of class inheriting from \code{"block.splsda"}. 
#' @param ncomp Which component to plot calculated from each data set. Has to
#' be lower than the minimum of \code{object$ncomp}.
#' @param col.per.group A named character of colours for each group class
#'   representation. Its names must match the levels of object$Y.
#' @param legend Logical. Whether the legend should be added. Default is TRUE.
#' @param legend.ncol Number of columns for the legend. Default to
#' \code{min(5,nlevels(x$Y))}.
#' @param \dots not used
#' @return none
#' @author Amrit Singh, Florian Rohart, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{block.splsda}} and http://www.mixOmics.org/mixDIABLO
#' for more details.
#' @references 
#' Singh A., Shannon C., Gautier B., Rohart F., Vacher M., Tebbutt S.
#' and Lê Cao K.A. (2019), DIABLO: an integrative approach for identifying key 
#' molecular drivers from multi-omics assays, Bioinformatics, 
#' Volume 35, Issue 17, 1 September 2019, Pages 3055–3062.
#' @keywords regression multivariate
#' @export
#' @examples
#' data('breast.TCGA')
#' Y = breast.TCGA$data.train$subtype
#' 
#' data = list(mrna =  breast.TCGA$data.train$mrna,
#' mirna =  breast.TCGA$data.train$mirna, prot =  breast.TCGA$data.train$protein)
#' 
#' # set number of component per data set
#' ncomp = 3
#' # set number of variables to select, per component and per data set (arbitrarily set)
#' list.keepX = list(mrna = rep(20, 3), mirna = rep(10,3), prot = rep(10,3))
#' 
#' # DIABLO using a full design where every block is connected
#' BC.diablo = block.splsda(X = data, Y = Y, ncomp = ncomp, keepX = list.keepX, design = 'full')
#' ## default col.per.group
#' plotDiablo(BC.diablo, ncomp = 1, legend = TRUE, col.per.group = NULL)
#' ## custom col.per.group
#' col.per.group <- color.mixo(1:3)
#' names(col.per.group) <- levels(Y)
#' plotDiablo(BC.diablo, ncomp = 1, legend = TRUE, col.per.group = col.per.group)
plotDiablo <- function(object,
                       ncomp = 1,
                       legend = TRUE,
                       legend.ncol,
                       col.per.group = NULL,
                       ...)
{
    if (!is.null(list(...)$x))
        .stop("use of 'x' has been deprecated. Use 'object' instead.")
    
    Y <- object$Y
    ## col.per.group
    col.per.group <- .change_if_null(col.per.group, default =  color.mixo(1:nlevels(Y)))
    col.per.group <- .get.cols.and.group(col.per.group = col.per.group, group = Y, object = object, n_ind = length(object$names$sample))
    col.per.group <- col.per.group$col.per.group
    
    #need to reorder variates and loadings to put 'Y' in last
    opar = par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
    
    indY=object$indY
    object$variates=c(object$variates[-indY],object$variates[indY])
    object$loadings=c(object$loadings[-indY],object$loadings[indY])
    
    VarX = do.call(cbind, lapply(object$variates, function(i) i[, ncomp]))
    datNames = colnames(VarX)
    
    if(ncol(VarX)<=2)
        stop("This function is only available when there are more than 3 blocks") # so 2 blocks + the outcome Y
    
    # check input parameters
    
    if (length(ncomp) != 1 | ncomp > min(object$ncomp))
        stop(paste0("'ncomp' must be a numeric value lower than ", min(object$ncomp),", which is min(object$ncomp)"))
    # end check parameters
    
    if(missing(legend.ncol))
        legend.ncol = min(5, nlevels(Y))
    
    numberOfCols = ncol(VarX)-1
    numberOfRows = numberOfCols #- 1
    
    mat = matrix(0, nrow = numberOfRows, ncol = numberOfRows)
    for(i in 1:nrow(mat))
    {
        for(j in 1:ncol(mat))
            mat[i,j] = paste(i,j, sep="_")
    }
    plotType = list(cor=mat[lower.tri(mat)], scatter=mat[upper.tri(mat)],
                    lab=diag(mat))#,
    #bar=paste(1:(numberOfRows-1), numberOfCols, sep="_"),
    #stackedbar=paste(numberOfRows, numberOfCols, sep="_"))
    
    par(mfrow = c(numberOfRows+1, numberOfCols), mar = rep.int(1/2, 4), oma = c(2,2,2,2))
    layout(matrix(c(1:(numberOfCols)^2, rep((numberOfCols)^2+1,numberOfCols)),numberOfRows+1,numberOfCols, byrow=TRUE),
           heights = c(rep(1,numberOfRows), 0.25 * floor(nlevels(Y)/legend.ncol)))
    for(i in 1:numberOfRows)
    {
        for(j in 1:numberOfCols)
        {
            ptype = unlist(lapply(plotType, function(x)
            {
                intersect(paste(i,j,sep="_"), x)
            }))
            splotMatPlot(x=VarX[, j], y=VarX[, i], datNames, Y, ptype, col.per.group)
            
            if(i == 1 & j %in% seq(2, numberOfRows, 1))
                Axis(side = 3, x=VarX[, i])
            
            if(j == numberOfRows & i %in% seq(1, numberOfRows-1, 1))
                Axis(side = 4, x=VarX[, i])
        }
    }
    #add legend
    plot(1:3,1:3,type="n",axes=FALSE,xlab="",ylab="")
    if(legend)
        legend("center",legend=levels(Y), col = col.per.group, pch = 19, ncol = legend.ncol, cex = 1.5)
    
    par(opar)
}

#' @rdname plotDiablo
#' @method plot sgccda
#' @export
plot.sgccda <- function(x, ...) plotDiablo(object = x, ...)

splotMatPlot = function(x, y, datNames, Y, ptype, col.per.group)
{
    if(names(ptype) == "cor")
    {
        plot(1, type = "n", axes = FALSE)
        r = round(cor(x, y), 2)
        text(1, 1, labels=r, cex = 0.6/strwidth(abs(r))*abs(r))
        box()
    }
    if(names(ptype) == "scatter")
        panel.ellipses(x=x, y=y, Y = Y, col.per.group = col.per.group)
    
    if(names(ptype) == "lab")
    {
        plot(1, type = "n", axes = FALSE)
        ind = as.numeric(unlist(lapply(strsplit(ptype, "_"), unique)))
        text(x=1, y=1, labels=datNames[ind], cex = 2)
        box()
    }
    # if(FALSE)
    # {
    #     if(names(ptype) == "bar")
    #     {
    #         Y2 = factor(as.character(Y), levels = groupOrder)
    #         boxplot(x ~ Y2, horizontal=TRUE, axes = FALSE, ylim = c(min(x)-3, max(x)),
    #                 col= col.per.group[groupOrder])
    #         axis(4, at=1:nlevels(Y2), labels=levels(Y2))
    #     }
    #     if(names(ptype) == "stackedbar")
    #     {
    #         Y2 = factor(as.character(Y), levels = groupOrder)
    #         bars = table(Y2)
    #         barplot(bars, col= color.mixo(match(levels(Y2), levels(Y))),
    #                 axes = FALSE)
    #         axis(4, at=seq(0,max(bars),length.out=5), labels=seq(0,max(bars),length.out=5))
    #     }
    # }
}

#' @importFrom ellipse ellipse
panel.ellipses = function(x, y, Y = Y, pch = par("pch"), col.lm = "red", axes = FALSE, col.per.group, ...)
{
    ind.gp = matrice = cdg = variance = list()
    for(i in 1:nlevels(Y))
        ind.gp[[i]] = which(as.numeric(Y)==i)
    
    matrice = lapply(ind.gp, function(z){matrix(c(x[z], y[z]), ncol = 2)})
    cdg = lapply(matrice, colMeans)
    variance = lapply(matrice, var)
    
    coord.ellipse = lapply(1:nlevels(Y), function(x){ellipse(variance[[x]], centre = cdg[[x]], ellipse.level = 0.95)})
    max.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, max)})
    min.ellipse = sapply(coord.ellipse, function(x){apply(x, 2, min)})
    ind.names = names(Y)
    cex = 0.5
    plot(x, y, xlab = "X.label", ylab = "Y.label", col = col.per.group[Y], pch=20, axes=axes,
         xlim = c(min(x, min.ellipse[1, ]), max(x, max.ellipse[1, ])), ylim = c(min(y, min.ellipse[2, ]), max(y, max.ellipse[2, ])))
    #text(x, y, ind.names, col = col, cex = cex)
    box()
    for (z in 1:nlevels(Y))
        points(coord.ellipse[[z]], type = "l", col = col.per.group[z])
    
}


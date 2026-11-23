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

# -------------------------------------------- #
# CIM for objects "pca","spca","ipca","sipca","mlsplsda","splsda","plsda","rcc",
#   "pls","spls","mlspls"
# -------------------------------------------- #
#' Clustered Image Maps (CIMs) ("heat maps")
#' 
#' This function generates color-coded Clustered Image Maps (CIMs) ("heat
#' maps") to represent "high-dimensional" data sets.
#' 
#' One matrix Clustered Image Map (default method) is a 2-dimensional
#' visualization of a real-valued matrix (basically
#' \code{\link{image}(t(mat))}) with rows and/or columns reordered according to
#' some hierarchical clustering method to identify interesting patterns.
#' Generated dendrograms from clustering are added to the left side and to the
#' top of the image. By default the used clustering method for rows and columns
#' is the \emph{complete linkage} method and the used distance measure is the
#' distance \emph{euclidean}.
#' 
#' In \code{"pca"}, \code{"spca"}, \code{"ipca"}, \code{"sipca"},
#' \code{"plsda"}, \code{"splsda"} and multilevel variants methods the
#' \code{mat} matrix is \code{object$X}.
#' 
#' For the remaining methods, if \code{mapping = "X"} or \code{mapping = "Y"}
#' the \code{mat} matrix is \code{object$X} or \code{object$Y} respectively. If
#' \code{mapping = "XY"}: \itemize{ \item in \code{rcc} method, the matrix
#' \code{mat} is created where element \eqn{(j,k)} is the scalar product value
#' between every pairs of vectors in dimension \code{length(comp)} representing
#' the variables \eqn{X_j} and \eqn{Y_k} on the axis defined by \eqn{Z_i} with
#' \eqn{i} in \code{comp}, where \eqn{Z_i} is the equiangular vector between
#' the \eqn{i}-th \eqn{X} and \eqn{Y} canonical variate.
#' 
#' \item in \code{pls}, \code{spls} and multilevel spls methods, if
#' \code{object$mode} is \code{"regression"}, the element \eqn{(j,k)} of the
#' matrix \code{mat} is given by the scalar product value between every pairs
#' of vectors in dimension \code{length(comp)} representing the variables
#' \eqn{X_j} and \eqn{Y_k} on the axis defined by \eqn{U_i} with \eqn{i} in
#' \code{comp}, where \eqn{U_i} is the \eqn{i}-th \eqn{X} variate. If
#' \code{object$mode} is \code{"canonical"} then \eqn{X_j} and \eqn{Y_k} are
#' represented on the axis defined by \eqn{U_i} and \eqn{V_i} respectively.}
#' 
#' The \code{blocks} parameter controls which blocks are to be included when 
#' \code{class(mat) == "block.pls" OR "block.spls"}. This can be a character or 
#' a integer vector. 
#' 
#' If using a multiblock object then \code{mapping} can be
#' set to \code{"multiblock"}. When done so, this will emulate the function of 
#' \code{cimDiablo()}, such that rows will denote each sample and all features
#' included in \code{blocks} will be shown as columns, coloured by which block
#' they inherit from. In this case, \code{blocks} can include any number of 
#' input blocks. If \code{mapping = "X", "Y" OR "XY"}, then it functions similarly
#' to if a \code{mixo_pls} object was being used. \code{blocks} has to be of length 2 
#' in this scenario.
#' 
#' By default four components will be displayed in the plot. At the top left is
#' the color key, top right is the column dendogram, bottom left is the row
#' dendogram, bottom right is the image plot. When \code{sideColors} are
#' provided, an additional row or column is inserted in the appropriate
#' location. This layout can be overriden by specifiying appropriate values for
#' \code{lwid} and \code{lhei}. \code{lwid} controls the column width, and
#' \code{lhei} controls the row height. See the help page for
#' \code{\link{layout}} for details on how to use these arguments.
#' 
#' For visualization of "high-dimensional" data sets, a nice zooming tool was
#' created. \code{zoom = TRUE} open a new device, one for CIM, one for zoom-out
#' region and define an interactive 'zoom' process: click two points at imagen
#' map region by pressing the first mouse button. It then draws a rectangle
#' around the selected region and zoom-out this at new device. The process can
#' be repeated to zoom-out other regions of interest.
#' 
#' The zoom process is terminated by clicking the second button and selecting
#' 'Stop' from the menu, or from the 'Stop' menu on the graphics window.
#' 
#' @param mat numeric matrix of values to be plotted. Alternatively, an object
#' of class inheriting from \code{"pca"}, \code{"spca"}, \code{"ipca"},
#' \code{"sipca"}, \code{"rcc"}, \code{"pls"}, \code{"spls"}, \code{"plsda"},
#' \code{"splsda"}, \code{"mlspls"}, \code{"mlsplsda"}, \code{"block.pls"} or
#' \code{"block.spls"} (where \code{"ml"} stands for multilevel).
#' @param color a character vector of colors such as that generated by
#' \code{\link{terrain.colors}}, \code{\link{topo.colors}},
#' \code{\link{rainbow}}, \code{\link{color.jet}} or similar functions.
#' @param row.names,col.names logical, should the name of rows and/or columns
#' of \code{mat} be shown? If \code{TRUE} (defaults) \code{rownames(mat)}
#' and/or \code{colnames(mat)} are used. Possible character vectors with row
#' and/or column labels can be used.
#' @param row.sideColors (optional) character vector of length \code{nrow(mat)}
#' containing the color names for a vertical side bar that may be used to
#' annotate the rows of \code{mat}.
#' @param col.sideColors (optional) character vector of length \code{ncol(mat)}
#' containing the color names for a horizontal side bar that may be used to
#' annotate the columns of \code{mat}.
#' @param row.cex,col.cex positive numbers, used as \code{cex.axis} in for the
#' row or column axis labeling. The defaults currently only use number of rows
#' or columns, respectively.
#' @param mapping character string indicating whether to map \code{"X"},
#' \code{"Y"} or \code{"XY"}-association matrix. Can also be \code{"multiblock"} 
#' when \code{class(mat) == "block.pls" OR "block.spls"}. See Details.
#' @param cluster character string indicating whether to cluster \code{"none"},
#' \code{"row"}, \code{"column"} or \code{"both"}. Defaults to \code{"both"}.
#' @param dist.method character vector of length two. The distance measure used
#' in clustering rows and columns. Possible values are \code{"correlation"} for
#' Pearson correlation and all the distances supported by \code{\link{dist}},
#' such as \code{"euclidean"}, etc.
#' @param clust.method character vector of length two. The agglomeration method
#' to be used for rows and columns. Accepts the same values as in
#' \code{\link{hclust}} such as \code{"ward"}, \code{"complete"}, etc.
#' @param cut.tree numeric vector of length two with components in [0,1]. The
#' height proportions where the trees should be cut for rows and columns, if
#' these are clustered.
#' @param comp atomic or vector of positive integers. The components to
#' adequately account for the data association. For a non sparse method, the
#' similarity matrix is computed based on the variates and loading vectors of
#' those specified components. For a sparse approach, the similarity matric is
#' computed based on the variables selected on those specified components. See
#' example. Defaults to \code{comp = 1:object$ncomp}.
#' @param blocks integer or character vector. Used when \code{class(mat) == 
#' "block.pls" OR "block.spls"}. Dictates which blocks will be visualised. See
#' Details.
#' @param transpose logical indicating if the matrix should be transposed for
#' plotting. Defaults to \code{FALSE}.
#' @param center either a logical value or a numeric vector of length equal to
#' the number of columns of \code{mat}. See \code{\link{scale}} function.
#' @param scale either a logical value or a numeric vector of length equal to
#' the number of columns of \code{mat}. See \code{\link{scale}} function.
#' @param cutoff numeric between 0 and 1. Variables with correlations below
#' this threshold in absolute value are not plotted. To use only when mapping
#' is "XY".
#' @param symkey Logical indicating whether the color key should be made
#' symmetric about 0. Defaults to \code{TRUE}.
#' @param keysize vector of length two, indicating the size of the color key.
#' @param keysize.label vector of length 1, indicating the size of the labels
#' and title of the color key.
#' @param zoom logical. Whether to use zoom for interactive zoom. See Details.
#' @param title,xlab,ylab title, \eqn{x}- and \eqn{y}-axis titles; default to
#' none.
#' @param margins numeric vector of length two containing the margins (see
#' \code{\link{par}(mar)}) for column and row names respectively.
#' @param lhei,lwid arguments passed to \code{layout} to divide the device up
#' into two (or three if a side color is drawn) rows and two columns, with the
#' row-heights \code{lhei} and the column-widths \code{lwid}.
#' @param legend A list indicating the legend for each group, the color vector,
#' title of the legend and cex.
#' @param save should the plot be saved? If so, argument to be set to either
#' \code{'jpeg'}, \code{'tiff'}, \code{'png'} or \code{'pdf'}.
#' @param name.save character string for the name of the file to be saved.
#' @return A list containing the following components: \item{M}{the mapped
#' matrix used by \code{cim}.} \item{rowInd, colInd}{row and column index
#' permutation vectors as returned by \code{\link{order.dendrogram}}.}
#' \item{ddr, ddc}{object of class \code{"dendrogram"} which describes the row
#' and column trees produced by \code{cim}.} \item{mat.cor}{the correlation
#' matrix used for the heatmap. Available only when mapping = "XY".}
#' \item{row.names, col.names}{character vectors with row and column labels
#' used.} \item{row.sideColors, col.sideColors}{character vector containing the
#' color names for vertical and horizontal side bars used to annotate the rows
#' and columns.}
#' @author Ignacio González, Francois Bartolo, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{heatmap}}, \code{\link{hclust}}, \code{\link{plotVar}},
#' \code{\link{network}} and
#' 
#' \url{https://mixomics.org/graphics/} for more details on all options
#' available.
#' @references Eisen, M. B., Spellman, P. T., Brown, P. O. and Botstein, D.
#' (1998). Cluster analysis and display of genome-wide expression patterns.
#' \emph{Proceeding of the National Academy of Sciences of the USA} \bold{95},
#' 14863-14868.
#' 
#' Weinstein, J. N., Myers, T. G., O'Connor, P. M., Friend, S. H., Fornace Jr.,
#' A. J., Kohn, K. W., Fojo, T., Bates, S. E., Rubinstein, L. V., Anderson, N.
#' L., Buolamwini, J. K., van Osdol, W. W., Monks, A. P., Scudiero, D. A.,
#' Sausville, E. A., Zaharevitz, D. W., Bunow, B., Viswanadhan, V. N., Johnson,
#' G. S., Wittes, R. E. and Paull, K. D. (1997). An information-intensive
#' approach to the molecular pharmacology of cancer. \emph{Science} \bold{275},
#' 343-349.
#' 
#' González I., Lê Cao K.A., Davis M.J., Déjean S. (2012). Visualising
#' associations between paired 'omics' data sets. \emph{BioData Mining};
#' \bold{5}(1).
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' @keywords multivariate iplot hplot graphs cluster
#' @example ./examples/cim-examples.R
#' @importFrom stats as.dendrogram as.dist
#' @export
cim <-
    function(mat = NULL,
             color = NULL,
             row.names = TRUE,
             col.names = TRUE,
             row.sideColors = NULL,
             col.sideColors = NULL,
             row.cex = NULL,
             col.cex = NULL,
             cutoff = 0,
             cluster = "both",
             dist.method = c("euclidean", "euclidean"),
             clust.method = c("complete", "complete"),
             cut.tree = c(0, 0),
             transpose = FALSE,
             symkey = TRUE,
             keysize = c(1, 1),
             keysize.label = 1,
             zoom = FALSE,
             title = NULL,
             xlab = NULL,
             ylab = NULL,
             margins = c(5, 5),
             lhei = NULL,
             lwid = NULL,
             comp = NULL,
             center = TRUE,
             scale = FALSE,
             mapping = "XY",
             legend = NULL,
             save = NULL,
             name.save = NULL,
             blocks = NULL)

    {
        class.object <- class(mat)
        
        object.single.omics = c("pca",
                                "spca",
                                "ipca",
                                "sipca",
                                "mixo_mlsplsda",
                                "mixo_splsda",
                                "mixo_plsda")
        object.rcc = c("rcc")
        object.pls = c("mixo_pls", 
                       "mixo_spls",
                       "mixo_mlspls")
        object.block.pls = c("block.pls",
                             "block.spls")
        object.all = c(object.single.omics,
                       object.rcc,
                       object.pls,
                       object.block.pls)
        
        #-- checking general input parameters -------------------------------------
        #--------------------------------------------------------------------------#
        
        #-- check that the user did not enter extra arguments
        arg.call = match.call()
        user.arg = names(arg.call)[-1]
        
        err = tryCatch(
            mget(names(formals()), sys.frame(sys.nframe())),
            error = function(e)
                e
        )
        
        if ("simpleError" %in% class(err))
            stop(err[[1]], ".", call. = FALSE)
        
        #-- color
        if (is.null(color))
            color = color.spectral(25)
        
        
        #-- cluster
        choices = c("both", "row", "column", "none")
        cluster = choices[pmatch(cluster, choices)]
        
        if (is.na(cluster))
            stop("'cluster' should be one of 'both', 'row', 'column' or 'none'.",
                 call. = FALSE)
        
        #-- cluster method
        if (!is.character(clust.method) |
            length(as.vector(clust.method)) != 2)
            stop("'clust.method' must be a character vector of length 2.",
                 call. = FALSE)
        
        choices = c("ward.D",
                    "single",
                    "complete",
                    "average",
                    "mcquitty",
                    "median",
                    "centroid")
        clust.method = choices[c(pmatch(clust.method[1], choices),
                                 pmatch(clust.method[2], choices))]
        
        if (any(is.na(clust.method)))
            stop("invalid clustering method.", call. = FALSE)
        
        #-- distance method
        if (!is.character(dist.method) |
            length(as.vector(dist.method)) != 2)
            stop("'dist.method' must be a character vector of length 2.", call. = FALSE)
        
        choices = c(
            "euclidean",
            "correlation",
            "maximum",
            "manhattan",
            "canberra",
            "binary",
            "minkowski"
        )
        dist.method = choices[c(pmatch(dist.method[1], choices),
                                pmatch(dist.method[2], choices))]
        
        if (any(is.na(dist.method)))
            stop("invalid distance method.", call. = FALSE)
        
        
        #-- checking general input arguments -------------------------------------
        #-------------------------------------------------------------------------#
        
        #-- color
        if (any(!sapply(color, function(color) {
            tryCatch(
                is.matrix(col2rgb(color)),
                error = function(e)
                    FALSE
            )
        })))
            stop("'color' must be a character vector of recognized colors.",
                 call. = FALSE)
        
        
        #-- row.sideColors
        if (any(!sapply(row.sideColors, function(row.sideColors) {
            tryCatch(
                is.matrix(col2rgb(row.sideColors)),
                error = function(e)
                    FALSE
            )
        })))
            stop(
                "color names for vertical side bar must be a character vector
    of recognized colors.",
                call. = FALSE
            )
        
        
        #-- col.sideColors
        if (any(!sapply(col.sideColors, function(col.sideColors) {
            tryCatch(
                is.matrix(col2rgb(col.sideColors)),
                error = function(e)
                    FALSE
            )
        })))
            stop(
                "color names for horizontal side bar must be a character vector
    of recognized colors.",
                call. = FALSE
            )
        
        
        #-- row.cex
        if (!is.null(row.cex)) {
            if (!is.numeric(row.cex) || length(row.cex) != 1)
                stop("'row.cex' must be a numerical value.", call. = FALSE)
        }
        
        #-- col.cex
        if (!is.null(col.cex)) {
            if (!is.numeric(col.cex) || length(col.cex) != 1)
                stop("'col.cex' must be a numerical value.", call. = FALSE)
        }
        
        #-- transpose
        if (!is.logical(transpose))
            stop("'transpose' must be a logical constant (TRUE or FALSE).",
                 call. = FALSE)
        
        #-- cut.tree
        if (!is.numeric(cut.tree) || length(cut.tree) != 2)
            stop("'cut.tree' must be a numeric vector of length 2.",
                 call. = FALSE)
        else {
            if (!(all(0 <= cut.tree & cut.tree <= 1)))
                stop("Components of 'cut.tree' must be between 0 and 1.",
                     call. = FALSE)
        }
        
        #-- keysize
        if (length(keysize) != 2 || any(!is.finite(keysize)))
            stop("'keysize' must be a numeric vector of length 2.",
                 call. = FALSE)
        
        #-- keysize.label
        if (length(keysize.label) != 1 || any(!is.finite(keysize)))
            stop("'keysize' must be a numeric vector of length 1.",
                 call. = FALSE)
        
        #-- zoom
        if (!is.logical(zoom))
            stop("'zoom' must be a logical constant (TRUE or FALSE).",
                 call. = FALSE)
        
        #-- margins
        if (!is.numeric(margins) || length(margins) != 2)
            stop("'margins' must be a numeric vector of length 2.",
                 call. = FALSE)
        
        #-- symkey
        if (!is.logical(symkey))
            stop("'symkey' must be a logical constant (TRUE or FALSE).",
                 call. = FALSE)
        
        #-- lhei
        if (!is.null(lhei)) {
            if (is.null(col.sideColors)) {
                if (length(lhei) != 2 | !is.numeric(lhei) | any(is.na(lhei)))
                    stop("'lhei' must be a numeric vector of length 2.",
                         call. = FALSE)
            }
            else {
                if (length(lhei) != 3 | !is.numeric(lhei) | any(is.na(lhei)))
                    stop("'lhei' must be a numeric vector of length 3.",
                         call. = FALSE)
            }
        }
        
        #-- lwid
        if (!is.null(lwid)) {
            if (is.null(row.sideColors)) {
                if (length(lwid) != 2 | !is.numeric(lwid) | any(is.na(lwid)))
                    stop("'lwid' must be a numeric vector of length 2.",
                         call. = FALSE)
            }
            else {
                if (length(lwid) != 3 | !is.numeric(lwid) | any(is.na(lwid)))
                    stop("'lwid' must be a numeric vector of length 3.",
                         call. = FALSE)
            }
        }
        
        #-- xlab
        xlab = as.graphicsAnnot(xlab)
        
        #-- ylab
        ylab = as.graphicsAnnot(ylab)
        
        #-- title
        title = as.graphicsAnnot(title)
        
        #-- threshold correlation
        if (!is.numeric(cutoff) | (cutoff > 1) | (cutoff < 0))
            stop("The value taken by 'cutoff' must be between 0 and 1",
                 call. = FALSE)
        
        #-- save
        if (!is.null(save)) {
            if (!save %in% c("jpeg", "tiff", "png", "pdf"))
                stop("'save' must be one of 'jpeg', 'png', 'tiff' or 'pdf'.",
                     call. = FALSE)
        }
        
        #-- name.save
        if (!is.null(name.save)) {
            if (!is.character(name.save) || length(name.save) > 1)
                stop("'name.save' must be a character.", call. = FALSE)
        } else {
            if (!is.null(save))
                name.save = paste0("cim_", gsub(".", "_", deparse(substitute(mat)),
                                                fixed = TRUE))
        }
        
        
        #-- mapping
        if (!any(class.object  %in%  object.single.omics)) {
          
          choices = c("XY", "X", "Y", "multiblock")
          mapping = choices[pmatch(mapping, choices)]
          
          if (is.na(mapping)) {
            stop("'mapping' should be one of 'XY', 'X', 'Y' or 'multiblock'.", call. = FALSE)
          }
          
          if (mapping == "multiblock") {
            if (!any(class.object %in% object.block.pls)) {
              stop("'mapping' can only equal 'multiblock' if 'mat' is of class 'block.pls' or 'block.spls'",
                   call. = FALSE)
            }
          }
        }
        
        #-- end checking --#
        #------------------#
        if (!is.null(save)) {
            while (dev.cur() > 2)
                dev.off()
            
            if (save == "jpeg")
                jpeg(
                    filename = paste0(name.save, ".jpeg"),
                    res = 600,
                    width = 4000,
                    height = 4000
                )
            if (save == "png")
                jpeg(
                    filename = paste0(name.save, ".png"),
                    res = 600,
                    width = 4000,
                    height = 4000
                )
            if (save == "tiff")
                tiff(
                    filename = paste0(name.save, ".tiff"),
                    res = 600,
                    width = 4000,
                    height = 4000
                )
            if (save == "pdf")
                pdf(file = paste0(name.save, ".pdf"))
            
        }
        
        if (any(class.object == "block.splsda"))
            stop("Please call the 'cimDiablo' function on your 'block.splsda' object",
                 call. = FALSE)
        
        
        if (!any(class.object %in% c(object.all, "matrix")))
            stop(
                "'mat' has to be a matrix or one of the following object: ",
                paste(object.all, collapse = ", "),
                ".",
                call. = FALSE
            )
        
        
        
        #-- checks for parameters when block object in use
        if (any(class.object %in% object.block.pls)) {
          
          if (mapping != "multiblock") {
            
            #-- blocks
            if (!is.null(blocks)) {
              if (!is.numeric(blocks) && !is.character(blocks)) {
                stop("'blocks' must be a numeric or character vector when `mapping != 'multiblock'",
                     call. = FALSE)
              }
              if (length(blocks) != 2) {
                stop("'blocks' must be a vector of length 2 when `mapping != 'multiblock'",
                     call. = FALSE)
              }
            } else {
              message(paste0("'blocks' defaulting to: '", mat$names$blocks[1], "' and '", mat$names$blocks[2], "'"))
              blocks <- mat$names$blocks[c(1,2)]
            }
            
            X.block.name <- blocks[[1]]
            Y.block.name <- blocks[[2]]
          } 
          else {
            if (is.null(blocks)) {
              blocks <- mat$names$blocks
            }
            
          }
          
          if (length(blocks) > length(mat$X)) {
            stop("'blocks' cannot be longer than 'mat$X'",
                 call. = FALSE)
          }
          
          if (is.numeric(blocks)) {
            if (any(blocks > length(mat$X))) {
              stop("All values of 'blocks' must be less than 'length(mat$X)'",
                   call. = FALSE)
            }
            blocks <- mat$names$blocks[blocks]
          }
          if (is.character(blocks)) {
            if (!all(blocks %in% mat$names$blocks)) {
              stop("All values of 'blocks' must be be found in 'mat$names$blocks'",
                   call. = FALSE)
            }
          }
          
          if (length(blocks) != length(unique(blocks))) {
            message("Adjusting 'blocks' to contain only unique values.")
            blocks <- unique(blocks)
          }
          
          rel.blocks.X <- mat$X[blocks]
          rel.blocks.variates <- mat$variates[blocks]
          rel.blocks.loadings <- mat$loadings[blocks]
        }
        
        #-- if mixOmics class
        if (any(class.object  %in%  object.all))
        {
            #-- general checks  -------------
            if (any(class.object %in% object.block.pls)) {
              if (length(unique(mat$ncomp)) != 1) {
                stop("'ncomp' across blocks need to be consistent.", call. = FALSE)
              }
              
              if (mapping != "multiblock") {
                ncomp = mat$ncomp[[X.block.name]]
                p = ncol(mat$X[[X.block.name]])
                q = ncol(mat$X[[Y.block.name]])
              } 
              else {
                # use first block as all blocks are checked to have equal ncomp
                ncomp = mat$ncomp[[1]]
                
                p = q = sum(sapply(mat$X[blocks], ncol))
              }
              
              if (length(unique(sapply(1:length(mat$X), FUN = function(x) { nrow(mat$X[[x]]) }))) != 1) {
                stop("number of rows of each block needs to be consistent.", call. = FALSE)
              }
              
              if (mapping != "multiblock") {
                n = nrow(mat$X[[X.block.name]])
              } else {
                # use first block as all blocks are checked to have equal nrow
                n = nrow(mat$X[[1]]) 
              }
              
            } else {
              ncomp = mat$ncomp
              
              p = ncol(mat$X)
              q = ncol(mat$Y)
              n = nrow(mat$X)
            }
            
            #-- comp
            if (is.null(comp)) {
                comp = 1:ncomp
            }
            if (length(comp) > 1) {
                comp = unique(comp)
                if (!is.numeric(comp) || any(comp < 1))
                    stop("invalid vector for 'comp'.", call. = FALSE)
                if (any(comp > ncomp))
                    stop("the elements of 'comp' must be <= ",
                         ncomp,
                         ".",
                         call. = FALSE)
            }
            
            if (length(comp) == 1) {
                if (is.null(comp) || !is.numeric(comp) || comp <= 0 || comp > ncomp)
                    stop("invalid value for 'comp'.", call. = FALSE)
                # comp <- c(comp, comp)
            }
            
            comp <- round(comp)
            
            # if object is a pls or spls with a univariate Y, or multivariate but
            # only 1 variable selected on all comp, we only plot a heatmap of X
            if (any(class.object %in%  object.pls))
            {
                temp = apply(mat$loadings$Y, 2,
                             function(x) {
                                 which(x != 0, arr.ind = TRUE)
                             })
                # gives which variables are selected
                num.variable.selected.Y = table(unlist(temp))
                
                if (length(num.variable.selected.Y) == 1)
                    #only one variable in Y to plot, will raise trouble so we
                    # switch from (s)pls to (s)plsda
                    class.object = "mixo_splsda"
            }
            #-- if c("mixo_pls","mixo_spls","mixo_mlspls") or pls with univarite Y ----
            ## or multivariate but only one Y kept in sparse model
            
            if (!any(class.object  %in%  object.single.omics)) {
              
                if (mapping == "multiblock") {
                  if (is.logical(row.names)) {
                    if (isTRUE(row.names)) {
                      row.names <- mat$names$sample
                    }
                  } else {
                    if (length(row.names) != n)
                      stop("'row.names' must be a character vector of length ",
                           n,
                           ".",
                           call. = FALSE)
                  }
                  
                  if (is.logical(col.names)) {
                    if (isTRUE(col.names)) {
                      col.names <- unname(unlist(mat$names$colnames[blocks]))
                    }
                  } else {
                    if (length(col.names) != p)
                      stop("'col.names' must be a character vector of length ",
                           p,
                           ".",
                           call. = FALSE)
                  }
                }
              
                if (!is.null(row.sideColors))
                {
                  row.sideColors = as.matrix(row.sideColors)
                  if (nrow(row.sideColors) != n)
                    stop(
                      "'row.sideColors' must be a colors character vector
                      (matrix) of length (nrow) ",
                      n,
                      ".",
                      call. = FALSE
                    )
                }
                if (!is.null(col.sideColors))
                {
                  col.sideColors = as.matrix(col.sideColors)
                  if (nrow(col.sideColors) != p)
                    stop(
                      "'col.sideColors' must be a colors character vector
                      (matrix) of length (nrow) ",
                      p,
                      ".",
                      call. = FALSE
                    )
                }
                
                if (mapping == "XY")
                {
                    if (is.logical(row.names)) {
                        if (isTRUE(row.names))
                            if (any(class.object %in% object.block.pls)) {
                              row.names = mat$names$colnames[[X.block.name]]
                            } else {
                              row.names = mat$names$colnames$X
                            }
                        else
                            row.names = rep("", p)
                    } else {
                        if (length(row.names) != p)
                            stop("'row.names' must be a character vector of length ",
                                 p,
                                 ".",
                                 call. = FALSE)
                    }
                    
                    if (is.logical(col.names))
                    {
                        if (isTRUE(col.names))
                            if (any(class.object %in% object.block.pls)) {
                              col.names = mat$names$colnames[[Y.block.name]]
                            } else {
                              col.names = mat$names$colnames$Y
                            }
                            
                        else
                            col.names = rep("", q)
                    } else {
                        if (length(col.names) != q)
                            stop("'col.names' must be a character vector of length ",
                                 q,
                                 ".",
                                 call. = FALSE)
                    }
                    
                    if (!is.null(row.sideColors))
                    {
                        row.sideColors = as.matrix(row.sideColors)
                        if (nrow(row.sideColors) != p)
                            stop(
                                "'row.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ",
                                p,
                                ".",
                                call. = FALSE
                            )
                    }
                    if (!is.null(col.sideColors))
                    {
                        col.sideColors = as.matrix(col.sideColors)
                        if (nrow(col.sideColors) != q)
                            stop(
                                "'col.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ",
                                q,
                                ".",
                                call. = FALSE
                            )
                    }
                }
                
                
                if (mapping == "X")
                {
                    if (is.logical(row.names))
                    {
                        if (isTRUE(row.names)) {
                            row.names = mat$names$sample
                        }
                        else
                            row.names = rep("", n)
                    } else {
                        if (length(row.names) != n)
                            stop("'row.names' must be a character vector of length ",
                                 n,
                                 ".",
                                 call. = FALSE)
                    }
                    if (is.logical(col.names))
                    {
                        if (isTRUE(col.names))
                          if (any(class.object %in% object.block.pls)) {
                            col.names = mat$names$colnames[[X.block.name]]
                          } else {
                            col.names = mat$names$colnames$Y
                          }
                        else
                            col.names = rep("", p)
                    } else {
                        if (length(col.names) != p)
                            stop("'col.names' must be a character vector of length ",
                                 p,
                                 ".",
                                 call. = FALSE)
                    }
                    
                    if (!is.null(row.sideColors))
                    {
                        row.sideColors = as.matrix(row.sideColors)
                        if (nrow(row.sideColors) != n)
                            stop(
                                "'row.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ",
                                n,
                                ".",
                                call. = FALSE
                            )
                    }
                    if (!is.null(col.sideColors))
                    {
                        col.sideColors = as.matrix(col.sideColors)
                        if (nrow(col.sideColors) != p)
                            stop(
                                "'col.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ",
                                p,
                                ".",
                                call. = FALSE
                            )
                    }
                }
                
                if (mapping == "Y")
                {
                    if (is.logical(row.names))
                    {
                        if (isTRUE(row.names))
                        {
                            if (any(class.object %in% object.rcc))
                                row.names = mat$names$sample
                            else
                                row.names = mat$names$sample
                        }
                        else
                            row.names = rep("", n)
                    } else {
                        if (length(row.names) != n)
                            stop("'row.names' must be a character vector of length ",
                                 n,
                                 ".",
                                 call. = FALSE)
                    }
                    if (is.logical(col.names))
                    {
                        if (isTRUE(col.names))
                          if (any(class.object %in% object.block.pls)) {
                            col.names = mat$names$colnames[[Y.block.name]]
                          } else {
                            col.names = mat$names$colnames$Y
                          }
                        else
                            col.names = rep("", q)
                    } else {
                        if (length(col.names) != q)
                            stop("'col.names' must be a character vector of length ",
                                 q,
                                 ".",
                                 call. = FALSE)
                    }
                    
                    if (!is.null(row.sideColors))
                    {
                        row.sideColors = as.matrix(row.sideColors)
                        if (nrow(row.sideColors) != n)
                            stop(
                                "'row.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ",
                                n,
                                ".",
                                call. = FALSE
                            )
                    }
                    if (!is.null(col.sideColors))
                    {
                        col.sideColors = as.matrix(col.sideColors)
                        if (nrow(col.sideColors) != q)
                            stop(
                                "'col.sideColors' must be a colors character vector
                    (matrix) of length (nrow) ",
                                q,
                                ".",
                                call. = FALSE
                            )
                    }
                  
                  
                }
            }
            
            #-- if NOT c("mixo_pls","mixo_spls","mixo_mlspls") or pls with univarite Y ----
            if (any(class.object %in%  object.single.omics))
            {
                #-- row.sideColors
                if (!is.null(row.sideColors))
                {
                    row.sideColors = as.matrix(row.sideColors)
                    if (nrow(row.sideColors) != n)
                        stop(
                            "'row.sideColors' must be a colors character vector
                (matrix) of length (nrow) ",
                            n,
                            ".",
                            call. = FALSE
                        )
                }
                
                #-- col.sideColors
                if (!is.null(col.sideColors))
                {
                    col.sideColors = as.matrix(col.sideColors)
                    if (nrow(col.sideColors) != p)
                        stop(
                            "'col.sideColors' must be a colors character vector
                (matrix) of length (nrow) ",
                            p,
                            ".",
                            call. = FALSE
                        )
                }
                sample.sideColors = row.sideColors
                
                ## ----- DA object
                if (any(class.object %in%  c("mixo_splsda", "mixo_plsda",
                                             'mixo_mlsplsda')))
                {
                    #-- row.names
                    if (is.logical(row.names))
                    {
                        if (isTRUE(row.names))
                            row.names = mat$names$sample
                    }
                    #-- col.names
                    if (is.logical(col.names))
                    {
                        if (isTRUE(col.names))
                            col.names = mat$names$colnames$X
                    }
                    if (any(class.object %in%  c("mixo_splsda", 'mixo_mlsplsda')))
                        keep.X = apply(abs(mat$loadings$X[, comp, drop = FALSE]), 1,
                                       sum) > 0
                    else
                        keep.X = apply(abs(mat$loadings$X), 1, sum) > 0
                    cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp],
                                 use = "pairwise")
                    X.mat = as.matrix(mat$variates$X[, comp])
                }
                else{
                    #-- row.names
                    if (is.logical(row.names))
                    {
                        if (isTRUE(row.names))
                            row.names = mat$names$sample
                    }
                    #-- col.names
                    if (is.logical(col.names))
                    {
                        if (isTRUE(col.names))
                            col.names = mat$names$X
                    }
                    if (any(class.object %in%  c("spca", "sipca")))
                        keep.X = apply(abs(mat$rotation[, comp]), 1, sum) > 0
                    else
                        keep.X = apply(abs(mat$rotation), 1, sum) > 0
                    cord.X = cor(mat$X[, keep.X], mat$x[, comp], use = "pairwise")
                    X.mat = as.matrix(mat$x[, comp])
                }
                
                
                #-- cheking center and scale
                if (!is.logical(center)) {
                    if (!is.numeric(center) || (length(center) != ncol(mat$X)))
                        stop(
                            "'center' should be either a logical value or a numeric
                vector of length equal to the number of columns of 'X'.",
                            call. = FALSE
                        )
                }
                if (!is.logical(scale)) {
                    if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
                        stop(
                            "'scale' should be either a logical value or a numeric
                vector of length equal to the number of columns of 'X'.",
                            call. = FALSE
                        )
                }
                ## ---- clustering ----------------------------------------------
                object = scale(mat$X[, keep.X], center = center, scale = scale)
                col.names = col.names[keep.X]
                
                if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[keep.X,])
                
                if ((cluster == "both") || (cluster == "row")) {
                    Rowv = rowMeans(X.mat)
                    
                    if (dist.method[1] == "correlation")
                        dist.mat = as.dist(1 - cor(t(as.matrix(
                            X.mat
                        )),
                        method = "pearson"))
                    else
                        dist.mat = dist(X.mat, method = dist.method[1])
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd,]
                    row.names = row.names[rowInd]
                    
                    if (!is.null(row.sideColors))
                        row.sideColors = as.matrix(row.sideColors[rowInd,])
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    Colv = rowMeans(cord.X)
                    
                    if (dist.method[2] == "correlation")
                        dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
                    else
                        dist.mat = dist(cord.X, method = dist.method[2])
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    if (!is.null(col.sideColors))
                        col.sideColors = as.matrix(col.sideColors[colInd,])
                }
                
                #-- calling the image.map function --------------------------------#
                #------------------------------------------------------------------#
                
                
                
                #-- output --------------------------------------------------------#
                #------------------------------------------------------------------#
                res = list(mat = object,
                           row.names = row.names,
                           col.names = col.names)
                
                if ((cluster == "both") || (cluster == "row")) {
                    res$rowInd = rowInd
                    res$ddr = ddr
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    res$colInd = colInd
                    res$ddc = ddc
                }
                
                class(res) = paste("cim", class.object[1], sep = "_")
                
            }
            else if (any(class.object %in%  object.rcc))
            {
                bisect = mat$variates$X[, comp] + mat$variates$Y[, comp]
                cord.X = cor(mat$X, bisect, use = "pairwise")
                cord.Y = cor(mat$Y, bisect, use = "pairwise")
                XY.mat = as.matrix(cord.X %*% t(cord.Y))
                
                #-- if mapping = "XY"
                if (mapping == "XY") {
                    object = XY.mat
                    
                    cut = list()
                    if (cutoff != 0) {
                        cut[[1]] = unlist(lapply(1:nrow(object),
                                                 function(x) {
                                                     any(abs(object[x, ]) > cutoff)
                                                 }))
                        object = object[cut[[1]], ]
                        if (dist.method[1] != "correlation")
                            cord.X = cord.X[cut[[1]], ]
                        
                        
                        if (is.null(nrow(object)) || nrow(object) == 0)
                            stop("cutoff value very high. No variable was selected.",
                                 call. = FALSE)
                        
                        
                        cut[[2]] = unlist(lapply(1:ncol(object),
                                                 function(x) {
                                                     any(abs(object[, x]) > cutoff)
                                                 }))
                        object = object[, cut[[2]]]
                        if (dist.method[2] != "correlation")
                            cord.Y = cord.Y[cut[[2]], ]
                        
                        
                        if (is.null(ncol(object)) || ncol(object) == 0)
                            stop("cutoff value very high. No variable was selected.",
                                 call. = FALSE)
                        
                    }
                    
                    if ((cluster == "both") || (cluster == "row")) {
                        #Rowv = rowMeans(XY.mat)
                        Rowv = rowMeans(cord.X)
                        
                        if (dist.method[1] == "correlation")
                            dist.mat = as.dist(1 - cor(t(as.matrix(
                                object
                            )),
                            method = "pearson"))
                        
                        else
                            dist.mat = dist(cord.X, method = dist.method[1])
                        
                        if (cutoff > 0) {
                            row.names = row.names[cut[[1]]]
                            if (!is.null(row.sideColors))
                                row.sideColors = as.matrix(row.sideColors[cut[[1]],])
                        }
                        
                        hcr = hclust(dist.mat, method = clust.method[1])
                        ddr = as.dendrogram(hcr)
                        ddr = reorder(ddr, Rowv)
                        rowInd = order.dendrogram(ddr)
                        object = object[rowInd,]
                        row.names = row.names[rowInd]
                        
                        if (!is.null(row.sideColors))
                            row.sideColors = as.matrix(row.sideColors[rowInd,])
                    }
                    
                    if ((cluster == "both") || (cluster == "column")) {
                        Colv = rowMeans(cord.Y)
                        
                        if (dist.method[2] == "correlation")
                            dist.mat = as.dist(1 - cor(as.matrix(object),
                                                       method = "pearson"))
                        else
                            dist.mat = dist(cord.Y, method = dist.method[2])
                        
                        if (cutoff > 0) {
                            col.names = col.names[cut[[2]]]
                            if (!is.null(col.sideColors))
                                col.sideColors = as.matrix(col.sideColors[cut[[2]],])
                        }
                        
                        hcc = hclust(dist.mat, method = clust.method[2])
                        ddc = as.dendrogram(hcc)
                        ddc = reorder(ddc, Colv)
                        colInd = order.dendrogram(ddc)
                        object = object[, colInd]
                        col.names = col.names[colInd]
                        
                        if (!is.null(col.sideColors))
                            col.sideColors = as.matrix(col.sideColors[colInd,])
                    }
                }
                
                #-- if mapping = "X"
                if (mapping == "X") {
                    #-- cheking center and scale
                    if (!is.logical(center)) {
                        if (!is.numeric(center) || (length(center) != ncol(mat$X)))
                            stop(
                                "'center' should be either a logical value or a
                    numeric vector of length equal to the number of columns of
                    'X'.",
                                call. = FALSE
                            )
                    }
                    if (!is.logical(scale)) {
                        if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
                            stop(
                                "'scale' should be either a logical value or a
                    numeric vector of length equal to the number of columns of
                    'X'.",
                                call. = FALSE
                            )
                    }
                    
                    object = scale(mat$X, center = center, scale = scale)
                    X.mat = as.matrix(mat$variates$X[, comp])
                    
                    if ((cluster == "both") || (cluster == "row")) {
                        Rowv = rowMeans(X.mat)
                        
                        if (dist.method[1] == "correlation")
                            dist.mat = as.dist(1 - cor(t(as.matrix(
                                X.mat
                            )),
                            method = "pearson"))
                        else
                            dist.mat = dist(X.mat, method = dist.method[1])
                        
                        hcr = hclust(dist.mat, method = clust.method[1])
                        ddr = as.dendrogram(hcr)
                        ddr = reorder(ddr, Rowv)
                        rowInd = order.dendrogram(ddr)
                        object = object[rowInd,]
                        row.names = row.names[rowInd]
                        
                        
                        if (!is.null(row.sideColors))
                            row.sideColors = as.matrix(row.sideColors[rowInd,])
                    }
                    
                    if ((cluster == "both") || (cluster == "column")) {
                        Colv = rowMeans(cord.X)
                        
                        if (dist.method[2] == "correlation")
                            dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
                        else
                            dist.mat = dist(cord.X, method = dist.method[2])
                        
                        hcc = hclust(dist.mat, method = clust.method[2])
                        ddc = as.dendrogram(hcc)
                        ddc = reorder(ddc, Colv)
                        colInd = order.dendrogram(ddc)
                        object = object[, colInd]
                        col.names = col.names[colInd]
                        
                        if (!is.null(col.sideColors))
                            col.sideColors = as.matrix(col.sideColors[colInd,])
                        
                    }
                    
                }
                
                #-- if mapping = "Y"
                if (mapping == "Y") {
                    #-- cheking center and scale
                    if (!is.logical(center)) {
                        if (!is.numeric(center) || (length(center) != ncol(mat$Y)))
                            stop(
                                "'center' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                                call. = FALSE
                            )
                    }
                    if (!is.logical(scale)) {
                        if (!is.numeric(scale) || (length(scale) != ncol(mat$Y)))
                            stop(
                                "'scale' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                                call. = FALSE
                            )
                    }
                    
                    object = scale(mat$Y, center = center, scale = scale)
                    Y.mat = as.matrix(mat$variates$Y[, comp])
                    
                    
                    if ((cluster == "both") || (cluster == "row")) {
                        Rowv = rowMeans(Y.mat)
                        
                        if (dist.method[1] == "correlation")
                            dist.mat = as.dist(1 - cor(t(as.matrix(
                                Y.mat
                            )),
                            method = "pearson"))
                        else
                            dist.mat = dist(Y.mat, method = dist.method[1])
                        
                        hcr = hclust(dist.mat, method = clust.method[1])
                        ddr = as.dendrogram(hcr)
                        ddr = reorder(ddr, Rowv)
                        rowInd = order.dendrogram(ddr)
                        object = object[rowInd,]
                        row.names = row.names[rowInd]
                        
                        if (!is.null(row.sideColors))
                            row.sideColors = as.matrix(row.sideColors[rowInd,])
                    }
                    
                    if ((cluster == "both") || (cluster == "column")) {
                        Colv = rowMeans(cord.Y)
                        
                        if (dist.method[2] == "correlation")
                            dist.mat = as.dist(1 - cor(t(cord.Y), method = "pearson"))
                        else
                            dist.mat = dist(cord.Y, method = dist.method[2])
                        
                        hcc = hclust(dist.mat, method = clust.method[2])
                        ddc = as.dendrogram(hcc)
                        ddc = reorder(ddc, Colv)
                        colInd = order.dendrogram(ddc)
                        object = object[, colInd]
                        col.names = col.names[colInd]
                        
                        if (!is.null(col.sideColors))
                            col.sideColors = as.matrix(col.sideColors[colInd,])
                        
                    }
                }
                
                #-- calling the image.map function --------------------------------#
                #------------------------------------------------------------------#
                
                
                
                #-- output --------------------------------------------------------#
                #------------------------------------------------------------------#
                res = list(mat = object,
                           row.names = row.names,
                           col.names = col.names)
                
                if ((cluster == "both") || (cluster == "row")) {
                    res$rowInd = rowInd
                    res$ddr = ddr
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    res$colInd = colInd
                    res$ddc = ddc
                }
                
                class(res) = "cim_rcc"
                
            }
            else if (any(class.object %in%  c(object.pls, object.block.pls)))
            {
                if (any(class.object %in% object.block.pls)) {
                  
                  if (any(class.object %in% c("block.spls")))
                  {
                    if (mapping == "multiblock") {
                      
                      num.keep.X <- sapply(rel.blocks.loadings, function(x) {
                        length(which(apply(abs(x), 1, sum) > 0))
                      })
                      
                      keep.X <- apply(abs(do.call("rbind", rel.blocks.loadings)), 1, sum) > 0
                    }
                    else {
                      keep.X <- apply(abs(mat$loadings[[X.block.name]][, comp, drop = FALSE]), 1, sum) > 0
                      keep.Y <- apply(abs(mat$loadings[[Y.block.name]][, comp, drop = FALSE]), 1, sum) > 0
                    }
                  } 
                  else {
                    if (mapping == "multiblock") {
                      num.keep.X <- sapply(rel.blocks.loadings, function(x) {
                        length(which(apply(abs(x), 1, sum) > 0))
                      })
                      keep.X <- apply(abs(do.call("rbind", rel.blocks.loadings)), 1, sum) > 0
                    } else {
                      keep.X <- apply(abs(mat$loadings[[X.block.name]]), 1, sum) > 0
                      keep.Y <- apply(abs(mat$loadings[[Y.block.name]]), 1, sum) > 0
                      
                      num.keep.X <- unname(sapply(mat$X, ncol))
                    }
                  }
                  
                  if (mapping  == "multiblock") {
                    
                    if (mat$mode == "canonical") {
                      bisect = rel.blocks.variates[[1]][, comp]
                      for (b in 2:length(rel.blocks.variates)) { bisect <- bisect + rel.blocks.variates[[b]][, comp] }
                      
                      cord.X = cor(do.call("cbind", rel.blocks.X)[, keep.X, drop = FALSE],
                                   bisect, use = "pairwise")
                    }
                    else {
                      cord.X = cor(do.call("cbind", rel.blocks.X)[, keep.X, drop = FALSE],
                                   do.call("cbind", rel.blocks.variates)[, comp], use = "pairwise")
                    }
                  } else {
                    if (mat$mode == "canonical") {
                      bisect = mat$variates[[X.block.name]][, comp] + mat$variates[[Y.block.name]][, comp]
                      cord.X = cor(mat$X[[X.block.name]][, keep.X, drop = FALSE],
                                   bisect, use = "pairwise")
                      cord.Y = cor(mat$X[[Y.block.name]][, keep.Y, drop = FALSE],
                                   bisect, use = "pairwise")
                    }
                    else {
                      cord.X = cor(mat$X[[X.block.name]][, keep.X, drop = FALSE],
                                   mat$variates[[X.block.name]][, comp], use = "pairwise")
                      cord.Y = cor(mat$X[[Y.block.name]][, keep.Y, drop = FALSE],
                                   mat$variates[[Y.block.name]][, comp], use = "pairwise")
                    }
                  }
                  
                  
                  
                  
                  
                } else {
                  if (any(class.object %in% c("mixo_spls", "mixo_mlspls")))
                  {
                      keep.X = apply(abs(mat$loadings$X[, comp, drop = FALSE]), 1,
                                     sum) > 0
                      keep.Y = apply(abs(mat$loadings$Y[, comp, drop = FALSE]), 1,
                                     sum) > 0
                  }
                  else
                  {
                      keep.X = apply(abs(mat$loadings$X), 1, sum) > 0
                      keep.Y = apply(abs(mat$loadings$Y), 1, sum) > 0
                  }
                  
                  
                  
                  if (mat$mode == "canonical") {
                      bisect = mat$variates$X[, comp] + mat$variates$Y[, comp]
                      cord.X = cor(mat$X[, keep.X, drop = FALSE],
                                   bisect, use = "pairwise")
                      cord.Y = cor(mat$Y[, keep.Y, drop = FALSE],
                                   bisect, use = "pairwise")
                  }
                  else {
                      cord.X = cor(mat$X[, keep.X, drop = FALSE],
                                   mat$variates$X[, comp], use = "pairwise")
                      cord.Y = cor(mat$Y[, keep.Y, drop = FALSE],
                                   mat$variates$X[, comp], use = "pairwise")
                  }
                }
                
                if (mapping != "multiblock") { XY.mat = as.matrix(cord.X %*% t(cord.Y)) }
                sample.sideColors = row.sideColors
                
                
                #-- if mapping = "multiblock"
                if (mapping == "multiblock") {
                  #-- cheking center and scale
                  if (!is.logical(center)) {
                    if (!is.numeric(center) || (length(center) != p))
                      stop(
                        "'center' should be either a logical value or a numeric
                    vector of length equal to the sum number of columns across all blocks.",
                    call. = FALSE
                      )
                  }
                  if (!is.logical(scale)) {
                    if (!is.numeric(scale) || (length(scale) != p))
                      stop(
                        "'scale' should be either a logical value or a numeric
                    vector of length equal to the sum number of columns across all blocks",
                    call. = FALSE
                      )
                  }
                  
                  
                  object = scale(do.call("cbind", rel.blocks.X), center = center, scale = scale)
                  object = object[, keep.X]
                  X.mat = as.matrix(do.call("cbind", rel.blocks.variates))
                  col.names = col.names[keep.X]
                  
                  if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(matrix(col.sideColors)[keep.X,])
                  else {
                    idx <- unlist(c(sapply(1:length(rel.blocks.X), function(x) { rep(x, num.keep.X[[x]]) })))
                    col.sideColors = brewer.pal(n = 12, name = 'Paired')[seq(2, 12, by = 2)]
                    col.sideColors = as.matrix(col.sideColors[idx])
                  }
                  
                  if ((cluster == "both") || (cluster == "row")) {
                    Rowv = rowMeans(X.mat)
                    
                    if (dist.method[1] == "correlation")
                      dist.mat = as.dist(1 - cor(t(as.matrix(X.mat)),
                      method = "pearson"))
                    else
                      dist.mat = dist(X.mat, method = dist.method[1])
                    
                    hcr = hclust(dist.mat, method = clust.method[1])
                    ddr = as.dendrogram(hcr)
                    ddr = reorder(ddr, Rowv)
                    rowInd = order.dendrogram(ddr)
                    object = object[rowInd,]
                    row.names = row.names[rowInd]
                    
                    if (!is.null(row.sideColors))
                      row.sideColors = as.matrix(matrix(row.sideColors)[rowInd,])
                  }
                  
                  if ((cluster == "both") || (cluster == "column")) {
                    Colv = rowMeans(cord.X)
                    
                    if (dist.method[2] == "correlation")
                      dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
                    else
                      dist.mat = dist(cord.X, method = dist.method[2])
                    
                    hcc = hclust(dist.mat, method = clust.method[2])
                    ddc = as.dendrogram(hcc)
                    ddc = reorder(ddc, Colv)
                    colInd = order.dendrogram(ddc)
                    object = object[, colInd]
                    col.names = col.names[colInd]
                    
                    col.sideColors = as.matrix(col.sideColors[colInd,])
                  }
                }
                
                #-- if mapping = "XY"
                if (mapping == "XY") {
                    object = XY.mat
                    
                    row.names = row.names[keep.X]
                    col.names = col.names[keep.Y]
                    
                    cut = list()
                    if (cutoff != 0) {
                        cut[[1]] = unlist(lapply(1:nrow(object),
                                                 function(x) {
                                                     any(abs(object[x, ]) > cutoff)
                                                 }))
                        object = object[cut[[1]], ]
                        if (dist.method[1] != "correlation")
                            cord.X = cord.X[cut[[1]], ]
                        
                        
                        if (is.null(nrow(object)) || nrow(object) == 0)
                            stop("cutoff value very high. No variable was selected.",
                                 call. = FALSE)
                        
                        
                        cut[[2]] = unlist(lapply(1:ncol(object),
                                                 function(x) {
                                                     any(abs(object[, x]) > cutoff)
                                                 }))
                        object = object[, cut[[2]]]
                        if (dist.method[2] != "correlation")
                            cord.Y = cord.Y[cut[[2]], ]
                        
                        
                        if (is.null(ncol(object)) || ncol(object) == 0)
                            stop("cutoff value very high. No variable was selected.",
                                 call. = FALSE)
                        
                    }
                    
                    
                    if (!is.null(row.sideColors))
                        row.sideColors = as.matrix(row.sideColors[keep.X,])
                    if (!is.null(col.sideColors))
                        col.sideColors = as.matrix(col.sideColors[keep.Y,])
                    
                    if ((cluster == "both") || (cluster == "row")) {
                        #Rowv = rowMeans(XY.mat)
                        Rowv = rowMeans(cord.X)
                        
                        if (dist.method[1] == "correlation")
                            dist.mat = as.dist(1 - cor(t(as.matrix(
                                object
                            )),
                            method = "pearson"))
                        else
                            #dist.mat = dist(mat, method = dist.method[1])
                            dist.mat = dist(cord.X, method = dist.method[1])
                        
                        if (cutoff > 0) {
                            row.names = row.names[cut[[1]]]
                            if (!is.null(row.sideColors))
                                row.sideColors = as.matrix(row.sideColors[cut[[1]],])
                        }
                        
                        hcr = hclust(dist.mat, method = clust.method[1])
                        ddr = as.dendrogram(hcr)
                        ddr = reorder(ddr, Rowv)
                        rowInd = order.dendrogram(ddr)
                        object = object[rowInd,]
                        row.names = row.names[rowInd]
                        
                        if (!is.null(row.sideColors))
                            row.sideColors = as.matrix(row.sideColors[rowInd,])
                    }
                    
                    if ((cluster == "both") || (cluster == "column")) {
                        #Colv = colMeans(mat)
                        Colv = rowMeans(cord.Y)
                        
                        if (dist.method[2] == "correlation")
                            dist.mat = as.dist(1 - cor(as.matrix(object),
                                                       method = "pearson"))
                        else
                            #dist.mat = dist(t(mat), method = dist.method[2])
                            dist.mat = dist(cord.Y, method = dist.method[2])
                        
                        if (cutoff > 0) {
                            col.names = col.names[cut[[2]]]
                            if (!is.null(col.sideColors))
                                col.sideColors = as.matrix(col.sideColors[cut[[2]],])
                        }
                        
                        hcc = hclust(dist.mat, method = clust.method[2])
                        ddc = as.dendrogram(hcc)
                        ddc = reorder(ddc, Colv)
                        colInd = order.dendrogram(ddc)
                        object = object[, colInd]
                        col.names = col.names[colInd]
                        
                        if (!is.null(col.sideColors))
                            col.sideColors = as.matrix(col.sideColors[colInd,])
                    }
                }
                
                #-- if mapping = "X"
                if (mapping == "X") {
                    #-- cheking center and scale
                    if (!is.logical(center)) {
                        if (!is.numeric(center) || (length(center) != ncol(mat$X)))
                            stop(
                                "'center' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'X'.",
                                call. = FALSE
                            )
                    }
                    if (!is.logical(scale)) {
                        if (!is.numeric(scale) || (length(scale) != ncol(mat$X)))
                            stop(
                                "'scale' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'X'.",
                                call. = FALSE
                            )
                    }
                    
                    if (any(class.object %in% object.block.pls)) {
                      tmp.df <- mat$X[[X.block.name]][, keep.X]
                      tmp.var <- mat$variates[[X.block.name]][, comp]
                    } else {
                      tmp.df <- mat$X[, keep.X]
                      tmp.var <- mat$variates$X[, comp]
                    }
                  
                    object = scale(tmp.df, center = center, scale = scale)
                    X.mat = as.matrix(tmp.var)
                    col.names = col.names[keep.X]
                    
                    rm(tmp.df)
                    rm(tmp.var)
                    
                    if (!is.null(col.sideColors))
                        col.sideColors = as.matrix(col.sideColors[keep.X,])
                    
                    
                    
                    if ((cluster == "both") || (cluster == "row")) {
                        Rowv = rowMeans(X.mat)
                        
                        if (dist.method[1] == "correlation")
                            dist.mat = as.dist(1 - cor(t(as.matrix(
                                X.mat
                            )),
                            method = "pearson"))
                        else
                            dist.mat = dist(X.mat, method = dist.method[1])
                        
                        hcr = hclust(dist.mat, method = clust.method[1])
                        ddr = as.dendrogram(hcr)
                        ddr = reorder(ddr, Rowv)
                        rowInd = order.dendrogram(ddr)
                        object = object[rowInd,]
                        row.names = row.names[rowInd]
                        
                        if (!is.null(row.sideColors))
                            row.sideColors = as.matrix(row.sideColors[rowInd,])
                    }
                    
                    if ((cluster == "both") || (cluster == "column")) {
                        Colv = rowMeans(cord.X)
                        
                        if (dist.method[2] == "correlation")
                            dist.mat = as.dist(1 - cor(t(cord.X), method = "pearson"))
                        else
                            dist.mat = dist(cord.X, method = dist.method[2])
                        
                        hcc = hclust(dist.mat, method = clust.method[2])
                        ddc = as.dendrogram(hcc)
                        ddc = reorder(ddc, Colv)
                        colInd = order.dendrogram(ddc)
                        object = object[, colInd]
                        col.names = col.names[colInd]
                        
                        if (!is.null(col.sideColors))
                            col.sideColors = as.matrix(col.sideColors[colInd,])
                        
                    }
                }
                
                #-- if mapping = "Y"
                if (mapping == "Y") {
                    #-- cheking center and scale
                    if (!is.logical(center)) {
                        if (!is.numeric(center) || (length(center) != ncol(mat$Y)))
                            stop(
                                "'center' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                                call. = FALSE
                            )
                    }
                    if (!is.logical(scale)) {
                        if (!is.numeric(scale) || (length(scale) != ncol(mat$Y)))
                            stop(
                                "'scale' should be either a logical value or a numeric
                    vector of length equal to the number of columns of 'Y'.",
                                call. = FALSE
                            )
                    }
                  
                    if (any(class.object %in% object.block.pls)) {
                      tmp.df <- mat$X[[Y.block.name]][, keep.Y]
                      tmp.var <- mat$variates[[Y.block.name]][, comp]
                    } else {
                      tmp.df <- mat$Y[, keep.Y]
                      tmp.var <- mat$variates$Y[, comp]
                    }
                  
                    object = scale(tmp.df, center = center, scale = scale)
                    Y.mat = as.matrix(tmp.var)
                    col.names = col.names[keep.Y]
                    
                    rm(tmp.df)
                    rm(tmp.var)
                    
                    if (!is.null(col.sideColors))
                        col.sideColors = as.matrix(col.sideColors[keep.Y,])
                    
                    
                    
                    if ((cluster == "both") || (cluster == "row")) {
                        Rowv = rowMeans(Y.mat)
                        
                        if (dist.method[1] == "correlation")
                            dist.mat = as.dist(1 - cor(t(as.matrix(
                                Y.mat
                            )),
                            method = "pearson"))
                        else
                            dist.mat = dist(Y.mat, method = dist.method[1])
                        
                        hcr = hclust(dist.mat, method = clust.method[1])
                        ddr = as.dendrogram(hcr)
                        ddr = reorder(ddr, Rowv)
                        rowInd = order.dendrogram(ddr)
                        object = object[rowInd,]
                        row.names = row.names[rowInd]
                        
                        if (!is.null(row.sideColors))
                            row.sideColors = as.matrix(row.sideColors[rowInd,])
                    }
                    
                    if ((cluster == "both") || (cluster == "column")) {
                        Colv = rowMeans(cord.Y)
                        
                        if (dist.method[2] == "correlation")
                            dist.mat = as.dist(1 - cor(t(cord.Y), method = "pearson"))
                        else
                            dist.mat = dist(cord.Y, method = dist.method[2])
                        
                        hcc = hclust(dist.mat, method = clust.method[2])
                        ddc = as.dendrogram(hcc)
                        ddc = reorder(ddc, Colv)
                        colInd = order.dendrogram(ddc)
                        object = object[, colInd]
                        col.names = col.names[colInd]
                        
                        if (!is.null(col.sideColors))
                            col.sideColors = as.matrix(col.sideColors[colInd,])
                    }
                }
                
                #-- calling the image.map function --------------------------------#
                #------------------------------------------------------------------#
                
                
                #-- output --------------------------------------------------------#
                #------------------------------------------------------------------#
                res = list(mat = object,
                           row.names = row.names,
                           col.names = col.names)
                
                if (!is.null(sample.sideColors)) {
                    res$sample.sideColors = sample.sideColors
                }
                
                if ((cluster == "both") || (cluster == "row")) {
                    res$rowInd = rowInd
                    res$ddr = ddr
                }
                
                if ((cluster == "both") || (cluster == "column")) {
                    res$colInd = colInd
                    res$ddc = ddc
                }
                class(res) = paste("cim", class.object[1], sep = "_")
            }
        } 
        else {
            #-- if matrix class  -------------------------------------
            
            #-- mat
            isMat = tryCatch(
                is.matrix(mat),
                error = function(e)
                    e
            )
            
            if ("simpleError" %in% class(isMat))
                stop(isMat[[1]], ".", call. = FALSE)
            
            if (!is.matrix(mat) || !is.numeric(mat))
                stop("'mat' must be a numeric matrix.", call. = FALSE)
            
            p = nrow(mat)
            q = ncol(mat)
            #-- row.names
            if (is.logical(row.names)) {
                if (isTRUE(row.names))
                    row.names = rownames(mat)
                else
                    row.names = rep("", p)
            }
            else {
                row.names = as.vector(row.names)
                if (length(row.names) != p)
                    stop("'row.names' must be a character vector of length ",
                         p,
                         ".",
                         call. = FALSE)
            }
            
            #-- col.names
            if (is.logical(col.names)) {
                if (isTRUE(col.names))
                    col.names = colnames(mat)
                else
                    col.names = rep("", q)
            }
            else {
                col.names = as.vector(col.names)
                if (length(col.names) != q)
                    stop("'col.names' must be a character vector of length ",
                         q,
                         ".",
                         call. = FALSE)
            }
            
            #-- row.sideColors
            if (!is.null(row.sideColors)) {
                row.sideColors = as.matrix(row.sideColors)
                if (nrow(row.sideColors) != p)
                    stop(
                        "'row.sideColors' must be a colors character vector (matrix)
            of length (nrow) ",
                        p,
                        ".",
                        call. = FALSE
                    )
            }
            
            #-- col.sideColors
            if (!is.null(col.sideColors)) {
                col.sideColors = as.matrix(col.sideColors)
                if (nrow(col.sideColors) != q)
                    stop(
                        "'col.sideColors' must be a colors character vector (matrix)
            of length (nrow) ",
                        q,
                        ".",
                        call. = FALSE
                    )
            }
            
            #-- clustering -------------------------------------------------------#
            #---------------------------------------------------------------------#
            object = mat
            if ((cluster == "both") || (cluster == "row")) {
                Rowv = rowMeans(mat)
                
                if (dist.method[1] == "correlation")
                    dist.mat = as.dist(1 - cor(t(as.matrix(mat)), method = "pearson"))
                else
                    dist.mat = dist(mat, method = dist.method[1])
                
                hcr = hclust(dist.mat, method = clust.method[1])
                ddr = as.dendrogram(hcr)
                ddr = reorder(ddr, Rowv)
                rowInd = order.dendrogram(ddr)
                object = mat[rowInd,]
                row.names = row.names[rowInd]
                
                if (!is.null(row.sideColors))
                    row.sideColors = as.matrix(row.sideColors[rowInd,])
            }
            
            if ((cluster == "both") || (cluster == "column")) {
                Colv = colMeans(mat)
                
                if (dist.method[2] == "correlation")
                    dist.mat = as.dist(1 - cor(as.matrix(mat), method = "pearson"))
                else
                    dist.mat = dist(t(mat), method = dist.method[2])
                hcc = hclust(dist.mat, method = clust.method[2])
                ddc = as.dendrogram(hcc)
                ddc = reorder(ddc, Colv)
                colInd = order.dendrogram(ddc)
                object = object[, colInd]
                col.names = col.names[colInd]
                
                if (!is.null(col.sideColors))
                    col.sideColors = as.matrix(col.sideColors[colInd,])
            }
            #-- calling the image.map function ------------------------------------#
            
            
            
            #-- output ------------------------------------------------------------#
            #----------------------------------------------------------------------#
            res = list(
                mat = object,
                row.names = row.names,
                col.names = col.names,
                row.sideColors = row.sideColors,
                col.sideColors = col.sideColors
            )
            
            if ((cluster == "both") || (cluster == "row")) {
                res$rowInd = rowInd
                res$ddr = ddr
            }
            
            if ((cluster == "both") || (cluster == "column")) {
                res$colInd = colInd
                res$ddc = ddc
            }
            
            class(res) = "cim_default"
            
        }
        
        #-- call imageMap  -------------------------------------
        opar = par(no.readonly = TRUE)
        
        try_plot <- tryCatch({
            imageMap(
                object,
                color = color,
                row.names = row.names,
                col.names = col.names,
                row.sideColors = row.sideColors,
                col.sideColors = col.sideColors,
                row.cex = row.cex,
                col.cex = col.cex,
                cluster = cluster,
                ddr = ddr,
                ddc = ddc,
                cut.tree = cut.tree,
                transpose = transpose,
                symkey = symkey,
                keysize = keysize,
                keysize.label = keysize.label,
                zoom = zoom,
                title = title,
                xlab = xlab,
                ylab = ylab,
                margins = margins,
                lhei = lhei,
                lwid = lwid
            )
        }, error = function(e)
            e)
        
        if (is(try_plot, "error")) {
            message(sprintf(
                "Error in cim plot: %s. See ?cim for help.",
                try_plot$message
            ))
        } else {
            
            #-- add to plot  -------------------------------------
            if (!is.null(legend) || mapping == "multiblock")
            {
                if (is.null(legend) && mapping == "multiblock") {
                  legend <- list()
                }
                if (is.null(legend$x))
                    legend$x = "topright"
                if (is.null(legend$bty))
                    legend$bty = "n"
                if (is.null(legend$cex))
                    legend$cex = 0.8
                
               # browser()
                if (any(class.object %in% object.block.pls)) {
                  if (is.null(legend$legend)) {
                    legend$legend <- blocks
                  }
                  if (is.null(legend$col)) {
                    legend$col <- unique(col.sideColors)
                  }
                }
                
                else if (any(class.object %in%  c("mixo_splsda", "mixo_plsda")))
                {
                    if (is.null(legend$legend))
                        legend$legend = mat$names$colnames$Y
                    
                    #-- col
                    if (is.null(legend$col)) {
                        if (!is.null(sample.sideColors))
                            legend$col = unique(as.matrix(sample.sideColors[order(map(mat$ind.mat)), 1]))
                    }
                    
                }
                else if (any(class.object %in%  c("mixo_mlsplsda")))
                {
                    if (is.null(legend$legend) && is.null(legend$col)) {
                        if (ncol(mat$multilevel) >= 2) {
                            df = data.frame(mat$multilevel[, 2], sample.sideColors[, 1])
                            df = unique(df)
                            legend$legend = as.character(df[, 1])
                            legend$col = as.character(df[, 2])
                        }
                        if (ncol(mat$multilevel) == 3) {
                            df = data.frame(mat$multilevel[, 3], sample.sideColors[, 2])
                            df = unique(df)
                            legend$legend = c(legend$legend, as.character(df[, 1]))
                            legend$col = c(legend$col, as.character(df[, 2]))
                        }
                    }
                }
                else if (any(class.object %in%  c("mixo_mlspls")))
                {
                    if (mapping != "XY") {
                        if (is.null(legend$legend) && is.null(legend$col)) {
                            if (ncol(mat$multilevel) >= 2) {
                                df = data.frame(mat$multilevel[, 2],
                                                sample.sideColors[, 1])
                                df = unique(df)
                                legend$legend = as.character(df[, 1])
                                legend$col = as.character(df[, 2])
                            }
                            if (ncol(mat$multilevel) == 3) {
                                df = data.frame(mat$multilevel[, 3],
                                                sample.sideColors[, 2])
                                df = unique(df)
                                legend$legend = c(legend$legend,
                                                  as.character(df[, 1]))
                                legend$col = c(legend$col, as.character(df[, 2]))
                            }
                        }
                    }
                }
                if (is.null(legend$legend))
                    stop("argument \"legend$legend\" is missing, with no default")
                
                #-- fill
                if (is.null(legend$fill))
                    legend$fill = legend$col
                
                
                par(mar = c(0, 0, 0, 0), new = TRUE)
                plot(
                    0,
                    0,
                    axes = FALSE,
                    type = "n",
                    xlab = "",
                    ylab = ""
                )
                
                if (!is.null(legend$title))
                {
                    legend(
                        x = legend$x,
                        y = legend$y,
                        legend = legend$legend,
                        col = legend$col,
                        fill = legend$fill,
                        bty = legend$bty,
                        title = legend$title,
                        cex = legend$cex
                    )
                } else{
                    legend(
                        x = legend$x,
                        y = legend$y,
                        legend = legend$legend,
                        col = legend$col,
                        fill = legend$fill,
                        bty = legend$bty,
                        cex = legend$cex
                    )
                }
                
                
            }
          else {
            
          }
            if (any(class.object %in% object.all) &
                !any(class.object %in%
                     object.single.omics) &
                mapping == "XY")
                res$mat.cor = object
            par(opar)
            
            if (!is.null(save))
                dev.off()
        }
        
        
        return(invisible(res))
    }

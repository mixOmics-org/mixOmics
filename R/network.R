#' Relevance Network for (r)CCA and (s)PLS regression
#' 
#' Display relevance associations network for (regularized) canonical
#' correlation analysis and (sparse) PLS regression. The function avoids the
#' intensive computation of Pearson correlation matrices on large data set by
#' calculating instead a pair-wise similarity matrix directly obtained from the
#' latent components of our integrative approaches (CCA, PLS, block.pls
#' methods). The similarity value between a pair of variables is obtained by
#' calculating the sum of the correlations between the original variables and
#' each of the latent components of the model. The values in the similarity
#' matrix can be seen as a robust approximation of the Pearson correlation (see
#' González et al. 2012 for a mathematical demonstration and exact formula).
#' The advantage of relevance networks is their ability to simultaneously
#' represent positive and negative correlations, which are missed by methods
#' based on Euclidean distances or mutual information. Those networks are
#' bipartite and thus only a link between two variables of different types can
#' be represented. The network can be saved in a .glm format using the
#' \code{igraph} package, the function \code{write.graph} and extracting the
#' output \code{object$gR}, see details.
#' 
#' \code{network} allows to infer large-scale association networks between the
#' \eqn{X} and \eqn{Y} datasets in \code{rcc} or \code{spls}. The output is a
#' graph where each \eqn{X}- and \eqn{Y}-variable corresponds to a node and the
#' edges included in the graph portray associations between them.
#' 
#' In \code{rcc}, to identify \eqn{X}-\eqn{Y} pairs showing relevant
#' associations, \code{network} calculate a similarity measure between \eqn{X}
#' and \eqn{Y} variables in a pair-wise manner: the scalar product value
#' between every pairs of vectors in dimension \code{length(comp)} representing
#' the variables \eqn{X} and \eqn{Y} on the axis defined by \eqn{Z_i} with
#' \eqn{i} in \code{comp}, where \eqn{Z_i} is the equiangular vector between
#' the \eqn{i}-th \eqn{X} and \eqn{Y} canonical variate.
#' 
#' In \code{spls}, if \code{object$mode} is \code{regression}, the similarity
#' measure between \eqn{X} and \eqn{Y} variables is given by the scalar product
#' value between every pairs of vectors in dimension \code{length(comp)}
#' representing the variables \eqn{X} and \eqn{Y} on the axis defined by
#' \eqn{U_i} with \eqn{i} in \code{comp}, where \eqn{U_i} is the \eqn{i}-th
#' \eqn{X} variate. If \code{object$mode} is \code{canonical} then \eqn{X} and
#' \eqn{Y} are represented on the axis defined by \eqn{U_i} and \eqn{V_i}
#' respectively.
#' 
#' Variable pairs with a high similarity measure (in absolute value) are
#' considered as relevant. By changing the cut-off, one can tune the relevance
#' of the associations to include or exclude relationships in the network.
#' 
#' \code{interactive=TRUE} open two device, one for association network, one
#' for scrollbar, and define an interactive process: by clicking either at each
#' end (\eqn{-} or \eqn{+}) of the scrollbar or at middle portion of this.
#' The position of the slider indicate which is the `cutoff' value associated
#' to the display network.
#' 
#' The network can be saved in a .glm format using the \pkg{igraph} package,
#' the function \code{write.graph} and extracting the output \code{obkect$gR}.
#' 
#' The interactive process is terminated by clicking the second button and
#' selecting \code{Stop} from the menu, or from the \code{Stop} menu on the graphics
#' window.
#' 
#' The \code{color.node} is a vector of length two, of any of the three kind of
#' \code{R} colors, i.e., either a color name (an element of \code{colors()}),
#' a hexadecimal string of the form \code{"#rrggbb"}, or an integer \code{i}
#' meaning \code{palette()[i]}. \code{color.node[1]} and \code{color.node[2]}
#' give the color for filled nodes of the \eqn{X}- and \eqn{Y}-variables
#' respectively. Defaults to \code{c("white", "white")}.
#' 
#' \code{color.edge} give the color to edges with colors corresponding to the
#' values in \code{mat}. Defaults to \code{color.GreenRed(100)} for negative
#' (green) and positive (red) correlations. We also propose other palettes of
#' colors, such as \code{color.jet} and \code{color.spectral}, see help on
#' those functions, and examples below. Other palette of colors from the stats
#' package can be used too.
#' 
#' \code{shape.node[1]} and \code{shape.node[2]} provide the shape of the nodes
#' associate to \eqn{X}- and \eqn{Y}-variables respectively. Current acceptable
#' values are \code{"circle"} and \code{"rectangle"}. Defaults to
#' \code{c("circle", "rectangle")}.
#' 
#' \code{lty.edge[1]} and \code{lty.egde[2]} give the line type to edges with
#' positive and negative weight respectively. Can be one of \code{"solid"},
#' \code{"dashed"}, \code{"dotted"}, \code{"dotdash"}, \code{"longdash"} and
#' \code{"twodash"}. Defaults to \code{c("solid", "solid")}.
#' 
#' \code{lwd.edge[1]} and \code{lwd.edge[2]} provide the line width to edges
#' with positive and negative weight respectively. This attribute is of type
#' double with a default of \code{c(1, 1)}.
#' 
#' @aliases network network.default network.rcc network.pls network.spls
#' @param mat numeric matrix of values to be represented. Alternatively,
#' an object from one of the following models: \code{mix_pls}, \code{plsda}, 
#' \code{mixo_spls}, \code{splsda}, \code{rcc}, \code{sgcca}, \code{rgcca}, 
#' \code{sgccda}.
#' @param comp atomic or vector of positive integers. The components to
#' adequately account for the data association. Defaults to \code{comp = 1:ncomp},
#' i.e. all components are used.
#' @param cutoff numeric value between \code{0} and \code{1}. The tuning
#' threshold for the relevant associations network (see Details).
#' @param row.names,col.names character vector containing the names of \eqn{X}-
#' and \eqn{Y}-variables.
#' @param graph.scale Numeric between 0 and 1 which alters the scale of the entire
#' plot. Increasing the value decreases the size of nodes and increases their distance
#' from one another. Defaults to 0.5.
#' @param size.node Numeric between 0 and 1 which determines the relative size of nodes.
#' Defaults to 0.5.
#' @param color.node vector of length two, the colors of the \eqn{X} and
#' \eqn{Y} nodes (see Details).
#' @param shape.node character vector of length two, the shape of the \eqn{X}
#' and \eqn{Y} nodes (see Details).
#' @param alpha.node Numeric between 0 and 1 which determines the opacity of nodes.
#' Only used in block objects.
#' @param color.edge vector of colors or character string specifying the colors
#' function to using to color the edges, set to default to
#' \code{color.GreenRed(100)} but other palettes can be chosen (see Details and
#' Examples).
#' @param lty.edge character vector of length two, the line type for the edges
#' (see Details).
#' @param lwd.edge vector of length two, the line width of the edges (see
#' Details).
#' @param show.edge.labels logical. If \code{TRUE}, plot association values as
#' edge labels (defaults to \code{FALSE}).
#' @param show.color.key Logical. If \code{TRUE} a color key should be plotted.
#' @param symkey Logical indicating whether the color key should be made
#' symmetric about 0. Defaults to \code{TRUE}.
#' @param keysize numeric value indicating the size of the color key.
#' @param keysize.label vector of length 1, indicating the size of the labels
#' and title of the color key.
#' @param breaks (optional) either a numeric vector indicating the splitting
#' points for binning \code{mat} into colors, or a integer number of break
#' points to be used, in which case the break points will be spaced equally
#' between \code{min(mat)} and \code{max(mat)}.
#' @param interactive logical. If \code{TRUE}, a scrollbar is created to change
#' the cutoff value interactively (defaults to \code{FALSE}). See Details.
#' @param save should the plot be saved ? If so, argument to be set either to
#' \code{'jpeg'}, \code{'tiff'}, \code{'png'} or \code{'pdf'}.
#' @param name.save character string giving the name of the saved file.
#' @param cex.edge.label the font size for the edge labels.
#' @param cex.node.name the font size for the node labels.
#' @param blocks a vector indicating the block variables to display.
#' @param block.var.names either a list of vector components for variable names
#' in each block or FALSE for no names. If TRUE, the columns names of the
#' blocks are used as names.
#' @param layout.fun a function. It specifies how the vertices will be placed
#' on the graph. See help(layout) in the igraph package. Defaults to
#' layout.fruchterman.reingold.
#' @param plot.graph logical. If \code{TRUE} (default), plotting window will be 
#' filled with network. If \code{FALSE}, then no graph will be plotted, though 
#' the return value of the function is the exact same.
#' @return \code{network} return a list containing the following components:
#' \item{M}{the correlation matrix used by \code{network}.} \item{gR}{a
#' \code{graph} object to save the graph for cytoscape use (requires to load
#' the \pkg{igraph} package).}
#' @section Warning: If the number of variables is high, the generation of the
#' network generation can take some time.
#' @author Ignacio González, Kim-Anh Lê Cao, AL J Abadi
#' @seealso \code{\link{plotVar}}, \code{\link{cim}},
#' \code{\link{color.GreenRed}}, \code{\link{color.jet}},
#' \code{\link{color.spectral}} and http: //www.mixOmics.org for more details.
#' @references Mathematical definition: González I., Lê Cao K-A., Davis, M.J.
#' and Déjean, S. (2012). Visualising associations between paired omics data
#' sets. J. Data Mining 5:19.
#' \url{http://www.biodatamining.org/content/5/1/19/abstract}
#' 
#' Examples and illustrations:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' 
#' Relevance networks:
#' 
#' Butte, A. J., Tamayo, P., Slonim, D., Golub, T. R. and Kohane, I. S. (2000).
#' Discovering functional relationships between RNA expression and
#' chemotherapeutic susceptibility using relevance networks. \emph{Proceedings
#' of the National Academy of Sciences of the USA} \bold{97}, 12182-12186.
#' 
#' Moriyama, M., Hoshida, Y., Otsuka, M., Nishimura, S., Kato, N., Goto, T.,
#' Taniguchi, H., Shiratori, Y., Seki, N. and Omata, M. (2003). Relevance
#' Network between Chemosensitivity and Transcriptome in Human Hepatoma Cells.
#' \emph{Molecular Cancer Therapeutics} \bold{2}, 199-205.
#' @keywords multivariate graphs dplot hplot iplot
#' @export
#' @example ./examples/network-examples.R
network <- function(mat,
                    comp = NULL,
                    blocks = c(1, 2),
                    cutoff = 0,
                    row.names = TRUE,
                    col.names = TRUE,
                    block.var.names = TRUE,
                    graph.scale = 0.5,
                    size.node = 0.5,
                    color.node = NULL,
                    shape.node = NULL,
                    alpha.node = 0.85,
                    cex.node.name = NULL,
                    color.edge = color.GreenRed(100),
                    lty.edge = "solid",
                    lwd.edge = 1,
                    show.edge.labels = FALSE,
                    cex.edge.label = 1,
                    show.color.key = TRUE,
                    symkey = TRUE,
                    keysize = c(1, 1),
                    keysize.label = 1,
                    breaks,
                    interactive = FALSE,
                    layout.fun = NULL,
                    save = NULL,
                    name.save = NULL,
                    plot.graph = TRUE
                    
)
{
    #-- checking general input parameters --------------------------------------#
    #---------------------------------------------------------------------------#
    
    #-- check that the user did not enter extra arguments
    arg.call = match.call()
    user.arg = names(arg.call)[-1]
    
    err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                   error = function(e) e)
    
    if ("simpleError" %in% class(err))
        stop(err[[1]], ".", call. = FALSE)
    
    function.arg = names(mget(names(formals()), sys.frame(sys.nframe())))
    not.arg = !(user.arg %in% function.arg)
    
    if (any(not.arg))
    {
        unused.arg = user.arg[not.arg]
        not.arg = which(not.arg) + 1
        output = rep("", length(not.arg))
        
        for (i in 1:length(not.arg))
        {
            output[i] = paste0(unused.arg[i], " = ", arg.call[[not.arg[i]]])
        }
        
        output = paste0("(", paste(output, collapse = ", "), ").")
        msg = "unused argument "
        if (length(not.arg) > 1) msg = "unused arguments "
        stop(msg, output, call. = FALSE)
    }
    
    #    #-- check blocks
    #    if(length(blocks) != 2)
    #    stop("We can only display 2 blocks",call.=FALSE)
    
    #-- save
    if (!is.null(save))
    {
        if (! save %in% c("jpeg","tiff","png","pdf"))
            stop("'save' must be one of 'jpeg', 'png', 'tiff' or 'pdf'.", call. = FALSE)
    }
    
    #-- name.save
    if (!is.null(name.save))
    {
        if (! is.character(name.save) || length(name.save) > 1)
            stop("'name.save' must be a character.", call. = FALSE)
    } else {
        if (!is.null(save))
            name.save = paste0("network_",gsub(".", "_", deparse(substitute(mat)) ,fixed = TRUE))
    }
    
    if (!is.null(save)){
        
        while (dev.cur()>1)
            dev.off()
        
        if (save == "jpeg")
            jpeg(filename = paste0(name.save,".jpeg"), res = 600, width = 4000, height = 4000)
        if (save == "png")
            jpeg(filename = paste0(name.save,".png"), res = 600, width = 4000, height = 4000)
        if (save == "tiff")
            tiff(filename = paste0(name.save,".tiff"), res = 600, width = 4000, height = 4000)
        if (save == "pdf")
            pdf(file = paste0(name.save,".pdf"))
        
    }
    
    if (!plot.graph & interactive) {
        stop("plot.graph cannot be FALSE if interactive = TRUE", call.=FALSE)
    }
    
    class.object = class(mat)
    object.pls=c("mixo_pls","mixo_spls","mixo_mlspls")
    object.rcc="rcc"
    object.blocks=c("sgcca","rgcca")
    
    if (! any(class.object %in% c(object.pls,object.rcc,object.blocks, "matrix")))
        stop( " 'network' is only implemented for the following objects: matrix, pls, plsda, spls, splsda, rcc, sgcca, rgcca, sgccda", call.=FALSE)
    
    
    if(any(class.object %in% c(object.rcc,object.pls)))
    {
        p = ncol(mat$X)
        if(any(class.object == "DA")) # object is DA
            mat$Y = mat$ind.mat
        
        q = ncol(mat$Y)
        n = nrow(mat$X)
        ncomp = mat$ncomp
        
        
        #-- comp
        if(is.null(comp))
            comp=1:mat$ncomp
        if (length(comp) == 1)
        {
            if(comp>ncomp)
            {
                stop("the elements of 'comp' must be smaller than or equal to ", ncomp, ".",
                     call. = FALSE)
            } else if ( !is.numeric(comp) || comp <= 0) {
                stop("invalid value for 'comp'.", call. = FALSE)
            }
        }
        
        if (length(comp) > 1)
        {
            if(length(comp) > ncomp)
                stop("the length of 'comp' must be smaller than or equal to ", ncomp, ".",
                     call. = FALSE)
            
            if (!is.numeric(comp) || any(comp < 1))
                stop("invalid vector for 'comp'.", call. = FALSE)
            
            if (any(comp > ncomp))
                stop("the elements of 'comp' must be <= ", ncomp, ".",
                     call. = FALSE)
        }
        
        comp = round(comp)
        
        
        #-- row.names
        row.names.plot = TRUE # whether to plot the label in the graph
        if (is.logical(row.names))
        {
            if(!isTRUE(row.names))
            {
                row.names.plot = FALSE
            }
            row.names = mat$names$colnames$X
        } else {
            row.names = as.vector(row.names)
            if (length(unique(row.names)) != p)
                stop("'row.names' must be a character vector of ", p, " unique entries.",
                     call. = FALSE)
        }
        
        if(row.names.plot == TRUE)
        {
            row.names.plot = row.names
        }else{
            row.names.plot = rep("",p)
        }
        
        #-- col.names
        col.names.plot = TRUE # whether to plot the label in the graph
        if (is.logical(col.names))
        {
            if(!isTRUE(col.names))
            {
                col.names.plot = FALSE
            }
            col.names = mat$names$colnames$Y
        } else {
            col.names = as.vector(col.names)
            if (length(col.names) != q)
                stop("'col.names' must be a character vector of ", q, " unique entries.",
                     call. = FALSE)
        }
        
        if(col.names.plot == TRUE)
        {
            col.names.plot = col.names
        }else{
            col.names.plot = rep("",q)
        }
        
        #-- end checking --#
        #------------------#
        
        #-- network ----------------------------------------------------------------#
        #---------------------------------------------------------------------------#
        if(any(class.object %in% object.rcc))
        {
            #-- similarity matrix --#
            bisect = mat$variates$X[, comp] + mat$variates$Y[, comp]
            cord.X = cor(mat$X, bisect, use = "pairwise")
            cord.Y = cor(mat$Y, bisect, use = "pairwise")
            mat = cord.X %*% t(cord.Y)
        } else if(any(class.object %in% object.pls)) {
            #-- variable selection --#
            if (all(class(mat) %in% "mixo_pls"))
            {
                keep.X = rep(TRUE,p)
                keep.Y = rep(TRUE,q)
            } else {
                keep.X = apply(abs(mat$loadings$X[, comp, drop=FALSE]), 1, sum) > 0
                keep.Y = apply(abs(mat$loadings$Y[, comp, drop=FALSE]), 1, sum) > 0
                
                row.names = row.names[keep.X]
                col.names = col.names[keep.Y]
            }
            
            #-- similarity matrix --#
            if (mat$mode == "canonical")
            {
                cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp], use = "pairwise")
                cord.Y = cor(mat$Y[, keep.Y], mat$variates$Y[, comp], use = "pairwise")
            } else {
                cord.X = cor(mat$X[, keep.X], mat$variates$X[, comp], use = "pairwise")
                cord.Y = cor(mat$Y[, keep.Y], mat$variates$X[, comp], use = "pairwise")
            }
            
            mat = cord.X %*% t(cord.Y)
        }
        
    } else if(any(class.object %in% object.blocks)) {
        
        mat$X <- mapply(mat$X, names(mat$X), FUN =  function(x, y)
        {
            colnames(x) <- paste0(colnames(x),"_", y)
            x
        }, SIMPLIFY = FALSE)
        
        # remove Y from the list of blocks for DA objects
        if(any(class.object == "DA"))
        {
            mat$names$blocks = mat$names$blocks [-mat$indY]
            mat$names$colnames = mat$names$colnames [-mat$indY]
            mat$ncomp = mat$ncomp [-mat$indY]
        }
        
        if (is.null(blocks))
        {
            if (any(mat$ncomp > 1))
            {
                blocks = mat$names$blocks[ which(mat$ncomp > 1)]
            } else {
                stop(("The number of components for each block is 1. The number of components must be superior or equal to 2."), call. = FALSE)
            }
        } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <= length(mat$names$blocks)) {
            blocks = mat$names$blocks[blocks]
        } else if (is.character(blocks)) {
            if (!all(blocks %in% mat$names$blocks))
                stop("One element of 'blocks' does not match with the names of the blocks")
        } else {
            stop("Incorrect value for 'blocks", call. = FALSE)
        }
        
        #-- comp
        if (is.null(comp))
        {
            comp = vector("list", length(blocks))
            names(comp) = blocks
            
            for (i in blocks)
                comp[[i]] = 1:mat$ncomp[i]
        }
        
        if (is.list(comp))
        {
            if (length(comp) != length(blocks))
                stop("'comp' must be either NULL a list of length ", length(blocks), ".",
                     call. = FALSE)
            
            if (!all(blocks %in% names(comp)))
                stop("names of 'comp' must be from {",
                     paste(blocks, collapse = ", "), "}.", call. = FALSE)
            
            for (i in blocks)
            {
                if (any(!is.finite(comp[[i]])))
                    stop("invalid value for 'comp' of the block '", i, "'.", call. = FALSE)
                
                if (any(comp[[i]] > mat$ncomp[i]))
                    stop("the elements of 'comp' for block '", i, "' must be <= ",
                         mat$ncomp[i], ".", call. = FALSE)
                
                if (any(comp[[i]] < 1))
                    stop("invalid value for 'comp' of the block '", i, "'.", call. = FALSE)
                
            }
            
        } else {
            stop("'comp' must be either NULL or a list of length ", length(blocks), ".",
                 call. = FALSE)
        }
        
        
        #-- block.var.names
        num.var = unlist(lapply(mat$X[blocks], ncol))
        
        
        if (is.logical(block.var.names))
        {
            if (length(block.var.names)==1)
                block.var.names=rep(block.var.names,length(blocks))
            if (length(block.var.names) != length(blocks))
                stop("'block.var.names' must be a logical vector of length 1 or ",  length(blocks),", or a list of length ",  length(blocks), ".",
                     call. = FALSE)
            
            vec=(which(block.var.names==FALSE))
            
            block.var.names = mat$names$colnames
            
            for (i in 1:length(blocks))
            {
                if (i %in% vec)
                    block.var.names[[blocks[i]]] = rep(" ",length(mat$names$colnames[[blocks[i]]]))
            }
        } else {
            if (is.list(block.var.names))
            {
                if (length(block.var.names) != length(blocks))
                {
                    stop("'block.var.names' must be a logical vector or a list of length ",  length(blocks), ".",
                         call. = FALSE)
                } else {
                    if (!all(unlist(lapply(block.var.names, is.vector))))
                        stop("each component of 'block.var.names' must be a vector.", call. = FALSE)
                    
                    block.var.names.length = unlist(lapply(block.var.names, length))
                    
                    if (any(block.var.names.length != num.var))
                        stop("components of 'block.var.names' must be vectors of length ",
                             paste(num.var, collapse = ", "), ".", call. = FALSE)
                    
                }
            } else {
                stop("'block.var.names' must be either a logical value or a list of length ",
                     length(blocks), ".", call. = FALSE)
            }
        }
        #-- network approach -------------------------------------------------------#
        #---------------------------------------------------------------------------#
        
        #-- calculation of the similarity matrix for each block --#
        coord = M_block = list()
        j = 1
        
        if (is(mat, "sgcca"))
        {
            for (k in blocks)
            {
                if (length(comp[[k]]) > 1)
                {
                    keep = (apply(abs(mat$loadings[[k]][, comp[[k]]]), 1, sum) > 0)
                } else {
                    keep = abs(mat$loadings[[k]][, comp[[k]]]) > 0
                }
                
                coord[[j]] = cor(mat$X[[k]][, keep], mat$variates[[k]][, comp[[k]]], use = "pairwise")
                j = j + 1
            }
        } else {
            for(k in blocks)
            {
                coord[[j]] = cor(mat$X[[k]], mat$variates[[k]][, comp[[k]]], use = "pairwise")
                j = j + 1
            }
        }
        
        node.X = node.Y = w = NULL
        l = 1
        
        for (j in 1:(length(blocks) - 1))
        {
            for (k in (j + 1):length(blocks))
            {
                if(!any(comp[[blocks[j]]] %in% comp[[blocks[k]]]))
                    stop("comp of block ",blocks[j], " is ",  comp[[blocks[j]]], " but comp of block ", blocks[k]," is ",comp[[blocks[k]]],
                         call. = FALSE)
                int.comp = intersect(comp[[blocks[j]]], comp[[blocks[k]]])
                
                object = coord[[j]][, comp[[blocks[j]]] %in% int.comp] %*% t(coord[[k]][, comp[[blocks[k]]] %in% int.comp])
                M_block[[l]] = object
                l = l + 1
                
                X = rownames(coord[[j]])
                Y = rownames(coord[[k]])
                
                rep.X = rep(X, each = length(Y))
                rep.Y = rep(Y, length(X))
                
                node.X = c(node.X, rep.X)
                node.Y = c(node.Y, rep.Y)
                
                w = c(w, as.vector(t(object)))
            }
        }
        
    } else {
        #-- mat
        if (!is.matrix(mat))
            stop("'mat' must be a numeric matrix.", call. = FALSE)
        
        if (length(dim(mat)) != 2)
            stop("'mat' must be a numeric matrix.")
        
        if (!is.numeric(mat))
            stop("'mat' must be a numeric matrix.")
        
        p = nrow(mat)
        q = ncol(mat)
        
        #-- row.names
        row.names.plot = TRUE # whether to plot the label in the graph
        if (is.logical(row.names))
        {
            if(!isTRUE(row.names))
            {
                row.names.plot = FALSE
            }
            row.names = rownames(mat)
        } else {
            row.names = as.vector(row.names)
            if (length(row.names) != p)
                stop("'row.names' must be a character vector of ", p, " unique entries.",
                     call. = FALSE)
        }
        
        if(row.names.plot == TRUE)
        {
            row.names.plot = row.names
        }else{
            row.names.plot = rep("",p)
        }
        
        #-- col.names
        col.names.plot = TRUE # whether to plot the label in the graph
        if (is.logical(col.names))
        {
            if(!isTRUE(col.names))
            {
                col.names.plot = FALSE
            }
            col.names = colnames(mat)
        } else {
            col.names = as.vector(col.names)
            if (length(col.names) != q)
                stop("'col.names' must be a character vector of ", q, " unique entries.",
                     call. = FALSE)
        }
        
        if(col.names.plot == TRUE)
        {
            col.names.plot = col.names
        }else{
            col.names.plot = rep("",q)
        }
    }
    
    #-- size.node
    if (!is.null(size.node)) {
      if (!is.finite(size.node) || size.node < 0 || size.node > 1) {
        stop("'size.node' must be a numerical value between 0 - 1.", call. = FALSE)
      }
      if (size.node == 0) {
        shape.node <- "none"
        size.node <- 1 # prevents error 
      }
    }
    
    #-- graph.scale
    if (!is.finite(graph.scale) || graph.scale < 0 || graph.scale > 1) {
      stop("'graph.scale' must be a numerical value between 0 - 1.", call. = FALSE)
    }
    
    #-- color.node
    if(any(class.object %in% object.blocks))
    {
        if (is.null(color.node))
        {
            ## see circosPlot
            color.blocks = brewer.pal(n = 12, name = 'Paired') #why 12?? ANS: bc max allowed n is 12 for this function
            if (length(blocks) > 6) {
                color.blocks <- colorRampPalette(color.blocks)(2*length(object$X))
            }
            color.node = color.blocks[seq(from = 1, to = 2*length(blocks), by = 2)]
            color.node = adjustcolor(color.node, alpha.f = alpha.node)
            names(color.node) = blocks
        } else {
            if (!is.vector(color.node) || length(color.node) != length(blocks))
                stop("'color.node' must be a vector of length ", length(blocks), ".", call. = FALSE)
        }
    } else {
        if(is.null(color.node))
            color.node = brewer.pal(n = 12, name = 'Paired')[c(1,3)]
            color.node = adjustcolor(color.node, alpha.f = alpha.node)
        
        if (!is.list(color.node))
        {
            if (!is.vector(color.node) || length(color.node) != 2)
                stop("'color.node' must be a vector of length 2.", call. = FALSE)
            
        } else {
            stop("'color.node' must be a vector of length 2.", call. = FALSE)
        }
    }
    
    if (any(!sapply(color.node, function(x){tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })))
        stop("'color.node' must be a character vector of recognized colors.",
             call. = FALSE)
    
    #-- shape.node
    if(any(class.object%in% object.blocks))
    {
        if (is.null(shape.node))
            shape.node = "circle"
        
        if (is.vector(shape.node))
        {
            if (length(shape.node) == 1)
                shape.node = rep(shape.node, length(blocks))
        }
        
        if (!is.list(shape.node))
        {
            if (!is.vector(shape.node) || length(shape.node) != length(blocks))
                stop("'shape.node' must be a character vector of length ", length(blocks), ".",
                     call. = FALSE)
        } else {
            stop("'shape.node' must be a numeric vector of length ", length(blocks), ".",
                 call. = FALSE)
        }
        
        if (!all(shape.node %in% c("none", "circle", "rectangle")))
            stop("elements of 'shape.node' must be from {'none', 'circle', 'rectangle'}.",
                 call. = FALSE)
        
        if (is.null(names(shape.node)))
            names(shape.node) = blocks
        
    } else {
        if(is.null(shape.node))
            shape.node=c("circle", "rectangle")
        if (length(shape.node)==1) {
          if(shape.node=="none") {
            shape.node=c("none", "none")
          }
        }
        
        if (!is.list(shape.node))
        {
            if (!is.vector(shape.node) || length(shape.node) != 2)
                stop("'shape.node' must be a vector of length 2.", call. = FALSE)
        } else {
            stop("'shape.node' must be a vector of length 2.", call. = FALSE)
        }
        
        if (!all(shape.node %in% c("none", "circle", "rectangle")))
            stop("elements of 'shape.node' must be from {'none', 'circle', 'rectangle'}.",
                 call. = FALSE)
        
    }
    
    #-- cex.node.name
    if (!is.null(cex.node.name)) {
      if (!is.finite(cex.node.name) || cex.node.name < 0 || length(cex.node.name)>1) {
        stop("'cex.node.name' must be a non-negative numerical value.", call. = FALSE)
      }
    }
    
        
    
    #-- color.edge
    if (length(color.edge) < 2 && (!is(color.edge, "function")))
        stop("'color.edge' must be a vector of length larger than or equal to 2.", call. = FALSE)
    
    if ((length(color.edge) %% 2) != 0 && (!is(color.edge, "function")) && isTRUE(symkey))
        stop("'color.edge' must be a vector of length an even number if 'symkey = TRUE'.", call. = FALSE)
    
    if (any(!sapply(color.edge, function(x) {tryCatch(is.matrix(col2rgb(x)), error = function(e) FALSE) })))
        stop("'color.edge' must be a character vector of recognized colors.", call. = FALSE)
    
    #-- lty.edge
    if(length(lty.edge) >2)
        stop("\"lty.edge\" must a character vector of up to 2 entries from
        'solid', 'dashed', 'dotted', 'dotdash', 'longdash', twodash' or 'blank'.
        see ?network",
             call.=FALSE)
    
    if (length(lty.edge)==1) lty.edge = c(lty.edge, lty.edge)
    
    choices = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash", "blank")
    lty.edge = choices[pmatch(lty.edge, choices,duplicates.ok=TRUE)]
    
    if (any(is.na(lty.edge)))
        stop("'lty.edge' should be from 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', twodash' or 'blank'.",
             call. = FALSE)
    
    
    #-- lwd.edge
    if(length(lwd.edge) >2)
        stop("'lwd.edge' must be a vector of up to 2 positive numbers.
        See ?network")
    if (length(lwd.edge)==1) 
        lwd.edge = c(lwd.edge, lwd.edge)
    
    if (length(lwd.edge) != 2 || any(!is.finite(lwd.edge)) || any(lwd.edge <= 0))
        stop("'lwd.edge' must be positive.")
    
    #-- show.edge.labels
    if (!is.logical(show.edge.labels))
        stop("'show.edge.labels' must be a logical constant (TRUE or FALSE).",
             call. = FALSE)
    
    #-- cex.edge.label
    if (!is.finite(cex.edge.label) || cex.edge.label < 0 || length(cex.edge.label)>1)
        stop("'cex.edge.label' must be a non-negative numerical value.", call. = FALSE)
    
    #-- show.color.key
    if (!is.logical(show.color.key))
        stop("'show.color.key' must be a logical constant (TRUE or FALSE).",
             call. = FALSE)
    
    #-- symkey
    if (!is.logical(symkey))
        stop("'symkey' must be a logical constant (TRUE or FALSE).",
             call. = FALSE)
    
    #-- keysize
    if (length(keysize) != 2 || any(!is.finite(keysize)))
        stop("'keysize' must be a numeric vector of length 2.",
             call. = FALSE)
    
    #-- keysize.label
    if (length(keysize.label) != 1 || any(!is.finite(keysize)))
        stop("'keysize' must be a numeric vector of length 1.",
             call. = FALSE)
    
    #-- interactive
    if (!is.logical(interactive))
        stop("'interactive' must be a logical constant (TRUE or FALSE).",
             call. = FALSE)
    
    #-- layout.fun
    if (!is.null(layout.fun) && !is(layout.fun, "function"))
        stop("'layout.fun' must be a valid layout function.", call. = FALSE)
    
    #-- end checking --#
    #------------------#
    
    
    #-- network approach -------------------------------------------------------#
    #---------------------------------------------------------------------------#
    if (!(any(class.object %in% object.blocks)))
        w = as.vector(t(mat))
    
    #-- check cutoff
    if (round(max(abs(w)), 2) == 0)
        stop("There is no correlation between these blocks whith these components.",
             "Try a different value of 'comp'.", call. = FALSE)
    
    if (!is.numeric(cutoff) || cutoff < 0 || cutoff > 1)
        stop("'cutoff' should be a numeric between 0 and 1", call. = FALSE)
    if(cutoff > max(abs(w)))
        stop("You have chosen a high cutoff value of ", cutoff, 
             " which is greaer than the max value in the similarity matrix which is ",
             round(max(abs(w)), 2), call. = FALSE)

        # Definition of nodes #
    #---------------------#
    #save(list=ls(),file="temp.Rdata")
    if(any(class.object %in% object.blocks))
    {
        group = NULL
        temp = lapply(mat$X, function(x) colnames(x))
        
        for (i in 1:length(temp))
        {
            group = c(group, rep(names(temp)[i], length(temp[[i]])))
        }
        
        nodes = data.frame(name = unlist(temp), group = group)
    } else if(any(class.object %in% object.pls)) {
        w = as.vector(t(mat))
        
        Xn=sum(keep.X) #number of non-zero parameters in X (over all comp)
        Yn=sum(keep.Y) #number of non-zero parameters in Y (over all comp)
        node.X = row.names#[keep.X]#paste0("X", 1:Xn)
        node.Y = col.names#[keep.Y]#paste0("Y", 1:Yn)
        
        row.names.plot = row.names.plot[keep.X]
        col.names.plot = col.names.plot[keep.Y]
        
        nodes = data.frame(name = c(node.X, node.Y),
                           group = c(rep("x", Xn), rep("y", Yn)))
        
        
        node.X = rep(node.X, each = Yn)
        node.Y = rep(node.Y, Xn)
    } else {
        node.X = row.names # paste0("X", 1:p)
        node.Y = col.names # paste0("Y", 1:q)
        
        nodes = data.frame(name = c(node.X, node.Y),
                           group = c(rep("x", p), rep("y", q)))
        
        
        node.X = rep(node.X, each = q)
        node.Y = rep(node.Y, p)
    }
    
    # Definition of edges #
    #---------------------#
    relations = data.frame(from = node.X, to = node.Y, weight = w)
    
    # edges colors #
    #--------------#
    id = bin.color(w, cutoff = cutoff, breaks = breaks,
                   col = color.edge, symkey = symkey)
    col.id = id$bin
    color.edge = id$col[col.id]
    
    # selection of the edges to incluir in the network #
    #--------------------------------------------------#
    idx = (abs(w) >= cutoff)
    relations = relations[idx, ]
    color.edge = color.edge[idx]
    
    # Generation of the graph with all the significant edges #
    #--------------------------------------------------------#
    gR = graph_from_data_frame(relations, directed = FALSE, vertices = nodes)
    
    # nodes attributes #
    #------------------#
    
    
    V(gR)$label.color = "black"
    
    V(gR)$label.family = "sans"
    
    if(any(class.object %in% object.blocks))
    {
        V(gR)$label = unlist(block.var.names)
        
        j = 1
        for (i in blocks)
        {
            V(gR)$color[V(gR)$group == i] = color.node[j]
            V(gR)$shape[V(gR)$group == i] = shape.node[j]
            if (shape.node[j] == "none") {
              V(gR)$label.color[V(gR)$group == i] = paste0(substr(color.node[j], 1, 7), "FF")
            }
            j = j + 1
        }
    } else {
        V(gR)$label = c(row.names.plot, col.names.plot)
        V(gR)$color = color.node[1]
        V(gR)$color[V(gR)$group == "y"] = color.node[2]
        
        V(gR)$shape = shape.node[1]
        if (shape.node[1] == "none") {
          V(gR)$label.color = paste0(substr(color.node[1], 1, 7), "FF")
        }
        V(gR)$shape[V(gR)$group == "y"] = shape.node[2]
        if (shape.node[2] == "none") {
          V(gR)$label.color[V(gR)$group == "y"] = paste0(substr(color.node[2], 1, 7), "FF")
        }
        
    }
    
    # edges attributes #
    #------------------#
    if (show.edge.labels)
        E(gR)$label = round(E(gR)$weight, 2)
    
    E(gR)$label.color = "black"
    
    E(gR)$color = color.edge
    
    #E(gR)$lty = lty.edge
    
    #E(gR)$width = lwd.edge
    #from 5.0.x: allow different edge/lwd for >0 and <0 vertices
    E(gR)$lty = lty.edge[1]
    E(gR)$lty[E(gR)$weight < 0] = lty.edge[2]
    
    E(gR)$width = lwd.edge[1]
    E(gR)$width[E(gR)$weight < 0] = lwd.edge[2]
    
    
    gR = delete_vertices(gR, which(degree(gR) == 0))
    
    # plot attributes #
    #-----------------#
    lwid = c(keysize[1], 4)
    lhei = c(keysize[2], 4)
    lmat = matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
    nc = length(id$col)
    x = seq(0, 1, length = nc + 2)
    z.mat = seq(0, 1, length = nc + 1)
    z.mat = matrix(z.mat, ncol = 1)
    
    if ((id$lim[1] < -cutoff) & (id$lim[2] < cutoff))
    {
        xv = c(0, x[nc + 1])
        lv = round(c(id$lim[1], -cutoff), 2)
        col = c(id$col, "white")
    }
    
    if ((id$lim[1] > -cutoff) & (id$lim[2] > cutoff))
    {
        xv = c(x[2], 1)
        lv = round(c(cutoff, id$lim[2]), 2)
        col = c("white", id$col)
    }
    
    if ((id$lim[1] < -cutoff) & (id$lim[2] > cutoff))
    {
        idn = max(which(id$breaks < 0))
        idp = min(which(id$breaks > 0))
        xv = c(0, x[idn + 1], x[idp], 1)
        lv = round(c(id$lim[1], -cutoff, cutoff, id$lim[2]), 2)
        col = c(id$col[1:idn], "white", id$col[(idn + 1):nc])
    }
    
    #-----------------------------------#
    # construction of the initial graph #
    #-----------------------------------#

    nn = vcount(gR)
    V(gR)$label.cex = ifelse(is.null(cex.node.name), 1, cex.node.name)
    E(gR)$label.cex = min(2.25 * cex.edge.label/log(nn), 1)
    cex0 = 2 * V(gR)$label.cex
    
    def.par = par(no.readonly = TRUE)
    dev.new()
    par(pty = "s", mar = c(0, 0, 0, 0),mfrow=c(1,1))
    upr.bnd <- 190*(graph.scale^2) + 10
    plot(1:upr.bnd, 1:upr.bnd, type = "n", axes = FALSE, xlab = "", ylab = "")
    cha = V(gR)$label
    cha = paste("", cha, "")
    xh = strwidth(cha, cex = ifelse(size.node==0.5, cex0, size.node*4)) * 1.5 * (2-graph.scale)
    yh = strheight(cha, cex = ifelse(size.node==0.5, cex0, size.node*4)) * 3 * (2-graph.scale)
  
    
    V(gR)$size = xh
    V(gR)$size2 = yh
    
    dev.off()
    
    if (is.null(layout.fun))
    {
        l = layout.fruchterman.reingold(gR, weights = (1 - abs(E(gR)$weight)))
    } else {
        l = layout.fun(gR)
    }
    
    if (isTRUE(!interactive))
    {
        if (isTRUE(show.color.key))
        {
            l = layout.fruchterman.reingold(gR, weights = (1 - abs(E(gR)$weight)))
        } else {
            l = layout.fun(gR)
        }
        
        if (isTRUE(!interactive))
        {
            if (isTRUE(show.color.key))
            {
                layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
                par(mar = c(5, 4, 2, 1), cex = 0.75)
                image(z.mat, col = col, xaxt = "n", yaxt = "n")
                box()
                par(usr = c(0, 1, 0, 1))
                axis(1, at = xv, labels = lv, cex.axis = keysize.label)
                title("Color key", font.main = 1, cex.main = keysize.label)
                par(def.par)
                par(new = TRUE)
            }
            
            par(pty = "s", mar = c(0, 0, 0, 0),mfrow=c(1,1))
            plot(gR, layout = l)
            par(def.par)
        }
    }
    
    #-----------------------#
    # procedure interactive #
    #-----------------------#
    gE.none = FALSE
    if (isTRUE(interactive))
    {
        
        # cutoff control bar #
        #-----------------------#
        min.cut = cutoff
        max.cut = max(abs(w))
        
        cutoff.old = cutoff
        
        dev.new("width" = 5, "height" = 2.7, xpos = -1)
        def.par = par(no.readonly = TRUE)
        
        cuts = seq(0, 1, length = 21)
        par(mai = c(0.25, 0.15, 0.3, 0.15), bg = gray(0.95))
        layout(matrix(c(0, 1, 0), ncol = 1, nrow = 3),
               widths = 1, heights = c(0.25, 1, 0.25))
        
        plot(cuts, type = "n", rep(0, 21), xlab = "", ylab = "",
             xlim = c(-0.10, 1.10), axes = FALSE)
        title("cutoff control", cex.main = 1.9, font.main = 1)
        text(0.5, -0.6, "value", cex = 1.5)
        text(0, -0.6, round(min.cut, 2), cex = 1.4)
        text(1, -0.6, round(max.cut, 2), cex = 1.4)
        mtext(min.cut, side = 1, line = -1, outer = FALSE, cex = 0.95)
        
        rect(-0.1, -0.3, -0.02, 0.3, col = "white")
        rect(1.02, -0.3, 1.1, 0.3, col = "white")
        points(1.06, 0, pch = 3, cex = 2.4)
        lines(c(-0.085, -0.035), c(0, 0))
        
        for (i in seq(0, 1, length = 21))
            lines(c(i, i), c(-0.22, 0.2))
        
        x = pos = 0
        rect(-0.01, -0.045, x, 0.04, col = "red")
        rect(x, -0.045, 1.01, 0.04, col = "white")
        
        bar.dev = dev.cur()
        
        # construction of the initial graph #
        #-----------------------------------#
        dev.new()
        net.dev = dev.cur()
        
        if (isTRUE(show.color.key))
        {
            layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
            par(mar = c(5, 4, 2, 1), cex = 0.75)
            image(z.mat, col = col, xaxt = "n", yaxt = "n")
            box()
            par(usr = c(0, 1, 0, 1))
            axis(1, at = xv, labels = lv, cex.axis = keysize.label)
            title("Color key", font.main = 1, cex.main = keysize.label)
            par(def.par)
            par(new = TRUE)
        }
        
        par(pty = "s", mar = c(0, 0, 0, 0))
        plot(gR, layout = l)
        par(def.par)
        
        old.pos = -1
        
        repeat {
            dev.set(bar.dev)
            
            z = locator(1, type = "n")
            x = z[[1]]
            y = z[[2]]
            
            if (is.null(z)) break
            
            if (0 <= x & x <= 1 & -0.22 <= y & y <= 0.22)
            {
                rect(0, -0.045, x, 0.04, col = "red")
                rect(x, -0.045, 1.01, 0.04, col = "white")
                pos = x
            }
            
            if (1.02 <= x & x <= 1.1 & -0.3 <= y & y <= 0.3)
            {
                x = pos + 0.05
                idx = which.min(abs(cuts - x))
                x = cuts[idx]
                pos = x
                rect(0, -0.045, x, 0.04, col = "red")
                rect(x, -0.045, 1.01, 0.04, col = "white")
            }
            
            if (-0.1 <= x & x <= -0.02 & -0.3 <= y & y <= 0.3)
            {
                x = pos - 0.05
                idx = which.min(abs(cuts - x))
                x = cuts[idx]
                pos = x
                rect(0, -0.045, x, 0.04, col = "red")
                rect(x, -0.045, 1.01, 0.04, col = "white")
            }
            
            if (old.pos != pos)
            {
                old.pos = pos
                rect(0.4, -0.8, 0.6, -1.5, col = gray(0.95), border = NA)
                cutoff = (max.cut - min.cut) * pos + min.cut
                mtext(round(cutoff, 3), side = 1, line = -1, cex = 0.9)
                
                
                # new graph plot #
                #----------------#
                dev.set(net.dev)
                
                if (cutoff >= cutoff.old)
                {
                    
                    # selection of the edges to remove of the network #
                    #-------------------------------------------------#
                    supp.edge = E(gR)[abs(E(gR)$weight) < cutoff]
                    
                    # Generation of the graph with all the significant edges #
                    #--------------------------------------------------------#
                    gE = delete.edges(gR, supp.edge)
                    gE = delete.vertices(gE, which(degree(gE) == 0))
                    
                    # graph plot #
                    #------------#
                    nn = vcount(gE)
                    V(gR)$label.cex = min(2.5 * cex.node.name/log(nn), 1)
                    E(gR)$label.cex = min(2.25 * cex.edge.label/log(nn), 1)
                    cex0 = 2 * V(gE)$label.cex
                    
                    def.par = par(no.readonly = TRUE)
                    
                    par(pty = "s", mar = c(0, 0, 0, 0))
                    plot(1:100, 1:100, type = "n", xaxt = "n")
                    cha = V(gE)$label
                    cha = paste("", cha, "")
                    xh = strwidth(cha, cex = cex0) * 1.5
                    yh = strheight(cha, cex = cex0) * 3
                    
                    V(gE)$size = xh
                    V(gE)$size2 = yh
                    
                    par(def.par)
                    
                    if (is.null(layout.fun))
                    {
                        l = layout.fruchterman.reingold(gE, weights = (1 - abs(E(gE)$weight)))
                    } else {
                        l = layout.fun(gE)
                    }
                    
                    if (isTRUE(show.color.key))
                    {
                        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
                        par(mar = c(5, 4, 2, 1), cex = 0.75)
                        image(z.mat, col = col, xaxt = "n", yaxt = "n")
                        box()
                        par(usr = c(0, 1, 0, 1))
                        axis(1, at = xv, labels = lv, cex.axis = keysize.label)
                        title("Color key", font.main = 1, cex.main = keysize.label)
                        par(def.par)
                        par(new = TRUE)
                    }
                    
                    par(pty = "s", mar = c(0, 0, 0, 0))
                    plot(gE, layout = l)
                    par(def.par)
                    
                    cutoff.old = cutoff
                } else {
                    # selection of the edges to incluir in the network #
                    #--------------------------------------------------#
                    supp.edge = E(gR)[abs(E(gR)$weight) < cutoff]
                    
                    # generation of the graph with all the significant edges #
                    #--------------------------------------------------------#
                    gE = delete.edges(gR, supp.edge)
                    gE = delete.vertices(gE, which(degree(gE) == 0))
                    
                    # graph plot #
                    #------------#
                    nn = vcount(gE)
                    V(gR)$label.cex = min(2.5 * cex.node.name/log(nn), 1)
                    E(gR)$label.cex = min(2.25 * cex.edge.label/log(nn), 1)
                    cex0 = 2 * V(gE)$label.cex
                    
                    def.par = par(no.readonly = TRUE)
                    
                    par(pty = "s", mar = c(0, 0, 0, 0))
                    plot(1:100, 1:100, type = "n", xaxt = "n")
                    cha = V(gE)$label
                    cha = paste("", cha, "")
                    xh = strwidth(cha, cex = cex0) * 1.5
                    yh = strheight(cha, cex = cex0) * 3
                    
                    V(gE)$size = xh
                    V(gE)$size2 = yh
                    
                    par(def.par)	
                    
                    if (is.null(layout.fun))
                    {
                        l = layout.fruchterman.reingold(gE, weights = (1 - abs(E(gE)$weight)))
                    } else {
                        l = layout.fun(gE)
                    }
                    
                    if (isTRUE(show.color.key))
                    {
                        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
                        par(mar = c(5, 4, 2, 1), cex = 0.75)
                        image(z.mat, col = col, xaxt = "n", yaxt = "n")
                        box()
                        par(usr = c(0, 1, 0, 1))						
                        axis(1, at = xv, labels = lv, cex.axis = keysize.label)
                        title("Color key", font.main = 1, cex.main = keysize.label)
                        par(def.par)
                        par(new = TRUE)
                    }
                    
                    par(pty = "s", mar = c(0, 0, 0, 0))
                    plot(gE, layout = l)
                    par(def.par)
                    
                    cutoff.old = cutoff
                }
                
                gE.none = TRUE
            }
            
        } # end loop
        
        if (gE.none != FALSE)
            gR = gE
    }
    res=list(gR = gR)
    
    
    if(any(class.object %in% object.blocks))
    {
        l = 1
        for (i in 1:(length(blocks)-1))
        {
            for (j in (i + 1):length(blocks))
            {
                M_block[[l]][abs(M_block[[l]]) < cutoff] = 0
                res[paste("M",blocks[i],blocks[j],sep="_")] = list(M_block[[l]])
                l = l + 1
            }
        }
    } else {
        mat[abs(mat) < cutoff] = 0
        res$M=mat
    }
    
    res$cutoff = cutoff
    
    if (!is.null(save))
        dev.off()
    
    return(invisible(res))
}

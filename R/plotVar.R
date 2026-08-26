#----------------------------------------------------------------------------------------------------------#
#-- Includes plotVar for PLS, sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#
#' Plot of Variables
#' 
#' This function provides variables representation for (regularized) CCA,
#' (sparse) PLS regression, PCA and (sparse) Regularized generalised CCA.
#' 
#' \code{plotVar} produce a "correlation circle", i.e. the correlations between
#' each variable and the selected components are plotted as scatter plot, with
#' concentric circles of radius one et radius given by \code{rad.in}. Each
#' point corresponds to a variable. For (regularized) CCA the components
#' correspond to the equiangular vector between \eqn{X}- and \eqn{Y}-variates.
#' For (sparse) PLS regression mode the components correspond to the
#' \eqn{X}-variates. If mode is canonical, the components for \eqn{X} and
#' \eqn{Y} variables correspond to the \eqn{X}- and \eqn{Y}-variates
#' respectively.
#' 
#' For \code{plsda} and \code{splsda} objects, only the \eqn{X} variables are
#' represented.
#' 
#' For \code{spls} and \code{splsda} objects, only the \eqn{X} and \eqn{Y}
#' variables selected on dimensions \code{comp} are represented.
#' 
#' The arguments \code{col}, \code{pch}, \code{cex} and \code{font} can be
#' either vectors of length two or a list with two vector components of length
#' \eqn{p} and \eqn{q} respectively, where \eqn{p} is the number of
#' \eqn{X}-variables and \eqn{q} is the number of \eqn{Y}-variables. In the
#' first case, the first and second component of the vector determine the
#' graphics attributes for the \eqn{X}- and \eqn{Y}-variables respectively.
#' Otherwise, multiple arguments values can be specified so that each point
#' (variable) can be given its own graphic attributes. In this case, the first
#' component of the list correspond to the \eqn{X} attributs and the second
#' component correspond to the \eqn{Y} attributs. Default values exist for this
#' arguments.
#' 
#' @aliases plotVar plotVar.rcc plotVar.pls plotVar.spls plotVar.plsda
#' plotVar.splsda plotVar.pca plotVar.spca plotVar.sgcca plotVar.rgcca
#' @param object object of class inheriting from \code{"rcc"}, \code{"pls"},
#' \code{"plsda"}, \code{"spls"}, \code{"splsda"}, \code{"pca"} or
#' \code{"spca"}.
#' @param comp integer vector of length two. The components that will be used
#' on the horizontal and the vertical axis respectively to project the
#' variables. By default, comp=c(1,2) except when style='3d', comp=c(1:3)
#' @param comp.select for the sparse versions, an input vector indicating the
#' components on which the variables were selected. Only those selected
#' variables are displayed. By default, comp.select=comp
#' @param plot if TRUE (the default) then a plot is produced. If not, the
#' summaries which the plots are based on are returned.
#' @param var.names either a character vector of names for the variables to be
#' plotted, or \code{FALSE} for no names. If \code{TRUE}, the col names of the
#' first (or second) data matrix is used as names.
#' @param blocks for an object of class \code{"rgcca"} or \code{"sgcca"}, a
#' numerical vector indicating the block variables to display.
#' @param X.label x axis titles.
#' @param Y.label y axis titles.
#' @param Z.label z axis titles (when style = '3d').
#' @param abline should the vertical and horizontal line through the center be
#' plotted? Default set to \code{FALSE}
#' @param col character or integer vector of colors for plotted character and
#' symbols, can be of length 2 (one for each data set) or of length (p+q) (i.e.
#' the total number of variables). See Details.
#' @param cex numeric vector of character expansion sizes for the plotted
#' character and symbols, can be of length 2 (one for each data set) or of
#' length (p+q) (i.e. the total number of variables).
#' @param pch plot character. A vector of single characters or integers, can be
#' of length 2 (one for each data set) or of length (p+q) (i.e. the total
#' number of variables). See \code{\link{points}} for all alternatives.
#' @param font numeric vector of font to be used, can be of length 2 (one for
#' each data set) or of length (p+q) (i.e. the total number of variables). See
#' \code{\link{par}} for details.
#' @param cutoff numeric between 0 and 1. Variables with correlations below
#' this cutoff in absolute value are not plotted (see Details).
#' @param rad.in numeric between 0 and 1, the radius of the inner circle.
#' Defaults to \code{0.5}.
#' @param title character indicating the title plot.
#' @param legend Logical when more than 3 blocks. Can be a character vector
#' when one or 2 blocks to customize the legend. See examples. Default is
#' FALSE.
#' @param legend.title title of the legend
#' @param style argument to be set to either \code{'graphics'},
#' \code{'lattice'}, \code{'ggplot2'} or \code{'3d'} for a style of plotting.
#' @param overlap Logical. Whether the variables should be plotted in one
#' single figure. Default is TRUE.
#' @param axes.box for style '3d', argument to be set to either \code{'axes'},
#' \code{'box'}, \code{'bbox'} or \code{'all'}, defining the shape of the box.
#' @param label.axes.box for style '3d', argument to be set to either
#' \code{'axes'}, \code{'box'}, \code{'both'}, indicating which labels to
#' print.
#' @return A list containing the following components: \item{x}{a vector of
#' coordinates of the variables on the x-axis.} \item{y}{a vector of
#' coordinates of the variables on the y-axis.} \item{Block}{the data block
#' name each variable belongs to.} \item{names}{the name of each variable,
#' matching their coordinates values.}
#' @author Ignacio González, Benoit Gautier, Francois Bartolo, Florian Rohart,
#' Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{cim}}, \code{\link{network}}, \code{\link{par}} and
#' http://www.mixOmics.org for more details.
#' @references González I., Lê Cao K-A., Davis, M.J. and Déjean, S. (2012).
#' Visualising associations between paired 'omics data sets. J. Data Mining
#' 5:19. \url{http://www.biodatamining.org/content/5/1/19/abstract}
#' @keywords multivariate hplot dplot
#' @export
#' @importFrom ellipse ellipse
#' @example ./examples/plotVar-examples.R
plotVar <-
    function(object,
             comp = NULL,
             comp.select = comp,
             plot = TRUE,
             var.names = NULL,
             blocks = NULL,
             # to choose which block data to plot, when using GCCA module
             X.label = NULL,
             Y.label = NULL,
             Z.label = NULL,
             abline = TRUE,
             col,
             cex,
             pch,
             font,
             cutoff = 0,
             rad.in = 0.5,
             title = "Correlation Circle Plot",
             legend = FALSE,
             legend.title = "Block",
             style = "ggplot2",
             # can choose between graphics,3d, lattice or ggplot2,
             overlap = TRUE,
             axes.box = "all",
             label.axes.box = "both")
    {
        class.object = class(object)
        object.pls=c("mixo_pls","mixo_spls","mixo_mlspls","mixo_mlsplsda","rcc")
        object.pca=c("ipca","sipca","pca","spca")
        object.blocks=c("sgcca","rgcca")
        
        #-- check that the user did not enter extra arguments
        arg.call = match.call()
        user.arg = names(arg.call)[-1]
        
        err = tryCatch(mget(names(formals()), sys.frame(sys.nframe())),
                       error = function(e) e)
        
        if ("simpleError" %in% class(err))
            stop(err[[1]], ".", call. = FALSE)
        
        #-- style
        if (!style %in% c("ggplot2", "lattice", "graphics","3d"))
            stop("'style' must be one of 'ggplot2', '3d' , lattice' or 'graphics'.", call. = FALSE)
        
        #-- plot
        if (length(plot) > 1)
            stop("'plot' must be single logical value.", call. = FALSE)
        else if (!is.logical(plot))
            stop("'plot' must be logical.", call. = FALSE)
        if(!plot)
        {
            style="N"}
        
        #-- axes.box
        if(style=="3d")
        {
            choices = c("axes", "box", "bbox", "all")
            axes.box = choices[pmatch(axes.box, choices)]
            
            if (is.na(axes.box))
                stop("'axes.box' should be a subset of {'axes', 'box', 'bbox', 'all'}.",
                     call. = FALSE)
            
            #-- label.axes.box
            choices = c("axes", "box", "both")
            label.axes.box = choices[pmatch(label.axes.box, choices)]
            
            if (is.na(label.axes.box))
                stop("'label.axes.box' should be one of 'axes', 'box' or 'both'.",
                     call. = FALSE)}
        
        
        
        ### Start: Validation of arguments
        ncomp = object$ncomp
        if (any(class.object %in% object.blocks))
        {
            #-- legend: logical only
            if (length(legend) != 1 || !is.logical(legend))
                stop("'legend' must be a logical value.", call. = FALSE)
            
            if (is.null(blocks))
            {
                blocks = names(object$X)#names$blocks
                
                if (any(class.object == "DA"))
                    blocks = names(object$X)#blocks[-object$indY]
                
            } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <= length(object$names$blocks)) {
                blocks = object$names$blocks[blocks]
            } else if (is.character(blocks)) {
                if (!any(blocks %in% object$names$blocks))
                    stop("One element of 'blocks' does not match with the names of the blocks")
            } else {
                stop("Incorrect value for 'blocks'", call. = FALSE)
            }
            object$variates = object$variates[names(object$variates) %in% blocks]
            object$names$colnames = object$names$colnames[names(object$names$colnames) %in% blocks]
            object$blocks = object$X[names(object$X) %in% blocks]
            
            if (any(object$ncomp[blocks] == 1))
            {
                stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "),"' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
            }
            ncomp = object$ncomp[blocks]
        } else if (any(class.object %in% c("rcc", "mixo_pls", "mixo_spls", "mixo_mlspls")) & all(class.object !="DA")) {
            #-- legend: logical or a name for X and Y
            if ( ! (length(legend) == 1 & is.logical(legend) || (length(legend)==2)))
                stop("'legend' must be a logical value or a vector of 2 names for X and Y.", call. = FALSE)
            
            if(length(legend)==2){
                blocks=legend
                legend=TRUE
            } else{
                blocks = c("X","Y")
            }
            
        } else {
            #-- legend: logical or a name for X
            if (length(legend) != 1 )
                stop("'legend' must be a logical value or a vector of 1 name for X.", call. = FALSE)
            
            if(is.logical(legend)){
                blocks = "X"
            }else{
                blocks = legend
                legend=TRUE
            }
        }
        #-- legend.title
        if (length(legend.title)>1)
            stop("'legend.title' needs to be a single value (length 1)")
        
        #-- ellipse.level
        if (!is.numeric(rad.in) | (rad.in > 1) | (rad.in < 0))
            stop("The value taken by 'rad.in' must be between 0 and 1", call. = FALSE)
        
        #-- cutoff correlation
        if (!is.numeric(cutoff) | (cutoff > 1) | (cutoff < 0))
            stop("The value taken by 'cutoff' must be between 0 and 1", call. = FALSE)
        
        #-- comp
        if(is.null(comp))
        {
            if (style=="3d")
            {
                comp = seq_len(3)
            } else {
                comp = seq_len(2)
            }
        }
        if (length(comp) != 2 && !(style=="3d"))
        {
            stop("'comp' must be a numeric vector of length 2.", call. = FALSE)
        } else if(length(comp) != 3 && (style=="3d")) {
            stop("'comp' must be a numeric vector of length 3.", call. = FALSE)
        }
        
        if (!is.numeric(comp))
            stop("Invalid vector for 'comp'.")
        
        if (any(ncomp < max(comp)) || min(comp) <= 0)
            stop("Each element of 'comp' must be positive <= ", min(object$ncomp), ".", call. = FALSE)
        
        comp1 = round(comp[1])
        comp2 = round(comp[2])
        if (style=="3d")
            comp3 = round(comp[3])
        
        #-- comp.select
        if (!is.null(comp.select))
        {
            if (!is.numeric(comp.select))
                stop("Invalid vector for 'comp'.", call. = FALSE)
            
            if (any(ncomp < max(comp.select)) || min(comp.select) <= 0)
                stop("Each element of 'comp.select' must be positive and <= ", max(object$ncomp), ".", call. = FALSE)
        } else {
            comp.select = comp
        }
        
        #-- abline
        if (length(abline) > 1)
        {
            stop("'abline' must be single logical value.", call. = FALSE)
        }else if (!is.logical(abline)) {
            stop("'abline' must be logical.", call. = FALSE)
        }
        
        #-- Start: Retrieve variates from object
        cord.X = sample.X = ind.var.sel = list()
        if(style=="3d")
        {
            if (any(class.object%in%  c(object.pls, object.blocks)))
            {
                if (any(class.object == "rcc"))
                {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2, comp3)] + object$variates$Y[, c(comp1, comp2, comp3)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y, object$variates$X[, c(comp1, comp2, comp3)] + object$variates$Y[, c(comp1, comp2, comp3)], use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                    
                } else if (any(class.object %in% "mixo_plsda")) {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2, comp3)], use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                    
                } else if (any(class.object %in%  "mixo_pls")) {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2, comp3)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y, if(object$mode ==  "canonical"){object$variates$Y[, c(comp1, comp2, comp3)]} else {object$variates$X[, c(comp1, comp2, comp3)]}, use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                    
                } else if (any(class.object %in%  c("mixo_splsda", "mixo_mlsplsda"))) {
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(unique(c(comp1, comp2, comp3, comp.select)), function(x){selectVar(object, comp = x)$name})))], # variables selected at least once on unique(comp1, comp2, comp3 and comp.select
                                      object$variates$X[, c(comp1, comp2, comp3, comp.select)], use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                    if (!is.null(comp.select))
                    {
                        cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$name}))), ,drop = FALSE]
                    }
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                    
                } else if (any(class.object %in%  c("mixo_spls", "mixo_mlspls"))) {
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(c(comp1, comp2, comp3), function(x){selectVar(object, comp = x)$X$name})))],
                                      object$variates$X[, c(comp1, comp2, comp3)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y[, colnames(object$Y) %in% unique(unlist(lapply(c(comp1, comp2, comp3), function(x){selectVar(object, comp = x)$Y$name})))],
                                      if(object$mode ==  "canonical"){object$variates$Y[, c(comp1, comp2, comp3)]} else {object$variates$X[, c(comp1, comp2, comp3)]}, use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                    ind.var.sel[[2]] = sample.X[[2]] = seq_len(length(colnames(object$Y)))
                    if (!is.null(comp.select)) {
                        cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$X$name}))), ,drop = FALSE]
                        cord.X[[2]] = cord.X[[2]][row.names(cord.X[[2]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$Y$name}))), , drop = FALSE]
                    }
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                    ind.var.sel[[2]] = which(colnames(object$Y) %in% rownames(cord.X[[2]]))
                } else {
                    cord.X = lapply(blocks, function(x){cor(object$blocks[[x]], object$variates[[x]][, c(comp1, comp2, comp3)], use = "pairwise")})
                    ind.var.sel = sample.X = lapply(object$blocks, function(x){seq_len(ncol(x))})
                    if (!is.null(comp.select)) {
                        cord.X = lapply(seq_len(length(cord.X)), function(z){cord.X[[z]][row.names(cord.X[[z]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, block = z, comp = x)[[1]]$name}))), ,drop = FALSE]})
                    }
                    for (i in seq_len(length(cord.X))){
                        ind.var.sel[[i]] = which(colnames(object$X) %in% rownames(cord.X[[i]]))
                    }
                }
            } else if (any(class.object %in%  object.pca)) {
                if (any(class.object %in%  c("sipca", "spca"))){
                    
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(c(comp1, comp2, comp3), function(x){selectVar(object, comp = x)$name})))],
                                      object$variates$X[, c(comp1, comp2, comp3)], use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                    if (!is.null(comp.select)) {
                        cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$name}))), ,drop = FALSE]
                    }
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                } else {
                    
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2, comp3)], use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                }
            }
        } else {
            if (any(class.object %in%  c(object.pls, object.blocks)))
            {
                if (any(class.object == "rcc"))
                {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2)] + object$variates$Y[, c(comp1, comp2)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y, object$variates$X[, c(comp1, comp2)] + object$variates$Y[, c(comp1, comp2)], use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                    
                } else if (any(class.object == "mixo_plsda")) {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                    
                } else if (any(class.object ==  "mixo_pls")) {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y, if(object$mode ==  "canonical"){object$variates$Y[, c(comp1, comp2)]} else {object$variates$X[, c(comp1, comp2)]}, use = "pairwise")
                    sample.X = lapply(cord.X, function(x){seq_len(nrow(x))})
                    
                } else if (any(class.object %in%  c("mixo_splsda", "mixo_mlsplsda"))) {
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$name}))), drop = FALSE],
                                      object$variates$X[, unique(c(comp1, comp2))], use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                    #if (!is.null(comp.select)) {
                    #   cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$name}))), ,drop = FALSE]
                    #}
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                    
                } else if (any(class.object %in%  c("mixo_spls", "mixo_mlspls"))) {
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$X$name}))), drop = FALSE],
                                      object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    cord.X[[2]] = cor(object$Y[, colnames(object$Y) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$Y$name}))), drop = FALSE],
                                      if(object$mode ==  "canonical")
                                      {
                                          object$variates$Y[, c(comp1, comp2)]
                                      } else {
                                          object$variates$X[, c(comp1, comp2)]
                                      }, use = "pairwise")
                    #ind.var.sel[[1]] =
                    sample.X[[1]] = seq_len(length(colnames(object$X)))
                    #ind.var.sel[[2]] =
                    sample.X[[2]] = seq_len(length(colnames(object$Y)))
                    #if (!is.null(comp.select)) {
                    #   cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$X$name}))), ,drop = FALSE]
                    #   cord.X[[2]] = cord.X[[2]][row.names(cord.X[[2]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$Y$name}))), , drop = FALSE]
                    #}
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                    ind.var.sel[[2]] = which(colnames(object$Y) %in% rownames(cord.X[[2]]))
                    
                } else { #block object
                    cord.X = lapply(blocks, function(x){cor(object$blocks[[x]], object$variates[[x]][, c(comp1, comp2)], use = "pairwise")})
                    ind.var.sel = sample.X = lapply(object$blocks, function(x){seq_len(ncol(x))})
                    if (!is.null(comp.select))
                    {
                        cord.X = lapply(seq_len(length(cord.X)), function(z){cord.X[[z]][row.names(cord.X[[z]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, block = blocks[z], comp = x)[[1]]$name}))), ,drop = FALSE]})
                    }
                    for (i in seq_len(length(cord.X)))
                    {
                        ind.var.sel[[i]] = which(colnames(object$blocks[[i]]) %in% rownames(cord.X[[i]]))
                    }
                }
            } else if (any(class.object %in%  object.pca)) {
                if (any(class.object %in%  c("sipca", "spca"))){
                    
                    cord.X[[1]] = cor(object$X[, colnames(object$X) %in% unique(unlist(lapply(comp.select, function(x){selectVar(object, comp = x)$name}))), drop = FALSE],
                                      object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    #ind.var.sel[[1]] =
                    sample.X[[1]] = seq_len(length(colnames(object$X)))
                    #if (!is.null(comp.select)) {
                    #    cord.X[[1]] = cord.X[[1]][row.names(cord.X[[1]]) %in% unique(unlist(lapply(comp.select, function(x) {selectVar(object, comp = x)$name}))), ,drop = FALSE]
                    #}
                    ind.var.sel[[1]] = which(colnames(object$X) %in% rownames(cord.X[[1]]))
                } else {
                    cord.X[[1]] = cor(object$X, object$variates$X[, c(comp1, comp2)], use = "pairwise")
                    ind.var.sel[[1]] = sample.X[[1]] = seq_len(length(colnames(object$X)))
                }
            }}
        
        # output a message if some variates are anti correlated among blocks
        if (any(class.object %in%  object.blocks))
        {
            VarX = lapply(seq_len(2), function(j){do.call(cbind, lapply(object$variates, function(i) i[, comp[j]]))})
            corX = lapply(VarX, cor)
            if(any(sapply(corX, function(j){any(j < 0)})))
                warning("We detected negative correlation between the variates of some blocks, which means that some clusters of variables observed on the correlation circle plot are not necessarily positively correlated.")
        }
        
        if (any(sapply(cord.X, nrow) == 0))
            stop("No variable selected on at least one block")
        
        if (any(sapply(cord.X, nrow) == 1)) {
          cord.X <- cord.X[-which(sapply(cord.X, nrow) == 1)]
        }
          
        
        #-- End: Retrieve variates from object
        
        #-- Names of labels X and Y
        if (is.null(X.label)) X.label = paste("Component ", comp1)
        if (is.null(Y.label)) Y.label = paste("Component ", comp2)
        if (is.null(Z.label) && style=="3d") Z.label = paste("Component ", comp3)
        
        if (!is.character(X.label))
            stop("'X.label' must be a character.", call. = FALSE)
        if (!is.character(Y.label))
            stop("'Y.label' must be a character.", call. = FALSE)
        
        
        #-- pch argument
        missing.pch = FALSE
        if (missing(pch))
        {
            missing.pch = TRUE
            if(style=="3d")
            {
                pch = unlist(lapply(seq_len(length(cord.X)), function(x){rep(c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")[x], sum(sapply(cord.X[x], nrow)))}))
            } else {
                pch = unlist(lapply(seq_len(length(cord.X)), function(x){rep(seq_len(20)[x], sum(sapply(cord.X[x], nrow)))}))
            }
            
        } else if (((is.vector(pch, mode = "double") || is.vector(pch, mode = "integer")) && !(style=="3d"))
                   || (is.vector(pch, mode = "character") && style=="3d")) {
            if (length(pch) != length(sample.X))
                stop.message('pch', sample.X)
            pch = unlist(lapply(seq_len(length(cord.X)), function(x){rep(pch[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(pch)) {
            if (length(pch) != length(sample.X) || length(unlist(pch)) != sum(sapply(sample.X, length)))
                stop.message('pch', sample.X)
            if (length(ind.var.sel) != 0)
                pch = lapply(seq_len(length(pch)), function(x){pch[[x]][ind.var.sel[[x]]]})
            pch = unlist(pch)
        } else if (style=="3d") {
            if (!all(pch %in% c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")) && style=="3d")
                stop("pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.",
                     call. = FALSE)
        }
        else {
            stop.message('pch', sample.X)
        }
        
        #-- col argument
        if (missing(col)) {
            if (length(cord.X) < 10) {
                col = unlist(lapply(seq_len(length(cord.X)), function(x){rep(color.mixo(x), sum(sapply(cord.X[x], nrow)))}))
            } else {
                col = unlist(lapply(seq_len(length(cord.X)), function(x){rep(color.jet(length(cord.X))[x], sum(sapply(cord.X[x], nrow)))}))
            }
        } else if (is.vector(col, mode = "double") | is.vector(col, mode = "character")) {
            if (length(col) != length(sample.X))
                stop.message('col', sample.X)
            col = unlist(lapply(seq_len(length(cord.X)), function(x){rep(col[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(col)) {
            if (length(col) != length(sample.X) || length(unlist(col)) != sum(sapply(sample.X, length)))
                stop.message('col', sample.X)
            if (length(ind.var.sel) != 0)
                col = lapply(seq_len(length(col)), function(x){col[[x]][ind.var.sel[[x]]]})
            col = unlist(col)
        } else {
            stop.message('col', sample.X)
        }
        
        #-- cex argument
        if (missing(cex)){
            if (style == "ggplot2"){
                cex = rep(5, sum(sapply(cord.X, nrow)))
            } else {
                cex = rep(1, sum(sapply(cord.X, nrow)))
            }
        } else if (is.vector(cex, mode = "double")) {
            if (length(cex) != length(cord.X))
                stop.message('cex', sample.X)
            cex = unlist(lapply(seq_len(length(cord.X)), function(x){rep(cex[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(cex)) {
            if (length(cex) != length(sample.X) || length(unlist(cex)) != sum(sapply(sample.X, length)))
                stop.message('cex', sample.X)
            if (length(ind.var.sel) != 0)
                cex = lapply(seq_len(length(cex)), function(x){cex[[x]][ind.var.sel[[x]]]})
            cex = unlist(cex)
        } else {
            stop.message('cex', sample.X)
        }
        
        #-- font argument
        if (missing(font)) {
            font = rep(1, sum(sapply(cord.X, nrow)))
        } else if (is.vector(font, mode = "numeric")) {
            if (length(font) != length(cord.X))
                stop.message('font', sample.X)
            font = unlist(lapply(seq_len(length(cord.X)), function(x){rep(font[x], sum(sapply(cord.X[x], nrow)))}))
        } else if (is.list(font)) {
            if (length(font) != length(sample.X) || length(unlist(font)) != sum(sapply(sample.X, length)))
                stop.message('font', sample.X)
            if (length(ind.var.sel) != 0)
                font = lapply(seq_len(length(font)), function(x){font[[x]][ind.var.sel[[x]]]})
            font = unlist(font)
        } else {
            stop.message('font', sample.X)
        }
        
        #-- var.names
        ind.group = cumsum(c(0, sapply(cord.X, nrow)))
        if (is.null(var.names)){
            var.names.list = unlist(sapply(cord.X, rownames))
            if (!missing.pch) {
                var.names = rep(FALSE, length(cord.X))
            } else {
                var.names = rep(TRUE, length(cord.X))
            }
        } else if (is.vector(var.names, mode = "logical")) {
            if (length(var.names) == 1){
                var.names = rep(var.names,length(cord.X))}
            else if (length(var.names) != length(cord.X))
                stop.message('var.names', sample.X)
            var.names.list = unlist(lapply(seq_len(length(var.names)), function(x){if(var.names[x]){rownames(cord.X[[x]])}
                else {pch[(ind.group[x] + 1) : ind.group[x + 1]]}}))
        } else if (is.list(var.names)) {
            if (length(var.names) != length(cord.X))
                stop.message('var.names', sample.X)
            length.var.names <- sum(sapply(seq_len(length(var.names)), function(x){if (!lapply(var.names, is.logical)[[x]]){
                if (length(ind.var.sel) == 0){
                    length(var.names[[x]])
                } else {
                    length(var.names[[x]][ind.var.sel[[x]]])
                }
            } else {0}}))
            
            length.cord.X.names <- sum(sapply(seq_len(length(var.names)), function(x){
                if(!lapply(var.names, is.logical)[[x]]){nrow(cord.X[[x]])}else {0}
            }))
            if (length.var.names != length.cord.X.names){
                stop.message('var.names', sample.X)
            }
            
            var.names.list = unlist(sapply(seq_len(length(var.names)), function(x){if(lapply(var.names, is.logical)[[x]]){
                if (var.names[[x]]) {
                    row.names(cord.X[[x]])
                } else {
                    pch[(ind.group[x] + 1) : ind.group[x + 1]]
                }
            } else {
                if (length(ind.var.sel) == 0){
                    as.character(var.names[[x]])
                } else {
                    as.character(var.names[[x]])[ind.var.sel[[x]]]
                }
            }
            }))
            var.names = sapply(var.names, function(x){if(is.logical(x)){x}else{TRUE}})
        } else {
            stop.message('var.names', sample.X)
        }
        
        #-- Start: Computation ellipse
        circle = list()
        circle[[1]] = ellipse(0, levels = 1, t = 1)
        circle[[2]] = ellipse(0, levels = 1, t = rad.in)
        circle = data.frame(do.call("rbind", circle), "Circle" = c(rep("Main circle", 100), rep("Inner circle", 100)))
        #-- End: Computation ellipse
        
        #-- Start: data set
        df = data.frame(do.call(rbind, cord.X), "Block" = paste0("Block: ", unlist(lapply(seq_len(length(cord.X)), function(z){rep(blocks[z], nrow(cord.X[[z]]))}))))
        if (style=="3d")
            names(df)[seq_len(3)] = c("x", "y","z")
        else
            names(df)[seq_len(2)] = c("x", "y")
        
        df$names = as.vector(var.names.list)
        
        df$pch = pch; df$cex = cex; df$col = col; df$font = font
        
        if(missing.pch)
            df$pch=1
        
        if (overlap)
        {
            df$Overlap = title
            df$Block = factor(unlist(lapply(seq_len(length(cord.X)), function(z){rep(blocks[z], nrow(cord.X[[z]]))})))
            if(style %in%c("ggplot2","lattice"))
                title=NULL # to avoid double title
        } else {
            df$Overlap = df$Block
            if(style %in%c("ggplot2","lattice"))
                df$Block = factor(unlist(lapply(seq_len(length(cord.X)), function(z){rep(blocks[z], nrow(cord.X[[z]]))})))
        }
        
        if (cutoff != 0){
            if(style=="3d")
                df = df[abs(df$x) > cutoff | abs(df$y) > cutoff | abs(df$z) > cutoff, ,drop = FALSE]
            else
                df = df[abs(df$x) > cutoff | abs(df$y) > cutoff, ,drop = FALSE]
            ind.group = c(0, cumsum(table(df$Block)[unique(df$Block)])) # add unique to have names of cumsum matching the order of the blocks in df
        }
        
        if (nrow(df) == 0)
            stop("Cutoff value very high for the components ", comp1, " and ", comp2, ".No variable was selected.")
        
        
        #-- End: data set
        #save(list=ls(),file="temp.Rdata")
        #-- Start: ggplot2
        if (style == "ggplot2" &  plot)
        {
            Block = NULL# R check
            # visible variable issues for x, y and Circle
            # according to http://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when
            # one hack is to set to NULL first.
            x = y = Circle = NULL
            
            #-- Initialise ggplot2
            p = ggplot(df, aes(x = x, y = y, color = Block)) +
                labs(title = title, x = X.label, y = Y.label) + theme_bw()
            
            for (i in levels(df$Block))
                p = p + geom_point(data = subset(df, df$Block == i), size = 0, shape = 0)
            
            #-- Display sample or var.names
            for (i in seq_len(length(var.names))){
                if (var.names[i]) {
                    p = p + geom_text(data = df[c((ind.group[i] + 1) : ind.group[i + 1]), ],
                                      label = df[c((ind.group[i] + 1) : ind.group[i + 1]), "names"],
                                      size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                                      color = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                      fontface = df[c((ind.group[i] + 1) : ind.group[i + 1]), "font"])
                } else {
                    p = p + geom_point(data = df[c((ind.group[i] + 1) : ind.group[i + 1]), ],
                                       size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                                       shape = df[c((ind.group[i] + 1) : ind.group[i + 1]), "pch"],
                                       color = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"])
                }
            }
            
            #-- Modify scale colour - Change X/Ylabel - split plots into Blocks
            cols <- unique(df[,c('Block', "col")])
            col.per.group <- cols[,'col']
            names(col.per.group) <- cols[,'Block']
            
            p = p + scale_colour_manual(values = col.per.group)
            p = p + scale_x_continuous(limits = c(-1, 1)) + scale_y_continuous(limits = c(-1, 1))
            p = p + facet_wrap(~ Overlap, ncol = 2, as.table = TRUE)
            
            #-- Legend
            if (!legend)
            {
                p = p + theme(legend.position="none")
            } else {
                p = p + guides(colour = guide_legend(override.aes = list(shape = 19,
                                                                         size = unique(df$cex))))
            }
            
            
            
            #-- abline
            if (abline)
                p = p + geom_vline(aes(xintercept = 0), linetype = 2,
                                   colour = "darkgrey") + geom_hline(aes(yintercept = 0),linetype = 2,
                                                                     colour = "darkgrey")
            
            #-- circle correlation
            for (i in c("Main circle", "Inner circle")){
                p = p + geom_path(data = subset(circle, Circle == i),
                                  aes(x = .data[["x"]], y = .data[["y"]]), color = "Black")
            }
            
            print(p)
        }
        #-- End: ggplot2
        
        #-- Start: Lattice
        if(style == "lattice" )
        {
            legend.lattice = list(space = "right", title = legend.title, cex.title = 1.25,
                                  points=list(col=unique(df$col),cex = unique(df$cex),pch = unique(df$pch)),
                                  text = list(blocks))
            
            if (overlap) {
                p = xyplot(y ~ x | Overlap, data = df, xlab = X.label, ylab = Y.label, main = title,
                           scales = list(x = list(relation = "free", limits = c(-1, 1)),
                                         y = list(relation = "free", limits = c(-1, 1))),
                           key=if (legend) {legend.lattice} else {NULL},
                           panel = function(x, y, ...) {
                               
                               #-- Abline
                               if (abline) {panel.abline(v = 0, lty = 2, col = "darkgrey")
                                   panel.abline(h = 0, lty = 2, col = "darkgrey")}
                               
                               #-- Display sample or row.names
                               for (i in seq_len(length(var.names))){
                                   if (var.names[i]) {
                                       panel.text(x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                  y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                  df[c((ind.group[i] + 1) : ind.group[i + 1]), "names"],
                                                  col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                  cex = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                                                  font = df[c((ind.group[i] + 1) : ind.group[i + 1]), "font"])
                                   } else {
                                       panel.points(x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                    y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                    col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                    cex = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                                                    pch = df[c((ind.group[i] + 1) : ind.group[i + 1]), "pch"])
                                   }
                               }
                           })
                print(p)
                
                panels = trellis.currentLayout(which = "panel")
                ind = which(panels == 1, arr.ind = TRUE)
                trellis.focus("panel",ind[2], ind[1],highlight = FALSE)
                for (i in seq_len(length(c("Main circle", "Inner circle")))){
                    panel.lines(x = circle[circle$Circle %in% c("Main circle", "Inner circle")[i], "x"],
                                y = circle[circle$Circle %in% c("Main circle", "Inner circle")[i], "y"],
                                col = "black")
                }
                trellis.unfocus()
            } else {
                p = xyplot(y ~ x | Block, data = df, xlab = X.label, ylab = Y.label, main = title, as.table = TRUE,
                           scales = list(x = list(relation = "free", limits = c(-1, 1)),
                                         y = list(relation = "free", limits = c(-1, 1))),
                           col = "white",
                           key=if (legend) {legend.lattice} else {NULL},
                )
                print(p)
                
                panels = trellis.currentLayout(which = "panel")
                for (k in seq_len(length(cord.X))) {
                    ind = which(panels == k, arr.ind = TRUE)
                    trellis.focus("panel",ind[2], ind[1],highlight = FALSE)
                    
                    if (var.names[k]){
                        panel.text(x = df[c((ind.group[k] + 1) : ind.group[k + 1]), "x"],
                                   y = df[c((ind.group[k] + 1) : ind.group[k + 1]), "y"],
                                   df[c((ind.group[k] + 1) : ind.group[k + 1]), "names"],
                                   col = df[c((ind.group[k] + 1) : ind.group[k + 1]), "col"],
                                   cex = df[c((ind.group[k] + 1) : ind.group[k + 1]), "cex"],
                                   font = df[c((ind.group[k] + 1) : ind.group[k + 1]), "font"])
                    } else {
                        panel.points(x = df[c((ind.group[k] + 1) : ind.group[k + 1]), "x"],
                                     y = df[c((ind.group[k] + 1) : ind.group[k + 1]), "y"],
                                     col = df[c((ind.group[k] + 1) : ind.group[k + 1]), "col"],
                                     cex = df[c((ind.group[k] + 1) : ind.group[k + 1]), "cex"],
                                     pch = df[c((ind.group[k] + 1) : ind.group[k + 1]), "pch"])
                    }
                    
                    for (i in seq_len(length(c("Main circle", "Inner circle")))){
                        panel.lines(x = circle[circle$Circle %in% c("Main circle", "Inner circle")[i], "x"],
                                    y = circle[circle$Circle %in% c("Main circle", "Inner circle")[i], "y"],
                                    col = "black")
                    }
                }
                trellis.unfocus()
            }
            
        }
        #-- End: Lattice
        
        #-- Start: graphics
        if(style=="graphics" )
        {
            
            
            if (overlap)
            {
                
                if(legend){
                    opar = par(no.readonly = TRUE)
                    par(mai=c( 1.360000, 1.093333, 1.093333,max(strwidth("Legend","inches"),max(strwidth(blocks,"inches"))+0.3)+0.2),xpd=TRUE)
                }
                
                plot(df$x, df$y, type = "n", xlab = X.label, ylab = Y.label, main = "", xlim = c(-1, 1), ylim = c(-1, 1))
                
                #-- Display sample or row.names
                for (i in seq_len(length(var.names))){
                    if (var.names[i]) {
                        text(x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                             y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                             labels = df[c((ind.group[i] + 1) : ind.group[i + 1]), "names"],
                             col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                             cex = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                             font = df[c((ind.group[i] + 1) : ind.group[i + 1]), "font"])
                    } else {
                        points(x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                               y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                               col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                               cex = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                               pch = df[c((ind.group[i] + 1) : ind.group[i + 1]), "pch"])
                    }
                }
                
                #-- legend
                if (legend)
                    legend(x = 1.09, y=0.2,
                           legend = blocks,
                           title=legend.title,
                           col = unique(df$col),
                           pch = unique(df$pch),
                           pt.cex = unique(df$cex),
                           bty = "n")
                
                #-- Abline
                if (abline)
                    abline(v = 0, h = 0, lty = 2, xpd = FALSE)
                
                #-- Ellipse
                for (i in c("Main circle", "Inner circle")){
                    lines(x = circle[circle$Circle == i, "x"], y = circle[circle$Circle == i, "y"], col = "black")
                }
                
                title(title)#, outer = TRUE, line = -1)
                
                if (legend) par(mai = opar$mai, xpd = opar$xpd)
                
            } else {
                opar <- par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
                #-- Define layout
                mat = matrix(seq_len((ceiling(length(cord.X)/2) * 2)), ceiling(length(cord.X)/2), min(length(cord.X), 2), byrow = TRUE)
                if (legend){
                    mat = matrix(rep(mat,each=2),nrow=nrow(mat),byrow=TRUE)
                    mat = cbind(mat,rep(max(mat) + 1, nrow(mat)))
                }
                
                layout(mat)
                for (k in seq_len(length(cord.X))){
                    #-- initialise plot
                    plot(df[df$Block %in% paste0("Block: ", blocks[k]), "x" ],
                         df[df$Block %in% paste0("Block: ", blocks[k]), "y" ],
                         type = "n", xlab = X.label, ylab = Y.label, main = paste0("Block: ", blocks[k]),
                         xlim = c(-1, 1), ylim = c(-1, 1))
                    
                    #-- Display sample or row.names
                    if (var.names[k]) {
                        text(x = df[df$Block %in% paste0("Block: ", blocks[k]), "x"],
                             y = df[df$Block %in% paste0("Block: ", blocks[k]), "y"],
                             labels = df[df$Block %in% paste0("Block: ", blocks[k]), "names"],
                             col = df[df$Block %in% paste0("Block: ", blocks[k]), "col"],
                             cex = df[df$Block %in% paste0("Block: ", blocks[k]), "cex"],
                             font = df[df$Block %in% paste0("Block: ", blocks[k]), "font"])
                    } else {
                        points(x = df[df$Block %in% paste0("Block: ", blocks[k]), "x"],
                               y = df[df$Block %in% paste0("Block: ", blocks[k]), "y"],
                               col = df[df$Block %in% paste0("Block: ", blocks[k]), "col"],
                               cex = df[df$Block %in% paste0("Block: ", blocks[k]), "cex"],
                               pch = df[df$Block %in% paste0("Block: ", blocks[k]), "pch"])
                    }
                    
                    #-- Abline
                    if (abline)
                        abline(v = 0, h = 0, lty = 2, xpd = FALSE)
                    
                    #-- Ellipse
                    for (i in c("Main circle", "Inner circle")){
                        lines(x = circle[circle$Circle == i, "x"], y = circle[circle$Circle == i, "y"], col = "black")
                    }
                }
                
                
                title(title, outer = TRUE, line = -1)
                if (length(cord.X) != max(mat) & length(cord.X) != 1){
                    for (i in seq_len((max(mat)-length(cord.X)))){
                        plot(1,1, type = "n", axes = FALSE, ann = FALSE)
                    }
                }
                if (legend)
                    legend("center",
                           legend = blocks,
                           title=legend.title,
                           col = unique(df$col),
                           pch = unique(df$pch),
                           cex = unique(df$cex),
                           bty = "n")
                
                par(opar)
            }
            
        }
        #-- End: graphics
        
        #-- Start: 3d
        if(style=="3d") {
            if(requireNamespace("rgl") == FALSE)
                stop("the rgl package is required for 3d plots")
            
            rgl::open3d()
            rgl::par3d(windowRect = c(500, 30, 1100, 630))
            Sys.sleep(0.5)
            
            if (!is.null(title)) {
                mat = matrix(seq_len(2), 2)
                rgl::layout3d(mat, heights = c(1, 10), model = "inherit")
                rgl::next3d()
                rgl::text3d(0, 0, 0, title)
                rgl::next3d()
            }
            
            rgl::par3d(userMatrix = rgl::rotationMatrix(pi/80, 1, -1/(100*pi), 0))
            
            
            
            
            
            if (legend) {
                rgl::legend3d(x="right",
                              legend = blocks,
                              col = unique(col),
                              pch = rep(16,length(unique(pch))),
                              pt.cex = unique(cex),
                              bty="n")
            }
            
            if (any(axes.box == "axes") || any(axes.box == "all"))
                rgl::axes3d(c('x','y','z'), pos = c(0, 0, 0), nticks = 2, at = c(-1.2, 1.2),
                            tick = FALSE, labels = "")
            
            for (i in seq_len(length(var.names))){
                if (var.names[i]) {
                    rgl::text3d(x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                z=df[c((ind.group[i] + 1) : ind.group[i + 1]), "z"],
                                texts = df[c((ind.group[i] + 1) : ind.group[i + 1]), "names"],
                                color = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                cex = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"],
                                font = df[c((ind.group[i] + 1) : ind.group[i + 1]), "font"])
                } else {
                    switch(unique(df[c((ind.group[i] + 1) : ind.group[i + 1]), "pch"]),
                           sphere = rgl::plot3d(x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                z=df[c((ind.group[i] + 1) : ind.group[i + 1]), "z"], type = "s",
                                                col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"], radius = cex/20, add = TRUE),
                           tetra = rgl::shapelist3d(rgl::tetrahedron3d(), x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                    y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                    z=df[c((ind.group[i] + 1) : ind.group[i + 1]), "z"],
                                                    col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                    size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"]/25),
                           cube = rgl::shapelist3d(rgl::cube3d(), x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                   y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                   z=df[c((ind.group[i] + 1) : ind.group[i + 1]), "z"],
                                                   col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                   size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"]/30),
                           octa = rgl::shapelist3d(rgl::octahedron3d(), x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                   y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                   z=df[c((ind.group[i] + 1) : ind.group[i + 1]), "z"],
                                                   col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                   size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"]/17),
                           icosa = rgl::shapelist3d(rgl::icosahedron3d(), x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                    y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                    z=df[c((ind.group[i] + 1) : ind.group[i + 1]), "z"],
                                                    col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                    size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"]/20),
                           dodeca = rgl::shapelist3d(rgl::dodecahedron3d(), x = df[c((ind.group[i] + 1) : ind.group[i + 1]), "x"],
                                                     y = df[c((ind.group[i] + 1) : ind.group[i + 1]), "y"],
                                                     z=df[c((ind.group[i] + 1) : ind.group[i + 1]), "z"],
                                                     col = df[c((ind.group[i] + 1) : ind.group[i + 1]), "col"],
                                                     size = df[c((ind.group[i] + 1) : ind.group[i + 1]), "cex"]/20))
                }
            }
            
            rgl::par3d(cex = 0.8)
            
            #-- draws axes --#
            if (any(axes.box == "axes") || any(axes.box == "all")) {
                if (any(label.axes.box == "axes") || any(label.axes.box == "both")) {
                    rgl::text3d(1.2, -0.05, 0, texts = X.label, cex = 0.8, color = "black")
                    rgl::text3d(0, 1.27, 0, texts = Y.label, cex = 0.8, color = "black")
                    rgl::text3d(0, -0.05, 1.2, texts = Z.label, cex = 0.8, color = "black")
                }
                X =  c(1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09, 1.09, 1.2, 1.09,  1.09,
                       0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
                       0.0, 0.0,  0.0, 0.0, 0.0,  0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4))
                
                Y = c(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.035, -0.035, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
                      1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09,
                      0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4))
                
                Z = c(0.0, 0.035, -0.035, 0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4),
                      0.0, 0.035, -0.035, 0.0, 0.0,  0.0, 0.0, 0.035*sin(pi/4), -0.035*sin(pi/4), 0.0, -0.035*sin(pi/4), 0.035*sin(pi/4),
                      1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09, 1.2, 1.09,  1.09)
                rgl::triangles3d(x = X, y = Y, z = Z, col = "black")
                
            }
            
            rgl::points3d(1.2, 0, 0, size = 0.1, alpha = 0)
            rgl::points3d(0, 1.2, 0, size = 0.1, alpha = 0)
            rgl::points3d(0, 0, 1.2, size = 0.1, alpha = 0)
            rgl::points3d(-1.2, 0, 0, size = 0.1, alpha = 0)
            rgl::points3d(0, -1.2, 0, size = 0.1, alpha = 0)
            rgl::points3d(0, 0, -1.2, size = 0.1, alpha = 0)
            
            #-- draws sphere --#
            rgl::spheres3d(0, 0, 0, radius = rad.in, front = "fill", back = "fill", emission = gray(0.9), alpha = 0.4)
            rgl::spheres3d(0, 0, 0, radius = rad.in, front = "line", back = "line", emission = gray(0.9))
            
            #-- draws axes/box and add axes labels --#
            if (any(axes.box == "box") || any(axes.box == "all")) {
                rgl::axes3d(marklen = 25)
                rgl::box3d()
                if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
                    rgl::mtext3d(X.label, "x-+", line = 1)
                    rgl::mtext3d(Y.label, "y-+", line = 1.5)
                    rgl::mtext3d(Z.label, "z+-", line = 1)
                }
            }
            
            if (any(axes.box == "bbox") || any(axes.box == "all")) {
                rgl::bbox3d(color = c("#333377", "black"), emission = gray(0.5),
                            specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)
                if (any(label.axes.box == "box") || any(label.axes.box == "both")) {
                    rgl::mtext3d(X.label, "x-+", line = 1)
                    rgl::mtext3d(Y.label, "y-+", line = 1.5)
                    rgl::mtext3d(Z.label, "z+-", line = 1)
                }
            }
            
            
        }
        #-- End: graphics
        if(plot){
            return(invisible(df))}
        else
            return(df)
        
    }

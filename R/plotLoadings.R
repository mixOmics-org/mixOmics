## ------------------------------------------------------------------------ ##
####                               Generic                                ####
## ------------------------------------------------------------------------ ##

#' Plot of Loading vectors
#' 
#' This function provides a horizontal bar plot to visualise loading vectors.
#' For discriminant analysis, it provides visualisation of highest or lowest
#' mean/median value of the variables with color code corresponding to the
#' outcome of interest.
#' 
#' 
#' The contribution of each variable for each component (depending on the
#' object) is represented in a barplot where each bar length corresponds to the
#' loading weight (importance) of the feature. The loading weight can be
#' positive or negative.
#' 
#' For discriminant analysis, the color corresponds to the group in which the
#' feature is most 'abundant'.  Note that this type of graphical output is
#' particularly insightful for count microbial data - in that latter case using
#' the \code{method = 'median'} is advised. Note also that if the parameter
#' \code{contrib} is not provided, plots are white.
#' 
#' For MINT analysis, \code{study="global"} plots the global loadings while
#' partial loadings are plotted when \code{study} is a level of
#' \code{object$study}. Since variable selection in MINT is performed at the
#' global level, only the selected variables are plotted for the partial
#' loadings even if the partial loadings are not sparse. See references.
#' Importantly for multi plots, the legend accounts for one subplot in the
#' layout design.
#' 
#' The \code{block} argument behavior varies depending on the object type. 
#' For \code{mixo_pls}, \code{mixo_spls}, \code{rcc}, \code{rgcca}, \code{sgcca},
#' \code{block} can be any number of blocks from \code{object$names$blocks}.
#' For \code{mint.pls}, \code{mint.spls}: when \code{study="all.partial"},
#' can only be one block from \code{object$names$blocks}.
#' 
#' @aliases plotLoadings.pls plotLoadings.spls
#' @param object object
#' @param comp integer value, which component to plot. Default is 1.
#' @param ndisplay integer indicating how many of the most important variables
#' are to be plotted (ranked by decreasing weights in each component).
#' Useful to lighten a graph.
#' @param style argument to be set to either \code{'graphics'} or \code{'ggplot2'} 
#' to indicate the style of the plot. Default is \code{'graphics'}.
#' @param layout Vector of two values (rows,cols) that indicates the layout of
#' the plot. If \code{layout} is provided, the remaining empty subplots are
#' still active. See details and examples. 
#' @param xlim Limit of the x-axis, vector of length 2.
#' @param col color of barplot, only for non-DA methods.
#' @param border Argument from \code{\link{barplot}}: indicates whether to draw
#' a border on the barplot and in which colour. Note that for -DA methods, 
#' when \code{contrib} is not provided, the colour is set to white and border black.
#' @param name.var A character vector or list indicating the names of the variables.
#' For single block analysis or when \code{study="all.partial"}, a vector of length equal 
#' to the number of variables in the block. For multi-block analysis, a list where each element 
#' is a vector of length equal to the number of variables in the corresponding block.
#' The names of the vector should match the names of the input data, see example.
#' @param size.name A numerical value giving the amount by which plotting the
#' variable name text should be magnified or reduced relative to the default.
#' @param title Title of the plot. Default is NULL.
#' @param size.title size of the title
#' @param subtitle subtitle for each plot, only used when several \code{block}
#' or \code{study} are plotted.
#' @param size.subtitle size of the subtitle
#' @param size.axis size of text on the X axis
#' @param X.label X axis label. Default is NULL.
#' @param Y.label Y axis label. Default is NULL.
#' @param size.labs size of the axis labels. 
#' @param block A single value or vector indicating which block to plot. 
#' See details for behavior depending on object type.
#' @param contrib a character set to 'max' or 'min' indicating if the color of
#' the bar should correspond to the group with the maximal or minimal
#' expression levels / abundance.
#' @param method a character set to 'mean' or 'median' indicating the criterion
#' to assess the contribution. We recommend using median in the case of count
#' or skewed data.
#' @param show.ties Logical. If TRUE then tie groups appear in the color set by
#' \code{col.ties}, which will appear in the legend. Ties can happen when
#' dealing with count data type. By default set to TRUE.
#' @param col.ties Color corresponding to ties, only used if
#' \code{show.ties=TRUE} and ties are present.
#' @param legend Logical indicating if the legend indicating the group outcomes
#' should be added to the plot. Default value is TRUE.
#' @param legend.color A color vector of length the number of group outcomes.
#' See examples.
#' @param legend.title A set of characters to indicate the title of the legend.
#' Default value is NULL.
#' @param size.legend A numerical value giving the amount by which plotting the
#' legend text should be magnified or reduced relative to the default.
#' @param study Indicates which study are to be plotted.  A character vector
#' containing some levels of \code{object$study}, "all.partial" to plot all
#' studies or "global" is expected.
#' @param \dots not used.
#' @return Invisibly returns a \code{data.frame} containing the contribution of 
#' features on each component. For supervised models the contributions for
#'  each class is also specified. See details.
#' @author Florian Rohart, Kim-Anh Lê Cao, Benoit Gautier, Al J Abadi
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{mint.pls}}, \code{\link{mint.spls}},
#' \code{\link{mint.plsda}}, \code{\link{mint.splsda}},
#' \code{\link{block.pls}}, \code{\link{block.spls}},
#' \code{\link{block.plsda}}, \code{\link{block.splsda}},
#' \code{\link{mint.block.pls}}, \code{\link{mint.block.spls}},
#' \code{\link{mint.block.plsda}}, \code{\link{mint.block.splsda}}
#' @references Rohart F. et al (2016, submitted). MINT: A multivariate
#' integrative approach to identify a reproducible biomarker signature across
#' multiple experiments and platforms.
#' 
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013). Multi-group
#' PLS Regression: Application to Epidemiology. In New Perspectives in Partial
#' Least Squares and Related Methods, pages 243-255. Springer.
#' 
#' Singh A., Shannon C., Gautier B., Rohart F., Vacher M., Tebbutt S.
#' and Lê Cao K.A. (2019), DIABLO: an integrative approach for identifying key 
#' molecular drivers from multi-omics assays, Bioinformatics, 
#' Volume 35, Issue 17, 1 September 2019, Pages 3055–3062.
#' 
#' Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009).
#' Sparse canonical methods for biological data integration: application to a
#' cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
#' 
#' Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris:
#' Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' @keywords multivariate
#' @export
#' @example ./examples/plotLoadings-examples.R
plotLoadings  <-
    function(object, ...)
        UseMethod("plotLoadings")

## ------------------------------------------------------------------------ ##
####                               Helpers                                ####
## ------------------------------------------------------------------------ ##

check.input.plotLoadings <- function(object,
                                     block,
                                     study,
                                     subtitle,
                                     size.name,
                                     size.legend,
                                     title,
                                     col,
                                     contrib,
                                     name.var,
                                     xlim)
{
    
    if (is.null(object$loadings))
        stop("'plotLoadings' should be used on object for which object$loadings is present.")
    
    # block
    # --
    if (missing(block))
    {
        if (!inherits(object, "DA"))
        {
            block = object$names$blocks
        } else  if (inherits(object, c("mixo_plsda", "mixo_splsda"))) {
            block = "X"
        } else {
            if (!is.null(object$indY))
            {
                block = object$names$blocks[-object$indY]
            } else {
                block = object$names$blocks
            }
        }
    }
    
    if (inherits(object, c("mixo_plsda", "mixo_splsda")) & (!all(block %in% c(1,"X")) | length(block) > 1 ))
        stop("'block' can only be 'X' or '1' for plsda and splsda object")
    
    if (inherits(object, c("mixo_plsda", "mixo_splsda","pca")))
    {
        object$indY = 2
    } else if (inherits(object, c("mixo_pls", "mixo_spls"))) {
        object$indY = 3 # we don't want to remove anything in that case, and 3 is higher than the number of blocks which is 2
    }
    
    if(!inherits(object, "DA"))
        object$indY = length(object$names$blocks)+1  # we don't want to remove anything in that case, and 3 is higher than the number of blocks which is 2
    
    if(is.numeric(block))
    {
        if(any(block>length(object$names$blocks[-object$indY])))
            stop("'block' needs to be lower than the number of blocks in the fitted model, which is ",length(object$names$blocks)-1)
        
    }else if(is.character(block) & any(is.na(match(block,object$names$blocks[-object$indY])))) {
        stop("Incorrect value for 'block', 'block' should be among the blocks used in your object: ", paste(object$names$blocks[-object$indY],collapse=", "), call. = FALSE)
    }
    
    
    if (!missing(subtitle))
    {
        if (length(subtitle)!=length(block))
            stop("'subtitle' indicates the subtitle of the plot for each block and it needs to be the same length as 'block'.")
    }
    
    if(!missing(study))
    {
        #study needs to be either: from levels(object$study), numbers from 1:nlevels(study) or "global"
        if (any(!study%in%c(levels(object$study), "global", "all.partial")))
            stop("'study' must be one of 'object$study' or 'all.partial' or 'global.")
        
        if (length(study)!=length(unique(study)))
            stop("Duplicate in 'study' not allowed")
    }
    
    # cex
    # --
    if (size.name <= 0)
        size.name = 0.7
    
    if (!missing(size.legend))
    {
        if(size.legend <= 0)
            size.legend = 0.8
    } else {
        size.legend = NULL
    }
    
    # contrib
    # --
    if(!missing(contrib))
    {
        if(length(contrib) > 1 | !all(contrib %in% c("min", "max")))
            stop("'contrib' must be either 'min' or 'max'")
        
    }
    
    # xlim
    #---
    if(!missing(xlim))
    {
        # check xlim, has to be a matrix with number of rows=number of blocks, or a vector of two values
        if(length(block) == 1 & !is.null(xlim))
        {
            if(length(xlim) !=2)
                stop("'xlim' must be a vector of length 2")
            
            xlim = matrix(xlim, nrow = 1)
        }
        
        if(length(block)>1 & !is.null(xlim))
        {
            if(is.matrix(xlim) && ( !nrow(xlim) %in%c(1, length(block))  | ncol(xlim) != 2 ))
                stop("'xlim' must be a matrix with ",length(block)," rows (length(block)) and 2 columns")
            
            if(is.vector(xlim))
            {
                if(length(xlim) !=2)
                    stop("'xlim' must be a matrix with ",length(block)," rows (length(block)) and 2 columns")
                
                xlim = matrix(xlim, nrow = 1)
            }
            
            if(nrow(xlim) != length(block)) # we complete xlim to have one xlim per block
                xlim = matrix(rep(xlim, length(block)), nrow = length(block), byrow=TRUE)
        }
        
    } else {
        xlim = NULL
    }
    
    #names.var
    #-----
    if(!is.null(name.var))
    {
        # For all.partial case, name.var should be a vector
        if (!missing(study) && any(study == "all.partial")) {
          if (length(block)>1){stop("When study = 'all.partial', block must be NULL or length of 1")}
            if (!is.vector(name.var)) {
                stop("When study = 'all.partial', 'name.var' should be a vector of length equal to the number of variables in the selected block")
            }
            # Check length matches the number of variables in the selected block
            if (length(name.var) != nrow(object$loadings[[block]])) {
                stop("When study = 'all.partial', 'name.var' should be a vector of length ", nrow(object$loadings[[block]]))
            }
        } else {
            # For other cases, name.var should be a list
            if (length(block) > 1 && length(block) != length(name.var)) {
                stop("'name.var' has to be a list of length the number of block to plot: ", length(block))
            }
            
            if (length(block) > 1) {
                for (block_i in block) {
                    if(length(name.var[[block_i]])!= nrow(object$loadings[[block_i]]))
                        stop("For block '", block_i,"', 'name.var' should be a vector of length ", nrow(object$loadings[[block_i]]))
                }
            } else {
                if(length(name.var)!= nrow(object$loadings[[block]]))
                    stop("For block '", block,"', 'name.var' should be a vector of length ", nrow(object$loadings[[block]]))
            }
        }
    }
    
    #title
    #-----
    if (!is.null(title) & !is.character(title))
        warning('title needs to be of type character')
    
    #col
    #-----
    if (!is.null(col) & (length(col) !=  1))
    {
        warning('col must be the of length 1, by default set to default colors')
        col = color.mixo(1)  # by default set to the colors in color.mixo (10 colors)
    }
    if (is.null(col))
        col = color.mixo(1) # by default set to the colors in color.mixo (10 colors)
    
    
    return(list(col = col, size.name = size.name, size.legend = size.legend, block = block, xlim = xlim))
}

layout.plotLoadings <- function(layout, plot, legend, block)
{
    # layout
    # --
        opar = par(no.readonly = TRUE)
        reset.mfrow = FALSE # if set to TRUE, the algorithm ends up with  par(mfrow=reset.mfrow)
        nResp = length(block) + length(block) * legend  #number of blocks *2 if legend is plotted
        
        if (is.null(layout))
        {
            # check if there are enough plots in mfrow
            omfrow = par("mfrow")
            available.plots = prod(omfrow)
            if (available.plots<nResp) # if not enough plots available, we create our new plot
            {
                
                if (legend)
                {
                    nRows = min(c(2, ceiling(nResp / 4)))
                    nCols = min(c(4, nResp))
                    layout(matrix(1 : (nCols * nRows), nRows, nCols, byrow=TRUE),rep(c(0.7,0.7 -0.4*legend),nCols/(1+legend)))
                } else {
                    nRows = min(c(3, ceiling(nResp/3)))
                    nCols = min(c(3, ceiling(nResp / nRows)))
                    
                    layout(matrix(1 : (nCols * nRows), nRows, nCols, byrow=TRUE))
                    
                }
                if (nRows * nCols < nResp)
                    devAskNewPage(TRUE)
                
                reset.mfrow=TRUE # we changed mfrow to suits our needs, so we reset it at the end
            }
        } else {
            if (length(layout) != 2 || !is.numeric(layout) || any(is.na(layout)))
                stop("'layout' must be a numeric vector of length 2.")
            
            nRows = layout[1]
            nCols = layout[2]
            par(mfrow = layout)
            
            if (nRows * nCols < nResp)
                devAskNewPage(TRUE)
        }
        
    
    return(list(reset.mfrow = reset.mfrow, opar = opar))
    
}


get.loadings.ndisplay <- function(object,
                                  comp,
                                  block,
                                  name.var,
                                  name.var.complete,
                                  ndisplay)
{
    ##selectvar
    selected.var = selectVar(object, comp = comp, block = block) # gives name and values of the blocks in 'block'
    name.selected.var = selected.var[[1]]$name
    value.selected.var = selected.var[[1]]$value
    
    # ndisplay
    # ------
    # if null set by default to all variables from selectVar
    if (is.null(ndisplay))
    {
        ndisplay.temp = length(name.selected.var)
    } else if (ndisplay > length(name.selected.var)) {
        message("'ndisplay' value is larger than the number of selected variables! It has been reseted to ", length(name.selected.var), " for block ", block)
        ndisplay.temp = length(name.selected.var)
    } else {
        ndisplay.temp = ndisplay
    }
    
    name.selected.var = name.selected.var[1:ndisplay.temp]
    value.selected.var = value.selected.var[1:ndisplay.temp,]
    
    #comp
    # ----
    if (inherits(object, c("mixo_pls","mixo_spls", "rcc")))# cause pls methods just have 1 ncomp, block approaches have different ncomp per block
    {
        ncomp = object$ncomp
        object$X = list(X = object$X, Y = object$Y) # so that the data is in object$X, either it's a pls or block approach
    } else {
        ncomp = object$ncomp[block]
    }
    
    if (any(max(comp) > ncomp))
        stop(paste("Argument 'comp' should be less or equal to ", ncomp))
    
    names.block = as.character(names(selected.var)[1]) #it should be one block and ncomp, so we take the first one
    
    X = object$X[names.block][[1]]
    
    #name.var
    ind.match = match(name.selected.var, colnames(X)) # look at the position of the selected variables in the original data X
    if(!is.null(name.var))
    {
        if(length(name.var)!= ncol(X))
            stop("For block '", names.block,"', 'name.var' should be a vector of length ", ncol(X))
        
        colnames.X = as.character(name.var[ind.match]) # get the
    }else{
        colnames.X = as.character(colnames(X))[ind.match]
    }
    X = X[, name.selected.var, drop = FALSE] #reduce the problem to ndisplay
    
    #completing colnames.X by the original names of the variables when missing
    if (name.var.complete == TRUE)
    {
        ind = which(colnames.X == "")
        if (length(ind) > 0)
            colnames.X[ind] = colnames(X)[ind]
    }
    
    
    return(list(X = X, names.block = names.block, colnames.X = colnames.X, name.selected.var = name.selected.var, value.selected.var = value.selected.var))
}



get.contrib.df <- function(Y,
                           X,
                           method,
                           contrib,
                           value.selected.var,
                           colnames.X,
                           name.selected.var,
                           legend.color,
                           col.ties)
{
    # Start: Initialisation
    which.comp = method.group = list()
    which.contrib = data.frame(matrix(FALSE, ncol = nlevels(Y) + 2, nrow = length(colnames.X),
                                      dimnames = list(name.selected.var, c(paste0("Contrib.", levels(Y)), "Contrib", "GroupContrib"))))
    # End: Initialisation
    
    # calculate the max.method per group for each variable, and identifies which group has the max max.method
    for(k in 1:ncol(X))
    {
        method.group[[k]] = tapply(X[, k], Y, method, na.rm=TRUE) #method is either mean or median
        # determine which group has the highest mean/median
        which.contrib[k, 1:nlevels(Y)] = (method.group[[k]]) == get(contrib)((method.group[[k]]), na.rm=TRUE) # contrib is either min or max
        
    }
    
    # we also add an output column indicating the group that is max
    # if ties, we set the color to white
    which.contrib$color = apply(which.contrib, 1, function(x)
    {
        if (length(which(x)) > 1)
        {
            return(col.ties)
        } else { # otherwise we use legend color provided
            return(legend.color[1 : nlevels(Y)][which(x)])
        }
    })
    
    which.contrib$GroupContrib = apply(which.contrib[, 1:(nlevels(Y))], 1, function(x)
    {
        if (length(which(x)) > 1)
        {
            return("tie")
        } else {
            return(levels(Y)[which(x)])
        }
    })
    
    method.group = do.call(rbind, method.group)
    df = data.frame(method.group, which.contrib, importance = value.selected.var)
    return(df)
}

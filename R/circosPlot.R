#' circosPlot for DIABLO
#' 
#' Displays variable correlation among different blocks
#' 
#' \code{circosPlot} function depicts correlations of variables selected with
#' \code{block.splsda} or \code{block.spls} among different blocks, using a
#' generalisation of the method presented in González et al 2012. If
#' \code{ncomp} is specified, then only the variables selected on that component
#' are displayed.
#' 
#' The \code{linkWidth} argument specifies the width of the links drawn.
#' If a vector of length 2 is provided, the smaller value will correspond to
#' a similarity values designated by \code{cutoff} argument, while the 
#' larger value will be used for a link with perfect similarity (1), if any.
#' @param object An object of class inheriting from \code{"block.plsda"},
#' \code{"block.splsda"}, \code{"block.pls"} or \code{"blocks.spls"}.
#' @param comp Numeric vector indicating which component to plot. Default to
#' all
#' @param cutoff Only shows links with a correlation higher than \code{cutoff}
#' @param blocks Character or integer vector indicating which blocks to show.
#' Default to all
#' @param color.Y a character vector of colors to be used for the levels of the
#' outcome
#' @param color.blocks a character vector of colors to be used for the blocks
#' @param color.cor a character vector of two colors. First one is for the
#' negative correlation, second one is for the positive correlation
#' @param var.names Optional parameter. A list of length the number of blocks
#' in \code{object$X}, containing the names of the variables of each block. If
#' \code{NULL}, the colnames of the data matrix are used.
#' @param showIntraLinks if TRUE, shows the correlation higher than the
#' threshold inside each block.
#' @param line if TRUE, shows the overall expression of the selected variables.
#' see examples.
#' @param size.legend size of the legend
#' @param ncol.legend number of columns for the legend
#' @param size.variables size of the variable labels
#' @param size.labels size of the block labels
#' @param legend Logical. Whether the legend should be added. Default is TRUE.
#' @param legend.title String. Name of the legend. Defaults to "Expression".
#' @param linkWidth Numeric. Specifies the range of sizes used for lines linking
#' the correlated variables (see details). Must be of length 2 or 1. Default to c(1). See details.
#' @param ... For object of class \code{block.splsda}, advanced plot parameters:
#' \itemize{
#'  \item \bold{var.adj} Numeric. Adjusts the radial location of variable names in 
#'  units of the arc band width. Positive values push feature names radially 
#'  outward. Default to -0.33. See examples.
#'  \item \bold{block.labels.adj} Numeric between -1 and 1. Adjusts the radial
#'  location of block names radially inward or outward. Default to 0. See examples.
#'  \item \bold{blocks.link} Character vector of blocks. If provided, only correlations
#'  from features of these blocks are shown using links. See examples.
#' }
#' For object of class \code{block.spls}, all listed and advanced arguments
#' passed to the \code{block.splsda} method.
#' @return If saved in an object, the circos plot will output the similarity
#' matrix and the names of the variables displayed on the plot (see
#' \code{attributes(object)}).
#' @author Michael Vacher, Amrit Singh, Florian Rohart, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{block.splsda}}, references and
#' http://www.mixOmics.org/mixDIABLO for more details.
#' @references Singh A., Gautier B., Shannon C., Vacher M., Rohart F., Tebbutt
#' S. and Lê Cao K.A. (2016). DIABLO: multi omics integration for biomarker
#' discovery. BioRxiv available here:
#' \url{http://biorxiv.org/content/early/2016/08/03/067611}
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' 
#' González I., Lê Cao K.A., Davis M.J., Déjean S. (2012). Visualising
#' associations between paired 'omics' data sets. \emph{BioData Mining};
#' \bold{5}(1).
#' @keywords regression multivariate
#' @example ./examples/circosPlot-examples.R
#' @importFrom reshape2 dcast
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate summarise arrange 
#' @name circosPlot
NULL
#' @export
#' @noRd
circosPlot <- function(object, ...) UseMethod('circosPlot')

#' @keywords Internal
.circosPlot <- 
    function(object,
             comp = 1:min(object$ncomp),
             cutoff,
             color.Y,
             blocks = NULL,
             color.blocks,
             color.cor,
             var.names = NULL,
             showIntraLinks = FALSE,
             line = FALSE,
             size.legend = 0.8,
             ncol.legend = 1,
             size.variables = 0.25,
             size.labels = 1,
             legend = TRUE,
             legend.title = "Expression",
             linkWidth = 1,
             ...)
    {
        
        if (inherits(object, c("block.plsda", "block.splsda"))) {
            indY <- object$indY
            object$variates[indY] <- NULL
            object$loadings[indY] <-  NULL
        }
        
        # to satisfy R CMD check that doesn't recognise x, y and group (in aes)
        Features = Exp = Dataset = Mean = NULL
        
        
        figSize = 800
        segmentWidth = 25
        linePlotWidth = 90
        
        # check input object
        if (!(inherits(object, c("block.plsda", "block.splsda", "block.pls", "block.spls"))))
            stop("circosPlot is only available for 'block.(s)pls(da)' objects")
        
        if (length(object$X) < 2)
            stop("This function is only available when there are more than 3 blocks
    (2 in object$X + an outcome object$Y)") # so 2 blocks in X + the outcome Y
        if (!is.null(blocks))
        { ## allow selection of blocks to show
            # TODO checks for valid blocks
            
            oblocks <- names(object$X)[-object$indY]
            
            if (is.numeric(blocks))
            {
                blocks <- oblocks[blocks]
            }
            
            drop.blocks <- which(!oblocks %in% blocks)
            if (length(drop.blocks) > 0) {
                object$X <- object$X[-drop.blocks]
                object$variates <- object$variates[-drop.blocks]
                object$loadings <- object$loadings[-drop.blocks]
            }
        }
        
        if (missing(cutoff))
            stop("'cutoff' is missing", call.=FALSE) # so 2 blocks in X + the outcome Y
        
        if(missing(color.Y))
        {
            color.Y = color.mixo(1:nlevels(object$Y))
        } else {
            if(length(color.Y) != nlevels(object$Y))
                stop("'color.Y' must be of length ", nlevels(object$Y))
            
        }
        if (missing(color.blocks))
        {
            color.blocks = brewer.pal(n = 12, name = 'Paired') #why 12?? ANS: bc max allowed n is 12 for this function
            if (length(object$X) > 6) {
                ## if more than 6 X, extend the palette
                ## 2*length(object$X) because we need dark and clear cols for each block
                color.blocks <- colorRampPalette(color.blocks)(2*length(object$X))
            }
            
        } else {
            if(length(color.blocks) != length(object$X))
                stop("'color.blocks' must be of length ", length(object$X))
            
            color.blocks.adj = adjustcolor(color.blocks, alpha.f = 0.5)
            #to get two shades of the same color per block
            
            color.blocks = c(rbind(color.blocks, color.blocks.adj))
            # to put the color next to its shaded color
        }
        
        if(missing(color.cor))
        {
            color.cor = c(colors()[134],  # blue, negative correlation
                          colors()[128])  # pale red, positive correlation
        } else {
            if(length(color.cor) != 2)
                stop("'color.cor' must be of length 2")
        }
        
        ## create unique feature names for circosPlot
        object$X <- mapply(object$X, names(object$X), FUN = function(mat,blockname)
        {
            colnames(mat) <- sprintf("%s_%s", colnames(mat), blockname)
            mat
        }, SIMPLIFY = FALSE)
        loadings.X <- object$loadings
        original.var.names <- lapply(loadings.X, rownames)
        loadings.X <- mapply(loadings.X, names(loadings.X), FUN = function(mat,blockname)
        {
            rownames(mat) <- sprintf("%s_%s", rownames(mat), blockname)
            mat
        }, SIMPLIFY = FALSE)
        object$loadings <- loadings.X
        ## DONE - create unique feature names for circosPlot
        Y = object$Y
        X = object$X
        
        #check var.names
        sample.X = lapply(object$loadings,
                          function(x){1 : nrow(x)})
        if (is.null(var.names))
        {
            var.names.list = original.var.names
        } else if (is.list(var.names)) {
            if (length(var.names) != length(object$loadings))
                stop.message('var.names', sample.X)
            
            if(sum(sapply(1 : length(var.names), function(x){
                length(var.names[[x]]) == length(sample.X[[x]])})) !=
               length(var.names))
                stop.message('var.names', sample.X)
            
            var.names.list = var.names
        } else {
            stop.message('var.names', sample.X)
        }
        
        
        if(any(comp > min(object$ncomp)))
        {
            warning("Limitation to ",min(object$ncomp),
                    " components, as determined by min(object$ncomp)")
            comp[which(comp > min(object$ncomp))] = min(object$ncomp)
        }
        comp = unique(sort(comp))
        
        ## check linkWidth
        invalid.linkWidth <- FALSE
        if (mode(linkWidth) == 'numeric')
        {
            if (length(linkWidth) == 1)
                linkWidth <- rep(linkWidth, 2)
            else if (length(linkWidth) != 2)
                invalid.linkWidth <- TRUE
            
            linkWidth <- sort(linkWidth)
        } else {
            invalid.linkWidth <- TRUE
        }
        
        if (invalid.linkWidth)
            stop("'linkWidth' must be a numeric of length 2 (or 1) specifying ",
                 "the range of widths used for link lines based on similarity measures.", call. = FALSE)
        
        keepA = lapply(object$loadings, function(i)
            apply(abs(i)[, comp, drop = FALSE], 1, sum) > 0)
        cord = mapply(function(x, y, keep){
            cor(x[, keep, drop = FALSE], y[, comp, drop = FALSE], use = "pairwise")
        }, x=X, y=object$variates,
        keep = keepA, SIMPLIFY = FALSE)
        
        simMatList = vector("list", length(X))
        for(i in 1:length(cord))
        {
            for(j in 1:length(cord))
            {
                simMatList[[i]][[j]] = cord[[i]] %*% t(cord[[j]])
            }
        }
        simMatList <- lapply(simMatList, function(i) do.call(cbind, i))
        simMat = do.call(rbind, simMatList)
        
        ## merge all similarity measures in one data.frame where rows are samples
        Xdat = as.data.frame(do.call(cbind, X)[, colnames(simMat)])
        ## add Y as another feature
        AvgFeatExp0 <- mutate(.data = Xdat, Y = Y)
        ## create a long data.frame with columns: CELL FEATURE FEATURE_VALUE -Y
        AvgFeatExp0 <- gather(data = AvgFeatExp0, Features, Exp, -Y)
        ## group by Y and feature for averaging and calcuate the mean and SD to 
        ## get a data.frame with columns: Y FEATURES MEAN SD
        AvgFeatExp0 <- group_by(.data = AvgFeatExp0, Y, Features)
        AvgFeatExp0 <- summarise(.data = AvgFeatExp0, Mean = mean(Exp, na.rm = TRUE), SD = sd(Exp,na.rm = TRUE))
        ## add a column as the dataset from which the feature came from:
        AvgFeatExp0$Dataset <- factor(rep(names(X), unlist(lapply(cord, nrow))),
                                      levels = names(X))[match(AvgFeatExp0$Features,colnames(Xdat))]
        # to match Xdat that is reordered in AvgFeatExp0
        featExp <- group_by(.data = AvgFeatExp0 , Dataset, Y)
        featExp <- arrange(.data = featExp , Mean)
        #Generate a circular plot (circos like) from a correlation matrix (pairwise)
        #
        # Args:
        #   simMat: the main correlation matrix.
        #         -> colnames == rownames (pairwise)  values = correlations
        #   featExp: data.frame holding the expression data.
        #   cutoff: minimum value for correlations (<threshold will be ignored)
        #   figSize: figure size
        #   segmentWidth: thickness of the segment (main circle)
        #   linePlotWidth: thickness of the line plot (showing expression data)
        #   showIntraLinks = display links intra segments
        
        
        
        # 1) Generate karyotype data
        chr = genChr(featExp, color.blocks = color.blocks)
        ## ------ re-order blocks based on object$X names
        ## define an ordered factor to re-order blocks so with given X, block orders
        ## are always the same
        Xblocks <- unique(paste0('chr',names(object$X)))
        chr$block <- factor(chr$chrom, levels = Xblocks, ordered = TRUE)
        chr <- chr[order(chr$block),]
        chr$block <- NULL
        
        chr.names = unique(chr$chrom) # paste("chr", 1:seg.num, sep="") 
        
        # 2) Generate Links
        blocks.link <- list(...)$blocks.link
        links = genLinks(chr, simMat, threshold=cutoff, blocks.link = blocks.link)
        if (nrow(links) < 1)
            warning("Choose a lower correlation threshold to highlight
    links between datasets")
        
        # 3) Plot
        # Calculate parameters
        circleR = (figSize / 2.0) -  segmentWidth - linePlotWidth
        linksR = circleR - segmentWidth
        linePlotR = circleR + segmentWidth
        chrLabelsR = (figSize / 2.0)
        
        # replace chr$name by the ones in var.names (matching)
        # matching var.names.list with object$loadings
        ind.match = match(chr$name, unlist(sapply(object$loadings,rownames)))
        chr$name.user = unlist(var.names.list)[ind.match]
        
        # Calculate the angles of the blocks and bands
        layout = .circos_layout(chr, blocks = chr.names, centre = figSize / 2)
        
        opar1=par("mar")
        par(mar=c(2, 2, 2, 2))
        
        plot(c(1,figSize), c(1,figSize), type="n", axes=FALSE, xlab="",
             ylab="", main="")
        
        # Plot blocks, bands and their labels
        .circos_draw_blocks(layout, r = circleR, width = segmentWidth,
                            label.r = chrLabelsR, size.variables = size.variables,
                            size.labels = size.labels, line = line, ...)
        # Plot links
        if(nrow(links)>0)
            .circos_draw_links(layout, links, r = linksR, linkWidth = linkWidth,
                               showIntraLinks = showIntraLinks, color.cor = color.cor)
        
        # Plot expression values
        lineCols = color.Y
        if(line==TRUE)
        {
            ## mean expression of each feature (rows) in each level of Y (columns)
            expr = dcast(as.data.frame(featExp), formula = Features ~ Y,
                         value.var = "Mean")
            .circos_draw_expression(layout, expr, r = linePlotR,
                                    width = linePlotWidth, col = lineCols,
                                    levels = levels(Y))
        }
        opar=par("xpd")
        par(xpd=TRUE) # to authorise the legend to be written outside the margin,
        #       otherwise it's too small
        # Plot legend
        if(legend == TRUE)
        {
            # First legeng bottom left corner
            legend(x=5, y = (circleR/4), title="Correlations",
                   c("Positive Correlation", "Negative Correlation"),
                   col = color.cor, pch = 19, cex=size.legend, bty = "n")
            # Second legend bottom righ corner
            if(line==TRUE)
                legend(x=figSize-(circleR/3), y = (circleR/3), title=legend.title,
                       legend=levels(Y),  ## changed PAM50 to Y
                       col = lineCols, pch = 19, cex=size.legend, bty = "n",ncol=ncol.legend)
            # third legend top left corner
            legend(x=figSize-(circleR/2), y = figSize, title="Correlation cut-off",
                   legend=paste("r", cutoff, sep = "="),
                   col = "black", cex=size.legend, bty = "n")
            
            legend(x=-circleR/4, y = figSize, legend=paste("Comp",
                                                           paste(comp,collapse="-")),
                   col = "black", cex=size.legend, bty = "n")
        }
        par(xpd=opar,mar=opar1)# put the previous default parameter for xpd
        
        # use original feature names in similarity matrix if possible
        feat_block <- colnames(simMat)
        # remove anything from last underscore to end of the feature name
        feat <- gsub('(.*)_\\w+', '\\1', feat_block)
        if(!any(duplicated(feat)))
        {
            dimnames(simMat) <- list(feat, feat)
        } else {
            cat("\n")
            cat("adding block name to feature names in the output similarity ")
            cat("matrix as there are similar feature names across blocks.")
            cat("\n")
        }
        return(invisible(simMat))
    }
        
#' @method circosPlot block.splsda
#' @rdname circosPlot
#' @export
circosPlot.block.splsda <- .circosPlot
        
#' @method circosPlot block.plsda
#' @rdname circosPlot
#' @export
circosPlot.block.plsda <- .circosPlot

#' @method circosPlot block.spls
#' @rdname circosPlot
#' @param group The grouping factor used when \code{line = TRUE}
#' @param Y.name Character, the name of the \code{Y} block
#' @export
circosPlot.block.spls <- function(object, ..., group = NULL, Y.name = 'Y')
{
    # when a block.spls object is supplied that uses the indY parameter, the name
    # of the object$X component for this dataframe is its proper name
    # for below checks, change it back to "Y"
    if (is.null(object$X$Y)) {
        names(object$X)[object$indY] <- "Y"
    }
    
    if (length(group) != nrow(object$X$Y))
        stop("group must be a factor of length: nrow(object$X$Y) = ", nrow(object$X$Y), "\n")
    object$Y <- factor(group)
    object$keepX <- c(object$keepX, list(Y = object$keepY))
    nX <- length(object$variates)
    
    block.names <- names(object$X)
    block.names[which(block.names == 'Y')] <- Y.name
    
    names(object$X) <-
        names(object$variates) <- 
        names(object$loadings) <- 
        block.names
    .circosPlot(object, ...)
}

#' @method circosPlot block.pls
#' @rdname circosPlot
#' @export
circosPlot.block.pls <- circosPlot.block.spls

## ------------- circosPlot utils
## Positions around the circle are given as angles in degrees which, like on a
## clock, increase clockwise: 0 points to 3 o'clock and 270 to 12 o'clock.

## x and y coordinates of the points at the given angles and distances (r)
## from the centre of the plot
.circos_xy <- function(angle, r, centre)
{
    rad <- angle * pi / 180
    list(x = centre + r * cos(rad), y = centre - r * sin(rad))
}

## Angular layout of the plot. Blocks are placed clockwise from 12 o'clock in
## the order of 'blocks', each followed by a gap of 'gap' degrees, and the rest
## of the circle is shared equally between the bands (features) of all blocks.
## Returns the bands of 'chr' (output of genChr) with their start, end and mid
## angles, ordered by block and position, and the start and end angles of
## each block.
.circos_layout <- function(chr, blocks, centre, gap = 2)
{
    band.start <- as.numeric(chr$chromStart)
    band.end <- as.numeric(chr$chromEnd)
    block.size <- vapply(blocks, function(b) max(band.end[chr$chrom == b]),
                         numeric(1), USE.NAMES = FALSE)
    deg.per.unit <- (360 - gap * length(blocks)) / sum(block.size)
    block.start <- 270 + deg.per.unit * (cumsum(block.size) - block.size) +
        gap * (seq_along(blocks) - 1)
    
    block.id <- match(chr$chrom, blocks)
    bands <- chr
    bands$angle.start <- block.start[block.id] + deg.per.unit * band.start
    bands$angle.end <- block.start[block.id] + deg.per.unit * band.end
    bands$angle.mid <- (bands$angle.start + bands$angle.end) / 2
    
    list(centre = centre,
         blocks = data.frame(block = blocks,
                             angle.start = block.start,
                             angle.end = block.start + deg.per.unit * block.size),
         bands = bands[order(block.id, band.start), ])
}

## Draw a band as a thick arc of the circle
.circos_arc <- function(angle.start, angle.end, r, centre, col, lwd)
{
    ## about 5 vertices per degree give a smooth arc
    n <- max(2, ceiling(5 * abs(angle.end - angle.start)))
    p <- .circos_xy(seq(angle.end, angle.start, length.out = n), r, centre)
    lines(p$x, p$y, col = col, lwd = lwd, lend = "butt")
}

## Write a label at the given angle and distance (r) from the centre. The label
## follows the radius (or the circle if tangent = TRUE) and is flipped on the
## left half of the plot so that it never reads upside down. The anchor point
## is first moved horizontally by 'shift' towards the centre, and the text is
## then written away from it.
.circos_label <- function(angle, r, centre, label, cex, shift, tangent = FALSE)
{
    ## direction of the radius, counter-clockwise from 3 o'clock, in (0, 360]
    direction <- 360 - angle %% 360
    right <- direction <= 90 || direction > 270
    if (tangent)
        srt <- if (direction <= 180) direction - 90 else direction - 270
    else
        srt <- if (right) direction else direction - 180
    
    p <- .circos_xy(angle, r, centre)
    text(p$x + if (right) -shift else shift, p$y, labels = label, srt = srt,
         pos = if (right) 4 else 2, offset = 1, cex = cex, col = "black")
}

## Draw the bands of each block with the feature names, then the block name
.circos_draw_blocks <- function(layout, r, width, label.r,
                                size.variables, size.labels, line,
                                var.adj = -0.33, ## radial adjustment of feature names
                                block.labels.adj = 0, ## adjustment of block label names
                                ...)
{
    bands <- layout$bands
    ## block names are written closer to the circle when there are no lines
    block.shift <- if (line) width else 55 + ceiling(block.labels.adj * 55)
    
    for (b in seq_len(nrow(layout$blocks)))
    {
        block <- layout$blocks[b, ]
        for (i in which(bands$chrom == block$block))
        {
            .circos_arc(bands$angle.start[i], bands$angle.end[i], r = r,
                        centre = layout$centre, col = bands$color[i], lwd = width)
            .circos_label(bands$angle.mid[i], r = r + width * var.adj,
                          centre = layout$centre,
                          label = as.character(bands$name.user[i]),
                          cex = size.variables, shift = width)
        }
        .circos_label((block$angle.start + block$angle.end) / 2, r = label.r,
                      centre = layout$centre, label = sub("^chr", "", block$block),
                      cex = size.labels, shift = block.shift, tangent = TRUE)
    }
}

## Draw the links between correlated features (output of genLinks). Each link
## is a cubic Bezier curve whose two control points are the centre of the
## circle. Its width grows linearly with the absolute similarity, from
## linkWidth[1] for the weakest link to linkWidth[2] for a similarity of 1.
.circos_draw_links <- function(layout, links, r, linkWidth, showIntraLinks,
                               color.cor)
{
    centre <- layout$centre
    bands <- layout$bands
    from <- .circos_xy(bands$angle.mid[match(links$feat1, bands$name)], r, centre)
    to <- .circos_xy(bands$angle.mid[match(links$feat2, bands$name)], r, centre)
    
    strength <- abs(links$value)
    weakest <- min(strength)
    if (weakest < 1) {
        lwd <- linkWidth[1] +
            (strength - weakest) / (1 - weakest) * (linkWidth[2] - linkWidth[1])
    } else {
        lwd <- rep(linkWidth[2], length(strength))
    }
    ## first colour for positive, second for negative similarities
    col <- add.alpha(ifelse(links$value < 0, color.cor[2], color.cor[1]),
                     alpha = 0.4)
    
    t <- seq(0, 1, length.out = 60)
    for (i in seq_len(nrow(links)))
    {
        if (!showIntraLinks && links$chr1[i] == links$chr2[i])
            next
        ## B(t) = (1-t)^3 from + 3t(1-t) centre + t^3 to
        x <- (1 - t)^3 * from$x[i] + 3 * t * (1 - t) * centre + t^3 * to$x[i]
        y <- (1 - t)^3 * from$y[i] + 3 * t * (1 - t) * centre + t^3 * to$y[i]
        lines(x, y, col = col[i], lwd = lwd[i], lend = "butt")
    }
}

## Draw around the circle the mean expression of the features of each block,
## with one line per level of Y. 'expr' has a 'Features' column and one column
## of means per level. Each block is drawn on its own scale.
.circos_draw_expression <- function(layout, expr, r, width, col, levels)
{
    bands <- layout$bands
    lvls <- setdiff(colnames(expr), "Features")
    ## leave a fifth of the track empty on each side
    r.range <- c(r + width / 5, r + width - width / 5)
    
    for (block in layout$blocks$block)
    {
        b <- bands[bands$chrom == block, ]
        means <- as.matrix(expr[match(b$name, expr$Features), lvls, drop = FALSE])
        if (all(is.na(means)))
            next
        m.range <- range(means, na.rm = TRUE)
        for (j in seq_along(lvls))
        {
            keep <- !is.na(means[, j])
            if (sum(keep) < 2)
                next
            if (diff(m.range) > 0) {
                radius <- r.range[1] + (means[keep, j] - m.range[1]) /
                    diff(m.range) * diff(r.range)
            } else {
                radius <- rep(mean(r.range), sum(keep))
            }
            p <- .circos_xy(b$angle.mid[keep], radius, layout$centre)
            n <- length(p$x)
            segments(p$x[-n], p$y[-n], p$x[-1], p$y[-1],
                     col = col[match(lvls[j], levels)], lwd = 1, lend = "butt")
        }
    }
}

genChr =function (expr, bandWidth = 1.0, color.blocks)
{
    # Generate the segments and calculate the
    # unique positions of the bands
    #
    # Args:
    #   expr : dataframe containing the features expression
    # example: colnames(concatFeatExp) "PAM50"    "Features" "Mean"     "SD"
    # "Dataset"
    #   bandWidth: thickness of each band
    #
    # Return:
    #   a data.frame that can be used with .circos_layout
    
    # expr can contains expression data for multiple diseases
    # here, we only use the Chrom and Dataset column and remove duplicates
    keeps = c("Features","Dataset") 
    expr = expr[keeps] 
    expr = unique(expr) 
    chrLengths = data.frame(table(expr$Dataset))  
    rownames(chrLengths) = chrLengths[,1] 
    chrLengths[,1] = NULL 
    colnames(chrLengths) = c( "Freq") 
    chrLengths[, "Count"] = rep(0, nrow(chrLengths)) 
    
    # Last column contains the bands' color
    # Create a color scheme
    #dark = c("brown3","darkgoldenrod","antiquewhite3","steelblue3")
    #clear = c("brown1","darkgoldenrod1","antiquewhite1","steelblue1")
    dark.clear = color.blocks#brewer.pal(n = 12, name = 'Paired')
    dark = dark.clear[seq(2, length(dark.clear), by = 2)]
    clear = dark.clear[seq(1, length(dark.clear), by = 2)]
    chrColScheme = data.frame(dark, clear)
    n_datasets = length(unique(expr$Dataset))
    chrColScheme = chrColScheme[c(1:n_datasets),]
    rownames(chrColScheme) = levels(factor(expr$Dataset))#alphabetical order
    
    seg.out = c() 
    for (i in 1:nrow(expr)){
        chrName    = paste("chr", as.character(expr[i,'Dataset'][[1]]), sep="") 
        dType = as.character(expr[i,'Dataset'][[1]]) 
        pStart = chrLengths[dType,'Count'] * bandWidth 
        pStop = chrLengths[dType,'Count'] * bandWidth + bandWidth  
        chrLengths[dType,'Count'] = chrLengths[dType,'Count'] + 1 
        fName = as.character(as.matrix(expr[i,'Features']))
        # added as.character() by amrit
        # Assign colors
        if (chrLengths[dType,'Count'] %% 2 == 0){
            chrCol = chrColScheme[dType,]$clear 
        } else{
            chrCol = chrColScheme[dType,]$dark 
        }
        seg.out = rbind(seg.out, c(chrName, pStart, pStop,  fName, chrCol)) 
    }

    colnames(seg.out) = c("chrom", "chromStart", "chromEnd", "name", "color")
    seg.out = as.data.frame(seg.out) 
    
    return(seg.out) 
}

#' @importFrom reshape2 melt 
genLinks = function(chr, simMat, threshold, blocks.link = NULL)
{
    
    # to satisfy R CMD check that doesn't recognise x, y and group (in aes)
    Var1=Var2=chrom=NULL
    
    
    # Generates the links corresponding to pairwise correlations
    #
    # Args:
    #   chr: ideogram structure (generated from genChr)
    #   simMat: main correlation matrix
    #   threshold: minimum correlation value
    #
    # Return:
    #   a data.frame with one row per link: chr1, po1, feat1, chr2, po2,
    #   feat2 and value (the similarity)
    #
    linkList = c() 
    # Remove matrix diagonal and the the mat in a list
    linkList = subset(melt(simMat), Var1!=Var2 ) 
    # Remove links below the threshold
    linkList = subset(linkList, abs(linkList$value) >= threshold) 
    
    ## Al: if we want links only to certain block
    if (!is.null(blocks.link)) {
        only.chr <- paste0('chr', blocks.link)
        only.names <- chr[chr$chrom %in% only.chr,]$name
        linkList <- linkList[linkList$Var2 %in% only.names,]
    }
    
    #First merge
    #' @importFrom dplyr rename
    linkList = rename(linkList, feat1=Var1, feat2=Var2)
    # CHANGED BY AMRIT
    linkList = merge(linkList, chr, by.x="feat1", by.y="name") 
    # Set the position in the middle of the band
    linkList$po1 = (as.numeric(linkList$chromStart)
                    + as.numeric(linkList$chromEnd)) / 2.0
    #' @importFrom dplyr rename
    linkList = rename(linkList, chr1=chrom)   # CHANGED BY AMRIT
    keeps = c("feat1","feat2","value","chr1","po1") 
    linkList = linkList[keeps] 
    
    #Second merge
    linkList = merge(linkList, chr, by.x="feat2", by.y="name") 
    linkList$po2 = (as.numeric(linkList$chromStart)
                    + as.numeric(linkList$chromEnd)) / 2.0
    #' @importFrom dplyr rename
    linkList = rename(linkList, chr2=chrom)   # CHANGED BY AMRIT
    keeps = c("chr1","po1","feat1","chr2","po2","feat2","value") 
    linkList = linkList[keeps] 
    
    return(linkList) 
}

## Add an alpha value to a colour
add.alpha = function(col, alpha=1){
    if(missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb)/255, 2, 
          function(x) 
              rgb(x[1], x[2], x[3], alpha=alpha))  
}
###-------------------------------------------------------------------------

#' Print Methods for CCA, (s)PLS, PCA and Summary objects
#' 
#' Produce \code{print} methods for class \code{"rcc"}, \code{"pls"},
#' \code{"spls"}, \code{"pca"}, \code{"rgcca"}, \code{"sgcca"} and
#' \code{"summary"}.
#' 
#' \code{print} method for \code{"rcc"}, \code{"pls"}, \code{"spls"}
#' \code{"pca"}, \code{"rgcca"}, \code{"sgcca"} class, returns a description of
#' the \code{x} object including: the function used, the regularization
#' parameters (if \code{x} of class \code{"rcc"}), the (s)PLS algorithm used
#' (if \code{x} of class \code{"pls"} or \code{"spls"}), the samples size, the
#' number of variables selected on each of the sPLS components (if \code{x} of
#' class \code{"spls"}) and the available components of the object.
#' 
#' \code{print} method for \code{"summary"} class, gives the (s)PLS algorithm
#' used (if \code{x} of class \code{"pls"} or \code{"spls"}), the number of
#' variates considered, the canonical correlations (if \code{x} of class
#' \code{"rcc"}), the number of variables selected on each of the sPLS
#' components (if \code{x} of class \code{"spls"}) and the available components
#' for Communalities Analysis, Redundancy Analysis and Variable Importance in
#' the Projection (VIP).
#' 
#' @aliases print print.rcc print.mixo_pls print.mixo_spls print.summary
#' print.pca print.spca print.rgcca print.sgcca
#' 
#' @param x object of class inherited from \code{"rcc"}, \code{"pls"},
#' \code{"spls"}, \code{"pca"}, \code{"spca"}, \code{"rgcca"}, \code{"sgcca"} or
#' \code{"summary"}.
#' @param ... not used currently.
#' 
#' @return none
#' @author Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao, Fangzhou Yao, Jeff Coquery, Al J Abadi.
#' @seealso \code{\link{rcc}}, \code{\link{pls}}, \code{\link{spls}},
#' \code{\link{vip}}.
#' @keywords regression multivariate
#' @example ./examples/print-examples.R
#------------------ print method for pls ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mixo_pls <-
    function(x, ...)
    {
        
        mode = paste("'", x$mode, "'", sep = "")
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" PLS with a", mode, "mode with", x$ncomp, "PLS components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
        
        cat(" No variable selection. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        cat(" loading vectors: see object$loadings \n")
        cat(" variates: see object$variates \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
        
    }
#------------------ print method for mint.pls ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mint.pls <-
    function(x, ...)
    {
        
        mode = paste("'", x$mode, "'", sep = "")
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" MINT PLS with a", mode, "mode with", x$ncomp, "PLS components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
        cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")
        
        cat(" No variable selection. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" loading vectors: see object$loadings \n")
        cat(" loading vectors per study: see object$loadings.partial \n")
        cat(" variates: see object$variates \n")
        cat(" variates per study: see object$variates.partial \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
    }

#------------------ print method for plsda ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mixo_plsda <-
    function(x, ...)
    {
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" PLS-DA (regression mode) with", x$ncomp, "PLS-DA components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")
        
        cat(" No variable selection. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" loading vectors: see object$loadings \n")
        cat(" variates: see object$variates \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow, cim \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
        cat("\n")
        cat(" Other functions: \n", "-------------------- \n")
        cat(" auc \n")
        
    }

#------------------ print method for mint.plsda ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mint.plsda <-
    function(x, ...)
    {
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" MINT PLS-DA (regression mode) with", x$ncomp, "PLS-DA components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")
        cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")
        
        cat(" No variable selection. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" loading vectors: see object$loadings \n")
        cat(" loading vectors per study: see object$loadings.partial \n")
        cat(" variates: see object$variates \n")
        cat(" variates per study: see object$variates.partial \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow, cim \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
        cat("\n")
        cat(" Other functions: \n", "-------------------- \n")
        cat(" perf, auc\n")
    }

#----------------- print method for spls ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mixo_spls <-
    function(x, ...)
    {
        
        mode = paste("'", x$mode, "'", sep = "")
        keepX = paste("[", x$keepX, "]", sep = "")
        keepY = paste("[", x$keepY, "]", sep = "")
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" sPLS with a", mode, "mode with", x$ncomp, "sPLS components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
        
        cat(" Selection of", keepX, "variables on each of the sPLS components on the X data set. \n")
        cat(" Selection of", keepY, "variables on each of the sPLS components on the Y data set. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" loading vectors: see object$loadings \n")
        cat(" variates: see object$variates \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
    }

#----------------- print method for mint.spls ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mint.spls <-
    function(x, ...)
    {
        
        mode = paste("'", x$mode, "'", sep = "")
        keepX = paste("[", x$keepX, "]", sep = "")
        keepY = paste("[", x$keepY, "]", sep = "")
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" MINT sPLS with a", mode, "mode with", x$ncomp, "sPLS components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y of dimensions:", nrow(x$Y), ncol(x$Y), "\n\n")
        cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")
        
        cat(" Selection of", keepX, "variables on each of the sPLS components on the X data set. \n")
        cat(" Selection of", keepY, "variables on each of the sPLS components on the Y data set. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" loading vectors: see object$loadings \n")
        cat(" loading vectors per study: see object$loadings.partial \n")
        cat(" variates: see object$variates \n")
        cat(" variates per study: see object$variates.partial \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
        cat("\n")
        cat(" Other functions: \n", "-------------------- \n")
        cat(" selectVar\n")
    }


#----------------- print method for splsda ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mixo_splsda <-
    function(x, ...)
    {
        
        keepX = paste("[", x$keepX, "]", sep = "")
        keepY = paste("[", x$keepY, "]", sep = "")
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        cat(" sPLS-DA (regression mode) with", x$ncomp, "sPLS-DA components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")
        
        cat(" Selection of", keepX, "variables on each of the sPLS-DA components on the X data set. \n")
        cat(" No Y variables can be selected. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" loading vectors: see object$loadings \n")
        cat(" variates: see object$variates \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow, cim \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
        cat("\n")
        cat(" Other functions: \n", "-------------------- \n")
        cat(" selectVar, tune, perf, auc \n")
    }

#----------------- print method for mint.splsda ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.mint.splsda <-
    function(x, ...)
    {
        
        keepX = paste("[", x$keepX, "]", sep = "")
        keepY = paste("[", x$keepY, "]", sep = "")
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        cat(" MINT sPLS-DA (regression mode) with", x$ncomp, "sPLS-DA components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y with", ncol(x$ind.mat) , "classes. \n\n")
        cat(" You entered a grouping factor with", nlevels(x$study), "studies. \n")
        
        cat(" Selection of", keepX, "variables on each of the sPLS-DA components on the X data set. \n")
        cat(" No Y variables can be selected. \n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" loading vectors: see object$loadings \n")
        cat(" loading vectors per study: see object$loadings.partial \n")
        cat(" variates: see object$variates \n")
        cat(" variates per study: see object$variates.partial \n")
        cat(" variable names: see object$names \n")
        
        cat("\n")
        cat(" Functions to visualise samples: \n", "-------------------- \n")
        cat(" plotIndiv, plotArrow, cim \n")
        cat("\n")
        cat(" Functions to visualise variables: \n", "-------------------- \n")
        cat(" plotVar, plotLoadings, network, cim \n")
        cat("\n")
        cat(" Other functions: \n", "-------------------- \n")
        cat(" selectVar, tune, perf, auc \n")
        
    }

#------------------ print method for rcc ------------------#
#' @name print
#' @rdname S3methods-print
#' @export
print.rcc <-
    function(x, ...)
    {
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" rCCA with", x$ncomp, "components and regularization parameters", x$lambda[1], "and", x$lambda[2], "for the X and Y data. \n")
        cat(" You entered data X of dimensions :", nrow(x$X), ncol(x$X), "\n")
        cat(" You entered data Y of dimensions :", nrow(x$Y), ncol(x$Y), "\n\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" canonical correlations: see object$cor \n")
        cat(" loading vectors: see object$loadings \n")
        cat(" variates: see object$variates \n")
        cat(" variable names: see object$names \n")
    }

# ------------------------ print for (s)pca --------------------------------
# #' @name print
#' @rdname S3methods-print
#' @export
print.pca <- function(x, ...)
{
    
    ind.show = min(10, x$ncomp)
    
    if (is(x, 'spca'))
    {
        cat("sparse PCA with", x$ncomp, "principal components. \n")
        cat("  Input data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        cat("  Number of selected variables on each prinicipal components:\n")
        print(x$keepX, print.gap = 3)
    }
    else
    {
        x$sdev=as.vector(x$sdev)
        names(x$sdev) = paste("PC", 1:length(x$sdev), sep = "")
        cat("  Eigenvalues for the first", ind.show, "principal components, see object$sdev^2:", "\n")
        print((x$sdev[1:ind.show])^2)
        cat("  \n") 
    }

    per.var = x$prop_expl_var$X
    cum.var = x$cum.var

    names(per.var) = paste("PC", 1:length(per.var), sep = "")
    names(cum.var) = paste("PC", 1:length(cum.var), sep = "")
    
    var.type <- ifelse(is(x, 'spca'), 'adjusted', '')
    
    cat("  Proportion of", var.type, "explained variance for the first", ind.show,
        "principal components, see object$prop_expl_var:", "\n")
    print(per.var[1:ind.show], print.gap = 6)
    cat("  \n")
    
    if (!is(x, 'spca'))
    {
        cat("  Cumulative proportion of", var.type, "explained variance for the first", ind.show, "principal components, see object$cum.var:", "\n")
        print(cum.var[1:ind.show], print.gap = 6)
        cat("  \n")
        
    }
    cat("  Other available components: \n", "-------------------- \n")
    cat("  loading vectors: see object$rotation \n")
    
    cat("  Other functions: \n", "-------------------- \n")
    if (is(x, 'spca'))
    {
        cat("  tune.spca, plotIndiv, plot, plotVar, selectVar, biplot\n")
    }
    else
    {
        cat("  plotIndiv, plot, plotVar, selectVar, biplot\n")
    }
    
}

# ------------------------ print for ipca -------------------------
# #' @name print
#' @rdname S3methods-print
#' @export
print.ipca <-
    function(x, ...)
    {
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" IPCA with", x$ncomp, "independent components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" unmixing matrix: see object$unmixing \n")
        cat(" independent principal components: see object$variates$X \n")
        cat(" mixing matrix: see object$mixing \n")
        cat(" kurtosis: see object$kurtosis \n")
        cat(" variable names: see object$names \n")
        cat(" independent loading vectors: see object$loadings \n")
    }

# ------------------------ print for sipca -------------------------
# #' @name print
#' @rdname S3methods-print
#' @export
print.sipca <-
    function(x, ...)
    {
        
        cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
        
        cat(" Sparse IPCA with", x$ncomp, "independent components. \n")
        cat(" You entered data X of dimensions:", nrow(x$X), ncol(x$X), "\n")
        
        cat(" Selection of", x$keepX, "variables on each of the principal components on the X data set. \n")
        
        cat(" Main numerical outputs: \n",
            "-------------------- \n")
        
        cat(" unmixing matrix: see object$unmixing \n")
        cat(" independent principal components: see object$variates$X \n")
        cat(" mxing matrix: see object$mixing \n")
        cat(" kurtosis: see object$kurtosis \n")
        cat(" variable names: see object$names \n")
        cat(" independent loading vectors: see object$loadings \n")
        
    }

# ------------------------ print for rgcca -------------------------
# #' @name print
#' @rdname S3methods-print
#' @export
print.rgcca <- function(x, ...)
{
    
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    
    # components
    for(k in 1:length(x$X)){
        cat(" rGCCA with", x$ncomp[k], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat("\n")
    
    # dimension
    for(k in 1 : length(x$X)){
        cat(" Dimension of block", k, 'is ', dim(x$X[[k]]), "\n")
    }
    cat("\n")
    cat(" Main numerical outputs: \n", "-------------------- \n")
    
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
    
    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar\n")
    
}


# ------------------------ print for sgcca -------------------------
# #' @name print
#' @rdname S3methods-print
#' @export
print.sgcca<- function(x, ...)
{
    
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    
    # components
    for(k in 1 : length(x$X)){
        cat(" sGCCA with", x$ncomp[[k]], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat("\n")
    
    # dimension
    for(k in 1 : length(x$X)){
        cat(" Dimension of block", k, 'is ', dim(x$X[[k]]), "\n")
    }
    cat("\n")
    
    # selected variables
    list.select = list()
    for(k in 1:length(x$X)){
        list.select[[k]] = apply(x$loadings[[k]], 2, function(x){sum(x!=0)})
        cat(" Selection of", list.select[[k]], "variables on each of the sGCCA components on the block", k, "\n")
    }
    cat("\n")
    cat(" Main numerical outputs: \n", "-------------------- \n")
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
    
    
    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network\n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar \n")
    
}


# ------------------------ print for sgcca -------------------------
# #' @name print
#' @rdname S3methods-print
#' @export
print.sgccda<- function(x, ...)
{
    
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    
    # components
    for(k in 1 : length(x$X)){
        cat(" sGCCA with", x$ncomp[[k]], "components on block", k, "named", x$names$blocks[k], "\n")
    }
    cat(" sGCCA with", x$ncomp[x$indY], "components on the outcome Y\n")
    cat("\n")
    
    
    # dimension
    for(k in 1 : length(x$X)){
        cat(" Dimension of block", k, 'is ', dim(x$X[[k]]), "\n")
    }
    cat(" Outcome Y has", nlevels(x$Y), "levels \n")
    cat("\n")
    
    # selected variables
    list.select = list()
    for(k in 1:length(x$X)){
        list.select[[k]] = apply(x$loadings[[k]], 2, function(x){sum(x!=0)})
        cat(" Selection of", list.select[[k]], "variables on each of the sGCCA components on the block", k, "\n")
    }
    cat("\n")
    cat(" Main numerical outputs: \n", "-------------------- \n")
    cat(" loading vectors: see object$loadings \n")
    cat(" variates: see object$variates \n")
    cat(" variable names: see object$names \n")
    
    
    cat("\n")
    cat(" Functions to visualise samples: \n", "-------------------- \n")
    cat(" plotIndiv, plotArrow, cimDiablo, plotDiablo \n")
    cat("\n")
    cat(" Functions to visualise variables: \n", "-------------------- \n")
    cat(" plotVar, plotLoadings, network, circosPlot \n")
    cat("\n")
    cat(" Other functions: \n", "-------------------- \n")
    cat(" selectVar, perf, auc \n")
    
}


#------- print for summary with (s)PLS object or rcc ---------#
##' @name print
#' @rdname S3methods-print
#' @export
print.summary <-
    function(x, ...)
    {
        
        print.gap = 4
        what = x$what
        digits = x$digits
        method = x$method
        
        #--------------------- output pls/spls ---------------------#
        if(isTRUE(method %in% c("pls", "spls"))){
            
            if (identical(method, "pls"))
            {
                cat(" PLS mode:", x$mode)
                cat("\n Number of variates considered:", x$ncomp, "\n")
            } else {
                cat(" sPLS mode:", x$mode)
                cat("\n Number of variates considered:", x$ncomp)
                cat("\n Number of X-variables selected on each of the sPLS components:", x$keepX)
                cat("\n Number of Y-variables selected on each of the sPLS components:", x$keepY, "\n")
            }
            
            #---------- affichage communaute ----------#
            if (any(what == "all") || any(what == "communalities"))
            {
                cat("\n\n Communalities Analysis:\n",
                    "----------------------")
                
                cat("\n X-Variables vs their own Variates: see object$CM.X$own \n")
                cat("\n X-Variables vs the opposite Variates: see object$CM.X$opp \n")
                cat("\n Y-Variables vs their own Variates: see object$CM.Y$opp \n")
                cat("\n Y-Variables vs the opposite Variates: see object$CM.Y$opp \n")
            }
            
            #--------- affichage redondance -----------#
            if (any(what == "all") || any(what == "redundancy"))
            {
                cat("\n\n Redundancy Analysis:\n",
                    "-------------------\n")
                
                cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
                cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
                cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
                cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
            }
            
            #---------- tableau VIP ---------#
            if (any(what == "all") || any(what == "VIP"))
            {
                cat("\n\n", "Variable Importance in the Projection (VIP): see object$VIP \n",
                    "------------------------------------------- \n\n")
            }
            
        }  #end if pls
        
        # ---------------------- output rcc ------------------------#
        if(identical(method, "rcc"))
        {
            print.gap = 4
            if (any(what == "all"))
            {
                cat(" Number of canonical variates considered:", x$ncomp, "\n")
                cat("\n Canonical correlations:",
                    "\n ----------------------\n")
                print(round(x$can.cor, digits = digits), print.gap = print.gap)
            }
            
            #-- affichage communaute --#
            if (any(what == "all") || any(what == "communalities"))
            {
                cat("\n\n Canonical Communalities Analysis:\n",
                    "--------------------------------")
                
                cat("\n X-Variables vs their own Canonical Variates: see object$Cm.X$own \n")
                cat("\n X-Variables vs the opposite Canonical Variates: see object$Cm.X$opp \n")
                cat("\n Y-Variables vs their own Canonical Variates: see object$Cm.Y$own \n")
                cat("\n Y-Variables vs the opposite Canonical Variates: see object$Cm.Y$opp \n")
            }
            
            #--------- affichage redondance -----------#
            if (any(what == "all") || any(what == "redundancy"))
            {
                cat("\n\n Redundancy Analysis:\n",
                    "-------------------\n")
                
                cat("\n X-Variables vs their own Variates: see object$Rd.X$own \n")
                cat("\n X-Variables vs the opposite Variates: see object$Rd.X$opp \n")
                cat("\n Y-Variables vs their own Variates: see object$Rd.Y$opp \n")
                cat("\n Y-Variables vs the opposite Variates: see object$Rd.Y$opp \n")
            }
            
        }  #end rcc
    }


# perf.diablo / sgccda.mthd
# perf.splsda = perf.plsda / splsda.mthd plsda.mthd
# perf.spls  = perf.pls / spls.mthd pls.mthd
#' @name print
#' @rdname S3methods-print
#' @export
print.perf.pls.mthd <- function(x, ...)
{
    spls.model <- length(x$features) > 0
    cat("\n")
    cat("Call:\n", deparse(x$call, width.cutoff = 500))
    cat("\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Cross-validated model performance measures:\n")
    cat(sprintf(" %s\n\n", paste0(names(x$measures), collapse = ', ')))
    if ('features' %in% names(x))
    if (spls.model)
    {
        cat(" Stable features selected on each component: see object$features \n\n")
    }
    cat(" See the help file ?perf \n")
    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
}

#' @name print
#' @rdname S3methods-print
#' @export
print.perf.plsda.mthd = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Error rate (overall or BER) for each component and for each distance: see object$error.rate \n")
    cat(" Error rate per class, for each component and for each distance: see object$error.rate.class \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" AUC values: see object$auc if auc = TRUE \n\n")
    
    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
    
}
#' @name print
#' @rdname S3methods-print
#' @export
print.perf.splsda.mthd = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Error rate (overall or BER) for each component and for each distance (averaged over the nrepeat): see object$error.rate \n")
    cat(" Error rate (overall or BER) for each component, for each distance and for each repeat: see object$error.rate.all \n")
    cat(" Error rate per class, for each component and for each distance: see object$error.rate.class \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" Stable features on each component: see object$features$stable \n")
    cat(" AUC values: see object$auc if auc = TRUE \n\n")
    
    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
    
}
#' @name print
#' @rdname S3methods-print
#' @export
print.perf.mint.splsda.mthd = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Study-specific error rate (overall, BER and error rate per class) for each component and for each distance: see object$study.specific.error \n")
    cat(" Global error rate (overall, BER and error rate per class) for each component and for each distance: see object$global.error \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" AUC values: see object$auc \n")
    cat(" AUC values per study: see object$auc.study if auc = TRUE \n\n")
    
    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
    
}
#' @name print
#' @rdname S3methods-print
#' @export
print.perf.sgccda.mthd = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Error rate (overall or BER) for each component and for each distance: see object$error.rate \n")
    cat(" Error rate per class, for each component and for each distance: see object$error.rate.class \n")
    cat(" Prediction values for each component: see object$predict \n")
    cat(" Classification of each sample, for each component and for each distance: see object$class \n")
    cat(" Stable features on each component: see object$features$stable \n")
    cat(" AUC values: see object$auc if auc = TRUE \n\n")
    
    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
    
}


#' @name print
#' @rdname S3methods-print
#' @export
print.tune.pca = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" for all principal components, see object$sdev, object$prop_expl_var and object$cum.var\n")
}

#' @name print
#' @rdname S3methods-print
#' @export
print.tune.spca = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Optimal keepX for each component, see object$choice.keepX \n\n")
    cat(" Visualisation functions: \n", "-------------------- \n")
    cat(" plot \n")
    
    
}

#' @name print
#' @rdname S3methods-print
#' @export
print.tune.rcc = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat("  lambda1 = ", x$opt.lambda1, ", see object$opt.lambda1\n", " lambda2 = ", x$opt.lambda2, ",  see object$opt.lambda2\n",
        "CV-score = ", x$opt.score, ", see object$opt.score\n")
}
#' @name print
#' @rdname S3methods-print
#' @export
print.tune.splsda = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Optimal keepX for each component, see object$choice.keepX \n")
    cat(" Optimal number of components: see object$choice.ncomp \n")
    cat(" Error rate for each tested keepX and for each component (averaged over the nrepeat, mean and standard deviation): see object$error.rate and object$error.rate.sd  \n")
    cat(" Error rate for each tested keepX, for each component and for each repeat: see object$error.rate.all \n")
    cat(" Error rate per class obtained with the optimal keepX, for each component and for each nrepeat: see object$error.rate.class \n")
    cat(" AUC when applicable, (note that those results may differ with sPLS-DA prediction, see details): see object$AUC \n\n")
    
    cat(" Other outputs are available, and details on those outputs in ?tune.splsda.  \n\n")
    
    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
}

#' @method print tune.pls
#' @rdname S3methods-print
#' @export
print.tune.pls = function(x, ...)
{
    cat("Call:\n", deparse(x$call, width.cutoff = 500), "\n")
    cat("\n  Main numerical outputs: \n",
        "-------------------- \n")
    cat("\n Optimal number of components based on Q2.total (regard with care): see object$choice.ncomp")
    meas <- switch(x$call$measure, cor = 'correlations', RSS = 'Residual Sum of Squares')
    if (is(x, 'tune.spls'))
        cat(sprintf("\n Optimal number of features selected based on %s from cross-validated components: see object$choice.keepX & object$choice.keepY", x$call$measure))
    cat("\n See object$measure.pred for prediction measures for X (value.t) and Y (value.u)\n")
    
    cat("\n  Visualisation Functions: \n", "--------------------")
    cat("\n  plot \n")
}

#' @method print tune.spls1
#' @rdname S3methods-print
#' @export
print.tune.spls1 = function(x, ...)
{
    cat("Call:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat("\n  Main numerical outputs: \n",
        "-------------------- \n")
    
    if (is(x, 'tune.spls1') & x$call$nrepeat > 2)
    {
        cat(sprintf("\n  Optimal keep%s for each component based on the provided tune measure, see object$choice.keep%s", c('X', 'Y'), c('X', 'Y')))
    }
    cat("\n  Optimal number of components: see object$choice.ncomp \n")
    
    cat("\n  Visualisation Functions: \n", "--------------------")
    cat("\n  plot \n")
}


#' @name print
#' @rdname S3methods-print
#' @export
print.tune.mint.splsda = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Optimal keepX for each component, see object$choice.keepX \n")
    cat(" Error rate for each tested keepX and for each component: see object$error.rate \n")
    cat(" Error rate per class obtained with the optimal keepX, for each component: see object$error.rate.class \n\n")
    
    cat(" Other outputs available, see ?tune.mint.splsda \n\n")
    
    cat(" Visualisation Functions: \n", "-------------------- \n")
    cat(" plot \n")
}
#' @name print
#' @rdname S3methods-print
#' @export
print.tune.block.splsda = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    cat(" Optimal keepX for each component, see object$choice.keepX \n")
    cat(" Error rate for each tested keepX and for each component (averaged over the nrepeat): see object$error.rate \n")
    cat(" Error rate for each tested keepX, for each component and for each repeat: see object$error.rate.all \n")
    cat(" Error rate per class obtained with the optimal keepX, for each component and for each nrepeat: see object$error.rate.class \n\n")
    
    cat(" Other outputs available, see ?tune.splsda.  \n\n")
    
}
#' @name print
#' @rdname S3methods-print
#' @export
print.predict = function(x, ...)
{
    cat("\nCall:\n", deparse(x$call, width.cutoff = 500), "\n\n")
    cat(" Main numerical outputs: \n",
        "-------------------- \n")
    if(is.list(x$predict)) #block analysis
    {
        cat(" Prediction values of the test samples for each block and each component: see object$predict \n")
        cat(" variates of the test samples for each block and each component: see object$variates \n")
        
        
        if(!is.null(x$dist)) #DA object, more outputs
        {
            cat(" Classification of the test samples for each distance, for each block and each component: see object$class \n")
            cat(" Majority vote of the test samples for each distance and each component: see object$vote \n")
        }
    }else{
        cat(" Prediction values of the test samples for each component: see object$predict \n")
        cat(" variates of the test samples: see object$variates \n")
        
        if(!is.null(x$dist)) #DA object, more outputs
            cat(" Classification of the test samples for each distance and for each component: see object$class \n")
    }
    
}

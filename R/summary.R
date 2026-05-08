#--------------------------------------------------------#
#-- Includes summary.pls, summary.spls and summary.rcc --#
#--------------------------------------------------------#

#--------------------- PLS and sPLS ---------------------#
#' Summary Methods for CCA and PLS objects
#' 
#' Produce \code{summary} methods for class \code{"rcc"}, \code{"pls"} and
#' \code{"spls"}.
#' 
#' The information in the \code{rcc}, \code{pls} or \code{spls} object is
#' summarised, it includes: the dimensions of \code{X} and \code{Y} data, the
#' number of variates considered, the canonical correlations (if \code{object}
#' of class \code{"rcc"}) and the (s)PLS algorithm used (if \code{object} of
#' class \code{"pls"} or \code{"spls"}) and the number of variables selected on
#' each of the sPLS components (if \code{x} of class \code{"spls"}).
#' 
#' \code{"communalities"} in \code{what} gives Communalities Analysis.
#' \code{"redundancy"} display Redundancy Analysis. \code{"VIP"} gives the
#' Variable Importance in the Projection (VIP) coefficients fit by \code{pls}
#' or \code{spls}. If \code{what} is \code{"all"}, all are given.
#' 
#' For class \code{"rcc"}, when a value to \code{cutoff} is specified, the
#' correlations between each variable and the equiangular vector between
#' \eqn{X}- and \eqn{Y}-variates are computed. Variables with at least one
#' correlation componente bigger than \code{cutoff} are showed. The defaults is
#' \code{cutoff=NULL} all the variables are given.
#' 
#' @param object object of class inherited from \code{"rcc"}, \code{"pls"} or
#' \code{"spls"}.
#' @param cutoff real between 0 and 1. Variables with all correlations
#' components below this cut-off in absolute value are not showed (see Details).
#' @param digits integer, the number of significant digits to use when
#' printing. Defaults to \code{4}.
#' @param what character string or vector. Should be a subset of
#' \code{c("all"}, \code{"summarised"}, \code{"communalities"},
#' \code{"redundancy"}, \code{"VIP"}). \code{"VIP"} is only available for
#' (s)PLS. See Details.
#' @param keep.var Logical. If \code{TRUE} only the variables with loadings not
#' zero (as selected by \code{spls}) are showed. Defaults to \code{FALSE}.
#' @param ... not used currently.
#' @return The function \code{summary} returns a list with components:
#' \item{ncomp}{the number of components in the model.} \item{cor}{the
#' canonical correlations.} \item{cutoff}{the cutoff used.}
#' \item{keep.var}{list containing the name of the variables selected.}
#' \item{mode}{the algorithm used in \code{pls} or \code{spls}.} \item{Cm}{list
#' containing the communalities.} \item{Rd}{list containing the redundancy.}
#' \item{VIP}{matrix of VIP coefficients.} \item{what}{subset of
#' \code{c("all"}, \code{"communalities"}, \code{"redundancy"}, \code{"VIP"}).}
#' \item{digits}{the number of significant digits to use when printing.}
#' \item{method}{method used: \code{rcc}, \code{pls} or \code{spls}.}
#' @author Sébastien Déjean, Ignacio González, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{rcc}}, \code{\link{pls}}, \code{\link{spls}},
#' \code{\link{vip}}.
#' @keywords regression multivariate
#' @examples
#' ## summary for objects of class 'rcc'
#' data(nutrimouse)
#' X <- nutrimouse$lipid
#' Y <- nutrimouse$gene
#' nutri.res <- rcc(X, Y, ncomp = 3, lambda1 = 0.064, lambda2 = 0.008)
#' more <- summary(nutri.res, cutoff = 0.65)
#' 
#' \dontrun{
#' ## summary for objects of class 'pls'
#' data(linnerud)
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.pls <- pls(X, Y)
#' more <- summary(linn.pls)
#' 
#' ## summary for objects of class 'spls'
#' data(liver.toxicity)
#' X <- liver.toxicity$gene
#' Y <- liver.toxicity$clinic
#' toxicity.spls <- spls(X, Y, ncomp = 3, keepX = c(50, 50, 50),
#' keepY = c(10, 10, 10))
#' more <- summary(toxicity.spls, what = "redundancy", keep.var = TRUE)
#' }
#' 
#' @name summary
NULL
## -------------------------------- (s)PLS -------------------------------- ##
#' @rdname summary
#' @method summary mixo_pls
#' @export
summary.mixo_pls <-
    function(object,
             what = c("all", "communalities", "redundancy", "VIP"),
             digits = 4,
             keep.var = FALSE,
             ...)
    {
        
        #-- initialisation des matrices --#
        ncomp = object$ncomp
        n = nrow(object$X)
        p = ncol(object$X)
        q = ncol(object$Y)
        result = list( )
        
        if (keep.var) {
            gp.X = apply(object$loadings$X, 1, sum) != 0
            gp.Y = apply(object$loadings$Y, 1, sum) != 0
        }
        else {
            gp.X = 1:p
            gp.Y = 1:q
        }
        
        #-- communalities --#
        if (any(what == "all") || any(what == "communalities")) {
            # communalities inter-group #
            #-------------------------#
            Cm.XvsV = cor(object$X[, gp.X], object$variate$Y[, 1:ncomp], use = "pairwise")
            Cm.XvsV = t(apply(Cm.XvsV^2, 1, cumsum))
            colnames(Cm.XvsV) = paste("comp", 1:ncomp)
            
            Cm.YvsU = cor(object$Y[, gp.Y], object$variate$X[, 1:ncomp], use = "pairwise")
            Cm.YvsU = t(apply(Cm.YvsU^2, 1, cumsum))
            colnames(Cm.YvsU) = paste("comp", 1:ncomp)
            
            # communalities intra-group #
            #-------------------------#
            Cm.XvsU = cor(object$X[, gp.X], object$variates$X[, 1:ncomp], use = "pairwise")
            Cm.XvsU = t(apply(Cm.XvsU^2, 1, cumsum))
            colnames(Cm.XvsU) = paste("comp", 1:ncomp)
            
            Cm.YvsV = cor(object$Y[, gp.Y], object$variates$Y[, 1:ncomp], use = "pairwise")
            Cm.YvsV = t(apply(Cm.YvsV^2, 1, cumsum))
            colnames(Cm.YvsV) = paste("comp", 1:ncomp)
            
            result$Cm.X = list(own = Cm.XvsU, opp = Cm.XvsV)
            result$Cm.Y = list(own = Cm.YvsV, opp = Cm.YvsU)
        }
        
        #-- redundancy --#
        if (any(what == "all") || any(what == "redundancy")) {
            Rd.XvsU = cor(object$X[, gp.X], object$variates$X[, 1:ncomp], use = "pairwise")
            Rd.XvsU = apply(Rd.XvsU^2, 2, sum)/p
            Rd.XvsV = cor(object$X[, gp.X], object$variates$Y[, 1:ncomp], use = "pairwise")
            Rd.XvsV = apply(Rd.XvsV^2, 2, sum)/p
            
            Rd.YvsU = cor(object$Y[, gp.Y], object$variates$X[, 1:ncomp], use = "pairwise")
            Rd.YvsU = apply(Rd.YvsU^2, 2, sum)/q
            Rd.YvsV = cor(object$Y[, gp.Y], object$variates$Y[, 1:ncomp], use = "pairwise")
            Rd.YvsV = apply(Rd.YvsV^2, 2, sum)/q
            
            own = cbind(Rd.XvsU, cumsum(Rd.XvsU))
            colnames(own) = c("Proportion", "Cumulative")
            opp = cbind(Rd.XvsV, cumsum(Rd.XvsV))
            colnames(opp) = c("Proportion", "Cumulative")
            result$Rd.X = list(own = own, opp = opp)
            
            own = cbind(Rd.YvsV, cumsum(Rd.YvsV))
            colnames(own) = c("Proportion", "Cumulative")
            opp = cbind(Rd.YvsU, cumsum(Rd.YvsU))
            colnames(opp) = c("Proportion", "Cumulative")
            result$Rd.Y = list(own = own, opp = opp)
        }
        
        #-- affichage --#
        if (any(what == "all") || any(what == "communalities") || 
            any(what == "redundancy") || any(what == "VIP")) {
            result$ncomp = ncomp
            result$mode = object$mode
            result$keep.var = list(X = colnames(object$X)[gp.X], Y = colnames(object$Y)[gp.Y])
        }
        
        #-- tableau VIP --#
        if (any(what == "all") || any(what == "VIP")) {
            VIP = vip(object)[gp.X, ]
            result$VIP = VIP
        }
        
        #-- valeurs sortantes --#
        result$what = what
        result$digits = digits
        if(is(object, "mixo_pls")) {
            result$method = 'pls'
        } 
        else {
            result$method = 'spls'
            result$keepX = object$keepX
            result$keepY = object$keepY
        } 
        
        class(result) = "summary"
        return(invisible(result))
    }

#' @rdname summary
#' @method summary mixo_spls
#' @export
summary.mixo_spls <- summary.mixo_pls

## --------------------------------- RCC ---------------------------------- ##
#' @rdname summary
#' @method summary rcc
#' @export
summary.rcc <-
    function(object,
             what = c("all", "communalities", "redundancy"),
             cutoff = NULL,
             digits = 4,
             ...)
        
    {
        
        #-- initialisation des matrices --#
        ncomp = object$ncomp
        p = ncol(object$X)
        q = ncol(object$Y)
        n = nrow(object$X)
        
        bisect = object$variates$X[, 1:ncomp] + object$variates$Y[, 1:ncomp]
        cord.X = cor(object$X, bisect, use = "pairwise")
        cord.Y = cor(object$Y, bisect, use = "pairwise")
        gp.X = 1:p
        gp.Y = 1:q
        
        if (!is.null(cutoff)) {
            # choix des variables avec au moins une coordonn?e # 
            # sup?rieur au cutoff                              #
            #--------------------------------------------------#
            gp.X = vector(mode = "numeric")
            gp.Y = vector(mode = "numeric")
            
            k = 1
            for (i in 1:p) {  
                if (any(abs(cord.X[i, ]) > cutoff)) { 
                    gp.X[k] = i
                    k = k + 1
                }
            }
            
            k = 1
            for (i in 1:q) {  
                if (any(abs(cord.Y[i, ]) > cutoff)) { 
                    gp.Y[k] = i
                    k = k + 1
                }
            }
        }
        
        if (length(gp.X) == 0 || length(gp.Y) == 0) 
            stop("Cutoff value very high for the components 1:ncomp.
             No variable was selected.")
        
        result = list( )
        
        #-- communaut? --#
        if (any(what == "all") || any(what == "communalities")) {
            # communaut? inter-groupe #
            #-------------------------#
            Cm.XvsV = cor(object$X[, gp.X], object$variate$Y[, 1:ncomp], use = "pairwise")
            Cm.XvsV = t(apply(Cm.XvsV^2, 1, cumsum))
            colnames(Cm.XvsV) = paste("comp", 1:ncomp)
            
            Cm.YvsU = cor(object$Y[, gp.Y], object$variate$X[, 1:ncomp], use = "pairwise")
            Cm.YvsU = t(apply(Cm.YvsU^2, 1, cumsum))
            colnames(Cm.YvsU) = paste("comp", 1:ncomp)
            
            # communaut? intra-groupe #
            #-------------------------#
            Cm.XvsU = cor(object$X[, gp.X], object$variates$X[, 1:ncomp], use = "pairwise")
            Cm.XvsU = scale(Cm.XvsU^2, center = FALSE, scale = (object$cor[1:ncomp])^2)
            Cm.XvsU = t(apply(Cm.XvsU, 1, cumsum))
            colnames(Cm.XvsU) = paste("comp", 1:ncomp)
            
            Cm.YvsV = cor(object$Y[, gp.Y], object$variates$Y[, 1:ncomp], use = "pairwise")
            Cm.YvsV = scale(Cm.YvsV^2, center = FALSE, scale = (object$cor[1:ncomp])^2)
            Cm.YvsV = t(apply(Cm.YvsV, 1, cumsum))
            colnames(Cm.YvsV) = paste("comp", 1:ncomp)
            
            result$Cm.X = list(own = Cm.XvsU, opp = Cm.XvsV)
            result$Cm.Y = list(own = Cm.YvsV, opp = Cm.YvsU)
        }
        
        #-- redondance --#
        if (any(what == "all") || any(what == "redundancy")) {
            Rd.XvsU = cor(object$X, object$variates$X[, 1:ncomp], use = "pairwise")
            Rd.XvsU = apply(Rd.XvsU^2, 2, sum)/p
            Rd.XvsV = cor(object$X, object$variates$Y[, 1:ncomp], use = "pairwise")
            Rd.XvsV = apply(Rd.XvsV^2, 2, sum)/p
            
            Rd.YvsU = cor(object$Y, object$variates$X[, 1:ncomp], use = "pairwise")
            Rd.YvsU = apply(Rd.YvsU^2, 2, sum)/q
            Rd.YvsV = cor(object$Y, object$variates$Y[, 1:ncomp], use = "pairwise")
            Rd.YvsV = apply(Rd.YvsV^2, 2, sum)/q
            
            own = cbind(Rd.XvsU, cumsum(Rd.XvsU))
            colnames(own) = c("Proportion", "Cumulative")
            rownames(own) = paste("comp", 1:ncomp)
            opp = cbind(Rd.XvsV, cumsum(Rd.XvsV))
            colnames(opp) = c("Proportion", "Cumulative")
            rownames(opp) = paste("comp", 1:ncomp)
            result$Rd.X = list(own = own, opp = opp)
            
            own = cbind(Rd.YvsV, cumsum(Rd.YvsV))
            colnames(own) = c("Proportion", "Cumulative")
            rownames(own) = paste("comp", 1:ncomp)
            opp = cbind(Rd.YvsU, cumsum(Rd.YvsU))
            colnames(opp) = c("Proportion", "Cumulative")
            rownames(opp) = paste("comp", 1:ncomp)
            result$Rd.Y = list(own = own, opp = opp)
        }
        
        #-- corr?lations canoniques --#
        can.cor = object$cor[1:ncomp]
        names(can.cor) = paste(1:ncomp, "th", sep = "")
        
        #-- valeurs sortantes --#
        result$can.cor = can.cor
        result$ncomp = ncomp
        result$cor = can.cor
        result$cutoff = cutoff
        result$keep = list(X = colnames(object$X)[gp.X], Y = colnames(object$Y)[gp.Y])
        result$what = what
        result$digits = digits
        result$method = 'rcc'
        
        class(result) = "summary"
        return(invisible(result))
    }

## --------------------------------- PCA ---------------------------------- ##
# from summary.prcomp, adapted for mixOmics
#' @rdname summary
#' @method summary pca
#' @export
summary.pca <-
    function (object, ...)
    {
        chkDots(...)
        vars <- object$prop_expl_var
        importance <- rbind(`Standard deviation` = object$sdev, `Proportion of Variance` = round(vars,
                                                                                                 5), `Cumulative Proportion` = round(cumsum(vars), 5))
        k <- ncol(object$rotation)
        colnames(importance) <- c(colnames(object$rotation), rep("",
                                                                 length(vars) - k))
        object$importance <- importance
        class(object) <- "summary.prcomp"
        object
    }


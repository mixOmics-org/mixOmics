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

#' Variable Importance in the Projection (VIP)
#' 
#' The function \code{vip} computes the influence on the \eqn{Y}-responses of
#' every predictor \eqn{X} in the model.
#' 
#' Variable importance in projection (VIP) coefficients reflect the relative
#' importance of each \eqn{X} variable for each \eqn{X} variate in the
#' prediction model. VIP coefficients thus represent the importance of each
#' \eqn{X} variable in fitting both the \eqn{X}- and \eqn{Y}-variates, since
#' the \eqn{Y}-variates are predicted from the \eqn{X}-variates.
#' 
#' VIP allows to classify the \eqn{X}-variables according to their explanatory
#' power of \eqn{Y}. Predictors with large VIP, larger than 1, are the most
#' relevant for explaining \eqn{Y}.
#' 
#' @param object object of class inheriting from \code{"pls"}, \code{"plsda"},
#' \code{"spls"} or \code{"splsda"}.
#' @return \code{vip} produces a matrix of VIP coefficients for each \eqn{X}
#' variable (rows) on each variate component (columns).
#' @author Sébastien Déjean, Ignacio Gonzalez, Florian Rohart, Al J Abadi
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{summary}}.
#' @references Tenenhaus, M. (1998). \emph{La regression PLS: theorie et
#' pratique}. Paris: Editions Technic.
#' @keywords regression multivariate
#' @export
#' @examples
#' data(linnerud)
#' X <- linnerud$exercise
#' Y <- linnerud$physiological
#' linn.pls <- pls(X, Y)
#' 
#' linn.vip <- vip(linn.pls)
#' 
#' barplot(linn.vip,
#' beside = TRUE, col = c("lightblue", "mistyrose", "lightcyan"),
#' ylim = c(0, 1.7), legend = rownames(linn.vip),
#' main = "Variable Importance in the Projection", font.main = 4)
vip <-
    function(object)
    {
        if (any(class(object) %in% c("mixo_plsda","mixo_splsda")))
        {
            object$Y = object$ind.mat
        } else if (any(class(object) %in% c("mixo_pls","mixo_spls"))) {
            #nothing
        } else {
            stop( " 'vip' is only implemented for the following objects: pls, plsda, spls, splsda", call.=FALSE)
        }
        #-- initialisation des matrices --#
        W = object$loadings$X
        H = object$ncomp
        q = ncol(object$Y)
        p = ncol(object$X)
        VIP = matrix(0, nrow = p, ncol = H)
        
        cor2 = cor(object$Y, object$variates$X, use = "pairwise")^2
        cor2 = as.matrix(cor2, nrow = q)
        
        VIP[, 1] = W[, 1]^2
        
        if (H > 1)
        {
            for (h in 2:H)
            {
                if (q == 1)
                {
                    Rd = cor2[, 1:h] 
                    VIP[, h] = Rd %*% t(W[, 1:h]^2) / sum(Rd)
                } else {
                    Rd = apply(cor2[, 1:h], 2, sum)
                    VIP[, h] = Rd %*% t(W[, 1:h]^2) / sum(Rd)
                }
            }
        }
        
        #-- valeurs sortantes --#
        VIP = sqrt(p * VIP)
        rownames(VIP) = rownames(W)
        colnames(VIP)= paste0("comp", 1:H)
        
        return(invisible(VIP))
    }

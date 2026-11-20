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

# -----------------------------------------------------------------------------
# calculation of the predicted area depending on the model.
# output to be passed to plotIndiv
# -----------------------------------------------------------------------------

# object: plsda or splsda object
# comp.predicted: prediction based on either component 1 or component 1:2
# dist: distance used in the predict function
# xlim: limit on the x-axis of the simulated variates
# ylim: limit on the y-axis of the simulated variates
# resolution: a total of resolution*resolution variates are simulated


# can only do a 2D prediction: cannot project 3D surface on 2D because we can
#     have multiple prediction for same point
# ex: variateXi_1 = variateXj_1 and variateXi_3 = variateXj_3 but
#     variateXi_2 !=variate Xj_3
# projection on comp1 and comp3 gives the same point,
#     but depending on variate_2, the prediction can be different

#' Calculate prediction areas
#' 
#' Calculate prediction areas that can be used in plotIndiv to shade the
#' background.
#' 
#' \code{background.predict} simulates \code{resolution*resolution} points
#' within the rectangle defined by xlim on the x-axis and ylim on the y-axis,
#' and then predicts the class of each point (defined by two coordinates). The
#' algorithm estimates the predicted area for each class, defined as the 2D
#' surface where all points are predicted to be of the same class. A polygon is
#' returned and should be passed to \code{\link{plotIndiv}} for plotting the
#' actual background.
#' 
#' Note that by default xlim and ylim will create a rectangle of simulated data
#' that will cover the plotted area of \code{plotIndiv}. However, if you use
#' \code{plotIndiv} with \code{ellipse=TRUE} or if you set \code{xlim} and
#' \code{ylim}, then you will need to adapt \code{xlim} and \code{ylim} in
#' \code{background.predict}.
#' 
#' Also note that the white frontier that defines the predicted areas when
#' plotting with \code{plotIndiv} can be reduced by increasing
#' \code{resolution}.
#' 
#' More details about the prediction distances in \code{?predict} and the
#' supplemental material of the mixOmics article (Rohart et al. 2017).
#' 
#' @param object A list of data sets (called 'blocks') measured on the same
#' samples. Data in the list should be arranged in matrices, samples x
#' variables, with samples order matching in all data sets.
#' @param comp.predicted Matrix response for a multivariate regression
#' framework. Data should be continuous variables (see block.splsda for
#' supervised classification and factor reponse)
#' @param dist distance to use to predict the class of new data, should be a
#' subset of \code{"centroids.dist"}, \code{"mahalanobis.dist"} or
#' \code{"max.dist"} (see \code{\link{predict}}).
#' @param xlim,ylim numeric list of vectors of length 2, giving the x and y
#' coordinates ranges for the simulated data. By default will be \eqn{1.2*} the
#' range of object$variates$X[,i]
#' @param resolution A total of \code{resolution*resolution} data are simulated
#' between xlim[1], xlim[2], ylim[1] and ylim[2].
#' @return \code{background.predict} returns a list of coordinates to be used
#' with \code{\link{polygon}} to draw the predicted area for each class.
#' @author Florian Rohart, Al J Abadi
#' @seealso \code{\link{plotIndiv}}, \code{\link{predict}},
#' \code{\link{polygon}}.
#' @references Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package
#' for 'omics feature selection and multiple data integration. PLoS Comput Biol
#' 13(11): e1005752
#' @example ./examples/background.predict-examples.R
#' @export
background.predict <- 
    function(object, comp.predicted = 1, dist = "max.dist",
             xlim = NULL, ylim = NULL, resolution = 100)
    {
        
        
        if (!any(class(object) %in% c("mixo_plsda", "mixo_splsda")))
            stop("'background.predict' can only be calculated for 'plsda'
        and 'splsda' objects")
        
        if (!any(dist %in% c("max.dist", "centroids.dist", "mahalanobis.dist")))
            stop(
                "Choose one of the three following distances: 'max.dist',
        'centroids.dist' or 'mahalanobis.dist'"
            )
        
        if (!comp.predicted %in% c(1, 2))
            stop("Can only show predicted background for 1 or 2 components")
        
        if (!is.null(xlim) && length(xlim) != 2)
            stop(
                "'xlim' must be a vector of two values, indicating the min
        and max of the simulated data on variates 1 (x-axis)"
            )
        
        if (!is.null(ylim) && length(ylim) != 2)
            stop(
                "'ylim' must be a vector of two values, indicating the min
        and max of the simulated data on variates 2 (y-axis)"
            )
        
        if (resolution <= 0)
            stop("'resolution' must be a positive value")
        
        # ... = arg to pass to plotIndiv
        #plotIndiv(object, style = "graphics", ...)
        #plot(-10:10,-10:10,type="n")
        
        
        ####################################
        # ---- simulating new data
        ####################################
        X = object$X
        Y = object$Y
        
        # we only need to simulate variates
        lim = apply(object$variates$X, 2, range) * 1.2
        if (is.null(xlim))
            xlim = lim[, 1]
        if (is.null(ylim))
            ylim = lim[, 2]
        
        lim = cbind(xlim, ylim)
        
        increment = apply(lim, 2, function(x) {
            sum(abs(x)) / resolution
        })
        incrementx = increment[1]
        incrementy = increment[2]#(abs(ylim[1]) + abs(ylim[2]))/resolution
        #incrementy = (abs(zlim[1]) + abs(zlim[2]))/resolution
        
        
        list.grid = lapply(1:2, function(x) {
            seq(lim[1, x], lim[2, x], increment[x])
        })
        grid = as.matrix(expand.grid(list.grid))
        
        t.pred = list(grid)
        ncomp = comp.predicted
        J = 1
        q = nlevels(Y)
        variatesX = list(object$variates$X)
        Y.prim = unmap(Y)
        
        
        
        ####################################
        # ---- estimate polygon
        ####################################
        poly.save = vector("list", length = nlevels(Y))
        G = cls = list()
        
        if (dist == "max.dist")
        {
            variatesX = list(X = object$variates [-2][[1]][, 1:comp.predicted,
                                                           drop = FALSE])
            Y = object$ind.mat
            means.Y = matrix(
                attr(Y, "scaled:center"),
                nrow = nrow(t.pred[[1]]),
                ncol = q,
                byrow = TRUE
            )
            sigma.Y = matrix(
                attr(Y, "scaled:scale"),
                nrow = nrow(t.pred[[1]]),
                ncol = q,
                byrow = TRUE
            )
            
            Cmat = crossprod(Y, variatesX[[1]])
            
            Y = object$Y
            
            #print(variatesX)
            
            Y.hat.temp = Y.hat = list()
            for (j in 1:ncomp)
            {
                A = matrix(
                    apply(variatesX[[1]][, 1:j, drop = FALSE], 2,
                          function(y) {
                              (norm(y, type = "2")) ^ 2
                          }),
                    nrow = nrow(t.pred[[1]]),
                    ncol = j,
                    byrow = TRUE
                )
                Y.hat.temp[[j]] = ((as.matrix(t.pred[[1]][, 1:j, drop = FALSE]) /
                                        A) %*%
                                       t(Cmat)[1:j, , drop = FALSE])
                #                *sigma.Y+means.Y
            }
            Ypred = sapply(Y.hat.temp, function(x) {
                x * sigma.Y + means.Y
            },
            simplify = "array")
            Y.hat[[1]] = Ypred
            
            cls$max.dist = lapply(1:J, function(x) {
                matrix(sapply(1:ncomp[x],
                              # List level
                              function(y) {
                                  apply(Y.hat[[x]][, , y, drop = FALSE], 1,
                                        # component level
                                        function(z) {
                                            paste(levels(Y)[which(z == max(z))], collapse = "/")
                                        }) # matrix level
                              }),
                       nrow = nrow(t.pred[[x]]),
                       ncol = ncomp[x])
            })
            cls$max.dist = lapply(1:J, function(x) {
                colnames(cls$max.dist[[x]]) =
                    paste0(rep("comp", ncomp[x]), 1:ncomp[[x]])
                
                rownames(cls$max.dist[[x]]) = rownames(t.pred[[x]])
                
                return(cls$max.dist[[x]])
            })
            names(cls$max.dist) = names(X)
            
        }
        
        if (dist == "mahalanobis.dist" | dist == "centroids.dist")
        {
            for (i in 1:J)
            {
                G[[i]] = sapply(1:q, function(x) {
                    apply(as.matrix(variatesX[[i]][Y.prim[, x] == 1, 1:ncomp[i] ,
                                                   drop = FALSE]), 2, mean)
                })
                if (ncomp[i] == 1)
                    G[[i]] = t(t(G[[i]]))
                else
                    G[[i]] = t(G[[i]])
                colnames(G[[i]]) = paste0("dim", c(1:ncomp[i]))
                rownames(G[[i]]) = levels(Y)
                
            }
            names(G) = names(X)
            
            
            # predicting class of simulated data
            if (dist == "centroids.dist")
            {
                ###Start: Centroids distance
                cl = list()
                centroids.fun = function(x, G, h, i) {
                    q = nrow(G[[i]])
                    x = matrix(x,
                               nrow = q,
                               ncol = h,
                               byrow = TRUE)
                    
                    if (h > 1) {
                        d = apply((x - G[[i]][, 1:h]) ^ 2, 1, sum)
                    }
                    else {
                        d = (x - G[[i]][, 1]) ^ 2
                    }
                    cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
                }
                
                for (i in 1:J)
                {
                    cl[[i]] = matrix(nrow = nrow(t.pred[[i]]), ncol = ncomp[i])
                    
                    for (h in 1:ncomp[[i]])
                    {
                        cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1,
                                      function(x) {
                                          centroids.fun(
                                              x = x,
                                              G = G,
                                              h = h,
                                              i = i
                                          )
                                      })
                        cl[[i]][, h] = cl.id
                    }
                }
                
                cls$centroids.dist = lapply(1:J, function(x) {
                    colnames(cl[[x]]) =
                        paste0(rep("comp", ncomp[x]), 1:ncomp[[x]])
                    
                    return(cl[[x]])
                })
                
            } else if (dist == "mahalanobis.dist") {
                ### Start: Mahalanobis distance
                cl = list()
                Sr.fun = function(x, G, Yprim, h, i) {
                    q = nrow(G[[i]])
                    Xe = Yprim %*% G[[i]][, 1:h]
                    #Xr = object$variates$X[, 1:h] - Xe
                    Xr = variatesX[[i]][, 1:h] - Xe
                    Sr = t(Xr) %*% Xr / nrow(Yprim)
                    Sr.inv = solve(Sr)
                    x = matrix(x,
                               nrow = q,
                               ncol = h,
                               byrow = TRUE)
                    if (h > 1) {
                        mat = (x - G[[i]][, 1:h]) %*% Sr.inv %*% t(x - G[[i]][, 1:h])
                        d = apply(mat ^ 2, 1, sum)
                    } else {
                        d = drop(Sr.inv) * (x - G[[i]][, 1]) ^ 2
                    }
                    cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
                }
                
                for (i in 1:J) {
                    cl[[i]] = matrix(nrow = nrow(t.pred[[1]]), ncol = ncomp[i])
                    
                    for (h in 1:ncomp[[i]]) {
                        cl.id = apply(
                            matrix(t.pred[[i]][, 1:h], ncol = h),
                            1,
                            Sr.fun,
                            G = G,
                            Yprim = Y.prim,
                            h = h,
                            i = i
                        )
                        cl[[i]][, h] = cl.id
                    }
                }
                
                cls$mahalanobis.dist = lapply(1:J, function(x) {
                    colnames(cl[[x]]) =
                        paste0(rep("comp", ncomp[x]), 1:ncomp[[x]])
                    
                    return(cl[[x]])
                })
            }
        }
        
        for (ind.area in 1:nlevels(Y))
        {
            ind1 = which(cls[[dist]][[1]][, comp.predicted] == levels(Y)[ind.area])
            
            if (length(ind1) > 0)
            {
                # if less than 8 direct neighbours, we keep the point => contour
                # from one point from the contour, we can only test the direct
                # neighbours to speed up
                area = t.pred[[1]][ind1, 1:2, drop = FALSE]#can only do it in 2d
                contour = NULL
                for (i in 1:nrow(area))
                {
                    areax = area[, 1]#as.numeric(as.character(area[,1]))
                    areay = area[, 2]#as.numeric(as.character(area[,2]))
                    
                    a = areax[i]
                    b = areay[i]
                    
                    
                    res = 0
                    for (x in c(a - incrementx, a, a + incrementx))
                    {
                        for (y in c(b - incrementy, b, b + incrementy))
                        {
                            temp =  intersect(which(areax ==
                                                        as.numeric(as.character(
                                                            x
                                                        ))),
                                              which(areay ==
                                                        as.numeric(as.character(
                                                            y
                                                        ))))
                            if (length(temp) > 0)
                                res = res + 1
                        }
                    }
                    
                    if (res != 9)
                        contour = c(contour, i)
                    
                    if (length(contour) == 2)
                        break
                }
                
                # now that we have two point of the contour,
                # we look for others in the direct neighbours.
                
                added = TRUE
                while (added)
                {
                    # as long as we're adding a point in contour,
                    # we keep looking for another one
                    added = FALSE
                    i = length(contour)
                    
                    point = contour[i]
                    
                    areax = area[, 1]#round(as.numeric(as.character(area[,1])),7)
                    areay = area[, 2]#round(as.numeric(as.character(area[,2])),7)
                    
                    a = areax[point]#round(areax[point],7)
                    b = areay[point]#round(areay[point],7)
                    
                    # we want to add the point (x,y) that has the lowest number of
                    # neighbour (the more extreme on the edge)
                    
                    neighbour = contour.temp = NULL
                    # around the point that is in the contour
                    for (x in c(a - incrementx, a, a + incrementx))
                    {
                        # around the point that is in the contour
                        for (y in c(b - incrementy, b, b + incrementy))
                        {
                            # (x,y) is a neighbour of (a,b) and I want to see
                            # whether it has 8+1 neighbours or less
                            res = 0
                            for (xx in c(x - incrementx, x, x + incrementx))
                            {
                                for (yy in c(y - incrementy, y, y + incrementy))
                                {
                                    # (xx,yy) is a neighbour of (x,y)
                                    temp =  intersect(which(abs(areax - xx) < 1e-5),
                                                      which(abs(areay - yy) < 1e-5))
                                    #which(area[,1]==as.numeric(x) & area[,2]==as.numeric(y))
                                    #print(xx)
                                    #print(yy)
                                    #print(temp)
                                    
                                    if (length(temp) > 0)
                                        res = res + 1
                                }
                            }
                            #print(res)
                            if (res != 9)
                                # if (x,y) has less than 8+1 neighbour,
                                # then it's on the edge and I want it,
                                # only if it's not already in contour
                            {
                                # recover which indice in area the point is
                                ind = intersect(which(abs(areax - x) < 1e-5),
                                                which(abs(areay - y) < 1e-5))
                                
                                # check whether it is already in contour
                                if (length(ind) > 0 &&
                                    sum(contour == ind) == 0)
                                {
                                    contour.temp = c(contour.temp, ind)
                                    neighbour = c(neighbour, res)
                                    #contour = c(contour, ind)
                                    added = TRUE
                                }
                            }
                            
                        }
                    }
                    
                    if (length(contour.temp) > 0)
                    {
                        contour = c(contour, contour.temp[which.min(neighbour)])
                    } else {
                        added = FALSE
                    }
                    
                }
                
                poly = area[contour, ]
                poly.save[[ind.area]] = poly
                
            }
            
        }
        names(poly.save) = levels(Y)#adjustcolor(color.mixo(ind.area), alpha.f=0.1)
        
        
        class(poly.save) = "background.predict"
        return(poly.save)
    }



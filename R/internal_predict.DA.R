################################################################################
# Author :
#   Florian Rohart,
#
# created: 24-05-2015
# last modified: 04-10-2017
#
# Copyright (C) 2015
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################

# ==============================================================================
# internal_predict.DA: prediction module for Discriminant Analysis.
#   Used in 'predict.mint.block.pls.R'
# ==============================================================================


#' @importFrom dplyr group_by row_number filter n
internal_predict.DA <- 
    function(object, out, q, dist, weights)
    {
        # a DA analysis (mint).(block).(s)plsda
        if (!is(object, "DA"))
            stop("'Object' is not from a Discriminant Analysis", call.=FALSE)
        
        out.DA = list()
        J = length(object$X) #at this stage we have a list of blocks
        p = lapply(object$X, ncol)
        t.pred = out$variates
        Y.hat = out$predict
        newdata = out$newdata #actually concat.newdata
        variatesX = object$variates[-(J + 1)];
        ncomp = object$ncomp
        
        Y = object$Y
        Y.prim = unmap(object$Y)
        G = cls = list()
        for (i in 1 : J)
        {
            G[[i]] = sapply(1:q, function(x)
            {apply(as.matrix(variatesX[[i]][Y.prim[, x] == 1,
                                            ,drop=FALSE]), 2, mean)})
            if (ncomp[i] == 1)
                G[[i]] = t(t(G[[i]]))
            else
                G[[i]] = t(G[[i]])
            colnames(G[[i]]) = paste0("dim", c(1:ncomp[i]))
            rownames(G[[i]]) = colnames(object$ind.mat)
            
        }
        names(G)=names(object$X)
        
        ### Start: Maximum distance
        if (any(dist == "all") || any(dist == "max.dist"))
        {
            cls$max.dist = lapply(1:J, function(x){matrix(sapply(1:ncomp[x],
                                                                 ### List level
                                                                 function(y){apply(Y.hat[[x]][, , y, drop = FALSE], 1,
                                                                                   ### component level
                                                                                   function(z){
                                                                                       paste(levels(Y)[which(z == max(z))], collapse = "/")
                                                                                   }) ### matrix level
                                                                 }), nrow = nrow(newdata[[x]]), ncol = ncomp[x])
            })
            cls$max.dist = lapply(1:J, function(x){colnames(cls$max.dist[[x]]) =
                paste0(rep("comp", ncomp[x]), 1 : ncomp[[x]]);
            rownames(cls$max.dist[[x]]) = rownames(newdata[[x]]);
            return(cls$max.dist[[x]])})
            names(cls$max.dist)=names(object$X)
        }
        
        
        ### Start: Centroids distance
        if (any(dist == "all") || any(dist == "centroids.dist"))
        {
            cl = list()
            centroids.fun = function(x, G, h, i) {
                q = nrow(G[[i]])
                x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
                
                if (h > 1) {
                    d = apply((x - G[[i]][, 1:h])^2, 1, sum)
                }
                else {
                    d = (x - G[[i]][, 1])^2
                }
                cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
            }
            
            for (i in 1 : J)
            {
                cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
                
                for (h in 1 : ncomp[[i]])
                {
                    cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1,
                                  function(x) {centroids.fun(x = x, G = G, h = h, i = i)})
                    cl[[i]][, h] = cl.id
                }
            }
            
            cls$centroids.dist = lapply(1:J, function(x){colnames(cl[[x]]) =
                paste0(rep("comp", ncomp[x]), 1 : ncomp[[x]]);
            rownames(cl[[x]]) = rownames(newdata[[x]]); return(cl[[x]])})
            names(cls$centroids.dist)=names(object$X)
        }### End: Centroids distance
        
        
        ### Start: Mahalanobis distance
        if (any(dist == "all") || any(dist == "mahalanobis.dist"))
        {
            cl = list()
            Sr.fun = function(x, G, Yprim, h, i) {
                q = nrow(G[[i]])
                Xe = Yprim %*% G[[i]][, 1:h]
                #Xr = object$variates$X[, 1:h] - Xe
                Xr = variatesX[[i]][, 1:h] - Xe
                Sr = t(Xr) %*% Xr/nrow(Yprim)
                Sr.inv = solve(Sr)
                x = matrix(x, nrow = q, ncol = h, byrow = TRUE)
                if (h > 1) {
                    mat = (x - G[[i]][, 1:h]) %*% Sr.inv %*% t(x - G[[i]][, 1:h])
                    d = apply(mat^2, 1, sum)
                } else {
                    d = drop(Sr.inv) * (x - G[[i]][, 1])^2
                }
                cl.id = paste(levels(Y)[which(d == min(d))], collapse = "/")
            }
            
            for (i in 1 : J){
                cl[[i]] = matrix(nrow = nrow(newdata[[1]]), ncol = ncomp[i])
                
                for (h in 1:ncomp[[i]]) {
                    cl.id = apply(matrix(t.pred[[i]][, 1:h], ncol = h), 1, Sr.fun,
                                  G = G, Yprim = Y.prim, h = h, i = i)
                    cl[[i]][, h] = cl.id
                }
            }
            
            cls$mahalanobis.dist = lapply(1:J, function(x){colnames(cl[[x]]) =
                paste0(rep("comp", ncomp[x]), 1 : ncomp[[x]]);
            rownames(cl[[x]]) = rownames(newdata[[x]]);return(cl[[x]])})
            names(cls$mahalanobis.dist)=names(object$X)
        } ### End: Mahalanobis distance
        
        out.DA$class = cls
        
        ### End if discriminant analysis is performed
        
        # at this stage, we have the classification of each sample for each dataset
        #   of object$X
        # now we need to combine the classification by vote (majority wins),
        #  only when more than one block, otherwise 'vote' is classic classification
        if (length(object$X)>1)
        {
            for (ijk in 1:length(out.DA$class))# loop on the dist
            {
                # create a temporary array to make computation on the lists easier
                temp=array(, c(nrow(newdata[[1]]), min(ncomp), J))
                for(i in 1:J)
                {
                    temp[, , i] = out.DA$class[[ijk]][[i]][, 1:min(ncomp)]
                    
                }
                # look at the majority vote for all dataset of object$X (with table)
                #   if more than a unique max, we put NA
                table.temp = apply(temp,c(1,2), function(x)
                {a=table(x); if (length(which(a==max(a)))==1)
                {b=names(which.max(a))}else{b=NA}; b})
                colnames(table.temp) =
                    colnames(out.DA$class[[ijk]][[i]])[1:min(ncomp)]
                rownames(table.temp) = rownames(out.DA$class[[ijk]][[i]])
                out.DA$MajorityVote[[ijk]] = table.temp
            }
            names(out.DA$MajorityVote) = names(out.DA$class)
            
            #save(list=ls(),file="temp.Rdata")
            #stop("")
            # weighted vote for each distance, each comp
            if(!is.null(weights))
            {
                out.WeightedVote = vector("list",length=length(out.DA$class))
                Group.2 = n = x =NULL #CRAN check
                for(i in 1:length(out.DA$class)){ #i distance
                    out = matrix(NA_real_,nrow=nrow(newdata[[1]]), ncol=min(ncomp))
                    rownames(out) = rownames(newdata[[1]])
                    colnames(out) = paste0("comp",1:min(ncomp))
                    
                    for(comp in 1:min(ncomp)){ #comp
                        data.temp=NULL
                        for(j in 1:J){ #block
                            data.temp = rbind(data.temp,
                                              out.DA$class[[i]][[j]][,comp,drop=FALSE])
                        }
                        colnames(data.temp)="pred"
                        temp=data.frame(data.temp,indiv=rownames(data.temp),
                                        weights=rep(weights,each=nrow(out.DA$class[[1]][[1]])))
                        ag = aggregate(temp$weights,
                                       by=list(temp$pred, temp$indiv), FUN=sum)
                        data_max <- group_by(.data = ag, Group.2)
                        data_max <- filter(.data = data_max, row_number(x)==n())
                        out.comp = as.matrix(data_max[,1])
                        rownames(out.comp) = as.matrix(data_max[,2])
                        colnames(out.comp) = paste0("comp",comp)
                        
                        out[,comp] = out.comp[match(rownames(out),
                                                    rownames(out.comp)),]
                    }
                    
                    out.WeightedVote[[i]] = out
                }
                names(out.WeightedVote) = names(out.DA$class)
                out.DA$WeightedVote = out.WeightedVote
                
                
                
                if(FALSE){
                    out.DA$WeightedVote = lapply(out.DA$class, function(x){
                        # x is a distance
                        class.per.comp = lapply(1:min(ncomp), function(y) {
                            matrix(sapply(x, function(z)
                                z[,y, drop = FALSE]),ncol=J)})
                        # combine the results per component
                        names(class.per.comp) = paste0("comp",1:min(ncomp))
                        class.weighted.per.comp = sapply(class.per.comp, function(y){
                            # for each component
                            apply(y,1,function(z){
                                # we aggregate the results of each individuals using the 'weights'
                                temp = aggregate(weights,list(as.character(z)),sum)
                                ind = which(temp[,2]== max (temp[,2]))
                                # if two max, then NA
                                if(length(ind) == 1)
                                {
                                    res = as.character(temp[ind, 1])
                                } else {
                                    res = NA
                                }
                                res
                                
                            })
                            
                        })
                        class.weighted.per.comp = matrix(class.weighted.per.comp,
                                                         nrow = nrow(class.per.comp[[1]]))
                        colnames(class.weighted.per.comp) = names(class.per.comp)
                        rownames(class.weighted.per.comp) =
                            rownames(out.DA$MajorityVote[[1]])
                        class.weighted.per.comp
                        
                    })
                }
                out.DA$weights = weights
                
            }
        }else{
            out.DA$MajorityVote = lapply(out.DA$class,function(x){x[[1]]})
        }
        
        block.object = c("block.pls", "block.spls", "block.plsda", "block.spsda")
        if (inherits(object, block.object) & J>1) # a block
        {
            out.DA$centroids = G
        }else{ #not a block
            out.DA$centroids = G[[1]]
            out.DA$class = out.DA$MajorityVote
        }
        if (any(dist == "all"))
            dist = "all"
        
        out.DA$dist = dist
        
        out.DA
        
    }






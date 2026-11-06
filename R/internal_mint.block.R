################################################################################
# Authors:
#   Florian Rohart,
#   Benoit Gautier,
#   Amrit Singh,
#   Kim-Anh Le Cao,
#
# created: 20-07-2014
# last modified: 04-10-2017
#
# Copyright (C) 2014
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Parts of this file come from the RGCCA package (GPL (>= 2)).
################################################################################


################################################################################
# Functions modified from RGCCA R-library
#   sgcca(), sgccak()
#
# Functions acquired from RGCCA R-library, see 'internal_mint.block_helpers.R'
#   cov2(), initsvd(), crossprod(),
#   defl.select(),
################################################################################

internal_mint.block <- 
    function (A, indY = NULL,  design = 1 - diag(length(A)),
              tau=NULL,#rep(1, length(A)),
              ncomp = rep(1, length(A)), scheme = "horst", scale = TRUE,
              init = "svd.single", tol = 1e-06,
              mode = "canonical", max.iter = 100,study = NULL, keepA = NULL,
              penalty = NULL, all.outputs = FALSE, misdata = NULL, is.na.A = NULL,
              ind.NA = NULL, ind.NA.col = NULL, remove.object=NULL)
    {
        # A: list of matrices
        # indY: integer, pointer to one of the matrices of A
        # design: design matrix, links between matrices. Diagonal must be 0
        # tau: numeric vector of length the number of blocks in \code{X}.
        #   Each regularization parameter will be applied on each block and takes
        #   the value between 0 (no regularisation) and 1.
        #   If tau = "optimal" the shrinkage paramaters are estimated for each block
        # ncomp: vector of ncomp, per matrix
        # scheme: a function "g", refer to the article (thanks Benoit)
        # scale: do you want to scale. mean is done by default and cannot be changed
        # init: one of "svd" or "random", initialisation of the algorithm
        # tol: nobody cares about this
        # mode: canonical, classic, invariant, regression
        # max.iter: nobody cares about this
        # study: factor for each matrix of A, must be a vector
        # keepA: keepX of spls for each matrix of A. must be a list.
        #   Each entry must be of the same length (max ncomp)
        # penalty: numeric vector of length the number of blocks in \code{X}.
        #   Each penalty parameter will be applied on each block and takes the value
        #   between 0 (no variable selected) and 1 (all variables included).
        # all.outputs: calculation of non-essential outputs
        #   (e.g. explained variance, loadings.Astar, etc)
        # misdata: optional. any missing values in the data? list,
        #   misdata[[q]] for each data set
        # is.na.A: optional. where are the missing values? list,
        #   is.na.A[[q]] for each data set (if misdata[[q]] == TRUE)
        # ind.NA: optional. which rows have missing values? list,
        #   ind.NA[[q]] for each data set.
        # ind.NA.col: optional. which col have missing values? list,
        #   ind.NA.col[[q]] for each data set.
        
        
        names(ncomp) = names(A)
        
        
        # center the data per study, per matrix of A, scale if scale=TRUE, option
        mean_centered = lapply(A, function(x)
        {mean_centering_per_study(x, study, scale)})
        
        A = lapply(mean_centered, function(x){as.matrix(x$concat.data)})
        
        #save rownames study
        mean_centered.rownames.study = vector("list", nlevels(study))
        for (m in 1:nlevels(study))
            mean_centered.rownames.study[[m]] = mean_centered[[1]]$rownames.study[[m]]
        
        rm(mean_centered) #free memory
        
        ni = table(study) #number of samples per study
        
        ### Start: Initialization parameters
        pjs = sapply(A, NCOL)
        nb_ind = NROW(A[[1]])
        J = length(A)
        R = A # R: residuals matrices, will be a list of length ncomp
        N = max(ncomp)
        AVE_inner = AVE_outer = rep(NA, max(ncomp))
        
        
        # keepA[[comp]] is a matrix where each row is all the keepX the test over
        #   the block (each block is a column)
        
        #number of models to be tested: either a keepA per component,
        #   or multiple (e.g. in tune functions)
        number.models.per.comp = sapply(keepA, nrow)
        one.model = !any( number.models.per.comp !=1)
        
        
        AVE_X = crit = loadings.partial.A = variates.partial.A = tau.rgcca = list()
        P = loadings.A = loadings.Astar = variates.A =  vector("list", J)
        
        
        if(one.model) # more outputs that what is needed for tune functions
        {
            for (k in 1:J)
                variates.A[[k]] = matrix(NA_real_, nb_ind, N)
            
            for (k in 1:J)
            {
                loadings.A[[k]] = matrix(NA_real_, pjs[[k]], N)
                if(all.outputs)
                    P[[k]] = loadings.Astar[[k]]= matrix(NA_real_, pjs[[k]], N)
            }
            
            for (k in 1:J)
            {
                loadings.partial.A[[k]] = variates.partial.A[[k]] = vector("list",
                                                                           length = nlevels(study))
                for(m in 1:nlevels(study))
                {
                    loadings.partial.A[[k]][[m]] = matrix(NA_real_,
                                                          nrow = NCOL(A[[k]]), ncol = N)
                    variates.partial.A[[k]][[m]] = matrix(NA_real_,
                                                          nrow = ni[m], ncol = N)
                }
            }
        } else {
            for (k in 1:J)
            {
                variates.A[[k]] = matrix(NA_real_, nb_ind,
                                         sum(number.models.per.comp))
                loadings.A[[k]] = matrix(NA_real_, pjs[[k]],
                                         sum(number.models.per.comp))
            }
            loadings.partial.A = variates.partial.A = NULL
            # not needed for tune functions
        }
        
        
        ndefl = ncomp - 1
        J2 = J-1
        
        if (is.vector(tau))
            tau = matrix(rep(tau, N), nrow = N, ncol = length(tau), byrow = TRUE)
        
        #save(list=ls(),file="temp.Rdata")
        
        # if missing values are not given as input (only when direct call to a
        #   (mint).(block).(s)pls(da)), we search for them here (takes time)
        if(is.null(misdata) &  is.null(is.na.A) & is.null(ind.NA) &
           is.null(ind.NA.col))
        {
            misdata = sapply(A, anyNA) # Detection of missing data per block
            misdata.all = any(misdata) # is there any missing data overall
            
            #save(list=ls(),file="temp.Rdata")
            if (misdata.all)
            {
                is.na.A = temp = vector("list",length=length(A))
                is.na.A[misdata] = lapply(A[misdata], is.na) # size n*p,
                #   which entry is na. might be none, but at least one in all the
                #   block will be a TRUE
                
                temp[misdata] = lapply(is.na.A[misdata], function(x)
                {which(x,arr.ind=TRUE)})
                ind.NA = lapply(temp, function(x){unique(x[,1])})
                ind.NA.col = lapply(temp, function(x){unique(x[,2])})
                
            }else {
                is.na.A = NULL
                ind.NA = ind.NA.col = NULL
            }
        } else{
            misdata.all = any(misdata)
        }
        
        if(all.outputs & J==2 & nlevels(study) == 1 & one.model)
            #(s)pls(da) models, we calculate mat.c
        {
            if(misdata.all)
            {
                p.ones = rep(1, ncol(A[[1]]))
                is.na.X = is.na.A[[1]]
            }
            mat.c = matrix(0, nrow = ncol(A[[1]]), ncol = N,
                           dimnames = list(colnames(A[[1]],  paste0("comp", 1:N))))
        } else {mat.c = NULL}
        
        ### End: Initialization parameters
        
        iter=NULL
        compteur = 0
        for (comp in 1 : N)
        {
            
            if(misdata.all)# replace NA in A[[q]] by 0
                for(j in c(1:J)[misdata])
                    R[[j]][is.na.A[[j]]]=0 # faster than using replace
                # if missing data, R is the one replace by 0 where NA are supposed to be
                
                # initialisation_by_svd, get the loadings.A
                loadings.A.init = initialisation_by_svd(R, indY, misdata, is.na.A, init = init)
                
                # loop on keepA[[comp]]: multiple values per block and we go through
                #   them. Need to have the same number of values per block.
                # we assume keepA[[comp]] is a grid here: columns are the blocks,
                #   rows are the different keepX
                for(ijk.keepA in 1:nrow(keepA[[comp]]))
                {
                    compteur = compteur +1
                    keepA.ijk = keepA[[comp]][ijk.keepA,]
                    
                    ### start - repeat/convergence
                    if (is.null(tau))
                    {
                        mint.block.result = sparse.mint.block_iteration(R, design,
                                                                        study = study, loadings.A = loadings.A.init,
                                                                        keepA = keepA.ijk, #keepA is one value per block
                                                                        scheme = scheme, max.iter = max.iter, tol = tol,
                                                                        penalty = penalty,
                                                                        misdata=misdata, is.na.A=is.na.A, ind.NA = ind.NA,
                                                                        all.outputs = all.outputs)
                    } else {
                        mint.block.result = sparse.rgcca_iteration(R, design,
                                                                   tau = if (is.matrix(tau)){tau[comp, ]} else {"optimal"},
                                                                   scheme = scheme, init = init, tol = tol,
                                                                   max.iter = max.iter, penalty = penalty,
                                                                   keepA = keepA.ijk, all.outputs = all.outputs)
                    }
                    ### end - repeat/convergence
                    
                    if(one.model)
                    {
                        # reshape outputs
                        for (k in 1 : J)
                        {
                            loadings.A[[k]][, comp] = mint.block.result$loadings.A[[k]]
                            variates.A[[k]][, comp] = mint.block.result$variates.A[, k]
                            
                            if(is.null(tau))
                            {
                                # recording loadings.partials, $Ai$study[,ncomp]
                                # recording variates.partials, $Ai[,ncomp]
                                for(k in 1:J)
                                {
                                    for(m in 1:nlevels(study))
                                    {
                                        loadings.partial.A[[k]][[m]][, comp] =
                                            matrix(mint.block.result$
                                                       loadings.partial.A.comp[[k]][[m]], ncol=1)
                                        variates.partial.A[[k]][[m]][, comp] =
                                            matrix(mint.block.result$
                                                       variates.partial.A.comp[[k]][[m]], ncol=1)
                                    }
                                }
                            }
                        }
                    } else {
                        # no record of partial component for multilple models, for gain of memory
                        for (k in 1 : J)
                        {
                            loadings.A[[k]][, compteur] =
                                mint.block.result$loadings.A[[k]]
                            variates.A[[k]][, compteur] =
                                mint.block.result$variates.A[, k]
                        }
                    }
                    
                    crit[[comp]] = mint.block.result$crit
                    tau.rgcca[[comp]] = mint.block.result$tau
                    if(all.outputs)
                        AVE_inner[comp] = mint.block.result$AVE_inner
                    
                    if(all.outputs & J==2 & nlevels(study) == 1 & one.model)
                        # mat.c, (s)pls(da)
                    {
                        if(misdata.all) #only one model, so misdata[1]=TRUE
                        {
                            R.temp = R[[1]]
                            R.temp[is.na.X] = 0
                            c = crossprod(R.temp, variates.A[[1]][,comp])
                            rm(R.temp) #free memory
                            
                            #save(list=ls(),file="temp.Rdata")
                            
                            t.norm = rep(crossprod(variates.A[[1]][,comp]), length(c))
                            
                            if(length(ind.NA.col[[1]])>0) # should always be true
                            {
                                temp = drop(variates.A[[1]][,comp]) %o% rep(1,
                                                                            length(ind.NA.col[[1]])) #p*n -> p * where there are NA
                                temp[is.na.X[,ind.NA.col[[1]],drop=FALSE]] = 0
                                t.norm[ind.NA.col[[1]]] = apply(temp,2, crossprod)
                            }
                            c = c / t.norm
                            mat.c[,comp] = c
                        } else {
                            mat.c[,comp] <- t(crossprod(variates.A[[1]][,comp],
                                                        R[[1]])) / drop(crossprod (variates.A[[1]][,comp]))
                        }
                    } else {
                        mat.c = NULL
                    }
                    
                    # deflation if there are more than 1 component and if we haven't
                    #   reached the max number of component (N)
                    if (N != 1 & comp != N)
                    {
                        
                        defla.result = defl.select(yy=mint.block.result$variates.A,
                                                   rr=R, nncomp=ndefl, nn=comp, nbloc = J, indY = indY,
                                                   mode = mode, aa = mint.block.result$loadings.A,
                                                   misdata=misdata, is.na.A=is.na.A, ind.NA = ind.NA.col)
                        
                        R = defla.result$resdefl
                        
                        if(!(all.outputs & one.model))
                            defla.result$resdefl=NULL
                        #free memory, only if not used in the loop below
                    }
                    
                    
                    if(all.outputs & one.model) #loadings.Astar
                    {
                        for (k in 1 : J)
                        {
                            if (N != 1)
                                P[[k]][, comp - 1] = defla.result$pdefl[[k]]
                        }
                        
                        if (comp == 1)
                        {
                            for (k in 1 : J)
                                loadings.Astar[[k]][, comp] = mint.block.result$loadings.A[[k]]
                        } else {
                            for (k in 1 : J)
                                loadings.Astar[[k]][, comp] =
                                    mint.block.result$loadings.A[[k]] - loadings.Astar[[k]][,
                                                                                            (1 : comp - 1), drop = FALSE] %*% drop(t(loadings.A[[k]][,
                                                                                                                                                     comp]) %*% P[[k]][, 1 : (comp - 1), drop = FALSE])
                        }
                    } else {
                        loadings.Astar = NULL
                    }
                    iter = c(iter, mint.block.result$iter)
                    
                } ### End loop on keepA
        } ### End loop on ncomp
        
        
        #### any model
        # loadings.A[[block]][1:p, all.keepA.tested]
        # variates.A[[block]][1:n, all.keepA.tested]
        
        #### a single model
        # loadings.partial.A[[block]][[study]][, 1:ncomp]
        # variates.partial.A[[block]][[study]][, 1:ncomp]
        # loadings.Astar[[block]][, 1:ncomp]
        
        if(one.model)
        {
            # only one model
            shave.matlist = function(mat_list, nb_cols) mapply(function(m, nbcomp)
                m[, 1:nbcomp, drop = FALSE], mat_list, nb_cols, SIMPLIFY = FALSE)
            shave.veclist = function(vec_list, nb_elts) mapply(function(m, nbcomp)
                m[1:nbcomp], vec_list, nb_elts, SIMPLIFY = FALSE)
            
            
            for (k in 1:J)
            {
                rownames(loadings.A[[k]]) = colnames(A[[k]])
                
                if(all.outputs)
                    rownames(loadings.Astar[[k]]) = colnames(A[[k]])
                
                rownames(variates.A[[k]]) = rownames(A[[k]])
                colnames(variates.A[[k]]) = colnames(loadings.A[[k]]) =
                    paste0("comp", 1:max(ncomp))
                if(all.outputs)
                    AVE_X[[k]] = apply(cor(A[[k]], variates.A[[k]])^2, 2, mean)
                
                if (is.null(tau))
                {
                    names(loadings.partial.A[[k]]) =
                        names(variates.partial.A[[k]]) = levels(study)
                    
                    for (m in 1:nlevels(study))
                    {
                        rownames(loadings.partial.A[[k]][[m]]) = colnames(A[[k]])
                        colnames(loadings.partial.A[[k]][[m]]) =
                            paste0("comp", 1:max(ncomp))
                        rownames(variates.partial.A[[k]][[m]]) = mean_centered.rownames.study[[m]]
                        colnames(variates.partial.A[[k]][[m]]) =
                            paste0("comp", 1:max(ncomp))
                    }
                }
            }
            
            variates.A = shave.matlist(variates.A, ncomp)
            
            if(all.outputs)
            {
                # AVE
                outer = matrix(unlist(AVE_X), nrow = max(ncomp))
                for (j in 1 : max(ncomp))
                    AVE_outer[j] = sum(pjs * outer[j, ])/sum(pjs)
                AVE_X = shave.veclist(AVE_X, ncomp)
                AVE = list(AVE_X = AVE_X, AVE_outer = AVE_outer,
                           AVE_inner = AVE_inner)
                names(AVE$AVE_X) = names(A)
                
                loadings.Astar = shave.matlist(loadings.Astar, ncomp)
                
                #calcul explained variance
                A_split=lapply(A, study_split, study) #split the data per study
                expl.A=lapply(1:length(A),function(x){
                    if (nlevels(study) == 1)
                    {
                        temp = explained_variance(A[[x]],
                                                  variates = variates.A[[x]], ncomp = ncomp[[x]])
                    } else {
                        temp = lapply(1:nlevels(study), function(y){
                            explained_variance(A_split[[x]][[y]],
                                               variates = variates.partial.A[[x]][[y]],
                                               ncomp = ncomp[[x]])})
                        temp[[length(temp)+1]] = explained_variance(A[[x]],
                                                                    variates = variates.A[[x]], ncomp = ncomp[[x]])
                        names(temp) = c(levels(study), "all data")
                    }
                    temp
                })
                names(expl.A) = names(A)
            } else {
                expl.A = NULL
                AVE = NULL
            }
            ### Start: output
            names(loadings.A) = names(variates.A) = names(A)
            
            if (is.null(tau))
                names(loadings.partial.A) = names(variates.partial.A) = names(A)
            
            names = lapply(1:J, function(x) {colnames(A[[x]])})
            names(names) = names(A)
            names[[length(names) + 1]] = row.names(A[[1]])
            names(names)[length(names)] = "indiv"
        } else {
            # multiple models (tune)
            
            #### any model
            # loadings.A[[block]][1:p, all.keepA.tested]
            # variates.A[[block]][1:n, all.keepA.tested]
            
            #### a single model
            # loadings.partial.A[[block]][[study]][, 1:ncomp]
            # variates.partial.A[[block]][[study]][, 1:ncomp]
            # loadings.Astar[[block]][, 1:ncomp]
            
            
            keepA.names = unlist(lapply(1:N, function(x){
                paste(paste0("comp",x),apply(keepA[[x]],1,function(x)
                    paste(x,collapse="_")), sep=":")
                
            }))
            
            for(k in 1:J)
                colnames(loadings.A[[k]]) = colnames(variates.A[[k]]) = keepA.names
            
            if (length(iter) == length(A)) {
                names(iter) = names(A)
            }
            names(loadings.A) =  names(variates.A) = names(A)
            
            expl.A = NULL
            AVE = NULL
            
        }
        
        out = list(A = A, indY = indY, ncomp = ncomp, mode = mode,
                   keepA = keepA,
                   variates = variates.A, loadings = loadings.A,
                   variates.partial= if(is.null(tau)) {variates.partial.A} ,
                   loadings.partial= if(is.null(tau)) {loadings.partial.A},
                   loadings.star = loadings.Astar,
                   names = list(sample = row.names(A[[1]]), colnames = lapply(A, colnames),
                                blocks = names(A)),
                   tol = tol, iter=iter, max.iter=max.iter,
                   design = design,
                   scheme = scheme,  crit = crit, AVE = AVE, mat.c = mat.c,
                   #defl.matrix = defl.matrix,
                   init = init,
                   scale = scale, tau = if(!is.null(tau)) tau.rgcca, study = study,
                   prop_expl_var = expl.A)
        ### End: Output
        
        return(out)
    }

# ------------------------------------------------------------------------------
# sgccak - Runs sgccak() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ------------------------------------------------------------------------------

sparse.mint.block_iteration = function (A, design, study = NULL, loadings.A,
                                        keepA = NULL,
                                        scheme = "horst", max.iter = 100, tol = 1e-06,
                                        misdata = NULL, is.na.A = NULL, ind.NA = NULL,
                                        penalty=NULL, all.outputs = FALSE)
{
    
    # keepA is a list of length the number of blocks. Each entry is a vector of
    #   numbers: variables to select for that block (component is fixed)
    # study is a vector
    # no check needed as this function is only used in
    #   internal_mint.block, in which the checks are conducted
    
    
    ### Start: Initialization parameters
    J = length(A)
    J2 = J-1
    pjs = sapply(A, NCOL)
    AVE_X = rep(0, J)
    if (!is.null(penalty))
        penalty = penalty * sqrt(pjs)
    
    iter = 1
    converg = crit = numeric()
    variates.A = Z = matrix(0, NROW(A[[1]]), J)
    
    g = function(x) switch(scheme, horst = x, factorial = x^2,
                           centroid = abs(x))
    
    # study split
    A_split = lapply(A, study_split, study)
    
    n = lapply(A_split, function(x){lapply(x,nrow)})
    p = lapply(A,ncol)
    
    nlevels_study = nlevels(study)
    ### End: Initialization parameters
    
    
    ### End: Initialisation "a" vector
    variates.partial.A.comp = NULL
    loadings.partial.A.comp = list()
    for (q in 1:J)
    {
        if(misdata[q])
        {
            loadings.temp = loadings.A[[q]]
            variates.A.temp = A[[q]] %*% loadings.temp
            
            # we only want the diagonal,
            #   which is the norm of each column of temp
            # loadings.A.norm = crossprod(temp)
            # variates.A[, q] = variates.A.temp / diag(loadings.A.norm)
            #only calculating the ones where there's a NA
            d.variates.A.norm = rep(crossprod(loadings.temp),
                                    length(variates.A.temp))
            
            if(length(ind.NA[[q]])>0) # should always be true
            {
                temp = drop(loadings.temp) %o% rep(1, length(ind.NA[[q]]))
                #p*n -> p * where there are NA
                temp[t(is.na.A[[q]][ind.NA[[q]],,drop=FALSE])] = 0
                d.variates.A.norm[ind.NA[[q]]] = apply(temp,2, crossprod)
            }
            
            variates.A[, q] = variates.A.temp / d.variates.A.norm
            
            # we can have 0/0, so we put 0
            a = is.na(variates.A[, q])
            if (any(a))
                variates.A[a, q] = 0
            
        }else{
            variates.A[, q] = A[[q]]%*%loadings.A[[q]]
        }
        loadings.A[[q]] = l2.norm(as.vector(loadings.A[[q]]))
        loadings.partial.A.comp[[q]] = list()
    }
    loadings.A_old = loadings.A
    
    ### Start Algorithm 1 Sparse generalized canonical analysis (See Variable
    # selection for generalized canonical correlation analysis (Tenenhaus))
    repeat {
        # variates.Aold = variates.A
        for (q in 1:J)
        {
            ### Start : !!! Impact of the diag of the design matrix !!! ###
            if (scheme == "horst")
                CbyCovq = design[q, ]
            
            if (scheme == "factorial")
                CbyCovq = design[q, ] * cov2(variates.A, variates.A[, q])
            
            if (scheme == "centroid")
                CbyCovq = design[q, ] * sign(cov2(variates.A, variates.A[, q]))
            ### End : !!! Impact of the diag of the design matrix !!! ###
            
            ### Step A start: Compute the inner components
            Z[, q] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
            Z_split = study_split(Z[,q,drop=FALSE],study)
            # split Z by the study factor
            
            ### Step A end: Compute the inner components
            
            
            ### Step B start: Computer the outer weight ###
            temp=0
            for (m in 1:nlevels_study)
            {
                loadings.partial.A.comp[[q]][[m]] =
                    crossprod(A_split[[q]][[m]],Z_split[[m]])
                temp=temp+loadings.partial.A.comp[[q]][[m]]
            }
            loadings.A[[q]] = temp
            
            # sparse using keepA / penalty
            if (!is.null(penalty))
            {
                loadings.A[[q]] = sparsity(loadings.A[[q]], keepA = NULL,
                                           penalty = penalty[q])
            }else{
                loadings.A[[q]] = sparsity(loadings.A[[q]], keepA[[q]],
                                           penalty = NULL)
            }
            
            loadings.A[[q]]=l2.norm(as.vector(loadings.A[[q]]))
            
            ### Step B end: Computer the outer weight ###
            if(misdata[q])
            {
                variates.A.temp = A[[q]] %*% loadings.A[[q]]
                d.variates.A.norm = rep(crossprod(loadings.A[[q]]),
                                        length(variates.A.temp))
                
                if(length(ind.NA[[q]])>0)
                {
                    temp = drop(loadings.A[[q]]) %o% rep(1, length(ind.NA[[q]]))
                    temp[t(is.na.A[[q]][ind.NA[[q]],,drop=FALSE])] = 0
                    d.variates.A.norm[ind.NA[[q]]] = apply(temp,2, crossprod)
                }
                variates.A[, q] = variates.A.temp / d.variates.A.norm
                
                # we can have 0/0, so we put 0
                a = is.na(variates.A[, q])
                if (any(a))
                    variates.A[a, q] = 0
                
            }else{
                variates.A[, q] =  A[[q]]%*%loadings.A[[q]]
            }
            
        }
        
        crit[iter] = sum(design * g(cov2(variates.A)))
        
        if (iter > max.iter)
            warning("The SGCCA algorithm did not converge", call. = FALSE)
        
        ### Start: Match algorithm with mixOmics algo (stopping point)
        diff.value <- max(sapply(1:J, function(x){crossprod(loadings.A[[x]] -
                                                                loadings.A_old[[x]])}))
        
        if (diff.value < tol | iter > max.iter)
            break
        ### End: Match algorithm with mixOmics algo (stopping point)
        
        loadings.A_old = loadings.A
        iter = iter + 1
    }
    ### End Algorithm 1 (See Variable selection for generalized canonical
    # correlation analysis (Tenenhaus))
    
    #calculation variates.partial.A.comp
    variates.partial.A.comp = apply(variates.A, 2, study_split, study)
    
    if(all.outputs){
        AVE_inner = sum(design * cor(variates.A)^2/2)/(sum(design)/2)
    } else{
        AVE_inner = NULL
    }
    
    names(loadings.A) = colnames(variates.A) =
        names(variates.partial.A.comp) = names(A)
    
    result = list(variates.A = variates.A, loadings.A = loadings.A, crit =
                      crit[which(crit != 0)],
                  AVE_inner = AVE_inner, loadings.partial.A.comp = loadings.partial.A.comp,
                  variates.partial.A.comp = variates.partial.A.comp, iter = iter)
    return(result)
}

# ------------------------------------------------------------------------------
# rgccak - Runs sgccak() modified from RGCCA
#   inputs: A - list of datasets each with the same number of rows (samples)
#           design - design matrix
#           ncomp - vector specifying number of components to keep per datasets
#   outputs:
# ------------------------------------------------------------------------------


sparse.rgcca_iteration = function (A, design, tau = "optimal", scheme = "horst",
                                   scale = FALSE, max.iter = 100, init = "svd.single", tol = .Machine$double.eps,
                                   keepA = NULL, penalty = NULL, all.outputs = FALSE)
{
    ### Start: Initialisation parameters
    A = lapply(A, as.matrix)
    J = length(A)
    n = NROW(A[[1]])
    pjs = sapply(A, NCOL)
    variates.A = matrix(0, n, J)
    if (!is.null(penalty))
        penalty = penalty * sqrt(pjs)
    ### End: Initialisation parameters
    
    if (!is.numeric(tau))
        tau = sapply(A, tau.estimate)
    
    loadings.A = alpha = M = Minv = K = list()
    which.primal = which((n >= pjs) == 1)
    which.dual = which((n < pjs) == 1)
    
    if (init == "svd.single")
    {
        for (j in which.primal)
            loadings.A[[j]] = initsvd(lapply(j, function(x)
            {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])
        
        for (j in which.dual)
        {
            alpha[[j]] = initsvd(lapply(j, function(x)
            {replace(A[[x]], is.na(A[[x]]), 0)})[[1]])
            K[[j]] = A[[j]] %*% t(A[[j]])
        }
    } else {
        stop("init should be 'svd.single'.")
    }
    
    N = n
    for (j in 1 : J)
    {
        if (j %in% which.primal)
        {
            M[[j]] = ginv(tau[j] * diag(pjs[j]) + (1 - tau[j]) * cov2(A[[j]]))
            loadings.A[[j]] = drop(1/sqrt(t(loadings.A[[j]]) %*% M[[j]] %*%
                                              loadings.A[[j]])) * M[[j]] %*% loadings.A[[j]]
        }
        
        if (j %in% which.dual)
        {
            M[[j]] = tau[j] * diag(n) + (1 - tau[j])/(N) * K[[j]]
            Minv[[j]] = ginv(M[[j]])
            alpha[[j]] = drop(1/sqrt(t(alpha[[j]]) %*% M[[j]] %*% K[[j]] %*%
                                         alpha[[j]])) * alpha[[j]]
            loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
        }
        variates.A[, j] = A[[j]] %*% loadings.A[[j]]
    }
    
    iter = 1
    converg = crit = numeric()
    Z = matrix(0, NROW(A[[1]]), J)
    loadings.A_old = loadings.A
    g = function(x)
        switch(scheme, horst = x, factorial = x^2, centroid = abs(x))
    
    repeat {
        variates.Aold = variates.A
        
        for (j in c(which.primal, which.dual))
        {
            
            if (scheme == "horst")
                CbyCovq = design[j, ]
            
            if (scheme == "factorial")
                CbyCovq = design[j, ] * cov2(variates.A, variates.A[, j])
            
            if (scheme == "centroid")
                CbyCovq = design[j, ] * sign(cov2(variates.A, variates.A[, j]))
            
            # Compute the inner components
            Z[, j] = rowSums(mapply("*", CbyCovq, as.data.frame(variates.A)))
            
            # Computer the outer weight
            if (j %in% which.primal)
                loadings.A[[j]] = drop(1/sqrt(t(Z[, j]) %*% A[[j]] %*% M[[j]] %*%
                                                  t(A[[j]]) %*% Z[, j])) * (M[[j]] %*% t(A[[j]]) %*% Z[, j])
            
            # Compute the outer weight
            if (j %in% which.dual)
            {
                alpha[[j]] = drop(1/sqrt(t(Z[, j]) %*% K[[j]] %*% Minv[[j]] %*%
                                             Z[, j])) * (Minv[[j]] %*% Z[, j])
                loadings.A[[j]] = t(A[[j]]) %*% alpha[[j]]
            }
            
            # sparse using keepA / penalty
            if (!is.null(keepA) || !is.null(penalty))
            {
                temp.norm = norm2(loadings.A[[j]])
                if (!is.null(keepA))
                {
                    loadings.A[[j]] = sparsity(loadings.A = loadings.A[[j]],
                                               keepA = keepA[[j]], penalty = NULL)
                } else if (!is.null(penalty)) {
                    loadings.A[[j]] = sparsity(loadings.A = loadings.A[[j]],
                                               keepA = NULL, penalty = penalty[j])
                }
                loadings.A[[j]] = (loadings.A[[j]]/norm2(loadings.A[[j]]))*
                    temp.norm
            }
            
            # Update variate
            variates.A[, j] = A[[j]] %*% loadings.A[[j]]
        }
        
        crit[iter] = sum(design * g(cov2(variates.A)))
        
        if (iter > max.iter)
            warning("The RGCCA algorithm did not converge")
        
        ### Start: Match algorithm with mixOmics algo (stopping point)
        if (max(sapply(1:J, function(x){crossprod(loadings.A[[x]] -
                                                  loadings.A_old[[x]])})) < tol | iter > max.iter)
            break
        ### End: Match algorithm with mixOmics algo (stopping point)
        
        loadings.A_old = loadings.A
        iter = iter + 1
    }
    
    #if (verbose)
    #plot(crit, xlab = "iteration", ylab = "criteria")
    
    if(all.outputs){
        AVE_inner = sum(design * cor(variates.A)^2/2)/(sum(design)/2)
    } else{
        AVE_inner = NULL
    }
    
    result = list(variates.A = variates.A, loadings.A = loadings.A,
                  crit = crit[which(crit != 0)], AVE_inner = AVE_inner, design = design,
                  tau = tau, scheme = scheme,iter=iter, keepA = keepA)
    return(result)
}

#############################################################################################################
# Authors:
#   Amrit Singh, University of British Columbia, Vancouver.
#   Florian Rohart, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 01-04-2015
# last modified: 27-05-2016
#
# Copyright (C) 2015
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# ----------------------------------------------------------------------------------------------------------
# perf.sgccda - Function to evaluate the performance of the fitted PLS (cross-validation)
#   inputs: object - object obtain from running sgccda or splsda
#           dist - to evaluate the classification performance
#           validation - type of validation
#           folds - number of folds if validation = "Mfold"
# ----------------------------------------------------------------------------------------------------------
#' @rdname perf
#' @importFrom utils relist
#' @method perf sgccda
#' @export
perf.sgccda <- 
  function (object,
            dist = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),
            validation = c("Mfold", "loo"),
            folds = 10,
            nrepeat = 1,
            auc = FALSE,
            progressBar = FALSE,
            signif.threshold = 0.01,
            BPPARAM = SerialParam(),
            seed = NULL,
            ...)
  {

    if (hasArg('cpus')) #defunct
    {
    stop("'cpus' has been replaced by BPPARAM. See documentation.")  
    }
    
    ### Start: Initialization parameters
    BPPARAM$RNGseed <- seed
    set.seed(seed)
    X = object$X; level.Y = object$names$colnames$Y;
    J = length(X)
    Y = object$Y#ind.mat; Y = map(Y); Y = factor(Y, labels = level.Y)
    n = nrow(X[[1]]);
    indY = object$indY
    max.iter =  object$max.iter
    ncomp <- max(object$ncomp[-indY])
    near.zero.var = !is.null(object$nzv) # if near.zero.var was used, we set it to TRUE. if not used, object$nzv is NULL
    
    if (any(dist == "all"))
    {
      dist.select = c("max.dist", "centroids.dist", "mahalanobis.dist")
    } else {
      dist.select = dist
    }
    ### End: Initialization parameters
    
    dist = match.arg(dist.select, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    
    ### Start: Check parameter validation / set up sample
    
    #-- check significance threshold
    signif.threshold <- .check_alpha(signif.threshold)
    
    if (length(validation) > 1 )
      validation = validation [1]
    
    if (!(validation %in% c("Mfold", "loo")))
      stop("Choose 'validation' among the two following possibilities: 'Mfold' or 'loo'")

    #-- folds supplied as a list --#
    # The list is validated once here, before the repeats are dispatched to
    # bplapply, so that errors and warnings reach the user directly.
    if (is.list(folds) && validation == "Mfold")
    {
      folds = .check_folds(folds, n = n, Y = Y)
      if (nrepeat > 1)
      {
        warning("'folds' is supplied as a list, so the cross-validation does not need to be repeated. ",
                "'nrepeat' is set to '1'.", call. = FALSE)
        nrepeat = 1
      }
    }
    
    #-- tells which variables are selected in the blocks --#
    
    keepX = object$keepX
    
    error.mat = error.mat.class = Y.all = predict.all = Y.predict = list.features = final.features = weights = crit = list()
    if (length(X) > 1)
      Y.mean = Y.mean.res = Y.weighted.vote = Y.weighted.vote.res = Y.vote = Y.vote.res = Y.WeightedPredict = Y.WeightedPredict.res = list()
    
    if (progressBar ==  TRUE) {
      cat(sprintf("\nPerforming repeated cross-validation with nrepeat = %s...\n", nrepeat))
      pb = txtProgressBar(style = 3)
    }
    
    ## ------------ function to perform repeat CV
    repeat_cv_perf.diablo <- function(nrep)
    {
      if (progressBar ==  TRUE) {
        setTxtProgressBar(pb = pb, value = nrep/nrepeat)
      }
      
      #-- define the folds --#
      if (validation ==  "Mfold")
      {
        if (is.list(folds))
        {
          # the list was validated in the calling function: use as is
          M = length(folds)
        } else if (is.null(folds) || !is.numeric(folds) || folds < 2 || folds > n)
        {
          stop("Invalid number of folds.")
        } else {
          M = round(folds)
          temp = stratified.subsampling(Y, folds = M)
          folds = temp$SAMPLE
          if(temp$stop > 0) # to show only once
            warning("At least one class is not represented in one fold, which may unbalance the error rate.\n  Consider a number of folds lower than the minimum in table(Y): ", min(table(Y)))
          
        }
      } else if (validation ==  "loo") {
        folds = split(1:n, rep(1:n, length = n))
        M = n
      } else {
        stop("validation can be only 'Mfold' or 'loo'")
      }
      M = length(folds)
      
      ### Start: Check parameter validation / set up sample
      
      ### Start: Training samples (X.training and Y.training) and Test samples (X.test / Y.test)
      X.training = lapply(folds, function(x){out=lapply(1:J, function(y) {X[[y]][-x, ]}); names(out) = names(X); out}) #need to name the block for prediction
      Y.training = lapply(folds, function(x) {Y[-x]});
      
      X.test = lapply(folds, function(x){out=lapply(1:J, function(y) {X[[y]][x, , drop = FALSE]}); names(out) = names(X); out})#need to name the block for prediction
      Y.test = lapply(folds, function(x) {Y[x]})
      ### End: Training samples (X.training and Y.training) and Test samples (X.test / Y.test)
      
      ### Estimation models
      
      ## helper function to return appropriate block.splsda args for CV
      block.splsda.args <- function(ind) {
        return(list(X = X.training[[ind]], Y = Y.training[[ind]], ncomp = ncomp, 
                    keepX = keepX,
                    design = object$design, max.iter =  max.iter, tol = object$tol,
                    near.zero.var = near.zero.var))
      }
      
      model = lapply(1:M, function(x) {suppressWarnings(do.call("block.splsda", block.splsda.args(x)))})
      ### CV AUC
      
      if (auc) {
        auc.rep.block.comp.fold <-
          lapply(1:M, function(x) {
            auroc.sgccda(
              object = model[[x]],
              plot = FALSE,
              print = FALSE
            )
          })
        
        
        blocks <- .name_list(names(auc.rep.block.comp.fold[[1]]))
        comps <- .name_list(names(auc.rep.block.comp.fold[[1]][[1]]))
        ## average AUC for each block and component across folds
        auc.rep.block.comp <- lapply(blocks, function(block) {
          lapply(comps, function(comp){
            Reduce("+", lapply(auc.rep.block.comp.fold, function(fold){
              fold[[block]][[comp]]
            }))/M
          })
        }
        )
        ## combined average AUC for all blocks and component for the repeat
        auc.rep.comp <- lapply(comps, function(comp){
          Reduce("+", lapply(blocks, function(block){
            auc.rep.block.comp[[block]][[comp]]
          }))/length(blocks)
        })
      }
      
      
      ### Retrieve convergence criterion
      crit = lapply(1:M, function(x){model[[x]]$crit})
      ### Retrieve weights
      weights = lapply(1:M, function(fold_i){
        fold_weights <- model[[fold_i]]$weights
        fold_weights$fold <- fold_i
        fold_weights$rep <- nrep
        fold_weights$block <- rownames(fold_weights)
        rownames(fold_weights) <- NULL
        fold_weights
        })
      weights <- Reduce(rbind , weights)
      
      ### Retrieve selected variables per component
      features = lapply(1 : J, function(x)
      {
        lapply(1 : object$ncomp[x], ### List level / ncomp level
               function(y)
               {
                 unlist(lapply(1:M,  ### Validation level
                               function(z)
                               {
                                 if (is.null(colnames(X[[x]])))
                                 {
                                   paste0("X", which(model[[z]]$loadings[[x]][, y] != 0))
                                 } else {
                                   if (!is.null(colnames(X[[x]])))
                                     colnames(X[[x]])[model[[z]]$loadings[[x]][, y] != 0]
                                 }
                               }
                 ))
               })
      })
      
      ### Start: Analysis feature selection
      # Statistics: stability
      list.features = lapply(1 : J, function(x){lapply(features[[x]], function(y){(sort(table(factor(y))/M, decreasing = TRUE))})})
      
      # Statistics: original model (object)
      final.features = lapply(1 : J, function(x)
      {
        lapply(1 : object$ncomp[x],function(y)
        {
          temp = as.data.frame(object$loadings[[x]][which(object$loadings[[x]][, y] != 0), y, drop = FALSE])
          if (is.null(colnames(X[[x]])))
          {
            row.names(temp) = paste0("X", which(object$loadings[[x]][, y] != 0))
          } else {
            if (!is.null(colnames(X[[x]])))
              row.names(temp) = colnames(X[[x]])[which(object$loadings[[x]][, y] != 0)]
          }
          names(temp) = "value.var"
          return(temp[sort(abs(temp[, 1]), index.return = TRUE, decreasing = TRUE)$ix, 1, drop = FALSE])
        })
      })
      
      list.features = lapply(1 : J, function(x){ names(list.features[[x]]) = paste0("comp", 1 : object$ncomp[x])
      return(list.features[[x]])})
      
      final.features = lapply(1 : J, function(x){ names(final.features[[x]]) = paste0("comp", 1 : object$ncomp[x])
      return(final.features[[x]])})
      ### End: Analysis feature selection
      
      ### Warning: no near.zero.var applies with sgcca
      
      ### Start: Prediction (score / class) sample test
      # Prediction model on test dataset
      predict.all = lapply(1:M, function(x) {predict.block.spls(model[[x]], X.test[[x]], dist = "all")})
      
      # Retrieve class prediction
      Y.predict = lapply(1:M, function(x) {predict.all[[x]]$class})
      
      ## Start: retrieve score for each component
      # Keep score values
      Y.all = lapply(1:M, function(x) {predict.all[[x]]$predict})
      
      # Reorganization list Y.all data / ncomp / folds
      Y.all = lapply(1 : J, function(x)
      {
        lapply(1 : object$ncomp[x], function(y)
        {
          lapply(1:M, function(z)
          {
            Y.all[[z]][[x]][, , y]
          })
        })
      })
      # Merge score
      Y.all = lapply(1 : J, function(x)
      {
        lapply(1 : object$ncomp[x], function(y)
        {
          do.call(rbind, Y.all[[x]][[y]])
        })
      })
      
      # Define row.names
      Y.all = lapply(1 : J, function(x)
      {
        lapply(1 : object$ncomp[x], function(y)
        {
          row.names(Y.all[[x]][[y]]) = row.names(X[[x]])[unlist(folds)]
          return(Y.all[[x]][[y]])
        })
      })
      
      # Define colnames
      Y.all = lapply(1 : J, function(x)
      {
        lapply(1 : object$ncomp[x], function(y)
        {
          colnames(Y.all[[x]][[y]]) = levels(Y)
          return(Y.all[[x]][[y]])
        })
        names(Y.all[[x]])=paste0("comp", 1:object$ncomp[x])
        return(Y.all[[x]])
      })
      ## End: retrieve score for each component
      
      #save(list=ls(),file="temp.Rdata")
      
      ## Start: retrieve class for each component
      # Reorganization input data / folds / dist.select
      Y.predict = lapply(1 : J, function(x)
      {
        lapply(1:M, function(y)
        {
          lapply(which(c("max.dist", "centroids.dist", "mahalanobis.dist") %in% dist.select), function(z)
          {
            Y.predict[[y]][[z]][[x]]
          })
        })
      })
      
      # Define row.names
      Y.predict = lapply(1 : J, function(x)
      {
        lapply(1:M, function(y)
        {
          lapply(1 : length(dist.select), function(z)
          {
            if (!is.null(row.names(X[[x]])[folds[[y]]]))
            {
              row.names(Y.predict[[x]][[y]][[z]]) = row.names(X[[x]])[folds[[y]]]
            } else {
              row.names(Y.predict[[x]][[y]][[z]]) = paste0("Ind", unlist(folds))
            }
            return(Y.predict[[x]][[y]][[z]])
          })
        })
      })
      
      # Reorganization list input / dist.select / folds
      Y.predict = lapply(1 : J, function(x)
      {
        lapply(1 : length(dist.select), function(y)
        {
          lapply(1:M, function(z)
          {
            Y.predict[[x]][[z]][[y]]
          })
        })
      })
      # Merge Score
      Y.predict = lapply(1 : J, function(x)
      {
        lapply(1 : length(dist.select), function(y)
        {
          do.call(rbind, Y.predict[[x]][[y]])
        })
      })
      
      # Define names
      Y.predict = lapply(1 : J, function(x)
      {
        names(Y.predict[[x]]) = dist.select
        return(Y.predict[[x]])
      })
      ## End: retrieve class for each component
      ### End: Prediction (score / class) sample test
      
      ### Start: Estimation error rate
      ## Start: Estimation overall error rate
      #Statistics overall error rate
      error.mat = lapply(1 : J, function(x)
      {
        lapply(dist.select , function(y)
        {
          apply(Y.predict[[x]][[y]], 2, function(z)
          {
            1 - sum(diag(table(factor(z, levels = levels(Y)), Y[unlist(folds)])))/length(Y)
          })
        })
      })
      
      # Merge error rate according to dist
      error.mat = lapply(1 : J, function(x)
      {
        do.call(cbind, error.mat[[x]])
      })
      
      # Define name
      error.mat = lapply(1 : J, function(x)
      {
        colnames(error.mat[[x]]) = dist.select
        return(error.mat[[x]])
      })
      ## End: Estimation overall error rate
      
      ## Start: Estimation error rate per class
      # Statistics error rate per class
      error.mat.class = lapply(1 : J, function(x)
      {
        lapply(dist.select , function(y)
        {
          apply(Y.predict[[x]][[y]], 2, function(z)
          {
            temp = diag(table(factor(z, levels = levels(Y)), Y[unlist(folds)]))
            1 - c(temp/summary(Y))#, sum(temp)/length(Y))
          })
        })
      })
      
      # Define names
      error.mat.class = lapply(1 : J, function(x)
      {
        names(error.mat.class[[x]]) = dist.select
        return(error.mat.class[[x]])
      })
      
      ## End: Estimation error rate per class
      ### End: Estimation error rate
      
      if (!is.null(names(X)))
      {
        names(error.mat) = names(X); names(error.mat.class) = names(X)
        names(list.features) = names(X); names(final.features) = names(X);
        names(Y.all) = names(X); names(Y.predict) = names(X)
      }
      
      
      ### Start: Supplementary analysis for sgcca
      if (length(X) > 1)
      {
        ###------------------------------------------------------------###
        ### Start: class of Average prediction
        # Reorganization dist.select / folds
        Y.mean =lapply(1:M, function(y)
        {
          predict.all[[y]][["AveragedPredict.class"]][[1]]
        })
        # Merge Score
        Y.mean = do.call(rbind, Y.mean)
        
        # Sort matrix
        Y.mean = Y.mean[sort(unlist(folds), index.return = TRUE)$ix, , drop = FALSE]
        
        
        # Estimation error.rate
        #Y.mean.res = sapply(1:max(object$ncomp[-(J + 1)]), function(x){temp = diag(table(factor(Y.mean[, x], levels = c(1:nlevels(Y))), Y))
        #                                                              c(temp/summary(Y), sum(temp)/length(Y))})
        Y.mean.res = sapply(1:ncomp, function(x)
        {
          mat = table(factor(Y.mean[, x], levels = levels(Y)), Y)
          mat2 <- mat
          diag(mat2) <- 0
          err = c(c(colSums(mat2)/summary(Y), sum(mat2)/length(Y)), mean(colSums(mat2)/colSums(mat)))
        })
        
        #Y.mean.res = t(Y.mean.res)
        colnames(Y.mean.res) = paste0("comp", 1:ncomp)
        row.names(Y.mean.res) = c(levels(Y), "Overall.ER", "Overall.BER")
        ### End: Average prediction
        ###------------------------------------------------------------###
        
        
        
        ###------------------------------------------------------------###
        ### Start: class of Weighted prediction
        # Reorganization dist.select / folds
        Y.WeightedPredict =lapply(1:M, function(y)
        {
          predict.all[[y]][["WeightedPredict.class"]][[1]]
        })
        # Merge Score
        Y.WeightedPredict = do.call(rbind, Y.WeightedPredict)
        
        # Sort matrix
        Y.WeightedPredict = Y.WeightedPredict[sort(unlist(folds), index.return = TRUE)$ix, , drop = FALSE]
        
        
        # Estimation error.rate
        #Y.mean.res = sapply(1:max(object$ncomp[-(J + 1)]), function(x){temp = diag(table(factor(Y.mean[, x], levels = c(1:nlevels(Y))), Y))
        #                                                              c(temp/summary(Y), sum(temp)/length(Y))})
        Y.WeightedPredict.res = sapply(1:ncomp, function(x)
        {
          mat = table(factor(Y.WeightedPredict[, x], levels = levels(Y)), Y)
          mat2 <- mat
          diag(mat2) <- 0
          err = c(c(colSums(mat2)/summary(Y), sum(mat2)/length(Y)), mean(colSums(mat2)/colSums(mat)))
        })
        
        #Y.mean.res = t(Y.mean.res)
        colnames(Y.WeightedPredict.res) = paste0("comp", 1:ncomp)
        row.names(Y.WeightedPredict.res) = c(levels(Y), "Overall.ER", "Overall.BER")
        ### End: Average prediction
        ###------------------------------------------------------------###
        
        
        ###------------------------------------------------------------###
        ## Start: retrieve (weighted) vote for each component
        # Reorganization dist.select / folds
        Y.weighted.vote = lapply(dist.select, function(x)
        {
          lapply(1:M, function(y)
          {
            predict.all[[y]][["WeightedVote"]][[x]]
          })
        })
        # Merge Score
        Y.weighted.vote = lapply(1 : length(dist.select), function(x)
        {
          do.call(rbind, Y.weighted.vote[[x]])
        })
        
        # Sort matrix
        Y.weighted.vote = lapply(1 : length(dist.select), function(x)
        {
          Y.weighted.vote[[x]][sort(unlist(folds), index.return = TRUE)$ix, , drop = FALSE]
        })
        
        names(Y.weighted.vote) = dist.select
        
        ## End: retrieve (weighted) vote for each component
        ### End: Prediction (score / class) sample test
        
        
        ## subjects with NA are considered false
        Y.weighted.vote.res = lapply(1 : length(dist.select), function(x)
        {
          apply(Y.weighted.vote[[x]], 2, function(y)
          {
            y[is.na(y)] <- nlevels(Y)+5   ## adding a new level for unsure subjects (replacing NA with this level)
            temp=table(factor(y, levels = c(levels(Y), nlevels(Y)+5)), Y)
            diag(temp) <- 0
            err = c(colSums(temp)/summary(Y), sum(temp)/length(Y), mean(colSums(temp)/summary(Y)))
            return(err=err)
          })
        })
        
        
        Y.weighted.vote.res = lapply(1 : length(dist.select), function(x)
        {
          colnames(Y.weighted.vote.res[[x]]) = paste0("comp", 1:max(object$ncomp[-(J + 1)]))
          row.names(Y.weighted.vote.res[[x]]) = c(levels(Y), "Overall.ER", "Overall.BER")
          return((Y.weighted.vote.res[[x]]))
        })
        names(Y.weighted.vote) = dist.select; names(Y.weighted.vote.res) = dist.select
        ###------------------------------------------------------------###
        
        
        
        ###------------------------------------------------------------###
        ## Start: retrieve Majority Vote (non weighted) for each component
        # Reorganization dist.select / folds
        Y.vote = lapply(dist.select, function(x)
        {
          lapply(1:M, function(y)
          {
            predict.all[[y]][["MajorityVote"]][[x]]
          })
        })
        # Merge Score
        Y.vote = lapply(1 : length(dist.select), function(x)
        {
          do.call(rbind, Y.vote[[x]])
        })
        
        # Sort matrix
        Y.vote = lapply(1 : length(dist.select), function(x)
        {
          Y.vote[[x]][sort(unlist(folds), index.return = TRUE)$ix, , drop = FALSE]
        })
        
        names(Y.vote) = dist.select
        
        ## End: retrieve (weighted) vote for each component
        ### End: Prediction (score / class) sample test
        
        
        ## subjects with NA are considered false
        Y.vote.res = lapply(1 : length(dist.select), function(x)
        {
          apply(Y.vote[[x]], 2, function(y)
          {
            y[is.na(y)] <- nlevels(Y)+5   ## adding a new level for unsure subjects (replacing NA with this level)
            temp=table(factor(y, levels = c(levels(Y), nlevels(Y)+5)), Y)
            diag(temp) <- 0
            err = c(colSums(temp)/summary(Y), sum(temp)/length(Y), mean(colSums(temp)/summary(Y)))
            return(err=err)
          })
        })
        
        
        Y.vote.res = lapply(1 : length(dist.select), function(x)
        {
          colnames(Y.vote.res[[x]]) = paste0("comp", 1:max(object$ncomp[-(J + 1)]))
          row.names(Y.vote.res[[x]]) = c(levels(Y), "Overall.ER", "Overall.BER")
          return((Y.vote.res[[x]]))
        })
        names(Y.vote) = dist.select; names(Y.vote.res) = dist.select
        
        ## End: retrieve non weighted vote for each component
        ###------------------------------------------------------------###
        repeat_cv_res <- list(
          error.mat = error.mat,
          error.mat.class = error.mat.class,
          Y.mean = Y.mean,
          Y.mean.res = Y.mean.res,
          Y.WeightedPredict.res = Y.WeightedPredict.res,
          Y.weighted.vote.res = Y.weighted.vote.res,
          Y.vote.res = Y.vote.res,
          Y.all = Y.all,
          predict.all = predict.all,
          Y.predict = Y.predict,
          list.features = list.features,
          final.features = final.features,
          crit = crit,
          weights = weights
        )
        
        if (auc) {
          repeat_cv_res$auc.rep.comp <- auc.rep.comp
        }
        
        return(repeat_cv_res)
      }
      ### End: Supplementary analysis for sgcca
    } ### end nrepeat
    
    ## a list of nreps for lapply
    nrep_list <- as.list(seq_len(nrepeat))
    names(nrep_list) <- paste0("nrep", nrep_list)
      
    # Execute using bplapply, whether in serial or parallel mode
    repeat_cv_perf.diablo_res <- BiocParallel::bplapply(nrep_list, function(nrep) {
      repeat_cv_perf.diablo(nrep)
    }, BPPARAM = BPPARAM)
    
    repeat_cv_perf.diablo_res  <- .unlist_repeat_cv_output(repeat_cv_perf.diablo_res)
    auc.rep.comp <- NULL ## R check pass
    list2env(repeat_cv_perf.diablo_res, envir = environment())
    rm(repeat_cv_perf.diablo_res)
    
    
    ###------------------------------------------------------------###
    ## we want to average the error per dataset over nrepeat for
    # error.rate, error.rate.per.class, AveragedPredict.error.rate, WeightedPredict.error.rate, MajorityVote.error.rate, WeightedVote.error.rate
    # with each time: error.rate, error.rate.sd, error.rate.all
    ###------------------------------------------------------------###
    
    if(nrepeat > 1)
    {
      #### error.rate and error.rate.per.class over nrepeat
      error.rate.all = error.mat
      error.rate = list()
      error.rate.sd = list()
      
      error.rate.per.class.all = error.mat.class
      error.rate.per.class = relist(0,skeleton = error.mat.class[[1]])
      error.rate.per.class.sd = relist(0,skeleton = error.mat.class[[1]])
      
      
      temp.error.rate = array(0, c(dim(error.rate.all[[1]][[1]]), nrepeat))
      temp.error.rate.per.class = array(0, c(dim(error.rate.per.class.all[[1]][[1]][[1]]), nrepeat))
      for(i in 1 : J)
      {
        for(nrep in 1:nrepeat)
          temp.error.rate[, , nrep] = error.rate.all[[nrep]][[i]]
        
        
        temp.error.rate.mean = apply(temp.error.rate, c(1,2), mean)
        temp.error.rate.sd = apply(temp.error.rate, c(1,2), sd)
        dimnames(temp.error.rate.mean) =  dimnames(temp.error.rate.sd) = dimnames(error.rate.all[[nrep]][[i]])
        
        error.rate[[i]] = temp.error.rate.mean
        error.rate.sd[[i]] = temp.error.rate.sd
        
        for(j in 1 : length(dist.select))
        {
          for(nrep in 1:nrepeat)
            temp.error.rate.per.class[, , nrep] = error.rate.per.class.all[[nrep]][[i]][[j]]
          
          temp.error.rate.per.class.mean = apply(temp.error.rate.per.class, c(1,2), mean)
          temp.error.rate.per.class.sd = apply(temp.error.rate.per.class, c(1,2), sd)
          dimnames(temp.error.rate.per.class.mean) =  dimnames(temp.error.rate.per.class.sd) = dimnames(error.rate.per.class.all[[nrep]][[i]][[j]])
          
          error.rate.per.class[[i]][[j]] = temp.error.rate.per.class.mean
          error.rate.per.class.sd[[i]][[j]] = temp.error.rate.per.class.sd
        }
        
      }
      names(error.rate) = names(error.rate.sd) = names(X)
      
      #### AveragePredict over nrepeat
      AveragedPredict.error.rate.all = Y.mean.res
      temp.AveragedPredict.error.rate = array(unlist(AveragedPredict.error.rate.all), c(dim(AveragedPredict.error.rate.all[[1]]),nrepeat))
      AveragedPredict.error.rate = apply(temp.AveragedPredict.error.rate, c(1,2), mean)
      AveragedPredict.error.rate.sd = apply(temp.AveragedPredict.error.rate, c(1,2), sd)
      dimnames(AveragedPredict.error.rate) = dimnames(AveragedPredict.error.rate.sd) = dimnames(AveragedPredict.error.rate.all[[1]])
      
      #### WeightedPredict over nrepeat
      WeightedPredict.error.rate.all = Y.WeightedPredict.res
      temp.WeightedPredict.error.rate = array(unlist(WeightedPredict.error.rate.all), c(dim(WeightedPredict.error.rate.all[[1]]),nrepeat))
      WeightedPredict.error.rate = apply(temp.WeightedPredict.error.rate, c(1,2), mean)
      WeightedPredict.error.rate.sd = apply(temp.WeightedPredict.error.rate, c(1,2), sd)
      dimnames(WeightedPredict.error.rate) = dimnames(WeightedPredict.error.rate.sd) = dimnames(WeightedPredict.error.rate.all[[1]])
      
      
      #### MajorityVote.error.rate and WeightedVote.error.rate over nrepeat
      MajorityVote.error.rate.all = Y.vote.res
      MajorityVote.error.rate = relist(0,skeleton = MajorityVote.error.rate.all[[1]])
      MajorityVote.error.rate.sd = relist(0,skeleton = MajorityVote.error.rate.all[[1]])
      
      WeightedVote.error.rate.all = Y.weighted.vote.res
      WeightedVote.error.rate = relist(0,skeleton = WeightedVote.error.rate.all[[1]])
      WeightedVote.error.rate.sd = relist(0,skeleton = WeightedVote.error.rate.all[[1]])
      
      temp.MajorityVote.error.rate = array(0, c(dim(MajorityVote.error.rate.all[[1]][[1]]), nrepeat))
      temp.WeightedVote.error.rate = array(0, c(dim(WeightedVote.error.rate.all[[1]][[1]]), nrepeat))
      for(j in 1 : length(dist.select))
      {
        for(nrep in 1:nrepeat)
        {
          temp.MajorityVote.error.rate[, , nrep] = MajorityVote.error.rate.all[[nrep]][[j]]
          temp.WeightedVote.error.rate[, , nrep] = WeightedVote.error.rate.all[[nrep]][[j]]
        }
        
        # MajorityVote.error.rate
        temp.MajorityVote.error.rate.mean = apply(temp.MajorityVote.error.rate, c(1,2), mean)
        temp.MajorityVote.error.rate.sd = apply(temp.MajorityVote.error.rate, c(1,2), sd)
        dimnames(temp.MajorityVote.error.rate.mean) =  dimnames(temp.MajorityVote.error.rate.sd) = dimnames(MajorityVote.error.rate.all[[nrep]][[j]])
        
        MajorityVote.error.rate[[j]] = temp.MajorityVote.error.rate.mean
        MajorityVote.error.rate.sd[[j]] = temp.MajorityVote.error.rate.sd
        
        # WeightedVote.error.rate
        temp.WeightedVote.error.rate.mean = apply(temp.WeightedVote.error.rate, c(1,2), mean)
        temp.WeightedVote.error.rate.sd = apply(temp.WeightedVote.error.rate, c(1,2), sd)
        dimnames(temp.WeightedVote.error.rate.mean) =  dimnames(temp.WeightedVote.error.rate.sd) = dimnames(WeightedVote.error.rate.all[[nrep]][[j]])
        
        WeightedVote.error.rate[[j]] = temp.WeightedVote.error.rate.mean
        WeightedVote.error.rate.sd[[j]] = temp.WeightedVote.error.rate.sd
        
      }
      
      
    } else {
      error.rate = error.mat[[1]]
      error.rate.per.class = error.mat.class[[1]]
      AveragedPredict.error.rate = Y.mean.res[[1]]
      WeightedPredict.error.rate = Y.WeightedPredict.res[[1]]
      MajorityVote.error.rate = Y.vote.res[[1]]
      WeightedVote.error.rate = Y.weighted.vote.res[[1]]
    }
    
    
    
    
    result = list()
    result$error.rate = error.rate
    if(nrepeat>1)
    {
      result$error.rate.sd = error.rate.sd
      result$error.rate.all = error.rate.all
    }
    
    result$error.rate.per.class = error.rate.per.class
    if(nrepeat>1)
    {
      result$error.rate.per.class.sd = error.rate.per.class.sd
      result$error.rate.per.class.all = error.rate.per.class.all
    }
    
    result$predict = Y.all
    result$class = Y.predict
    
    result$features$stable = list.features
    #result$features$final = final.features
    
    if (length(X) > 1)
    {
      result$AveragedPredict.class = Y.mean
      result$AveragedPredict.error.rate = AveragedPredict.error.rate
      if(nrepeat>1)
      {
        result$AveragedPredict.error.rate.sd = AveragedPredict.error.rate.sd
        result$AveragedPredict.error.rate.all = AveragedPredict.error.rate.all
      }
      
      result$WeightedPredict.class = Y.WeightedPredict
      result$WeightedPredict.error.rate = WeightedPredict.error.rate
      if(nrepeat>1)
      {
        result$WeightedPredict.error.rate.sd = WeightedPredict.error.rate.sd
        result$WeightedPredict.error.rate.all = WeightedPredict.error.rate.all
      }
      
      result$MajorityVote = Y.vote
      result$MajorityVote.error.rate = MajorityVote.error.rate
      if(nrepeat>1)
      {
        result$MajorityVote.error.rate.sd = MajorityVote.error.rate.sd
        result$MajorityVote.error.rate.all = MajorityVote.error.rate.all
      }
      
      result$WeightedVote = Y.weighted.vote
      result$WeightedVote.error.rate = WeightedVote.error.rate
      if(nrepeat>1)
      {
        result$WeightedVote.error.rate.sd = WeightedVote.error.rate.sd
        result$WeightedVote.error.rate.all = WeightedVote.error.rate.all
      }
      weights <- Reduce(rbind, weights)
      weights <- weights[order(weights$block),]
      result$weights <- weights
      
    }
    
    
    # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
    # one ncomp_opt for each dist, each overall/BER and each prediction framework (AveragePredict, WeightedPredict, MajorityVote, WeightedVote
    
    measure = c("Overall.ER","Overall.BER") # one of c("overall","BER")
    
    ncomp_opt = vector("list", length = 4)
    names(ncomp_opt) = c("AveragedPredict", "WeightedPredict", "MajorityVote", "WeightedVote")
    
    if(nrepeat > 2 & min(object$ncomp) >1)
    {
      for(prediction_framework in names(ncomp_opt))
      {
        if(prediction_framework %in% c("AveragedPredict", "WeightedPredict"))
        {
          ncomp_opt[[prediction_framework]] = matrix(NA, nrow = length(measure), ncol = 1,
                                                     dimnames = list(measure))
          
          for (measure_i in measure)
          {
            mat.error.rate = sapply(get(paste0(prediction_framework, ".error.rate.all")), function(x){x[measure_i,]})
            ncomp_opt[[prediction_framework]][measure_i,] = t.test.process(t(mat.error.rate), alpha = signif.threshold)
          }
          
        } else {
          ncomp_opt[[prediction_framework]] = matrix(NA, nrow = length(measure), ncol = length(dist.select),
                                                     dimnames = list(measure, dist.select))
          
          for (measure_i in measure)
          {
            for (ijk in dist.select)
            {
              mat.error.rate = sapply(get(paste0(prediction_framework, ".error.rate.all")), function(x){x[[ijk]][measure_i,]})
              ncomp_opt[[prediction_framework]][measure_i, ijk] = t.test.process(t(mat.error.rate), alpha = signif.threshold)
            }
          }
        }
      }
    }
    
    result$meth = "sgccda.mthd"
    class(result) = "perf.sgccda.mthd"
    result$call = match.call()
    result$crit = crit
    result$choice.ncomp = ncomp_opt
    
    if (auc) {
      auc.comp <- lapply(.name_list(names(auc.rep.comp[[1]])), function(comp) {
        Reduce("+", lapply(seq_len(nrepeat), function(rep){
          auc.rep.comp[[rep]][[comp]]
        })
        )/nrepeat
      })
      
      ## mean AUC for all blocks
      result$auc <- auc.comp
    }
    
    
    return(invisible(result))
  }

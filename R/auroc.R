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

#' Area Under the Curve (AUC) and Receiver Operating Characteristic (ROC)
#' curves for supervised classification
#' 
#' Calculates the AUC and plots ROC for supervised models from s/plsda,
#' mint.s/plsda and block.plsda, block.splsda or wrapper.sgccda functions.
#' 
#' For more than two classes in the categorical outcome Y, the AUC is
#' calculated as one class vs. the other and the ROC curves one class vs. the
#' others are output.
#' 
#' The ROC and AUC are calculated based on the predicted scores obtained from
#' the \code{predict} function applied to the multivariate methods
#' (\code{predict(object)$predict}). Our multivariate supervised methods
#' already use a prediction threshold based on distances (see \code{predict})
#' that optimally determine class membership of the samples tested. As such AUC
#' and ROC are not needed to estimate the performance of the model (see
#' \code{perf}, \code{tune} that report classification error rates). We provide
#' those outputs as complementary performance measures.
#' 
#' The pvalue is from a Wilcoxon test between the predicted scores between one
#' class vs the others.
#' 
#' External independent data set (\code{newdata}) and outcome
#' (\code{outcome.test}) can be input to calculate AUROC. The external data set
#' must have the same variables as the training data set (\code{object$X}).
#' 
#' If \code{object} is a named list of multiple \code{plsda} and \code{splsda}
#' objects, ensure that these models each have a response variable with the same
#' levels. Additionally, \code{newdata} and \code{outcome.test} cannot be passed
#' to this form of \code{auroc}. 
#' 
#' If \code{newdata} is not provided, AUROC is calculated from the training
#' data set, and may result in overfitting (too optimistic results).
#' 
#' Note that for mint.plsda and mint.splsda objects, if \code{roc.study} is
#' different from "global", then \code{newdata}), \code{outcome.test} and
#' \code{sstudy.test} are not used.
#' 
#' @aliases auroc auroc.mixo_plsda auroc.mixo_splsda auroc.list auroc.mint.plsda
#' auroc.mint.splsda auroc.sgccda
#' 
#' @param object Object of class inherited from one of the following supervised
#' analysis function: "plsda", "splsda", "mint.plsda", "mint.splsda",
#' "block.splsda" or "wrapper.sgccda". Alternatively, this can be a named list
#' of plsda and splsda objects if multiple models are to be compared. Note that 
#' these multiple models need to have used the same levels in the response variable.
#' @param newdata numeric matrix of predictors, by default set to the training
#' data set (see details).
#' @param outcome.test Either a factor or a class vector for the discrete
#' outcome, by default set to the outcome vector from the training set (see
#' details).
#' @param study.test For MINT objects, grouping factor indicating which samples
#' of \code{newdata} are from the same study. Overlap with \code{object$study}
#' are allowed.
#' @param multilevel Sample information when a newdata matrix is input and when
#' multilevel decomposition for repeated measurements is required. A numeric
#' matrix or data frame indicating the repeated measures on each individual,
#' i.e. the individuals ID. See examples in \code{splsda}.
#' @param plot Whether the ROC curves should be plotted, by default set to TRUE
#' (see details).
#' @param roc.comp Specify the component (integer) up to which the ROC
#'   will be calculated and plotted from the multivariate model, default to 1.
#' @param roc.block Specify the block number (integer) or the name of the block
#' (set of characters) for which the ROC will be plotted for a block.plsda or
#' block.splsda object, default to 1.
#' @param roc.study Specify the study for which the ROC will be plotted for a
#' mint.plsda or mint.splsda object, default to "global".
#' @param title Character, specifies the title of the plot.
#' @param print Logical, specifies whether the output should be printed.
#' @param ... external optional arguments for plotting - \code{line.col} for
#' custom colors and \code{legend.title} for custom legend title
#' @return Depending on the type of object used, a list that contains: The AUC
#' and Wilcoxon test pvalue for each 'one vs other' classes comparison
#' performed, either per component (splsda, plsda, mint.plsda, mint.splsda), or
#' per block and per component (wrapper.sgccda, block.plsda, blocksplsda).
#' @author Benoit Gautier, Francois Bartolo, Florian Rohart, Al J Abadi
#' @seealso \code{\link{tune}}, \code{\link{perf}}, and https://mixomics.org/
#' for more details.
#' @keywords regression multivariate
#' @example ./examples/auroc-examples.R
#' @export
auroc <- function(object, ...)
    UseMethod("auroc")


# PLSDA object
# ----------------------
#' @rdname auroc
#' @method auroc mixo_plsda
#' @export
auroc.mixo_plsda <- 
    function(
        object,
        newdata = object$input.X,
        outcome.test = as.factor(object$Y),
        multilevel = NULL,
        plot = TRUE,
        roc.comp = NULL,
        title = NULL,
        print=TRUE,
        ...)
    {
        if(dim(newdata)[[1]] != length(outcome.test))
            stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]],
                 " elements.",call. = FALSE)
        
        
        if (is.null(roc.comp))
        {
            roc.comp <- object$ncomp
        }
        
        if( length(roc.comp) != 1)
            stop("`roc.comp' must be a single integer")
        
        data = list()
        statauc.res = graph = list()
        data$outcome=factor(outcome.test)
        
        # note here: the dist does not matter as we used the predicted scores only
        res.predict = predict.mixo_spls(object, newdata = newdata,
                                        dist = "max.dist", multilevel = multilevel)$predict
        
        for (i in seq_len(object$ncomp))
        {
            title <- paste0("ROC Curve Using Comp(s): ",paste0(seq_len(i), collapse = ', '))
            data$data=res.predict[,,i]
            temp = statauc(data, plot = ifelse(i%in%roc.comp,plot,FALSE), title=title,...)
            statauc.res[[paste0("Comp", i, sep = "")]] = temp[[1]]
            graph[[paste0("Comp", i, sep = "")]] = temp$graph
        }
        if (isTRUE(print))
            print(statauc.res)
        
        return(invisible(c(statauc.res,graph=graph)))
    }

#' @rdname auroc
#' @export
auroc.mixo_splsda <- auroc.mixo_plsda


# Multiple (s)PLSDA objects
# ----------------------
#' @rdname auroc
#' @method auroc list
#' @export
auroc.list <- 
  function(
    object,
    plot = TRUE,
    roc.comp = NULL,
    title = NULL,
    print = TRUE,
    ...) 
  {
    
    # set baseline ncomp and response levels to check all objects against
    # these need to be constant so any deviation from these results in an error
    base.levels <- levels(object[[1]]$Y)
    base.ncomp <- object[[1]]$ncomp
    
    # for the sake of visual clutter as well as distinguishing them via linetypes
    # via ggplot, a maximum of 6 models can be handled
    if (length(object) > 6) {
      stop("Can take a maximum of SIX (s)plsda objects")
    }
    
    # apply checks on each model
    for (obj in object) {
      
      # check it is a plsda or splsda object
      if (!(any(class(obj) %in% c("mixo_plsda", "mixo_splsda")))) {
        stop("Combined auroc can only take 'plsda' and 'splsda' objects",
             call. = FALSE)
      }
      # check that the levels of the response variable is consistent
      if (length(setdiff(base.levels, levels(obj$Y))) != 0) {
        stop("Combined auroc must have models which utilise the same response variable",
             call. = FALSE)
      }
      # check the ncomp is consistent
      if (base.ncomp != obj$ncomp) {
        stop("Combined auroc must have models which have the same ncomp",
             call. = FALSE)
      }
    }
    # handle default roc..comp and ensure it is a single value
    if (is.null(roc.comp)) { roc.comp <- base.ncomp }
    if (length(roc.comp) != 1) { stop("`roc.comp' must be a single integer") }
    
    # initialise returned objects
    auc.list <- list()
    df <- data.frame(matrix(NA, nrow=0, ncol=4))
    Specificity <- Sensitivity <- Outcome <-  model <- NULL
    
    for (idx in seq_len(length(object))) { # for each model
      
      obj <- object[[idx]]
      
      data <- list()
      statauc.res <- list()
      
      # use all training data as testing data - note this likely results in overestimated efficacy
      newdata <- obj$input.X 
      #browser()
      data$outcome <- as.factor(obj$Y)
      
      # generate predictions
      res.predict = predict.mixo_spls(obj, newdata = newdata,
                                      dist = "max.dist")$predict
      data$data <- res.predict[,,roc.comp] # extract predictions on specified component 
      temp = statauc(data) # generate AUROC data
      auc.list[[names(object)[idx]]] <- temp[[1]] # extract auc values
      temp$df[, "model"] <- rep(names(object)[idx], nrow(temp$df)) # add model column to df
      df <- rbind(df, temp$df) # add all plot vertices to df
    }
    # output of statauc has the AUC values included in the Outcome column. This
    # results in them being unique across different models
    # this line homogenises the Outcome levels so a proper legend can be used
    df$Outcome <- substr(df$Outcome, start = 1, stop=regexpr(":", df$Outcome)-1)
    
    # default scenario for title
    if(is.null(title)) {title = paste0("ROC Curve Using Comp(s): ", paste0(seq_len(roc.comp), collapse = ', ')) }
    
    # set "best" order of line types. if only 2-3 models, these will be easiest to distinguish
    linetypes <- c("solid", "dotted", "dotdash", "twodash", "dashed", "longdash")
    
    # generate plot
    p = ggplot(df, aes(x=Specificity, y=Sensitivity)) + 
      geom_line(aes(linetype=model, color=Outcome), linewidth = 1.1) +
      xlab("100 - Specificity (%)") + 
      ylab("Sensitivity (%)") + 
      scale_x_continuous(breaks=seq(0, 100, by = 10)) + 
      scale_y_continuous(breaks=seq(0, 100, by = 10))+
      scale_linetype_manual(values=linetypes[1:length(object)])
    p = p + 
      geom_abline(intercept = 1) + 
      theme(legend.key.size = 
              unit(1.5, "cm"), 
            plot.title = element_text(lineheight=.8, face="bold"), 
            legend.title = element_text(size=14, face="bold")) + 
      ggtitle(title) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    if (plot) { plot(p) } # plot if desired
    if (print) { print(auc.list) } # print auc stats if desired
    
    return(invisible(list(auc=auc.list,
                          graph=p)))
    
  }

# MINT object
# ----------------------
#' @rdname auroc
#' @method auroc mint.plsda
#' @export
auroc.mint.plsda <- 
    function(
        object,
        newdata = object$X,
        outcome.test = as.factor(object$Y),
        study.test = object$study,
        multilevel = NULL,
        plot = TRUE,
        roc.comp = NULL,
        roc.study = "global",
        title=NULL,
        print=TRUE,
        ...)
    {
        if(length(roc.study) != 1)
            stop("`roc.study' must be a single entry,
    either `global' or one of levels(object$study)")
        
        if (is.null(roc.comp))
        {
            roc.comp <- object$ncomp
        }
        
        if( length(roc.comp) != 1)
            stop("`roc.comp' must be a single integer")
        
        if(roc.study == "global"){
            if(dim(newdata)[[1]] != length(outcome.test))
                stop("Factor outcome.test must be a factor with ",dim(newdata)[[1]],
                     " elements.",call. = FALSE)
            
            if(dim(newdata)[[1]]!=length(study.test))
                stop("Factor study.test must be a factor with ",dim(newdata)[[1]],
                     " elements.",call. = FALSE)
            study.test=factor(study.test)
            title.temp = NULL
            
        } else {
            
            # check study
            if (!roc.study%in%c(levels(object$study)))
                stop("'roc.study' must be one of 'levels(object$study)'")
            
            ind.study = object$study == roc.study
            newdata = object$X[ind.study, ]
            outcome.test = as.factor(object$Y[ind.study])
            study.test = factor(object$study[ind.study])
            title.temp = paste0(", Study ", roc.study)
            
        }
        
        data=list()
        statauc.res = graph = list()
        data$outcome=factor(outcome.test)
        
        # note here: the dist does not matter as we used the predicted scores only
        res.predict = predict.mixo_spls(object, newdata = newdata, dist = "max.dist",
                                        multilevel = multilevel, study.test = study.test)$predict
        
        data$data=res.predict[,,roc.comp]
        
        if (is.null(title)) {
            title=paste0("ROC Curve Using Comp(s): ", paste0(seq_len(roc.comp), collapse = ', '), title.temp)
        }
        
        temp = statauc(data, plot = ifelse(roc.comp%in% roc.comp,plot,FALSE), title=title,...)
        statauc.res[[paste0("Comp", roc.comp, sep = "")]] = temp[[1]]
        graph[[paste0("Comp", roc.comp, sep = "")]] = temp$graph
        if (isTRUE(print))
            print(statauc.res)
        
        return(invisible(c(statauc.res,graph=graph)))
        
    }

#' @rdname auroc
#' @method auroc mint.splsda
#' @export
auroc.mint.splsda <- auroc.mint.plsda

# block.splsda object
# ----------------------
#' @rdname auroc
#' @importFrom methods is
#' @method auroc sgccda
#' @export
auroc.sgccda <- function(
    object,
    newdata = object$X,
    outcome.test = as.factor(object$Y),
    multilevel = NULL,
    plot = TRUE,
    roc.block = 1L,
    roc.comp = NULL,
    title=NULL,
    print=TRUE,
    ...)
{
    
    data=list()
    auc.mean = graph=list()
    data$outcome=factor(outcome.test)
    
    if (is.null(roc.comp))
    {
        roc.comp <- object$ncomp
    }
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict.block.spls(object, newdata = newdata,
                                       dist = "max.dist", multilevel = multilevel)$predict
    block.all = names(res.predict)
    if (is(roc.block, "numeric")) {
        roc.block <- as.integer(roc.block)
        lb <- length(names(res.predict))
        if (roc.block > lb)
            stop(sprintf("roc.block cannot be greater than %s", lb ))
        block.temp = names(res.predict[roc.block])
    } else if (is(roc.block, "character")) {
        block.temp = roc.block 
    } else {
        stop("'roc.block' should be an integer or character")
    }
    title.temp = title
    for(j in seq_len(length(res.predict)))
    {
        for (i in seq_len(object$ncomp[j]))
        {
            data$data=res.predict[[j]][,,i]
            if (is.null(title.temp)) {
                title=paste0("ROC Curve\nBlock: ", names(res.predict)[j],
                            ", Using Comp(s): ",paste0(seq_len(i), collapse = ', '))
            }
            
            plot.temp =
                ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp,
                       plot, FALSE)
            temp = statauc(data, plot = plot.temp, title = title, ...)
            auc.mean[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp[[1]]
            graph[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp$graph
            
        }
        out = c(auc.mean,graph=graph)
    }
    if (isTRUE(print))
        print(auc.mean)
    
    return(invisible(out))
}

# mint.block.splsda object
# ----------------------
#' @rdname auroc
#' @method auroc mint.block.plsda
#' @export
auroc.mint.block.plsda <- function(
    object,
    newdata = object$X,
    study.test = object$study,
    outcome.test = as.factor(object$Y),
    multilevel = NULL,
    plot = TRUE,
    roc.block = 1,
    roc.comp = NULL,
    title=NULL,
    print=TRUE,
    ...)
{
    
    data=list()
    auc.mean = graph=list()
    data$outcome=factor(outcome.test)
    study.test=factor(study.test)
    
    if (is.null(roc.comp))
    {
        roc.comp <- object$ncomp
    }
    
    # note here: the dist does not matter as we used the predicted scores only
    res.predict  =  predict.mixo_spls(object, newdata = newdata,
                                      study.test=study.test,dist = "max.dist", multilevel = multilevel)$predict
    block.temp = names(res.predict[roc.block])
    
    for(j in seq_len(length(res.predict)))
    {
        for (i in seq_len(object$ncomp[j]))
        {
            data$data=res.predict[[j]][,,i]
            if (is.null(title)) {
                title=paste("ROC Curve\nBlock: ", names(res.predict)[j],
                            ", comp: ",i, sep="")
            }
            
            plot.temp =
                ifelse(i%in%roc.comp && names(res.predict)[j]%in%block.temp,
                       plot, FALSE)
            temp = statauc(data, plot = plot.temp, title = title, ...)
            auc.mean[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp[[1]]
            graph[[names(res.predict)[j]]][[paste0("comp",i,sep = "")]] =
                temp$graph
            
        }
        out = c(auc.mean,graph=graph)
    }
    if (isTRUE(print))
        print(auc.mean)
    
    return(invisible(out))
}

#' @rdname auroc
#' @method auroc mint.block.splsda
#' @export
auroc.mint.block.splsda <- auroc.mint.block.plsda

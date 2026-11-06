#############################################################################################################
# Authors:
#   Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France#   Kim-Anh Le Cao, University of Queensland Diamantina Institute, Brisbane, Australia
#   Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#
# created: 2016
# last modified: 2016
#
# Copyright (C) 2016
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
#############################################################################################################


statauc <- function(data = NULL, plot = FALSE, title = NULL, line.col = NULL, legend.title = NULL){
  res.predict = data$data; outcome = data$outcome
  
  
  ann_text = matrix(ncol=2,nrow=nlevels(outcome))
  colnames(ann_text) = c("AUC", "p-value")
  
  df = NULL; seauc = NULL; zauc = NULL
  for (i in 1 : nlevels(outcome)){
    tempout = outcome
    levels(tempout)[-i] = "Other(s)"
    tempout = factor(tempout, c(levels(outcome)[i], "Other(s)"))
    temp = roc.default(response = tempout, predictor = res.predict[, i,drop=FALSE])
    
    # print(coords(temp,x="local maximas" ,ret="threshold"))
    # print(coords(temp,x="best" ,ret="threshold",best.method = "youden"))
    # print(coords(temp,x="best" ,ret="threshold",best.method = "closest.topleft"))
    
    seauc = sqrt((0.25 + (sum(table(tempout)) -2) * 0.083333)/(prod(table(tempout))))
    zauc = (temp$auc-0.5)/seauc
    zauc[zauc < 0] = - zauc[zauc < 0]
    #for (k in unique(temp$specificities)){
    #    ind = which(temp$specificities == k)
    #  temp$sensitivities[ind] = temp$sensitivities[rev(ind)]
    #}
    
    #reorder sensitivities when several  identical specificities
    a=split(temp$sensitivities,
            factor(temp$specificities,levels=unique(temp$specificities)))
    b=lapply(a,function(x){rev(x)})
    temp$sensitivities = do.call("c",b)
    
    if (nlevels(outcome) == 2){
      ann_text = matrix(ncol=2,nrow=1)
      colnames(ann_text) = c("AUC", "p-value")
      df = rbind(df, cbind(temp$specificities, temp$sensitivities, paste(paste(levels(outcome)[1], levels(outcome)[2], sep = " vs "),":",signif(temp$auc, 4))))
      #ann_text[i , 1] =
      ann_text[i , 1] = signif(temp$auc, 4)
      ann_text[i , 2] = signif((1 - pnorm(zauc, 0, 1))*2, 4)
      break
    } else {
      df = rbind(df, cbind(temp$specificities, temp$sensitivities, paste(levels(outcome)[i], "vs Other(s):",signif(temp$auc, 4))))
      #ann_text[i , 1] =
      ann_text[i , 1] = as.numeric(signif(temp$auc, 4))
      ann_text[i , 2] = as.numeric(signif((1 - pnorm(zauc, 0, 1))*2, 4))
    }
  }
  #define rownames for ann_text
  if (nlevels(outcome) == 2){
    rownames(ann_text) = paste(levels(outcome)[1], levels(outcome)[2], sep = " vs ")
  }else{
    rownames(ann_text) = paste(levels(outcome), "vs Other(s)")
  }
  
  df = data.frame(df, stringsAsFactors = FALSE)
  names(df) = c("Specificity", "Sensitivity", "Outcome")
  df$Specificity = 100 - 100*as.numeric(df$Specificity)
  df$Sensitivity = 100*as.numeric(df$Sensitivity)
  
  if(plot)
  {
    Sensitivity = Specificity = Outcome = NULL #R check
    if(is.null(title))
      title = "ROC Curve"
    else
      title=title
    p = ggplot(df, aes(x=Specificity, y=Sensitivity, group = Outcome, colour = Outcome)) + xlab("100 - Specificity (%)") + ylab("Sensitivity (%)") + geom_line(linewidth = 1.5) + scale_x_continuous(breaks=seq(0, 100, by = 10)) + scale_y_continuous(breaks=seq(0, 100, by = 10))
    p = p + geom_abline(intercept = 1) + theme(legend.key.size = unit(1.5, "cm"), plot.title = element_text(lineheight=.8, face="bold"), legend.title = element_text(size=14, face="bold")) + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
    
    if (!is.null(line.col)) {
      lc <- length(unique(df$Outcome))
      if ( !identical(lc, length(try(line.col))) ) stop(sprintf("line.col should be character string of length %s", lc))
      p <- p + scale_color_manual(values = line.col)
    }
    
    if (!is.null(legend.title)) {
      if ( !is(legend.title, "character" )) stop("lgened.title should be character")
      p <- p + guides(col = guide_legend(title = legend.title))
    }
    
    plot(p)
  } else {
    p=NULL
  }
  return(list(ann_text,graph=p,df=df))
  
}


###ROC/AUC

roc.default <- function(response, predictor,
                        
                        auc=TRUE,
                        
                        levels=base::levels(as.factor(response))
                        
                        
) {
  
  
  # Response / Predictor
  original.predictor <- predictor # store a copy of the original predictor (before converting ordered to numeric and removing NA)
  original.response <- response # store a copy of the original predictor (before converting ordered to numeric)
  
  # remove NAs if requested
  nas <- is.na(response) | is.na(predictor)
  if (any(nas)) {
    na.action <- grep(TRUE, nas)
    class(na.action) <- "omit"
    response <- response[!nas]
    attr(response, "na.action") <- na.action
    predictor <- predictor[!nas]
    attr(predictor, "na.action") <- na.action
  }
  
  splitted <- split(predictor, response)
  controls <- splitted[[as.character(levels[1])]]
  cases <- splitted[[as.character(levels[2])]]
  
  # Remove patients not in levels
  patients.in.levels <- response %in% levels
  if (!all(patients.in.levels)) {
    response <- response[patients.in.levels]
    predictor <- predictor[patients.in.levels]
  }
  
  # update 13/01/17: first level as control to force directionality:
  # > : if the predictor values for the control group are higher than the values of
  #the case group (controls > t >= cases)
  
  #if (median(controls) <= median(cases))
  #direction <- "<"
  #else if (median(controls) > median(cases))
  direction <- ">"
  
  
  # create the roc object
  roc <- list()
  class(roc) <- "roc"
  
  # compute SE / SP
  so = sort(unique(c(controls, cases)))
  thresholds <-((c(-Inf, so) + c(so, +Inf))/2)
  
  #thresholds <-((c(-Inf, sort(unique(c(controls, cases)))) + c(sort(unique(c(controls, cases))), +Inf))/2)
  #perf.matrix <- sapply(thresholds, roc.utils.perfs, controls=controls, cases=cases, direction=direction)
  perf.matrix <- roc.utils.perfs.fast.all.threshold(thresholds, controls=controls, cases=cases, direction=direction)
  perfs <- list(se=perf.matrix[2,], sp=perf.matrix[1,])
  
  se <- perfs$se
  sp <- perfs$sp
  
  
  # store the computations in the roc object
  roc$sensitivities <- se
  roc$specificities <- sp
  roc$thresholds <- thresholds
  roc <- sort(roc)
  roc$direction <- direction
  roc$cases <- cases
  roc$controls <- controls
  
  # compute AUC
  if (auc)
    roc$auc <- auc_roc(roc)
  
  roc$call <- match.call()
  roc$original.predictor <- original.predictor
  roc$original.response <- original.response
  roc$predictor <- predictor
  roc$response <- response
  roc$levels <- levels
  
  
  # return roc
  return(roc)
}

# returns a vector with two elements, sensitivity and specificity, given the threshold at which to evaluate the performance, the values of controls and cases and the direction of the comparison, a character '>' or '<' as controls CMP cases

# only the same as roc.utils.perfs if thresholds starts low and ends large
# (added 0 and 1 value)
roc.utils.perfs.fast.all.threshold <- function(thresholds, controls, cases, direction) {
  a=cut(cases,thresholds)
  b=cut(controls, thresholds)
  if (direction == '>') {
    se=c(0,cumsum(table(a))/length(cases))
    sp=c(1, (length(controls)- cumsum(table(b)))/length(controls))
  }  else if (direction == '<') {
    se=c(1, (length(cases)- cumsum(table(a)))/length(cases))
    sp=c(0, cumsum(table(b))/length(controls))
    
  }
  return(rbind(sp=sp, se=se))
}


roc.utils.perfs <- function(threshold, controls, cases, direction) {
  if (direction == '>') {
    tp <- sum(cases <= threshold)
    tn <- sum(controls > threshold)
  }
  else if (direction == '<') {
    tp <- sum(cases >= threshold)
    tn <- sum(controls < threshold)
  }
  # return(c(sp, se))
  return(c(sp=tn/length(controls), se=tp/length(cases)))
}

# sort roc curve. Make sure specificities are increasing.
sort.roc <- function(roc) {
  if (is.unsorted(roc$specificities)) {
    roc$sensitivities <- rev(roc$sensitivities)
    roc$specificities <- rev(roc$specificities)
    roc$thresholds <- rev(roc$thresholds)
  }
  return(roc)
}

#' @importFrom stats coefficients
auc_roc <- function(roc,
                    # Partial auc definition
                    partial.auc=FALSE, # false (consider total area) or numeric length 2: boundaries of the AUC to consider, between 0 and 1, or 0 and 100 if percent is TRUE
                    partial.auc.focus=c("specificity", "sensitivity"), # if partial.auc is not FALSE: do the boundaries
                    partial.auc.correct=FALSE,
                    allow.invalid.partial.auc.correct = FALSE,
                    ... # unused required to allow roc passing arguments to plot or ci.
) {
  if (!identical(partial.auc, FALSE)) {
    partial.auc.focus <- match.arg(partial.auc.focus)
  }
  
  percent <- FALSE
  
  # Validate partial.auc
  if (! identical(partial.auc, FALSE) & !(is.numeric(partial.auc) && length(partial.auc)==2))
    stop("partial.auc must be either FALSE or a numeric vector of length 2")
  
  # Ensure partial.auc is sorted with partial.auc[1] >= partial.auc[2]
  partial.auc <- sort(partial.auc, decreasing=TRUE)
  # Get and sort the sensitivities and specificities
  roc <- sort(roc)
  se <- roc$se
  sp <- roc$sp
  
  # Full area if partial.auc is FALSE
  if (identical(partial.auc, FALSE)) {
    if (is(roc, "smooth.roc") && ! is.null(roc$smoothing.args) && roc$smoothing.args$method == "binormal") {
      coefs <- coefficients(roc$model)
      auc <- unname(pnorm(coefs[1] / sqrt(1+coefs[2]^2)) * ifelse(percent, 100^2, 1))
    }
    else {
      diffs.x <- sp[-1] - sp[-length(sp)]
      means.vert <- (se[-1] + se[-length(se)])/2
      auc <- sum(means.vert * diffs.x)
    }
  }
  # Partial area
  else {
    if (partial.auc.focus == "sensitivity") {
      # if we focus on SE, just swap and invert x and y and the computations for SP will work
      x <- rev(se)
      y <- rev(sp)
    }
    else {
      x <- sp
      y <- se
    }
    
    # find the SEs and SPs in the interval
    x.inc <- x[x <= partial.auc[1] & x >= partial.auc[2]]
    y.inc <- y[x <= partial.auc[1] & x >= partial.auc[2]]
    # compute the AUC strictly in the interval
    diffs.x <- x.inc[-1] - x.inc[-length(x.inc)]
    means.vert <- (y.inc[-1] + y.inc[-length(y.inc)])/2
    auc <- sum(means.vert * diffs.x)
    # add the borders:
    if (length(x.inc) == 0) { # special case: the whole AUC is between 2 se/sp points. Need to interpolate from both
      diff.horiz <- partial.auc[1] - partial.auc[2]
      # determine indices
      idx.hi <- match(FALSE, x < partial.auc[1])
      idx.lo <- idx.hi - 1
      # proportions
      proportion.hi <- (x[idx.hi] - partial.auc[1]) / (x[idx.hi] - x[idx.lo])
      proportion.lo <- (partial.auc[2] - x[idx.lo]) / (x[idx.hi] - x[idx.lo])
      # interpolated y's
      y.hi <- y[idx.hi] + proportion.hi * (y[idx.lo] - y[idx.hi])
      y.lo <- y[idx.lo] - proportion.lo * (y[idx.lo] - y[idx.hi])
      # compute AUC
      mean.vert <- (y.hi + y.lo)/2
      auc <- mean.vert*diff.horiz
    }
    else { # if the upper limit is not exactly present in SPs, interpolate
      if (!(partial.auc[1] %in% x.inc)) {
        # find the limit indices
        idx.out <- match(FALSE, x < partial.auc[1])
        idx.in <- idx.out - 1
        # interpolate y
        proportion <- (partial.auc[1] - x[idx.out]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.out] + proportion * (y[idx.in] - y[idx.out])
        # add to AUC
        auc <- auc + (partial.auc[1] - x[idx.in]) * (y[idx.in] + y.interpolated)/2
      }
      if (!(partial.auc[2] %in% x.inc)) { # if the lower limit is not exactly present in SPs, interpolate
        # find the limit indices in and out
        #idx.out <- length(x) - match(TRUE, rev(x) < partial.auc[2]) + 1
        idx.out <- match(TRUE, x > partial.auc[2]) - 1
        idx.in <- idx.out + 1
        # interpolate y
        proportion <- (x[idx.in] - partial.auc[2]) / (x[idx.in] - x[idx.out])
        y.interpolated <- y[idx.in] + proportion * (y[idx.out] - y[idx.in])
        # add to AUC
        auc <- auc + (x[idx.in] - partial.auc[2]) * (y[idx.in] + y.interpolated)/2
      }
    }
  }
  
  # In percent, we have 100*100 = 10,000 as maximum area, so we need to divide by a factor 100
  if (percent)
    auc <- auc/100
  
  # Correction according to McClish DC, 1989
  if (all(!identical(partial.auc, FALSE), partial.auc.correct)) { # only for pAUC
    min <- roc.utils.min.partial.auc(partial.auc, percent)
    max <- roc.utils.max.partial.auc(partial.auc, percent)
    # The correction is defined only when auc >= min
    if (!allow.invalid.partial.auc.correct && auc < min) {
      warning("Partial AUC correction not defined for ROC curves below the diagonal.")
      auc <- NA
    }
    else if (percent) {
      auc <- (100+((auc-min)*100/(max-min)))/2 # McClish formula adapted for %
    }
    else {
      auc <- (1+((auc-min)/(max-min)))/2 # original formula by McClish
    }
  }
  # Prepare the AUC to return with attributes
  auc <- as.vector(auc) # remove potential pre-existing attributes
  attr(auc, "partial.auc") <- partial.auc
  attr(auc, "percent") <- percent
  attr(auc, "roc") <- roc
  if (!identical(partial.auc, FALSE)) {
    attr(auc, "partial.auc.focus") <- partial.auc.focus
    attr(auc, "partial.auc.correct") <- partial.auc.correct
  }
  class(auc) <- c("auc", class(auc))
  return(auc)
}

coords <- function(roc, x, input=c("threshold", "specificity", "sensitivity"), ret=c("threshold", "specificity", "sensitivity"), as.list=FALSE, drop=TRUE, best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5), ...) {
  # make sure x was provided
  roc$percent=FALSE
  if (missing(x) || length(x) == 0)
    stop("'x' must be a numeric or character vector of positive length.")
  # match input 
  input <- match.arg(input)
  # match return
  ret <- roc.utils.match.coords.ret.args(ret)
  # make sure the sort of roc is correct
  roc <- sort(roc)
  
  if (is.character(x)) {
    x <- match.arg(x, c("all", "local maximas", "best"))
    partial.auc <- attr(roc$auc, "partial.auc")
    if (x == "all") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        se <- roc$se
        sp <- roc$sp
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$se[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          sp <- roc$sp[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
        }
        else {
          se <- roc$se[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          sp <- roc$sp[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      co <- coords(roc, x=thres, input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (is(co, "matrix"))
        colnames(co) <- rep(x, dim(co)[2])
      else if (is(co, "list") && is(co[[1]], "list"))
        names(co) <- rep(x, length(co))
      return(co)
    }
    else if (x == "local maximas") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        se <- roc$se
        sp <- roc$sp
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$se[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          sp <- roc$sp[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
        }
        else {
          se <- roc$se[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          sp <- roc$sp[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      lm.idx <- roc.utils.max.thresholds.idx(thres, sp=sp, se=se)
      co <- coords(roc, x=thres[lm.idx], input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (is(co,  "matrix"))
        colnames(co) <- rep(x, dim(co)[2])
      else if (is(co, "list") && is(co[[1]], "list"))
        names(co) <- rep(x, length(co))
      return(co)
    }
    else { # x == "best"
      # What kind of "best" do we want?
      # Compute weights
      if (is.numeric(best.weights) && length(best.weights) == 2)
        r <- (1 - best.weights[2]) / (best.weights[1] * best.weights[2]) # r should be 1 by default
      else
        stop("'best.weights' must be a numeric vector of length 2")
      # Compute optimality criterion and store it in the optim.crit vector
      best.method <- match.arg(best.method[1], c("youden", "closest.topleft", "topleft")) # cheat: allow the user to pass "topleft"
      if (best.method == "youden") {
        optim.crit <- roc$sensitivities + r * roc$specificities
      }
      else if (best.method == "closest.topleft" || best.method == "topleft") {
        fac.1 <- ifelse(roc$percent, 100, 1)
        optim.crit <- - ((fac.1 - roc$sensitivities)^2 + r * (fac.1 - roc$specificities)^2)
      }
      
      # Filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        thres <- roc$thresholds[optim.crit==max(optim.crit)]
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          optim.crit <- (optim.crit)[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
        else {
          optim.crit <- (optim.crit)[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      co <- coords(roc, x=thres, input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (is(co, "matrix"))
        colnames(co) <- rep(x, dim(co)[2])
      else if (is(co, "list") && is(co[[1]], "list"))
        names(co) <- rep(x, length(co))
      return(co)
    }
  }
  else if (is.numeric(x)) {
    if (length(x) > 1) { # make this function a vector function
      if (as.list) {
        res <- lapply(x, function(x) coords.roc(roc, x, input, ret, as.list))
        names(res) <- x
      }
      else {
        res <- sapply(x, function(x) coords.roc(roc, x, input, ret, as.list))
        if (length(ret) == 1) {# sapply returns a vector instead of a matrix
          res <- t(res)
          rownames(res) <- ret
        }
        colnames(res) <- x
      }
      return(res)
    }
    if (input == "threshold") {
      res <- c(x, as.vector(roc.utils.perfs(x, roc$controls, roc$cases, roc$direction)) * ifelse(roc$percent, 100, 1))
    }
    if (input == "specificity") {
      if (x < 0 || x > ifelse(roc$percent, 100, 1))
        stop("Input specificity not within the ROC space.")
      if (x %in% roc$sp) {
        idx <- match(x, roc$sp)
        res <- c(roc$thresholds[idx], roc$sp[idx], roc$se[idx])
      }
      else { # need to interpolate
        idx.next <- match(TRUE, roc$sp > x)
        proportion <-  (x - roc$sp[idx.next - 1]) / (roc$sp[idx.next] - roc$sp[idx.next - 1])
        int.se <- roc$se[idx.next - 1] - proportion * (roc$se[idx.next - 1] - roc$se[idx.next])
        res <- c(NA, x, int.se)
      }
    }
    if (input == "sensitivity") {
      if (x < 0 || x > ifelse(roc$percent, 100, 1))
        stop("Input sensitivity not within the ROC space.")
      if (x %in% roc$se) {
        idx <- length(roc$se) + 1 - match(TRUE, rev(roc$se) == x)
        res <- c(roc$thresholds[idx], roc$sp[idx], roc$se[idx])
      }
      else { # need to interpolate
        idx.next <- match(TRUE, roc$se < x)
        proportion <- (x - roc$se[idx.next]) / (roc$se[idx.next - 1] - roc$se[idx.next])
        int.sp <- roc$sp[idx.next] + proportion * (roc$sp[idx.next - 1] - roc$sp[idx.next])
        res <- c(NA, int.sp, x)
      }
    }
    # Deduce additional tn, tp, fn, fp, npv, ppv
    ncases <- ifelse(is(roc, "smooth.roc"), length(attr(roc, "roc")$cases), length(roc$cases))
    ncontrols <- ifelse(is(roc, "smooth.roc"), length(attr(roc, "roc")$controls), length(roc$controls))
    se <- res[3]
    sp <- res[2]
    if (roc$percent) {
      tp <- se * ncases / 100
      fn <- ncases - tp
      tn <- sp * ncontrols / 100
      fp <- ncontrols - tn
      substr.percent <- 100
      npv <- 100 * tn / (tn + fn)
      ppv <- 100 * tp / (tp + fp)
    }
    else {
      tp <- se * ncases
      fn <- ncases - tp
      tn <- sp * ncontrols
      fp <- ncontrols - tn
      npv <- tn / (tn + fn)
      ppv <- tp / (tp + fp)
      substr.percent <- 1
    }
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    if (as.list) {
      list <- list(threshold=res[1], specificity=sp, sensitivity=se, accuracy=accuracy, tn=tn, tp=tp, fn=fn, fp=fp, npv=npv, ppv=ppv, "1-specificity"=substr.percent-sp, "1-sensitivity"=substr.percent-se, "1-accuracy"=substr.percent-accuracy, "1-npv"=substr.percent-npv, "1-ppv"=substr.percent-ppv)
      list <- list[ret]
      if (drop == FALSE) {
        list <- list(list)
        names(list) <- x
      }
      return(list)
    }
    else {
      res <- as.matrix(res)
      res <- rbind(res, accuracy, tn, tp, fn, fp, npv, ppv, substr.percent-sp, substr.percent-se, substr.percent-accuracy, substr.percent-npv, substr.percent-ppv)
      rownames(res) <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
      colnames(res) <- x
      return(res[ret,, drop=drop])
    }
  }
  else {
    stop("'x' must be a numeric or character vector.")
  }
}


# Arguments which can be returned by coords
roc.utils.match.coords.ret.args <- function(x) {
  
  valid.ret.args <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
  x <- replace(x, x=="t", "threshold")
  x <- replace(x, x=="npe", "1-npv")
  x <- replace(x, x=="ppe", "1-ppv")
  match.arg(x, valid.ret.args, several.ok=TRUE)
}

# Find all the local maximas of the ROC curve. Returns a logical vector
roc.utils.max.thresholds.idx <- function(thresholds, sp, se) {
  reversed <- FALSE
  if (is.unsorted(sp)) {
    # make sure SP is sorted increasingly, and sort thresholds accordingly
    thresholds <- rev(thresholds)
    sp <- rev(sp)
    se <- rev(se)
    reversed <- TRUE
  }
  # TODO: find whether the duplicate check is still needed.
  # Should have been fixed by passing only c(controls, cases)
  # instead of whole 'predictor' to roc.utils.thresholds in roc.default
  # but are there other potential issues like that?
  dup <- duplicated(data.frame(sp, se))
  thresholds <- thresholds[!dup]
  sp <- sp[!dup]
  se <- se[!dup]
  # Find the local maximas
  if (length(thresholds) == 1) {
    local.maximas <- TRUE # let's consider that if there is only 1 point, we should print it.
  }
  else if (length(thresholds) == 2) {
    local.maximas <- c(se[1] > se[2], sp[2] > sp[1])
  }
  else {
    local.maximas <- se[1] > se[2]
    for (i in 2:(length(thresholds)-1)) {
      if (sp[i] > sp[i-1] & se[i] > se[i+1])
        local.maximas <- c(local.maximas, TRUE)
      else if (sp[i] > sp[i-1] & se[i] == 0)
        local.maximas <- c(local.maximas, TRUE)
      else if (se[i] > se[i-1] & sp[i] == 1)
        local.maximas <- c(local.maximas, TRUE)
      else
        local.maximas <- c(local.maximas, FALSE)
    }
    local.maximas <- c(local.maximas, sp[length(thresholds)] > sp[length(thresholds)-1])
  }
  if (any(dup)) {
    lms <- rep(FALSE, length(dup))
    lms[!dup] <- local.maximas
    local.maximas <- lms
  }
  if (reversed)
    rev(local.maximas)
  
  # Remove +-Inf at the limits of the curve
  #local.maximas <- local.maximas & is.finite(thresholds)
  # Question: why did I do that? It breaks coords.roc when partial.auc contains only the extreme point
  
  return(local.maximas)
}


# from pROC v1.8.0
# Compute the min/max for partial AUC
# ... with an auc
roc.utils.min.partial.auc.auc <- function(auc) {
  roc.utils.min.partial.auc(attr(auc, "partial.auc"), attr(auc, "percent"))
}

roc.utils.max.partial.auc.auc <- function(auc) {
  roc.utils.max.partial.auc(attr(auc, "partial.auc"), attr(auc, "percent"))
}

# ... with partial.auc/percent
roc.utils.min.partial.auc <- function(partial.auc, percent) {
  if (!identical(partial.auc, FALSE)) {
    min <- sum(ifelse(percent, 100, 1)-partial.auc)*abs(diff(partial.auc))/2/ifelse(percent, 100, 1)
  }
  else {
    min <- 0.5 * ifelse(percent, 100, 1)
  }
  return(min)
}

roc.utils.max.partial.auc <- function(partial.auc, percent) {
  if (!identical(partial.auc, FALSE)) {
    max <- abs(diff(partial.auc))
  }
  else {
    max <- 1 * ifelse(percent, 100, 1)
  }
  return(max)
}

# from pROC v1.8.0
coords.roc <- function(roc, x, input=c("threshold", "specificity", "sensitivity"), ret=c("threshold", "specificity", "sensitivity"), as.list=FALSE, drop=TRUE, best.method=c("youden", "closest.topleft"), best.weights=c(1, 0.5), ...) {
  # make sure x was provided
  if (missing(x) || length(x) == 0)
    stop("'x' must be a numeric or character vector of positive length.")
  # match input
  input <- match.arg(input)
  # match return
  ret <- roc.utils.match.coords.ret.args(ret)
  # make sure the sort of roc is correct
  roc <- sort(roc)
  
  if (is.character(x)) {
    x <- match.arg(x, c("all", "local maximas", "best"))
    partial.auc <- attr(roc$auc, "partial.auc")
    if (x == "all") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        se <- roc$se
        sp <- roc$sp
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$se[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          sp <- roc$sp[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
        }
        else {
          se <- roc$se[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          sp <- roc$sp[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      co <- coords(roc, x=thres, input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (is(co, "matrix"))
        colnames(co) <- rep(x, dim(co)[2])
      else if (is(co, "list") && is(co[[1]], "list"))
        names(co) <- rep(x, length(co))
      return(co)
    }
    else if (x == "local maximas") {
      # Pre-filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        se <- roc$se
        sp <- roc$sp
        thres <- roc$thresholds
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          se <- roc$se[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          sp <- roc$sp[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
        }
        else {
          se <- roc$se[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          sp <- roc$sp[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      lm.idx <- roc.utils.max.thresholds.idx(thres, sp=sp, se=se)
      co <- coords(roc, x=thres[lm.idx], input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (is(co, "matrix"))
        colnames(co) <- rep(x, dim(co)[2])
      else if (is(co, "list") && is(co[[1]], "list"))
        names(co) <- rep(x, length(co))
      return(co)
    }
    else { # x == "best"
      # What kind of "best" do we want?
      # Compute weights
      if (is.numeric(best.weights) && length(best.weights) == 2)
        r <- (1 - best.weights[2]) / (best.weights[1] * best.weights[2]) # r should be 1 by default
      else
        stop("'best.weights' must be a numeric vector of length 2")
      # Compute optimality criterion and store it in the optim.crit vector
      best.method <- match.arg(best.method[1], c("youden", "closest.topleft", "topleft")) # cheat: allow the user to pass "topleft"
      if (best.method == "youden") {
        optim.crit <- roc$sensitivities + r * roc$specificities
      }
      else if (best.method == "closest.topleft" || best.method == "topleft") {
        fac.1 <- ifelse(roc$percent, 100, 1)
        optim.crit <- - ((fac.1 - roc$sensitivities)^2 + r * (fac.1 - roc$specificities)^2)
      }
      
      # Filter thresholds based on partial.auc
      if (is.null(roc$auc) || identical(partial.auc, FALSE)) {
        thres <- roc$thresholds[optim.crit==max(optim.crit)]
      }
      else {
        if (attr(roc$auc, "partial.auc.focus") == "sensitivity") {
          optim.crit <- (optim.crit)[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]]
          thres <- roc$thresholds[roc$se <= partial.auc[1] & roc$se >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
        else {
          optim.crit <- (optim.crit)[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]]
          thres <- roc$thresholds[roc$sp <= partial.auc[1] & roc$sp >= partial.auc[2]][optim.crit==max(optim.crit)]
        }
      }
      if (length(thres) == 0)
        return(NULL)
      co <- coords(roc, x=thres, input="threshold", ret=ret, as.list=as.list, drop=drop)
      if (is(co, "matrix"))
        colnames(co) <- rep(x, dim(co)[2])
      else if (is(co, "list") && is(co[[1]], "list"))
        names(co) <- rep(x, length(co))
      return(co)
    }
  }
  else if (is.numeric(x)) {
    if (length(x) > 1) { # make this function a vector function
      if (as.list) {
        res <- lapply(x, function(x) coords.roc(roc, x, input, ret, as.list))
        names(res) <- x
      }
      else {
        res <- sapply(x, function(x) coords.roc(roc, x, input, ret, as.list))
        if (length(ret) == 1) {# sapply returns a vector instead of a matrix
          res <- t(res)
          rownames(res) <- ret
        }
        colnames(res) <- x
      }
      return(res)
    }
    if (input == "threshold") {
      res <- c(x, as.vector(roc.utils.perfs(x, roc$controls, roc$cases, roc$direction)) * ifelse(roc$percent, 100, 1))
    }
    if (input == "specificity") {
      if (x < 0 || x > ifelse(roc$percent, 100, 1))
        stop("Input specificity not within the ROC space.")
      if (x %in% roc$sp) {
        idx <- match(x, roc$sp)
        res <- c(roc$thresholds[idx], roc$sp[idx], roc$se[idx])
      }
      else { # need to interpolate
        idx.next <- match(TRUE, roc$sp > x)
        proportion <-  (x - roc$sp[idx.next - 1]) / (roc$sp[idx.next] - roc$sp[idx.next - 1])
        int.se <- roc$se[idx.next - 1] - proportion * (roc$se[idx.next - 1] - roc$se[idx.next])
        res <- c(NA, x, int.se)
      }
    }
    if (input == "sensitivity") {
      if (x < 0 || x > ifelse(roc$percent, 100, 1))
        stop("Input sensitivity not within the ROC space.")
      if (x %in% roc$se) {
        idx <- length(roc$se) + 1 - match(TRUE, rev(roc$se) == x)
        res <- c(roc$thresholds[idx], roc$sp[idx], roc$se[idx])
      }
      else { # need to interpolate
        idx.next <- match(TRUE, roc$se < x)
        proportion <- (x - roc$se[idx.next]) / (roc$se[idx.next - 1] - roc$se[idx.next])
        int.sp <- roc$sp[idx.next] + proportion * (roc$sp[idx.next - 1] - roc$sp[idx.next])
        res <- c(NA, int.sp, x)
      }
    }
    # Deduce additional tn, tp, fn, fp, npv, ppv
    ncases <- ifelse(is(roc, "smooth.roc"), length(attr(roc, "roc")$cases), length(roc$cases))
    ncontrols <- ifelse(is(roc, "smooth.roc"), length(attr(roc, "roc")$controls), length(roc$controls))
    se <- res[3]
    sp <- res[2]
    if (roc$percent) {
      tp <- se * ncases / 100
      fn <- ncases - tp
      tn <- sp * ncontrols / 100
      fp <- ncontrols - tn
      substr.percent <- 100
      npv <- 100 * tn / (tn + fn)
      ppv <- 100 * tp / (tp + fp)
    }
    else {
      tp <- se * ncases
      fn <- ncases - tp
      tn <- sp * ncontrols
      fp <- ncontrols - tn
      npv <- tn / (tn + fn)
      ppv <- tp / (tp + fp)
      substr.percent <- 1
    }
    accuracy <- (tp + tn) / (tp + tn + fp + fn)
    if (as.list) {
      list <- list(threshold=res[1], specificity=sp, sensitivity=se, accuracy=accuracy, tn=tn, tp=tp, fn=fn, fp=fp, npv=npv, ppv=ppv, "1-specificity"=substr.percent-sp, "1-sensitivity"=substr.percent-se, "1-accuracy"=substr.percent-accuracy, "1-npv"=substr.percent-npv, "1-ppv"=substr.percent-ppv)
      list <- list[ret]
      if (drop == FALSE) {
        list <- list(list)
        names(list) <- x
      }
      return(list)
    }
    else {
      res <- as.matrix(res)
      res <- rbind(res, accuracy, tn, tp, fn, fp, npv, ppv, substr.percent-sp, substr.percent-se, substr.percent-accuracy, substr.percent-npv, substr.percent-ppv)
      rownames(res) <- c("threshold", "specificity", "sensitivity", "accuracy", "tn", "tp", "fn", "fp", "npv", "ppv", "1-specificity", "1-sensitivity", "1-accuracy", "1-npv", "1-ppv")
      colnames(res) <- x
      return(res[ret,, drop=drop])
    }
  }
  else {
    stop("'x' must be a numeric or character vector.")
  }
}


#' Compute evaluation criteria for PLS, sPLS, PLS-DA, sPLS-DA, MINT and DIABLO
#' 
#' Function to evaluate the performance of the fitted PLS, sparse PLS, PLS-DA,
#' sparse PLS-DA, MINT (mint.splsda) and DIABLO (block.splsda) models using
#' various criteria.
#' 
#' Procedure. The process of evaluating the performance of a fitted model
#' \code{object} is similar for all PLS-derived methods; a cross-validation
#' approach is used to fit the method of \code{object} on \code{folds-1}
#' subsets of the data and then to predict on the subset left out. Different
#' measures of performance are available depending on the model. Parameters
#' such as \code{logratio}, \code{multilevel}, \code{keepX} or \code{keepY} are
#' retrieved from \code{object}.
#' 
#' Parameters. If \code{validation = "Mfold"}, M-fold cross-validation is
#' performed. \code{folds} specifies the number of folds to generate. The folds
#' also can be supplied as a list of vectors containing the indexes defining
#' each fold as produced by \code{split}; the folds are then fixed and
#' \code{nrepeat} is set to 1 for DIABLO objects. When using \code{validation =
#' "Mfold"} with a number of folds, make sure that you repeat the process
#' several times (as the results will be highly dependent on the random splits
#' and the sample size).
#' 
#' If \code{validation = "loo"}, leave-one-out cross-validation is performed
#' (in that case, there is no need to repeat the process).
#' 
#' Measures of performance. For fitted PLS and sPLS regression models,
#' \code{perf} estimates the mean squared error of prediction (MSEP),
#' \eqn{R^2}, and \eqn{Q^2} to assess the predictive perfity of the model using
#' M-fold or leave-one-out cross-validation. Note that only the \code{classic},
#' \code{regression} and \code{invariant} modes can be applied. For sPLS, the
#' MSEP, \eqn{R^2}, and \eqn{Q^2} criteria are averaged across all folds. Note
#' that for PLS and sPLS objects, perf is performed on the pre-processed data
#' after log ratio transform and multilevel analysis, if any.
#' 
#' Sparse methods. The sPLS, sPLS-DA and sgccda functions are run on several
#' and different subsets of data (the cross-folds) and will certainly lead to
#' different subset of selected features. Those are summarised in the output
#' \code{features$stable} (see output Value below) to assess how often the
#' variables are selected across all folds. Note that for PLS-DA and sPLS-DA
#' objects, perf is performed on the original data, i.e. before the
#' pre-processing step of the log ratio transform and multilevel analysis, if
#' any. In addition for these methods, the classification error rate is
#' averaged across all folds.
#' 
#' The mint.sPLS-DA function estimates errors based on Leave-one-group-out
#' cross validation (where each levels of object$study is left out (and
#' predicted) once) and provides study-specific outputs
#' (\code{study.specific.error}) as well as global outputs
#' (\code{global.error}). Note the mint perf methods do not use \code{seed}
#' or \code{BPPARAM} arguments. 
#' 
#' AUROC. For PLS-DA, sPLS-DA, mint.PLS-DA, mint.sPLS-DA, and block.splsda
#' methods: if \code{auc=TRUE}, Area Under the Curve (AUC) values are
#' calculated from the predicted scores obtained from the \code{predict}
#' function applied to the internal test sets in the cross-validation process,
#' either for all samples or for study-specific samples (for mint models).
#' Therefore we minimise the risk of overfitting. For block.splsda model, the
#' calculated AUC is simply the blocks-combined AUC for each component
#' calculated using \code{auroc.sgccda}.  See \code{\link{auroc}} for more
#' details. Our multivariate supervised methods already use a prediction
#' threshold based on distances (see \code{predict}) that optimally determine
#' class membership of the samples tested. As such AUC and ROC are not needed
#' to estimate the performance of the model. We provide those outputs as
#' complementary performance measures. See more details in our mixOmics
#' article.
#' 
#' Prediction distances. See details from \code{?predict}, and also our
#' supplemental material in the mixOmics article.
#' 
#' Repeats of the CV-folds. Repeated cross-validation implies that the whole CV
#' process is repeated a number of times (\code{nrepeat}) to reduce variability
#' across the different subset partitions. In the case of Leave-One-Out CV
#' (\code{validation = 'loo'}), each sample is left out once (\code{folds = N}
#' is set internally) and therefore nrepeat is by default 1.
#' 
#' BER is appropriate in case of an unbalanced number of samples per class as
#' it calculates the average proportion of wrongly classified samples in each
#' class, weighted by the number of samples in each class. BER is less biased
#' towards majority classes during the performance assessment.
#' 
#' For \code{sgccda} objects, we provide weighted measures (e.g. error rate) in
#' which the weights are simply the correlation of the derived components of a
#' given block with the outcome variable Y.
#' 
#' More details about the PLS modes in \code{?pls}.
#'
#' @param object object of class inherited from \code{"pls"}, \code{"plsda"},
#' \code{"spls"}, \code{"splsda"} or \code{"mint.splsda"}. The function will
#' retrieve some key parameters stored in that object.
#' @param dist only applies to an object inheriting from \code{"plsda"},
#' \code{"splsda"} or \code{"mint.splsda"} to evaluate the classification
#' performance of the model. Should be a subset of \code{"max.dist"},
#' \code{"centroids.dist"}, \code{"mahalanobis.dist"}. Default is \code{"all"}.
#' See \code{\link{predict}}.
#' @param validation character.  What kind of (internal) validation to use,
#' matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is
#' \code{"Mfold"}.
#' @param folds either a positive integer giving the number of folds in the
#' Mfold cross-validation, or a list of integer vectors where each element
#' gives the row indices of the samples in one test fold. A list is accepted
#' for \code{pls}, \code{spls} and DIABLO (\code{block.plsda},
#' \code{block.splsda}) objects. See Details.
#' @param nrepeat Number of times the Cross-Validation process is repeated.
#' This is an important argument to ensure the estimation of the performance to
#' be as accurate as possible.
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' performance of the model.
#' @param progressBar by default set to \code{FALSE} to output the progress bar
#' of the computation.
#' @param signif.threshold numeric between 0 and 1 indicating the significance
#' threshold required for improvement in error rate of the components. Default
#' to 0.01.
#' @template arg/BPPARAM 
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' Note 'seed' is not required or used in perf.mint.plsda as this method uses loo cross-validation
#' @param ... not used
#' 
#' @return For PLS and sPLS models, \code{perf} produces a list with the
#' following components for every repeat: 
#' \item{MSEP}{Mean Square Error Prediction for each \eqn{Y} variable, only 
#' applies to object inherited from \code{"pls"}, and \code{"spls"}. Only 
#' available when in regression (s)PLS.} 
#' \item{RMSEP}{Root Mean Square Error Prediction for each \eqn{Y} variable, only 
#' applies to object inherited from \code{"pls"}, and \code{"spls"}. Only 
#' available when in regression (s)PLS.} 
#' \item{R2}{a matrix of \eqn{R^2} values of the \eqn{Y}-variables for models 
#' with \eqn{1, \ldots ,}\code{ncomp} components, only applies to object
#' inherited from \code{"pls"}, and \code{"spls"}. Only available when in 
#' regression (s)PLS.}
#' \item{Q2}{if \eqn{Y} contains one variable, a vector of \eqn{Q^2} values
#' else a list with a matrix of \eqn{Q^2} values for each \eqn{Y}-variable.
#' Note that in the specific case of an sPLS model, it is better to have a look
#' at the Q2.total criterion, only applies to object inherited from
#' \code{"pls"}, and \code{"spls"}. Only available when in regression (s)PLS.} 
#' \item{Q2.total}{a vector of \eqn{Q^2}-total values for models with \eqn{1, 
#' \ldots ,}\code{ncomp} components, only applies to object inherited from 
#' \code{"pls"}, and \code{"spls"}. Available in both (s)PLS modes.}
#' \item{RSS}{Residual Sum of Squares across all selected features and the 
#' components.}
#' \item{PRESS}{Predicted Residual Error Sum of Squares across all selected 
#' features and the components.}
#' \item{features}{a list of features selected across the 
#' folds (\code{$stable.X} and \code{$stable.Y}) for the \code{keepX} and
#' \code{keepY} parameters from the input object. Note, this will be \code{NULL} 
#' if using standard (non-sparse) PLS.} 
#' \item{cor.tpred, cor.upred}{Correlation between the 
#' predicted and actual components for X (t) and Y (u)} 
#' \item{RSS.tpred, RSS.upred}{Residual Sum of Squares between the
#' predicted and actual components for X (t) and Y (u)}
#' 
#' 
#'
#' For PLS-DA and sPLS-DA models, \code{perf} produces a matrix of classification
#' error rate estimation. The dimensions correspond to the components in the
#' model and to the prediction method used, respectively. Note that error rates
#' reported in any component include the performance of the model in earlier
#' components for the specified \code{keepX} parameters (e.g. error rate
#' reported for component 3 for \code{keepX = 20} already includes the fitted
#' model on components 1 and 2 for \code{keepX = 20}).
#' \item{error.rate}{Prediction error rate for each dist and measure}
#' \item{auc}{AUC values per component averaged over the \code{nrepeat}}
#' \item{auc.all}{AUC values per component per repeat}
#' \item{predict}{A list of length ncomp that os predicted values of each sample for each class}
#' \item{features}{a list of features selected across the folds ($stable.X) for the keepX parameters from the input object.}
#' \item{choice.ncomp}{Otimal number of components for the model for each prediction distance using one-sided t-tests that test 
#' for a significant difference in the mean error rate (gain in prediction) when components are added to the model.}
#' \item{class}{A list which gives the predicted class of each sample for each dist and each of the ncomp components}
#' 
#' For mint.splsda models, \code{perf} produces the following outputs:
#' \item{study.specific.error}{A list that gives BER, overall error rate and
#' error rate per class, for each study} \item{global.error}{A list that gives
#' BER, overall error rate and error rate per class for all samples}
#' \item{predict}{A list of length \code{ncomp} that produces the predicted
#' values of each sample for each class} \item{class}{A list which gives the
#' predicted class of each sample for each \code{dist} and each of the
#' \code{ncomp} components. Directly obtained from the \code{predict} output.}
#' \item{auc}{AUC values} \item{auc.study}{AUC values for each study in mint
#' models}
#' 
#' For sgccda models (i.e. block (s)PLS-DA models), \code{perf} produces the following outputs:
#' \item{error.rate}{Prediction error rate for each block of \code{object$X}
#' and each \code{dist}} \item{error.rate.per.class}{Prediction error rate for
#' each block of \code{object$X}, each \code{dist} and each class}
#' \item{predict}{Predicted values of each sample for each class, each block
#' and each component} \item{class}{Predicted class of each sample for each
#' block, each \code{dist}, each component and each nrepeat} \item{features}{a
#' list of features selected across the folds (\code{$stable.X} and
#' \code{$stable.Y}) for the \code{keepX} and \code{keepY} parameters from the
#' input object.} \item{AveragedPredict.class}{if more than one block, returns
#' the average predicted class over the blocks (averaged of the \code{Predict}
#' output and prediction using the \code{max.dist} distance)}
#' \item{AveragedPredict.error.rate}{if more than one block, returns the
#' average predicted error rate over the blocks (using the
#' \code{AveragedPredict.class} output)} \item{WeightedPredict.class}{if more
#' than one block, returns the weighted predicted class over the blocks
#' (weighted average of the \code{Predict} output and prediction using the
#' \code{max.dist} distance). See details for more info on weights.}
#' \item{WeightedPredict.error.rate}{if more than one block, returns the
#' weighted average predicted error rate over the blocks (using the
#' \code{WeightedPredict.class} output.)} \item{MajorityVote}{if more than one
#' block, returns the majority class over the blocks. NA for a sample means that
#' there is no consensus on the predicted class for this particular sample over
#' the blocks.} \item{MajorityVote.error.rate}{if more than one block, returns
#' the error rate of the \code{MajorityVote} output}
#' \item{WeightedVote}{if more than one block, returns the weighted majority
#' class over the blocks. NA for a sample means that there is no consensus on
#' the predicted class for this particular sample over the blocks.}
#' \item{WeightedVote.error.rate}{if more than one block, returns the error
#' rate of the \code{WeightedVote} output} \item{weights}{Returns the weights
#' of each block used for the weighted predictions, for each nrepeat and each
#' fold} 
#' \item{choice.ncomp}{For supervised models; returns the optimal number
#' of components for the model for each prediction distance using one-sided
#' t-tests that test for a significant difference in the mean error rate (gain
#' in prediction) when components are added to the model. See more details in
#' Rohart et al 2017 Suppl. For more than one block, an optimal ncomp is
#' returned for each prediction framework.}
#' 
#' @author Ignacio González, Amrit Singh, Kim-Anh Lê Cao, Benoit Gautier,
#' Florian Rohart, Al J Abadi
#' @seealso \code{\link{predict}}, \code{\link{nipals}},
#' \code{\link{plot.perf}}, \code{\link{auroc}} and \url{www.mixOmics.org} for
#' more details.
#' @references 
#' Singh A., Shannon C., Gautier B., Rohart F., Vacher M., Tebbutt S.
#' and Lê Cao K.A. (2019), DIABLO: an integrative approach for identifying key 
#' molecular drivers from multi-omics assays, Bioinformatics, 
#' Volume 35, Issue 17, 1 September 2019, Pages 3055–3062.
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' 
#' MINT:
#' 
#' Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017). MINT: A
#' multivariate integrative approach to identify a reproducible biomarker
#' signature across multiple experiments and platforms. BMC Bioinformatics
#' 18:128.
#' 
#' PLS and PLS citeria for PLS regression: Tenenhaus, M. (1998). \emph{La
#' regression PLS: theorie et pratique}. Paris: Editions Technic.
#' 
#' Chavent, Marie and Patouille, Brigitte (2003). Calcul des coefficients de
#' regression et du PRESS en regression PLS1. \emph{Modulad n}, \bold{30} 1-11.
#' (this is the formula we use to calculate the Q2 in perf.pls and perf.spls)
#' 
#' Mevik, B.-H., Cederkvist, H. R. (2004). Mean Squared Error of Prediction
#' (MSEP) Estimates for Principal Component Regression (PCR) and Partial Least
#' Squares Regression (PLSR). \emph{Journal of Chemometrics} \bold{18}(9),
#' 422-429.
#' 
#' sparse PLS regression mode:
#' 
#' Lê Cao, K. A., Rossouw D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. \emph{Statistical
#' Applications in Genetics and Molecular Biology} \bold{7}, article 35.
#' 
#' One-sided t-tests (suppl material):
#' 
#' Rohart F, Mason EA, Matigian N, Mosbergen R, Korn O, Chen T, Butcher S,
#' Patel J, Atkinson K, Khosrotehrani K, Fisk NM, Lê Cao K-A&, Wells CA&
#' (2016). A Molecular Classification of Human Mesenchymal Stromal Cells. PeerJ
#' 4:e1845.
#' @keywords regression multivariate
#' @export
#' @example ./examples/perf-examples.R
## ------------------------------- Generic -------------------------------- ##
perf <- function(object, ...)
    UseMethod("perf")

## ------------------------------------------------------------------------ ##
####                              (s)PLS(DA)                              ####
## ------------------------------------------------------------------------ ##

## -------------------------------- (s)PLS -------------------------------- ##
#' @rdname perf
#' @export
perf.mixo_pls <- function(object,
                          validation = c("Mfold", "loo"),
                          folds,
                          progressBar = FALSE,
                          nrepeat = 1,
                          BPPARAM = SerialParam(),
                          seed = NULL,
                          ...)
{
    # checking args and initialize params
    ncomp = object$ncomp
    spls.model <- is(object, 'mixo_spls')
    progressBar <- .check_logical(progressBar)
    BPPARAM$RNGseed <- seed
    
    # run CV in parallel depending on BPPARAM
    repeat_names <- .name_list(char = seq_len(nrepeat))
    result <- bplapply(X = repeat_names, FUN = function(nrep) {
        ## progress bar
        if (progressBar == TRUE)
            .progressBar(nrep/nrepeat)
        ## CV
        .perf.mixo_pls_cv(object, validation = validation, folds = folds, nrep = nrep)
    }, BPPARAM = BPPARAM)
    
    measures <- lapply(result, function(x){
        x$measures
    })
    
    measures <- Reduce(rbind, measures)
    measures <- as.data.frame(measures)
    
    ## R CMD check stuff
    measure <- feature <- comp <- block <- stability <- value <- NULL
    lower <- upper <- keepX <- keepY <- NULL
    
    measure.names <- .name_list(unique(measures$measure))
    measures <- lapply(measure.names, function(meas) {
        ## ------ value of measures across repeats
        df <- measures %>% 
            filter(measure == meas) %>% 
            mutate(measure = NULL) %>% 
            as.data.frame()
        
        ## ------ summary of measures across repeats
        df.summ <- df %>%  
            group_by(feature, comp) %>% 
            summarise(mean = mean(value, na.rm = TRUE), 
                      sd = sd(value, na.rm = TRUE)) %>% 
            as.data.frame()
        
        list(values = df, summary = df.summ)
    })

    ## ------ feature stability
    if (spls.model)
    {
        features <- lapply(result, function(x){
            x$features
        })
        
        features <- Reduce(rbind, features) %>% 
            group_by(feature, comp, block) %>% 
            summarise(stability = mean(stability, na.rm = TRUE))
        
        features <- as.data.frame(features)
        features <- lapply(list(stability.X = 'X', stability.Y = 'Y'), function(z){
            lapply(.name_list(unique(features$comp)), function(n.comp){
                
                df <- features
                df <- filter(df, block == z & comp == n.comp)
                df <- df[,c('feature', 'stability')]
                vec <- df$stability
                names(vec) <- df$feature
                sort(vec, decreasing = TRUE)
            })
        })
    } else
    {
        features <- NULL
    }
    
    result <- list(measures = measures,
                   features = features)
    mc <- mget(names(formals())[-1], sys.frame(sys.nframe()))
    ## replace function, object with unevaluated call
    mc <- as.call(c(as.list(match.call())[1:2], mc))
    result <- c(list(call = mc), result)
    class(result) <- "perf.pls.mthd"
    
    return(result)
    
}

#' @rdname perf
#' @export
perf.mixo_spls  <- perf.mixo_pls

#' @noRd
#' @keywords Internal
.perf.mixo_pls_cv <- function(object,
                              validation = c("Mfold", "loo"),
                              folds,
                              nrep = 1,
                              ...)
{
    # changes to bypass the loop for the Q2
    
    ## -------- checks -------- ##
    if (object$mode == 'invariant')
        stop("'perf' is only available for (s)pls with modes: 'regression', 'canonical' or 'classic'.  Object has mode 'invariant'", call. = FALSE)
    
    validation = match.arg(validation)
    
    ## ---------- CV ---------- ##
    ## ------------- initialise
    # these are the centered and scaled matrices output from pls, we remove $nzv if needed
    if (length(object$nzv$Position)>0)
    {
        X = object$X[, -object$nzv$Position]
    } else {
        X = object$X
    }
    Y = object$Y
    
    scale = object$scale
    tol = object$tol
    max.iter = object$max.iter
    mode = object$mode
    ncomp = object$ncomp
    n = nrow(X)
    p = ncol(X)
    q = ncol(Y)
    
    if (any(is.na(X)) || any(is.na(Y)))
        stop("missing data in 'X' and/or 'Y'. Use 'nipals' for dealing with NAs.", call. = FALSE)
    
    
    #-- tells which variables are selected in X and in Y --#
    if (is(object, "mixo_spls"))
    {
        keepX = object$keepX
        keepY = object$keepY
    } else {
        keepX = rep(ncol(X), ncomp)
        keepY = rep(ncol(Y), ncomp)
    }
    
    #-- define the folds --#
    if (validation == "Mfold")
    {
        if (is.list(folds))
        {
            # folds supplied as a list: validated by the shared helper
            folds = .check_folds(folds, n = n)
            M = length(folds)
        } else {
            if (is.null(folds) || !is.finite(folds) || folds < 2 || folds > n)
            {
                stop("Invalid number of folds.", call. = FALSE)
            } else {
                M = round(folds)
                folds = split(sample(1:n), rep(1:M, length = n))
            }
        }
    } else {
        folds = split(1:n, rep(1:n, length = n))
        M = n
    }
    
    #-- initialize new objects --#
    if (mode == 'canonical'){
        RSS = rbind(rep(n - 1, p), matrix(nrow = ncomp, ncol = p))
        # RSS.indiv is the reconstructed matrix X
        #RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = p)})
        #RSS.indiv[[1]] = X
        press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = p)})
        PRESS.inside = Q2 = matrix(nrow = ncomp, ncol = p)
    }else{
        RSS = rbind(rep(n - 1, q), matrix(nrow = ncomp, ncol = q))
        # RSS.indiv is the reconstructed matrix Y
        #RSS.indiv = lapply(1 : (ncomp + 1), function(x){matrix(NA, nrow = n, ncol = q)})
        #RSS.indiv[[1]] = Y # KA changed
        press.mat = lapply(1 : ncomp, function(x){matrix(NA, nrow = n, ncol = q)})
        PRESS.inside = Q2 = matrix(nrow = ncomp, ncol = q)
    }
    
    MSEP.mat = Ypred = array(0, c(n, q, ncomp))
    MSEP = R2 = matrix(nrow = ncomp, ncol = q)
    
    # to store the predicted components
    t.pred.cv = matrix(nrow = nrow(X), ncol = ncomp)
    u.pred.cv = matrix(nrow = nrow(X), ncol = ncomp)
    
    # to record feature stability, a list of form
    # list(X = list(comp1 = c(feature1 = 0.99, ...), 
    #               comp2 = c(feature2 = 0.98, ...)), 
    #      Y = ...)
    features <-
        lapply(list(X = X, Y = Y), function(Z){
            features <- vector(mode = 'numeric', length = ncol(Z))
            names(features) <- colnames(Z)
            features <- lapply(seq_len(ncomp), function(x) features)
            names(features) <- paste0('comp', seq_len(ncomp))
            
            return(features)
        })
    
    
    # ====  loop on h = ncomp is only for the calculation of Q2 on each component
    for (h in 1:ncomp)
    {
        #-- initialising arguments --#
        tt = object$variates$X[, h]
        u = object$variates$Y[, h]
        b = object$loadings$Y[, h]
        #nx = p - keepX[h]
        #ny = q - keepY[h]
        
        # only used for matrices deflation across dimensions
        c = crossprod(X, tt)/drop(crossprod(tt))  #object$mat.c[, h]
        d = crossprod(Y, tt)/drop(crossprod(tt))  #object$mat.d[, h]
        e = crossprod(Y, u)/drop(crossprod(u))    
        
        # deflate matrices
        X = X - tt %*% t(c)
        
        #-- mode classic
        if (mode == "classic")
            Y = Y - tt %*% t(b)
        #-- mode regression
        if (mode == "regression")
            Y = Y - tt %*% t(d)
        #-- mode canonical 
        if (mode == "canonical")
            Y = Y - u %*% t(e)
        #-- mode invariant: Y is unchanged
        
        # update RSS for X/Y deflated
        if(mode == 'canonical'){  # based on X
            RSS[h + 1, ] =  colSums((X)^2)   # ==  colSums((X - tt %*% t(c))^2) if we had not deflated
        }else{ # regression, invariant, classic
            RSS[h + 1, ] = colSums((Y)^2)  # 
        }
        
    } # end h to calculate RSS   
    
    
    
    # ======== loop on i for cross validation ===================#
    for (i in 1:M)
    {
        # initialise the train / test datasets
        omit = folds[[i]]
        X.train = object$X[-omit, , drop = FALSE]
        Y.train = object$Y[-omit, , drop = FALSE]
        X.test = object$X[omit, , drop = FALSE]
        Y.test = object$Y[omit, , drop = FALSE]
        
        # New loop to calculate prediction
        for (h in 1:ncomp)
        { 
            #-- for MSEP and R2 criteria, no loop on the component as we do a spls with ncomp
            ##if (h == 1)
            #{
            #nzv = (apply(X.train, 2, var) > .Machine$double.eps) # removed in v6.0.0 so that MSEP, R2 and Q2 are obtained with the same data
            # re-added in >6.1.3 to remove constant variables
            nzv.X = (apply(X.train, 2, var) > .Machine$double.eps)
            nzv.Y = (apply(Y.train, 2, var) > .Machine$double.eps)
            
            # creating a keepX/Y.temp that can change for each fold, depending on nzv.X/Y
            keepX.temp = keepX
            keepY.temp = keepY
            if(any(keepX.temp > sum(nzv.X)))
                keepX.temp[which(keepX.temp>sum(nzv.X))] = sum(nzv.X)
            if(any(keepY.temp > sum(nzv.Y)))
                keepY.temp[which(keepY.temp>sum(nzv.Y))] = sum(nzv.Y)
            # TODO clarify the iterative nzv process in docs -- give it a better name (these are actually !nzv)
            # here h = 1 because we deflate at each step then extract the vectors for each h comp
            spls.res = spls(X.train[, nzv.X, drop = FALSE], Y.train[, nzv.Y, drop = FALSE], ncomp = 1, mode = mode, max.iter = max.iter, tol = tol, 
                            keepX = keepX.temp[h], keepY = keepY.temp[h], near.zero.var = FALSE, scale = scale)
            Y.hat = predict.mixo_spls(spls.res, X.test[, nzv.X, drop = FALSE])$predict
            
            # added the stop msg
            if(sum(is.na(Y.hat))>0) stop('Predicted Y values include NA')  
            
            # replaced h by 1; Y.hat is the prediction of the test samples for all q variable in comp h = 1
            Ypred[omit, nzv.Y, h] = Y.hat[, , 1]
            MSEP.mat[omit, nzv.Y, h] = (Y.test[, nzv.Y] - Y.hat[, , 1])^2
            
            
            # Q2 criterion: buidling directly from spls object
            u.cv = spls.res$variates$Y[, 1]
            t.cv = spls.res$variates$X[, 1]
            a.cv = spls.res$loadings$X[, 1]
            b.cv = spls.res$loadings$Y[, 1, drop = FALSE]
            
            # reg coefficients:
            c.cv = crossprod(X.train, u.cv) / drop(crossprod(u.cv)) 
            d.cv = crossprod(Y.train, t.cv) / drop(crossprod(t.cv)) # d.cv \neq to b.cv as d.cv is normed wrt to t.cv
            e.cv = crossprod(Y.train, u.cv) / drop(crossprod(u.cv)) 
            
            # calculate predicted components and store
            t.pred = c(X.test %*% a.cv)
            t.pred.cv[omit,h] = t.pred    # needed for tuning
            b.pred = crossprod(Y.test, t.pred)
            b.pred.cv = b.pred/ drop(sqrt(crossprod(b.pred)))
            u.pred.cv[omit,h] = Y.test[, nzv.Y] %*% b.cv  # needed for tuning, changed instead of b.pred.cv
            
            # predicted reg coeff, could be removed
            e.pred.cv = crossprod(as.matrix(Y.test), Y.test %*% b.pred.cv) / drop(crossprod(Y.test %*% b.pred))
            d.pred.cv = crossprod(as.matrix(Y.test), t.pred) / drop(crossprod(t.pred)) 
            
            # deflate matrices X
            X.train = X.train - t.cv %*% t(c.cv)
            X.test = X.test - t.pred %*% t(c.cv)
            # deflate matrices X      
            #-- mode classic
            if (mode == "classic"){
                Y.train[, nzv.Y] = Y.train[, nzv.Y] - t.cv %*% t(b.cv)  # could be pred on b
                Y.test[, nzv.Y] = Y.test[, nzv.Y] - t.pred %*% t(b.cv)
            }
            #-- mode regression
            if (mode == "regression"){
                Y.train = Y.train - t.cv %*% t(d.cv) # could be pred d.pred.cv? does not decrease enough
                Y.test[, nzv.Y] = Y.test[, nzv.Y] - Y.hat[, , 1]   # == Y.test - t.pred %*% t(d.cv) 
            }
            
            #-- mode canonical  ## KA added
            if (mode == "canonical"){
                Y.train = Y.train - u.cv %*% t(e.cv)  # could be pred on e
                Y.test = Y.test - (Y.test[, nzv.Y] %*% b.cv) %*% t(e.cv)  # here u.pred = Y.test %*% b.cv (b.pred.cv gives similar results)
            }
            #-- mode invariant: Y is unchanged
            
            # calculate predicted matrix X.hat or Y.hat based on X.test
            if(mode == 'canonical'){  # Xa c' = t c'
                #X.hat.cv = t.pred %*% t(c.cv), calculated earlier
                press.mat[[h]][omit, ] = X.test        # == X.test - X.hat.cv
            }else{ #  if(mode == 'regression'){  # Xa d' = t d'
                #Y.hat.cv = t.pred %*% t(d.cv), calculated earlier
                press.mat[[h]][omit, ] = Y.test        # == Y.test - Y.hat.cv
            }  
            
            
            # Record selected features in each set
            if (is(object,"mixo_spls"))
            {
                X.feature <- as.numeric(names(features$X[[h]]) %in% selectVar(spls.res, comp = 1)$X$name)
                Y.feature <- as.numeric(names(features$Y[[h]]) %in% selectVar(spls.res, comp = 1)$Y$name)
                # TODO using comp = 1 after deflation: this is problematic if, say, folds = 3, keepX = c(2, 100) (max 4 features (2 folds x 2 features) should be output for comp2 before calculating stability)
                features$X[[h]] <- features$X[[h]] + X.feature / length(folds)
                features$Y[[h]] <- features$Y[[h]] + Y.feature / length(folds)
            }
            
        } #  end loop on h ncomp
    } # end i (cross validation)
    
    
    # store results for each comp
    for (h in 1:ncomp){
        #-- compute the Q2 criterion --#
        # norm is equivalent to summing here the squared press values:
        PRESS.inside[h, ] = apply(press.mat[[h]], 2, function(x){norm(x, type = "2")^2})
        
        if(mode != 'canonical'){
            Q2[h, ] = 1 - PRESS.inside[h, ] / RSS[h, ]
            MSEP[h, ] = apply(as.matrix(MSEP.mat[, , h]), 2, mean)
            R2[h, ] = (diag(cor(object$Y, Ypred[, , h])))^2
        } # if mode == canonical, do not output
    }
    
    #-- output -----------------------------------------------------------------#
    #---------------------------------------------------------------------------#
    Q2.total = matrix(1 - rowSums(PRESS.inside) / rowSums(RSS[-(ncomp+1), , drop = FALSE]),
                      nrow = 1, ncol = ncomp,
                      dimnames = list("Q2.total", paste0("comp", seq_len(ncomp))))
    
    # set up dimnames and outputs
    result = list()
    
    if(mode != 'canonical'){
        rownames(MSEP) = rownames(R2) = rownames(Q2) = paste0("comp", seq_len(ncomp))
        colnames(MSEP) = colnames(R2) = colnames(Q2) = object$names$colnames$Y
        
        result$MSEP = t(MSEP)
        result$RMSEP = sqrt(t(MSEP))
        #result$MSEP.mat = MSEP.mat  
        result$R2 = t(R2)
        result$Q2 = t(Q2)  # remove this output when canonical mode?
    }
    
    result$Q2.total =  Q2.total
    RSS <- t(RSS) ## bc all others are transposed
    PRESS = t(PRESS.inside)
    result$RSS <- RSS[,-1, drop = FALSE] ## drop q/p
    result$PRESS <- PRESS
    if (ncol(object$Y) > 1)
    {
        # TODO ensure these are in fact no more necessary
        #result$d.cv = d.cv  # KA added  
        #result$b.cv = b.cv  # KA added 
        #result$c.cv = c.cv  # KA added 
        #result$u.cv = u.cv  # KA added 
        #result$a.cv = a.cv  # KA added 
        #result$t.pred.cv = t.pred.cv  # needed for tuning
        #result$u.pred.cv = u.pred.cv  # needed for tuning
        
        # extract the predicted components per dimension, take abs value
        result$cor.tpred = diag(abs(cor(t.pred.cv, object$variates$X)))
        result$cor.tpred = t(data.matrix(result$cor.tpred, rownames.force = TRUE))
        result$cor.upred = diag(abs(cor(u.pred.cv, object$variates$Y)))
        result$cor.upred = t(data.matrix(result$cor.upred, rownames.force = TRUE))
        
        # RSS: no abs values here
        result$RSS.tpred = apply((t.pred.cv - object$variates$X)^2, 2, sum)/(nrow(X) -1)
        result$RSS.tpred  = t(data.matrix(result$RSS.tpred, rownames.force = TRUE))
        result$RSS.upred = apply((u.pred.cv - object$variates$Y)^2, 2, sum)/(nrow(X) -1)
        result$RSS.upred  = t(data.matrix(result$RSS.upred, rownames.force = TRUE))
    }
    result <- mapply(result, names(result), FUN = function(arr, measure) {
        arr <- data.matrix(arr)
        col.names <- seq_len(ncomp)
        if (ncol(arr) == ncomp)
            colnames(arr) <- col.names
        else
            stop("unexpected dimension for entry in perf measures: ", measure)
        if (nrow(arr) == 1)
            rownames(arr) <- measure
        arr
    }, SIMPLIFY = FALSE)
    
    ## melt by comp
    result <- lapply(result, FUN = function(arr, nrep) {
        arr <- melt(arr)
        colnames(arr) <- c('feature', 'comp', 'value')
        if (nlevels(arr$feature) == 1) ## for Y-level measures (ass opossed to Y_feature level) such as Q2.total
            arr$feature <- factor('Y')
        arr$nrep <- nrep
        arr
    }, nrep = nrep)
    col.names <- names(result[[1]])
    #' @importFrom reshape2 melt
    result <- melt(result, id.vars = col.names)
    colnames(result) <- c(col.names, 'measure')
    
    result <- list(measures = result)
    #---- stability of features -----#
    if (is(object, "mixo_spls"))
    {
        features <- lapply(features, function(x)
        {
            x <- lapply(x, function(stab) round(stab, 2))
            df <- data.frame(x)
            df <- data.frame(feature = rownames(df), df)
            df
        })
        features <- melt(features, id.vars = 'feature', value.name = 'stability', variable.name = 'comp')
        ## add block name column instead of default 'L1'
        colnames(features) <- c(rev(rev(colnames(features))[-1]), 'block')
        features$nrep <- nrep
        
        result$features <- features
    }
    return(invisible(result))
}

## ------------------------------- (s)PLSDA ------------------------------- ##
#' @rdname perf
#' @importFrom methods hasArg
#' @export
perf.mixo_plsda <- function(object,
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
    
    #-- initialising arguments --#
    set.seed(seed)
    BPPARAM$RNGseed <- seed
    # these data are the centered and scaled X output or the unmapped(Y) scaled and centered
    X = object$input.X
    level.Y = object$names$colnames$Y  #to make sure the levels are ordered
    Y = object$Y
    ncomp = object$ncomp
    n = nrow(X)
    
    logratio = object$logratio
    if (is.null(logratio))
        logratio = "none"
    
    multilevel = object$multilevel # repeated measurement and Y
    near.zero.var = !is.null(object$nzv) # if near.zero.var was used, we set it to TRUE. if not used, object$nzv is NULL
    
    #-- tells which variables are selected in X and in Y --#
    
    if (is(object, "mixo_splsda"))
    {
        keepX = object$keepX
    } else {
        keepX = rep(ncol(X), ncomp)
    }
    
    tol = object$tol
    max.iter = object$max.iter
    scale = object$scale
    
    # initialize new objects:
    features = list()
    for(k in 1:ncomp)
        features[[k]] = NA
    
    # check input arguments
    
    if (hasArg(method.predict))
        stop("'method.predict' argument has been replaced by 'dist' to match the 'tune' function")
    method.predict = NULL # to pass R CMD check
    
    dist = match.arg(dist, choices = c("all", "max.dist", "centroids.dist", "mahalanobis.dist"), several.ok = TRUE)
    if (any(dist == "all"))
    {
        nmthdd = 3
        dist = c("max.dist", "centroids.dist", "mahalanobis.dist")
    } else {
        nmthdd = length(dist)
    }
    
    if (length(validation) > 1 )
        validation = validation [1]
    if (!(validation %in% c("Mfold", "loo")))
        stop("Choose 'validation' among the two following possibilities: 'Mfold' or 'loo'")
    
    if (validation == "loo")
    {
        if (nrepeat != 1)
            warning("Leave-One-Out validation does not need to be repeated: 'nrepeat' is set to '1'.")
        nrepeat = 1
    }

    if (any(table(object$Y) <= 1)) {
        stop(paste("Cannot evaluate performance when a class level ('", 
                   names(table(object$Y))[which(table(object$Y) == 1)],
                   "') has only a single associated sample.", sep = ""))
    }
    
    if (nrepeat < 3 && validation != "loo") {
        warning("Values in '$choice.ncomp' will reflect component count with the minimum error rate rather than the best based on a one-way t.test")
    }
    
    if (!is.logical(progressBar))
        stop("'progressBar' must be either TRUE or FALSE")
    
    measure = c("overall","BER") # one of c("overall","BER")
    
    
    if (!(logratio %in% c("none", "CLR")))
        stop("Choose one of the two following logratio transformation: 'none' or 'CLR'")
    #fold is checked in 'MCVfold'
    
    
    
    
    #-- check significance threshold
    signif.threshold <- .check_alpha(signif.threshold)
    
    #---------------------------------------------------------------------------#
    #-- logration + multilevel approach ----------------------------------------#
    # we can do logratio and multilevel on the whole data as these transformation are done per sample
    X = logratio.transfo(X = X, logratio = logratio)
    if (!is.null(multilevel))
    {
        Xw = withinVariation(X, design = multilevel)
        X = Xw
    }
    #-- logratio + multilevel approach -----------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    # -------------------------------------
    # added: first check for near zero var on the whole data set
    if (near.zero.var == TRUE)
    {
        nzv = nearZeroVar(X)
        if (length(nzv$Position > 0))
        {
            warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
            X = X[, -nzv$Position, drop=TRUE]
            
            if (ncol(X)==0)
                stop("No more predictors after Near Zero Var has been applied!")
            
            if (any(keepX > ncol(X)))
                keepX = ncol(X)
            
        }
    }
    # and then we start from the X data set with the nzv removed
    
    #---------------------------------------------------------------------------#
    #-- NA calculation      ----------------------------------------------------#
    
    misdata = c(X=anyNA(X), Y=FALSE) # Detection of missing data. we assume no missing values in the factor Y
    
    if (any(misdata))
    {
        is.na.A = is.na(X)
        
        #ind.NA = which(apply(is.na.A, 1, sum) > 0) # calculated only once
        #ind.NA.col = which(apply(is.na.A, 2, sum) > 0) # calculated only once
    } else {
        is.na.A = NULL
        #ind.NA = ind.NA.col = NULL
    }
    #-- NA calculation      ----------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    list.features = list()
    
    mat.error.rate = mat.sd.error = mat.mean.error = error.per.class.keepX.opt = error.per.class.keepX.opt.mean = list()
    error.per.class = list()
    final=list()
    
    for (measure_i in measure)
    {
        mat.sd.error[[measure_i]] = matrix(0,nrow = ncomp, ncol = length(dist),
                                           dimnames = list(c(paste0('comp', 1 : ncomp)), dist))
        mat.mean.error[[measure_i]] = matrix(0,nrow = ncomp, ncol = length(dist),
                                             dimnames = list(c(paste0('comp', 1 : ncomp)), dist))
        error.per.class.keepX.opt[[measure_i]] = list()
        error.per.class.keepX.opt.mean[[measure_i]] = list()
        mat.error.rate[[measure_i]]=list()
        for(ijk in dist)
        {
            mat.error.rate[[measure_i]][[ijk]] = matrix(0, nrow = ncomp, ncol = nrepeat,
                                                        dimnames = list(c(paste0('comp', 1 : ncomp)), c(paste0('nrep', 1 : nrepeat))))
            
            error.per.class.keepX.opt[[measure_i]][[ijk]] = array(0, c(nlevels(Y), nrepeat, ncomp),
                                                                  dimnames = list(c(levels(Y)), c(paste0('nrep', 1 : nrepeat)), c(paste0('comp', 1:ncomp, sep=''))))
            
            error.per.class.keepX.opt.mean[[measure_i]][[ijk]] = matrix(nrow = nlevels(Y), ncol = ncomp,
                                                                        dimnames = list(c(levels(Y)), c(paste0('comp', 1 : ncomp))))
        }
    }
    
    if(auc == TRUE)
    {
        auc.mean=list()
        auc.all=list()
    }
    
    prediction.all = class.all = auc.mean = auc.all = list()
    for(ijk in dist)
    {
        class.all[[ijk]] = array(0, c(nrow(X),  nrepeat ,ncomp),
                                 dimnames = list(rownames(X),c(paste0('nrep', 1 : nrepeat)),c(paste0('comp', 1 : ncomp))))
    }
    
    class.object=class(object)
    
    for (comp in 1 : ncomp)
    {
        if (progressBar == TRUE)
            cat("\ncomp",comp, "\n")
        
        
        if(comp > 1)
        {
            choice.keepX = keepX[1 : (comp - 1)]
        } else {
            choice.keepX = NULL
        }
        test.keepX = keepX[comp]
        names(test.keepX) = test.keepX
        #test.keepX is a value
        
        # estimate performance of the model for each component
        set.seed(seed)
        result = MCVfold.spls (X, Y, multilevel = multilevel, validation = validation, folds = folds, nrepeat = nrepeat, ncomp = comp,
                               choice.keepX = choice.keepX, test.keepX = test.keepX, test.keepY = nlevels(Y),
                               measure = measure, dist = dist, scale=scale,
                               near.zero.var = near.zero.var,
                               auc = auc, progressBar = progressBar, class.object = class.object,
                               BPPARAM = BPPARAM,
                               misdata = misdata, is.na.A = is.na.A)#, ind.NA = ind.NA, ind.NA.col = ind.NA.col)
        
        # ---- extract stability of features ----- # NEW
        if (is(object, "mixo_splsda"))
            list.features[[comp]] = result$features$stable
        
        for (ijk in dist)
        {
            for (measure_i in measure)
            {
                mat.error.rate[[measure_i]][[ijk]][comp,] = result[[measure_i]]$mat.error.rate[[ijk]][1,]
                mat.mean.error[[measure_i]][comp, ijk]=result[[measure_i]]$error.rate.mean[[ijk]]
                if (!is.null(result[[measure_i]]$error.rate.sd))
                {
                    mat.sd.error[[measure_i]][comp, ijk]=result[[measure_i]]$error.rate.sd[[ijk]]
                } else {
                    mat.sd.error= NULL
                }
                # confusion matrix for keepX.opt, for each nrep
                error.per.class.keepX.opt[[measure_i]][[ijk]][ , ,comp] = result[[measure_i]]$confusion[[ijk]]
                
                # confusion matrix for keepX.opt, averaged over all nrep
                error.per.class.keepX.opt.mean[[measure_i]][[ijk]][ ,comp] = apply(result[[measure_i]]$confusion[[ijk]],1 , mean)
            }
            
            #prediction of each samples for each fold and each repeat, on each comp
            class.all[[ijk]][, , comp] = result$class.comp[[ijk]][,,1]
        }
        prediction.all[[comp]] = array(unlist(result$prediction.comp),c(nrow(result$prediction.comp[[1]]), ncol(result$prediction.comp[[1]]), nrepeat),
                                       dimnames = c(dimnames(result$prediction.comp[[1]])[1:2], list(paste0("nrep",1:nrepeat))))#[[1]][, , 1] #take only one component [[1]] and one of test.keepX [,,1]
        
        if(auc == TRUE)
        {
            auc.all[[comp]] = lapply(result$auc.all, function(x) x[,,1])
            auc.mean[[comp]] = result$auc[, , 1]
        }
    }
    
    names(prediction.all) = paste0("comp", seq_len(ncomp))
    
    # calculating the number of optimal component based on t.tests and the error.rate.all, if more than 3 error.rates(repeat>3)
    ncomp_opt = matrix(NA, nrow = length(measure), ncol = length(dist),
                       dimnames = list(measure, dist))
    
    for (measure_i in measure)
    {
        for (ijk in dist) {
            if(nrepeat > 2 & ncomp >1)
            {
                ncomp_opt[measure, ijk] = t.test.process(t(mat.error.rate[[measure_i]][[ijk]]), alpha = signif.threshold)
            }
            else
            {
                if(length(min(t(rowMeans(mat.error.rate[[measure_i]][[ijk]])))) > 1){
                    # avoid multiple replacements - fix for issue #303
                    print("More than one component found with same minimum error rate, selecting smaller component...")
                    ncomp_opt[measure, ijk] = which(t(rowMeans(mat.error.rate[[measure_i]][[ijk]])) == min(t(rowMeans(mat.error.rate[[measure_i]][[ijk]]))))[1]
                } else{
                    ncomp_opt[measure, ijk] = which(t(rowMeans(mat.error.rate[[measure_i]][[ijk]])) == min(t(rowMeans(mat.error.rate[[measure_i]][[ijk]]))))
                }
            }
        }
    }
    
    
    
    
    result = list(error.rate = mat.mean.error,
                  error.rate.sd = mat.sd.error,
                  error.rate.all = mat.error.rate,
                  error.rate.class = error.per.class.keepX.opt.mean[[1]],
                  error.rate.class.all = error.per.class.keepX.opt[[1]],
                  predict = prediction.all,
                  class = class.all,
                  choice.ncomp = ncomp_opt)
    
    if(auc)
    {
        names(auc.mean) = c(paste0("comp", seq_len(ncomp)))
        result$auc = auc.mean
        
        names(auc.all) = c(paste0("comp", seq_len(ncomp)))
        result$auc.all =auc.all
    }
    
    if (is(object, "mixo_splsda"))
    {
        names(list.features) = paste0("comp", seq_len(ncomp))
        result$features$stable = list.features
    }
    
    if (progressBar == TRUE)
        cat('\n')
    
    # added
    if (near.zero.var == TRUE)
        result$nzvX = nzv$Position
    
    if (is(object, "mixo_splsda"))
    {
        method = "splsda.mthd"
    } else if (is(object, "mixo_plsda")) {
        method = "plsda.mthd"
    } else {
        warning("Something that should not happen happened. Please contact us.")
    }
    class(result) = c("perf",paste(c("perf", method), collapse ="."))
    result$call = match.call()
    
    
    #updated outputs
    return(invisible(result))
}

#' @rdname perf
#' @export
perf.mixo_splsda <- perf.mixo_plsda

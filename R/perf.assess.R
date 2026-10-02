#' Compute evaluation criteria for PLS, sPLS, PLS-DA, sPLS-DA, MINT and DIABLO
#' 
#' Function to evaluate the performance of the fitted PLS, sparse PLS, PLS-DA,
#' sparse PLS-DA, MINT (mint.splsda) and DIABLO (block.splsda) models using
#' various criteria.
#' 
#' This function is built upon 'perf()' but instead of assessing model performance 
#' across components 1:ncomp only assesses performance of the given model
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
#' calculated AUC is simply the blocks-combined AUC
#' calculated using \code{auroc.sgccda}.  See \code{\link{auroc}} for more
#' details. Our multivariate supervised methods already use a prediction
#' threshold based on distances (see \code{predict}) that optimally determine
#' class membership of the samples tested. As such AUC and ROC are not needed
#' to estimate the performance of the model. We provide those outputs as
#' complementary performance measures.
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
#' \code{"spls"}, \code{"splsda"}. \code{"sgccda"} or \code{"mint.splsda"}. The function will
#' retrieve some key parameters stored in that object.
#' @param validation a character string.  What kind of (internal) validation to use,
#' matching one of \code{"Mfold"} or \code{"loo"} (see below). Default is
#' \code{"Mfold"}. For MINT methods only \code{"loo"} will be used. 
#' @param folds either a positive integer giving the number of folds in the
#' Mfold cross-validation, or a list of integer vectors where each element
#' gives the row indices of the samples in one test fold. A list is accepted
#' for \code{pls}, \code{spls} and DIABLO (\code{block.plsda},
#' \code{block.splsda}) objects. See Details.
#' @param nrepeat numierc. Number of times the Cross-Validation process is repeated.
#' This is an important argument to ensure the estimation of the performance to
#' be as accurate as possible. Default it 1. 
#' @param dist only applies to an object inheriting from \code{"plsda"},
#' \code{"splsda"} or \code{"mint.splsda"} to evaluate the classification
#' performance of the model. Should be a subset of \code{"max.dist"},
#' \code{"centroids.dist"}, \code{"mahalanobis.dist"}. Default is \code{"all"}.
#' See \code{\link{predict}}.
#' @param auc if \code{TRUE} calculate the Area Under the Curve (AUC)
#' performance of the model.
#' @param progressBar by default set to \code{FALSE} to output the progress bar
#' of the computation.
#' @template arg/BPPARAM 
#' @param seed set a number here if you want the function to give reproducible outputs. 
#' Not recommended during exploratory analysis. Note if RNGseed is set in 'BPPARAM', this will be overwritten by 'seed'. 
#' Note 'seed' is not required or used in perf.mint.plsda as this method uses loo cross-validation
#' @param ... not used

#' @return For PLS and sPLS models:
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
#' \item{RSS}{Residual Sum of Squares across all selected features}
#' \item{PRESS}{Predicted Residual Error Sum of Squares across all selected features}
#' \item{cor.tpred, cor.upred}{Correlation between the 
#' predicted and actual components for X (t) and Y (u)} 
#' \item{RSS.tpred, RSS.upred}{Residual Sum of Squares between the
#' predicted and actual components for X (t) and Y (u)}
#' 
#' 
#'
#' For PLS-DA and sPLS-DA models:
#' \item{error.rate}{Prediction error rate for each dist and measure}
#' \item{auc}{AUC value averaged over the \code{nrepeat}}
#' \item{auc.all}{AUC values per repeat}
#' \item{predict}{Predicted values of each sample for each class}
#' \item{class}{A list which gives the predicted class of each sample for each dist and each of the ncomp components}
#' 
#' For mint.splsda models:
#' \item{study.specific.error}{A list that gives BER, overall error rate and
#' error rate per class, for each study} 
#' \item{global.error}{A list that gives
#' BER, overall error rate and error rate per class for all samples}
#' \item{predict}{A list of length \code{ncomp} that produces the predicted
#' values of each sample for each class} 
#' \item{class}{A list which gives the
#' predicted class of each sample for each \code{dist}.}
#' \item{auc}{AUC values} \item{auc.study}{AUC values for each study in mint models}
#' 
#' For sgccda models (i.e. block (s)PLS-DA models):
#' \item{error.rate}{Prediction error rate for each block of \code{object$X}
#' and each \code{dist}} 
#' \item{error.rate.per.class}{Prediction error rate for
#' each block of \code{object$X}, each \code{dist} and each class}
#' \item{predict}{Predicted values of each sample for each class and each block} 
#' \item{class}{Predicted class of each sample for each
#' block, each \code{dist}, and each nrepeat} 
#' \item{AveragedPredict.class}{if more than one block, returns
#' the average predicted class over the blocks (averaged of the \code{Predict}
#' output and prediction using the \code{max.dist} distance)}
#' \item{AveragedPredict.error.rate}{if more than one block, returns the
#' average predicted error rate over the blocks (using the
#' \code{AveragedPredict.class} output)} 
#' \item{WeightedPredict.class}{if more than one block, returns the weighted predicted class over the blocks
#' (weighted average of the \code{Predict} output and prediction using the
#' \code{max.dist} distance). See details for more info on weights.}
#' \item{WeightedPredict.error.rate}{if more than one block, returns the
#' weighted average predicted error rate over the blocks (using the
#' \code{WeightedPredict.class} output.)} 
#' \item{MajorityVote}{if more than one block, returns the majority class over the blocks. NA for a sample means that
#' there is no consensus on the predicted class for this particular sample over
#' the blocks.} 
#' \item{MajorityVote.error.rate}{if more than one block, returns the error rate of the \code{MajorityVote} output}
#' \item{WeightedVote}{if more than one block, returns the weighted majority
#' class over the blocks. NA for a sample means that there is no consensus on
#' the predicted class for this particular sample over the blocks.}
#' \item{WeightedVote.error.rate}{if more than one block, returns the error
#' rate of the \code{WeightedVote} output} 
#' \item{weights}{Returns the weights of each block used for the weighted predictions, for each nrepeat and each
#' fold} 

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
#' @example ./examples/perf.assess-examples.R
## ------------------------------- Generic -------------------------------- ##
perf.assess <- function(object, ...)
    UseMethod("perf.assess")
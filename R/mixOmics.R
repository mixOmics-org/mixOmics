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

# ========================================================================================================
# mixOmics: perform one of the package's function depending on the input data (list or matrix, vector or categerical data, etc)
# ========================================================================================================

#' PLS-derived methods: one function to rule them all!
#' 
#' \emph{This is the documentation for mixOmics function from mixOmics package.
#' For package documentation refer to \code{help(package='mixOmics')}}
#' 
#' This function performs one of the PLS derived methods included in the
#' mixOmics package that is the most appropriate for your input data, one of
#' (mint).(block).(s)pls(da) depending on your input data (single data, list of
#' data, discrete outcome, \dots{})
#' 
#' This function performs one of the PLS derived methods included in the
#' mixOmics package that is the most appropriate for your input data, one of
#' (mint).(block).(s)pls(da).
#' 
#' If your input data \code{X} is a matrix, then the algorithm is directed
#' towards one of (mint).(s)pls(da) depending on your input data \code{Y}
#' (factor for the discrete outcome directs the algorithm to DA analysis) and
#' whether you input a \code{study} parameter (MINT analysis) or a \code{keepX}
#' parameter (sparse analysis).
#' 
#' If your input data \code{X} is a list of matrices, then the algorithm is
#' directed towards one of (mint).block.(s)pls(da) depending on your input data
#' \code{Y} (factor for the discrete outcome directs the algorithm to DA
#' analysis) and whether you input a \code{study} parameter (MINT analysis) or
#' a \code{keepX} parameter (sparse analysis).
#' 
#' More details about the PLS modes in \code{?pls}.
#' 
#' @param X Input data. Either a matrix or a list of data sets (called
#' 'blocks') matching on the same samples. Data should be arranged in samples x
#' variables, with samples order matching in all data sets.
#' @param Y Outcome. Either a numeric matrix of responses or a factor or a
#' class vector for the discrete outcome.
#' @param indY To supply if Y is missing, indicates the position of the outcome
#' in the list X
#' @param study grouping factor indicating which samples are from the same
#' study
#' @param ncomp If \code{X} is a data matrix, \code{ncomp} is a single value.
#' If \code{X} is a list of data sets, \code{ncomp} is a numeric vector of
#' length the number of blocks in \code{X}. The number of components to include
#' in the model for each block (does not necessarily need to take the same
#' value for each block).
#' @param keepX Number of variables to keep in the \eqn{X}-loadings
#' @param keepY Number of variables to keep in the \eqn{Y}-loadings
#' @param design numeric matrix of size (number of blocks) x (number of blocks)
#' with only 0 or 1 values. A value of 1 (0) indicates a relationship (no
#' relationship) between the blocks to be modelled. If \code{Y} is provided
#' instead of \code{indY}, the \code{design} matrix is changed to include
#' relationships to \code{Y}.
#' @param tau numeric vector of length the number of blocks in \code{X}. Each
#' regularization parameter will be applied on each block and takes the value
#' between 0 (no regularisation) and 1. If tau = "optimal" the shrinkage
#' paramaters are estimated for each block and each dimension using the Schafer
#' and Strimmer (2005) analytical formula.
#' @param mode character string. What type of algorithm to use, (partially)
#' matching one of \code{"regression"}, \code{"canonical"}, \code{"invariant"}
#' or \code{"classic"}. See Details.
#' @param scale Logical. If scale = TRUE, each block is standardized to zero
#' means and unit variances (default: TRUE).
#' @param tol Convergence stopping value.
#' @param max.iter integer, the maximum number of iterations.
#' @param near.zero.var Logical, see the internal \code{\link{nearZeroVar}}
#' function (should be set to TRUE in particular for data with many zero
#' values). Setting this argument to FALSE (when appropriate) will speed up the
#' computations. Default value is FALSE
#' @return none
#' @author Florian Rohart, Kim-Anh Lê Cao, Al J Abadi
#' @seealso \code{\link{pls}}, \code{\link{spls}}, \code{\link{plsda}},
#' \code{\link{splsda}}, \code{\link{mint.pls}}, \code{\link{mint.spls}},
#' \code{\link{mint.plsda}}, \code{\link{mint.splsda}},
#' \code{\link{block.pls}}, \code{\link{block.spls}},
#' \code{\link{block.plsda}}, \code{\link{block.splsda}},
#' \code{\link{mint.block.pls}}, \code{\link{mint.block.spls}},
#' \code{\link{mint.block.plsda}}, \code{\link{mint.block.splsda}}
#' @references
#' 
#' mixOmics article:
#' 
#' Rohart F, Gautier B, Singh A, Lê Cao K-A. mixOmics: an R package for 'omics
#' feature selection and multiple data integration. PLoS Comput Biol 13(11):
#' e1005752
#' 
#' MINT models:
#' 
#' Rohart F, Eslami A, Matigian, N, Bougeard S, Lê Cao K-A (2017). MINT: A
#' multivariate integrative approach to identify a reproducible biomarker
#' signature across multiple experiments and platforms. BMC Bioinformatics
#' 18:128.
#' 
#' Eslami, A., Qannari, E. M., Kohler, A., and Bougeard, S. (2013). Multi-group
#' PLS Regression: Application to Epidemiology. In New Perspectives in Partial
#' Least Squares and Related Methods, pages 243-255. Springer.
#' 
#' Integration of omics data sets:
#' 
#' Singh A, Gautier B, Shannon C, Vacher M, Rohart F, Tebbutt S, Lê Cao K-A.
#' DIABLO: an integrative, multi-omics, multivariate method for multi-group
#' classification. \url{https://doi.org/10.1101/067611}
#' 
#' Lê Cao, K.-A., Martin, P.G.P., Robert-Granie, C. and Besse, P. (2009).
#' Sparse canonical methods for biological data integration: application to a
#' cross-platform study. \emph{BMC Bioinformatics} \bold{10}:34.
#' 
#' Lê Cao, K.-A., Rossouw, D., Robert-Granie, C. and Besse, P. (2008). A sparse
#' PLS for variable selection when integrating Omics data. \emph{Statistical
#' Applications in Genetics and Molecular Biology} \bold{7}, article 35.
#' 
#' Tenenhaus A., Phillipe C., Guillemot V., Lê Cao K-A. , Grill J. , Frouin V.
#' (2014), Variable selection for generalized canonical correlation analysis,
#' \emph{Biostatistics}, doi: 10.1093/biostatistics. PMID: 24550197.
#' 
#' Sparse SVD:
#' 
#' Shen, H. and Huang, J. Z. (2008). Sparse principal component analysis via
#' regularized low rank matrix approximation. \emph{Journal of Multivariate
#' Analysis} \bold{99}, 1015-1034.
#' 
#' PLS-DA:
#' 
#' Lê Cao K-A, Boitard S and Besse P (2011). Sparse PLS Discriminant Analysis:
#' biologically relevant feature selection and graphical displays for
#' multiclass problems. BMC Bioinformatics 12:253.
#' 
#' PLS:
#' 
#' Tenenhaus, M. (1998). \emph{La regression PLS: theorie et pratique}. Paris:
#' Editions Technic.
#' 
#' Wold H. (1966). Estimation of principal components and related models by
#' iterative least squares. In: Krishnaiah, P. R. (editors), \emph{Multivariate
#' Analysis}. Academic Press, N.Y., 391-420.
#' 
#' Abdi H (2010). Partial least squares regression and projection on latent
#' structure regression (PLS Regression). \emph{Wiley Interdisciplinary
#' Reviews: Computational Statistics}, 2(1), 97-106.
#' 
#' On multilevel analysis:
#' 
#' Liquet, B., Lê Cao, K.-A., Hocini, H. and Thiebaut, R. (2012) A novel
#' approach for biomarker selection and the integration of repeated measures
#' experiments from two platforms. \emph{BMC Bioinformatics} \bold{13}:325.
#' 
#' Westerhuis, J. A., van Velzen, E. J., Hoefsloot, H. C., and Smilde, A. K.
#' (2010). Multivariate paired data analysis: multilevel PLSDA versus OPLSDA.
#' \emph{Metabolomics}, \bold{6}(1), 119-128.
#' 
#' Visualisations:
#' 
#' González I., Lê Cao K.-A., Davis, M.J. and Déjean S. (2012) Visualising
#' associations between paired 'omics data sets.
#' BioData Mining 5:19.
#' @keywords multivariate hplot dplot
#' @export
#' @example ./examples/mixOmics-examples.R
mixOmics <- function(X,
                     Y,
                     indY, #only use if Y not provided
                     study, #mint
                     ncomp,
                     keepX, #sparse
                     keepY, #sparse
                     design, #block
                     tau = NULL, # rgcca, number between 0,1 or "optimal"
                     mode = c("regression", "canonical", "invariant", "classic"), #for DA only regression
                     scale,
                     tol =  1e-06,
                     max.iter = 100,
                     near.zero.var = FALSE)

{
    mode <- match.arg(mode)
    .check_regression_mode <- function(mode) {
        if (mode != 'regression') {
            warning("using `mode = 'regression'` for discriminant analysis\n")
        }
        return(NULL)
    }
    
    if (is.list(X) & !is.data.frame(X))# either rgcca, sgcca,sgcca-DA, mint.block, mint.block-DA
    {
        
        #need to check if Y or indY is a factor, unmap it and then do the checks (no other factors etc)
        if ((missing(indY)& missing(Y)) & is.null(tau))
            stop("Either 'Y', 'indY' or 'tau' is needed")
        
        if (is.null(tau)) # SGCCA/mint
        {
            
            isfactorY = FALSE
            
            
            if (!missing(Y))
            {
                if (is.list(Y) & !is.data.frame(X)) stop("Y must be a matrix or a factor")
                
                if (is.factor(Y)) {
                    #Y = as.factor(Y)
                    isfactorY = TRUE
                }
                
            }else if (!missing(indY)) {
                temp = X[[indY]] #not called Y to not be an input of the wrappers
                if (is.factor(temp)) {
                    #temp = as.factor(temp)
                    isfactorY = TRUE
                }
            }else if (missing(indY)) {
                stop("Either 'Y' or 'indY' is needed")
                
            }
            
            
            
            
            if (isfactorY)# either block.plsda/block.splsda/mint.block.plsda/mint.block.splsda
            {
                .check_regression_mode(mode)
                if (missing(keepX))
                {
                    if (missing(study)) #block.plsda
                    {
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a block Partial Least Squares - Discriminant Analysis is being performed (block.PLS-DA)")
                        res = block.plsda(X = X, Y = Y, indY = indY, ncomp = ncomp,design = design,
                                          scale = scale, tol = tol, max.iter = max.iter,near.zero.var = near.zero.var)
                        
                    } else {# mint.block.plsda
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a mint block Partial Least Squares - Discriminant Analysis is being performed (mint.block.PLS-DA)")
                        res = mint.block.plsda(X = X, Y = Y, indY = indY,study = study, ncomp = ncomp,design = design,
                                               scale = scale, tol = tol, max.iter = max.iter,near.zero.var = near.zero.var)
                    }
                    
                    
                } else {
                    if (missing(study))# block.splsda
                    {
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a block sparse Partial Least Squares - Discriminant Analysis is being performed (block.sPLS-DA)")
                        
                        res = mint.block.splsda(X = X, Y = Y, indY = indY, ncomp = ncomp,keepX = keepX,
                                                design = design, scale = scale, tol = tol,
                                                max.iter = max.iter,near.zero.var = near.zero.var)
                        
                        
                    } else {# mint.block.splsda
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a mint block sparse Partial Least Squares - Discriminant Analysis is being performed (mint.block.sPLS-DA)")
                        res = mint.block.splsda(X = X, Y = Y, indY = indY, ncomp = ncomp,study = study,keepX = keepX,
                                                design = design, scale = scale, tol = tol,
                                                max.iter = max.iter,near.zero.var = near.zero.var)
                    }
                    
                }
                
            } else { # either block.pls/block.spls/mint.block.pls/mint.block.spls
                
                
                if (missing(keepX) )
                {
                    if (missing(study)) #block.pls
                    {
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a block Partial Least Squares is being performed (block.PLS)")
                        res = block.pls(X = X, Y = Y, indY = indY, ncomp = ncomp,design = design,
                                        mode = mode,scale = scale, tol = tol, max.iter = max.iter,near.zero.var = near.zero.var)
                        
                    } else {# mint.block.pls
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a mint block Partial Least Squares is being performed (mint.block.PLS)")
                        res = mint.block.pls(X = X, Y = Y, indY = indY,study = study, ncomp = ncomp,design = design,
                                             mode = mode,scale = scale, tol = tol, max.iter = max.iter,near.zero.var = near.zero.var)
                    }
                    
                    
                } else {
                    if (missing(study))# block.spls
                    {
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a block sparse Partial Least Squares is being performed (block.sPLS)")
                        res = block.spls(X = X, Y = Y, indY = indY, ncomp = ncomp,keepX = keepX,
                                         design = design, mode =  mode,scale = scale, tol = tol,
                                         max.iter = max.iter,near.zero.var = near.zero.var)
                        
                        
                    } else {# mint.block.spls
                        if (missing(scale))
                            scale = FALSE
                        
                        message("a mint block sparse Partial Least Squares is being performed (mint.block.sPLS)")
                        res = mint.block.spls(X = X, Y = Y, indY = indY, ncomp = ncomp,study = study,keepX = keepX,
                                              design = design, mode =  mode,scale = scale, tol = tol,
                                              max.iter = max.iter,near.zero.var = near.zero.var)
                        
                    }
                    
                }
                
            }
            
        } else { # RGCCA
            
            
            if (!missing(study)) {message("'study' is not used")}
            
            if (missing(keepX) ) #RGCCA
            {
                message("A RGCCA analysis is being performed")
                if (missing(scale))
                    scale = FALSE
                
                res = wrapper.rgcca(X = X,design = design,tau = tau,ncomp = ncomp,
                                    max.iter = max.iter, scale = scale, tol = tol)
                
            } else { #sparse RGCCA
                if (missing(scale))
                    scale = FALSE
                
                message("A sparse RGCCA analysis is being performed")
                res = wrapper.rgcca(X = X,design = design,tau = tau,ncomp = ncomp,keepX = keepX,
                                    max.iter = max.iter, scale = scale, tol = tol)
                
                
            }
        }
        
        
        
        #end if (is.list(X))
    } else {#either pls,spls, plsda, splsda or mint. pls/spls/plsda/splsda
        if (missing(Y))
            stop("Y is missing")
        if (is.list(Y) & !is.data.frame(X))
            stop("Y must be a matrix or a factor")
        
        if (missing(mode)) mode = "regression"
        #check for unused inputs (scheme, etc etc) - now hardcoded scheme and init so dont check for these
        # if (!is.null(tau) | !missing(design) | !missing(init) | !missing(scheme))
        if (!is.null(tau) | !missing(design))
        {
            if (!is.null(tau))
                message("'tau' is not used")
            if (!missing(design))
                message("'design' is not used")
            # if (!missing(init))
            #     message("'init' is not used")
            # if (!missing(scheme))
            #     message("'scheme' is not used")
            
            stop("unused input parameters")
        }
        
        
        if (is.factor(Y))#either plsda, splsda
        {
            .check_regression_mode(mode)
            #Check.entry.pls.single(X, ncomp, keepX) # to have the warnings relative to X and Y, instead of blocks
            if (length(Y)!=nrow(X))
                stop("unequal number of rows in 'X' and 'Y'.")
            
            if (missing(keepX) & missing(keepY))  #plsda, mint.plsda
            {
                if (missing(study))
                {
                    if (missing(scale))
                        scale = TRUE
                    
                    message("a Partial Least Squares - Discriminant Analysis is being performed (PLS-DA)")
                    res = plsda(X = X, Y = Y, ncomp = ncomp, 
                                          max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                } else {# mint
                    if (missing(scale))
                        scale = FALSE
                    
                    message("a mint Partial Least Squares - Discriminant Analysis is being performed (mint.PLS-DA)")
                    res = mint.plsda(X = X, Y = Y, ncomp = ncomp, study = study,
                                     max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
                
            } else {#splsda, mint.splsda
                if (missing(study))
                {
                    if (missing(scale))
                        scale = TRUE
                    
                    message("a sparse Partial Least Squares - Discriminant Analysis is being performed (sPLS-DA)")
                    res = splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                } else {# mint
                    if (missing(scale))
                        scale = FALSE
                    
                    message("a mint sparse Partial Least Squares - Discriminant Analysis is being performed (mint.sPLS-DA)")
                    res = mint.splsda(X = X, Y = Y, ncomp = ncomp, study = study,keepX = keepX,
                                      max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            }
            
        } else { #pls or spls
            
            
            #Check.entry.pls(X, Y, ncomp, keepX, keepY) # to have the warnings relative to X and Y, instead of blocks
            
            if (missing(keepX) & missing(keepY))  #pls, mint.pls
            {
                if (missing(study))
                {
                    if (missing(scale))
                        scale = TRUE
                    
                    message("a Partial Least Squares is being performed (PLS)")
                    res = pls(X = X, Y = Y, ncomp = ncomp, mode = mode,
                                        max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                    
                } else { # mint
                    if (missing(scale))
                        scale = FALSE
                    
                    message("a mint Partial Least Squares is being performed (mint.PLS)")
                    res = mint.pls(X = X, Y = Y, ncomp = ncomp, mode = mode, study = study,
                                   max.iter = max.iter, tol = tol, near.zero.var = near.zero.var,scale = scale)
                }
                
            } else {
                if (missing(study))
                {
                    if (missing(scale))
                        scale = TRUE
                    
                    message("a sparse Partial Least Squares is being performed (sPLS)")
                    res = spls(X = X, Y = Y, ncomp = ncomp, mode = mode, keepX = keepX,keepY = keepY,
                                         max.iter = max.iter, tol = tol,
                                         near.zero.var = near.zero.var,scale = scale)
                } else {
                    if (missing(scale))
                        scale = FALSE
                    
                    message("a mint sparse Partial Least Squares is being performed (mint.sPLS)")
                    res = mint.spls(X = X, Y = Y, ncomp = ncomp, mode = mode, study = study,keepX = keepX,keepY = keepY,
                                    max.iter = max.iter, tol = tol,
                                    near.zero.var = near.zero.var,scale = scale)
                    
                }
            }
            
            
        }
    }
    cl = match.call()
    res$call = cl
    class(res) = c("mixOmics",class(res))
    return(invisible(res))
    
    
}

#############################################################################################################
# Author :
#   Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, The University of Queensland, Brisbane, QLD
#   Al J Abadi, Melbourne Integartive Genomics, The University of Melbourne, VIC
#
# created:  27-07-2017
# last modified: 2019
#
# Copyright (C) 2019
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
#############################################################################################################

#' Test whether adding a component improves the results
#'
#' determines the optimum number of components based on significance of
#' improvement in error rates
#'
#' @param mat.error.rate matrix of error rates, each column corresponding to
#' a component's classification error rate, in increasing order for components.
#' @param alpha significance threshold for t test. By default 0.01.
#' @param alternative Character, passed to stats::t.test
#'
#' @return integer, the optimal number of components
#'
#' @examples 
#' t.test.process(data.frame(comp1=100:104, comp2=20:24, comp3=10:14))
#' #> 3
#' t.test.process(data.frame(comp1=100:114, comp2=20:24, comp3=50:54))
#' #> 2
#' t.test.process(data.frame(comp1=10:14, comp2=20:24, comp3=50:54))
#' #> 1
#' @noRd
t.test.process <- function(mat.error.rate, alpha = 0.01, alternative = "greater")
{
    if (nrow(mat.error.rate) < 3)
    {
        return(NA)
    }
    ## ----- helper function to calculate pvalues for two columns of a data.frame:
    .calc_pval <- function(df, col1, col2) {
        x <- df[, col1]
        y <- df[, col2]
        pval <- tryCatch({t.test(x, y, alternative = alternative)$p.value},
                         error = function(e) e)
        
        if (!is.numeric(pval)) {
            ## if error rates constant - look at means only (super rare)
            if (is(pval, "error") && grepl("data are essentially constant", x = pval$message)) {
                ## significant if there's at least 5% improvement in error rates
                if (mean(y)/(mean(x) + 0.01) < 0.95) {
                    pval <- 0
                } else {
                    pval <- 1
                }
            } else { ## if unexpected condition
                .unexpected_err(trying_to = "choose the optimum number of components")
            }
        } else if (is.nan(pval)) { ## error rates constant and the same, not significant
            pval <- 1
        }
        return(pval)
    }
    
    max_comp <- ncol(mat.error.rate) ## number max of components included
    pval <- 1
    opt_comp <- 1 ## initialise the first optimal number of components
    next_comp <- 2 # initialise the first candidate component to compare
    while (opt_comp < max_comp & next_comp <= max_comp)
    {
        ## t.test of "is adding X comp improves the overall results"
        pval <- .calc_pval(mat.error.rate, opt_comp, next_comp)
        if (pval < alpha) {
            opt_comp <- next_comp
        }
        next_comp <- next_comp + 1
    }
    
    return(opt_comp)
}

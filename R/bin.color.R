################################################################################
# Authors:
#   Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
#
# created: 2009
# last modified:
#
# Copyright (C) 2009
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################

bin.color <-
    function(mat, cutoff, breaks, col, symkey) 
    {
        if (isTRUE(symkey)) {
            max.mat = max(abs(mat))
            min.mat = -max.mat
        }
        else {
            max.mat = max(mat)
            min.mat = min(mat)
        }
        
        if (missing(breaks) || is.null(breaks)) {
            if (is(col,"function")) breaks = 32
            else breaks = length(col)
        }
        
        if (length(breaks) == 1) {
            if (isTRUE(symkey)) {
                if ((breaks/2) - trunc(breaks/2) != 0)
                    stop("'breaks' must be a even number if 'symkey = TRUE'",
                         call. = FALSE)
                
                if (cutoff == 0) {
                    breaks = c(seq(min.mat, max.mat, length = breaks + 1))
                } else {
                    nb = breaks/2
                    breaks = c(seq(min.mat, -cutoff, length = nb + 1), 0,
                               seq(cutoff, max.mat, length = nb + 1))
                    id = which(breaks == 0)
                    breaks = breaks[-c(id - 1, id + 1)]
                }
            } else {
                breaks = breaks + 1
                
                if ((min.mat < -cutoff) & (max.mat < cutoff))
                    breaks = seq(min.mat, -cutoff, length = breaks)
                
                if ((min.mat > -cutoff) & (max.mat > cutoff))
                    breaks = seq(cutoff, max.mat, length = breaks)
                
                if ((min.mat < -cutoff) & (max.mat > cutoff)) {
                    if (cutoff == 0){
                        breaks = c(seq(min.mat, max.mat, length = breaks))
                    } else {
                        long = max.mat - min.mat - 2*cutoff
                        bin = long/breaks
                        breaks = seq(cutoff, -min.mat, by = bin)
                        o = order(breaks, decreasing = TRUE)
                        breaks = c(-breaks[o], 0, seq(cutoff, max.mat, by = bin))
                        id = which(breaks == 0)
                        breaks = breaks[-c(id - 1, id + 1)]
                    }
                }
            }
        }
        
        ncol = length(breaks) - 1
        
        if (is(col, "function"))
            col = col(ncol)
        
        if (length(breaks) != length(col) + 1)
            stop("must have one more break than colour", call. = FALSE)
        
        min.breaks = min(breaks)
        max.breaks = max(breaks)
        
        mat[mat < min.breaks] = min.breaks
        mat[mat > max.breaks] = max.breaks
        
        bin = .bincode(as.double(mat), as.double(breaks), TRUE, TRUE)
        
        return(invisible(list(bin = bin, col = col, breaks = breaks, lim =
                                  c(min.mat, max.mat))))
    }

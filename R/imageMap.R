###############################################################################
# Authors:
#   Ignacio Gonzalez,
#   Francois Bartolo,
#   Kim-Anh Le Cao,
#
# created: 2009
# last modified: 2015
#
# Copyright (C) 2009
#
# This file is part of mixOmics, licensed under AGPL (>= 3).
# SPDX-License-Identifier: AGPL-3.0-or-later
################################################################################


# --------------------------------------------
# imageMap
# --------------------------------------------


imageMap <-
    function (mat,
              color,
              row.names,
              col.names,
              row.sideColors,
              col.sideColors,
              row.cex = NULL,
              col.cex = NULL,
              cluster,
              ddr,
              ddc,
              cut.tree = c(0, 0),
              transpose = FALSE,
              symkey = TRUE,
              keysize = c(1, 1),
              keysize.label,
              zoom = FALSE,
              title = NULL,
              xlab = NULL,
              ylab = NULL,
              margins = c(5, 5),
              lhei = NULL,
              lwid = NULL)
    {
        #-- image map -------------------------------------------------------------#
        #----------
        
        if (isTRUE(symkey)) {
            max.mat = max(abs(mat), na.rm = TRUE)
            min.mat = -max.mat
        }
        else {
            max.mat = max(mat, na.rm = TRUE)
            min.mat = min(mat, na.rm = TRUE)
        }
        
        if (isTRUE(transpose)) {
            mat = t(mat)
            
            temp = col.sideColors
            col.sideColors = row.sideColors
            row.sideColors = temp
            
            temp = col.names
            col.names = row.names
            row.names = temp
            
            if (cluster == "both")
            {temp = ddc
            ddc = ddr
            ddr = temp}
            else if (cluster == "column")
                ddr=ddc
            else if (cluster == "row")
                ddc=ddr
            
            lhei = NULL
            lwid = NULL
        }
        
        nr = nrow(mat)
        nc = ncol(mat)
        
        #-- row.cex and col.cex
        if (is.null(row.cex)) row.cex = min(1, 0.2 + 1/log10(nr))
        if (is.null(col.cex)) col.cex = min(1, 0.2 + 1/log10(nc))
        
        #-- breaks
        breaks = length(color) + 1
        breaks = seq(min.mat, max.mat, length = breaks)
        
        nbr = length(breaks)
        ncol = nbr - 1
        
        min.breaks = min(breaks)
        max.breaks = max(breaks)
        mat[mat < min.breaks] = min.breaks
        mat[mat > max.breaks] = max.breaks
        mat = t(mat)
        #-- layout matrix
        lmat = matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE)
        csc = rsc = FALSE
        
        if (!is.null(col.sideColors)) {
            lmat = rbind(lmat[1, ], c(NA, 3), lmat[2, ] + 1)
            n.csc = ncol(col.sideColors)
            if (is.null(lhei)) {
                lhei = c(keysize[2], 0.15 + 0.1 * (n.csc - 1), 4)
            }
            csc = TRUE
        }
        
        if (!is.null(row.sideColors)) {
            lmat = cbind(lmat[, 1], c(rep(NA, nrow(lmat) - 1), nrow(lmat) + 2),
                         lmat[, 2] + c(rep(0, nrow(lmat) - 1), 1))
            if (is.null(lwid)) {
                n.rsc = ncol(row.sideColors)
                lwid = c(keysize[2], 0.15 + 0.1 * (n.rsc - 1), 4)
            }
            rsc = TRUE
        }
        
        lmat[is.na(lmat)] = 0
        
        if (is.null(lhei))
            lhei = c(keysize[2], 4)
        
        if (is.null(lwid))
            lwid = c(keysize[1], 4)
        
        if (isTRUE(zoom)) {
            graphics.off()
            dev.new(pos=-1)#x11(xpos = -1)
        }
        
        op = par(no.readonly = TRUE)
        on.exit(par(op))
        layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
        
        #-- layout 1 --#
        par(mar = c(5, 2, 2, 1), cex = 0.75)
        
        z = seq(0, 1, length = length(color))
        z = matrix(z, ncol = 1)
        image(z, col = color, xaxt = "n", yaxt = "n")
        box()
        par(usr = c(0, 1, 0, 1))
        lv = c(min.breaks, (3*min.breaks + max.breaks)/4,
               (min.breaks + max.breaks)/2,
               (3*max.breaks + min.breaks)/4, max.breaks)
        xv = (as.numeric(lv) - min.mat) / (max.mat - min.mat)
        axis(1, at = xv, labels = round(lv, 2), cex.axis = keysize.label)
        title("Color key", font.main = 1, cex.main = keysize.label)
        
        #-- layout 2 --#
        par(mar = c(ifelse(cut.tree[2] != 0, 0.5, 0), 0,
                    ifelse(!is.null(title), 5, 0), margins[2]))
        if ((cluster == "both") || (!transpose && cluster == "column") ||
            (transpose && cluster == "row" )) {
            h = attr(ddc, "height")
            plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none",
                 ylim = c(cut.tree[2] * h, h))
        }
        else {
            plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")
        }
        
        if (!is.null(title))
            title(title, cex.main = 1.5 * op[["cex.main"]])
        
        #-- layout 3 --#
        if (isTRUE(csc)) {
            par(mar = c(0.5, 0, 0, margins[2]))
            sideColors = as.vector(col.sideColors)
            img = matrix(c(1:(n.csc * nc)), ncol = n.csc, byrow = FALSE)
            
            image(1:nc, 1:n.csc, img, col = sideColors, axes = FALSE, xlab = "",
                  ylab = "")
            abline(h = 1:(n.csc - 1) + 0.5, lwd = 2,
                   col = ifelse(par("bg") == "transparent", "white", par("bg")))
        }
        
        #-- layout 4 --#
        par(mar = c(margins[1], 0, 0, ifelse(cut.tree[1] != 0, 0.5, 0)))
        if ((cluster == "both") || (cluster == "row" & !transpose) ||
            (cluster == "column" & transpose)) {
            h = attr(ddr, "height")
            plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none",
                 xlim = c(h, cut.tree[1] * h))
        }
        else {
            plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")
        }
        
        #-- layout 5 --#
        if (isTRUE(rsc)) {
            par(mar = c(margins[1], 0, 0, 0.5))
            n.rsc = ncol(row.sideColors)
            r.sideColors = row.sideColors[, n.rsc:1]
            sideColors = as.vector(r.sideColors)
            img = matrix(1:(n.rsc * nr), nrow = n.rsc, byrow = TRUE)
            
            image(1:n.rsc, 1:nr, img, col = sideColors, axes = FALSE, xlab = "",
                  ylab = "")
            abline(v = 1:(n.rsc - 1) + 0.5, lwd = 2,
                   col = ifelse(par("bg") == "transparent", "white", par("bg")))
        }
        
        #-- layout 6 --#
        par(mar = c(margins[1], 0, 0, margins[2]))
        image(1:nc, 1:nr, mat, xlim = 0.5 + c(0, nc), ylim = 0.5 +
                  c(0, nr), axes = FALSE, xlab = "", ylab = "", col = color,
              breaks = breaks)
        axis(1, 1:nc, labels = col.names, las = 2, line = -0.5, tick = 0,
             cex.axis = col.cex)
        
        if (!is.null(xlab))
            mtext(xlab, side = 1, line = margins[1] - 1.25)
        
        axis(4, 1:nr, labels = row.names, las = 2, line = -0.5, tick = 0,
             cex.axis = row.cex)
        
        if (!is.null(ylab))
            mtext(ylab, side = 4, line = margins[2] - 1.25)
        
        #-- ZOOM --#
        #----------#
        flag1 = flag2 = FALSE
        
        if (isTRUE(zoom)) {
            nD = dev.cur()
            zone = FALSE
            
            repeat {
                dev.set(nD)
                
                repeat {
                    loc = locator(1, type = "n")
                    
                    if (is.null(loc) & zone == TRUE) break
                    if (is.null(loc) & !isTRUE(flag1)) break
                    flag1 = TRUE
                    
                    x1 = round(loc[[1]] - 0.5) + 0.5
                    y1 = round(loc[[2]] - 0.5) + 0.5
                    
                    if (!(x1 < 0 | x1 > nc + 0.5 | y1 < 0 | y1 > nr + 0.5)) break
                }
                
                if (is.null(loc) & zone == TRUE) break
                
                if (!is.null(loc) & zone == TRUE) {
                    rect(xleft.old, ybottom.old, xright.old, ytop.old,
                         border = "white")
                    points(x1.old, y1.old, type = "p", pch = 3,
                           cex = 2, col = "white")
                }
                
                if (is.null(loc) & zone == FALSE) {
                    break
                }
                else {
                    x1.old = x1
                    y1.old = y1
                    
                    points(x1, y1, type = "p", pch = 3, cex = 2)
                    
                    repeat {
                        loc = locator(1, type = "n")
                        
                        if (is.null(loc) & zone == TRUE) break
                        if (is.null(loc) & !isTRUE(flag2)) break
                        flag2 = TRUE
                        
                        x2 = round(loc[[1]] - 0.5) + 0.5
                        y2 = round(loc[[2]] - 0.5) + 0.5
                        
                        if (!(x2 < 0 | x2 > nc + 0.5 | y2 < 0 | y2 > nr + 0.5)) {
                            zone = TRUE; break }
                    }
                    
                    if (is.null(loc) & zone == TRUE) break
                    if (is.null(loc) & !isTRUE(flag2)) break
                    
                    xleft.old = min(x1, x2)
                    xright.old = max(x1, x2)
                    ybottom.old = min(y1, y2)
                    ytop.old = max(y1, y2)
                    
                    rect(xleft.old, ybottom.old, xright.old, ytop.old)
                }
                
                dev.new()#x11()
                plot.par = par(no.readonly = TRUE)
                
                if (isTRUE(zone)) {
                    xleft = xleft.old + 0.5
                    ybottom = ybottom.old + 0.5
                    xright = xright.old - 0.5
                    ytop = ytop.old - 0.5
                    nr.zoom = length(xleft:xright)
                    nc.zoom = length(ybottom:ytop)
                    mat.zoom = matrix(mat[xleft:xright, ybottom:ytop],
                                      nrow = nr.zoom, ncol = nc.zoom)
                    rlab.zoom = col.names[xleft:xright]
                    clab.zoom = row.names[ybottom:ytop]
                    r.cex = min(1.2, 0.2 + 1/log10(nr.zoom))
                    c.cex = min(1.2, 0.2 + 1/log10(nc.zoom))
                    
                    layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
                    
                    # layout 1
                    par(mar = c(5, 2, 2, 1), cex = 0.75)
                    image(z, col = color, xaxt = "n", yaxt = "n")
                    box()
                    par(usr = c(0, 1, 0, 1))
                    lv = c(min.breaks, (3*min.breaks + max.breaks)/4,
                           (min.breaks + max.breaks)/2,
                           (3*max.breaks + min.breaks)/4, max.breaks)
                    xv = (as.numeric(lv) - min.mat) / (max.mat - min.mat)
                    axis(1, at = xv, labels = round(lv, 2))
                    title("Color key", font.main = 1)
                    
                    #-- layout 2 --#
                    par(mar = c(ifelse(cut.tree[2] != 0, 0.5, 0), 0,
                                ifelse(!is.null(title), 5, 0), margins[2]))
                    if ((cluster == "both") || (cluster == "column" & !transpose) ||
                        (cluster == "row" & transpose)) {
                        plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none",
                             xlim = c(xleft - 0.5, xright + 0.5))
                    }
                    else {
                        plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")
                    }
                    
                    if (!is.null(title))
                        title(title, cex.main = 1.5 * op[["cex.main"]])
                    
                    # layout 3
                    if (isTRUE(csc)) {
                        par(mar = c(0.5, 0, 0, margins[2]))
                        sideColors = as.vector(col.sideColors[xleft:xright, ])
                        img = matrix(c(1:(n.csc * nr.zoom)), ncol = n.csc,
                                     byrow = FALSE)
                        
                        image(1:nr.zoom, 1:n.csc, img, col = sideColors,
                              axes = FALSE, xlab = "", ylab = "")
                        abline(h = 1:(n.csc - 1) + 0.5, lwd = 2,
                               col = ifelse(par("bg") == "transparent", "white",
                                            par("bg")))
                    }
                    
                    #-- layout 4 --#
                    par(mar = c(margins[1], 0, 0, ifelse(cut.tree[1] != 0, 0.5, 0)))
                    if ((cluster == "both") || (cluster == "row" & !transpose) ||
                        (cluster == "column" & transpose)) {
                        plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i",
                             leaflab = "none", ylim = c(ybottom - 0.5, ytop + 0.5))
                    }
                    else {
                        plot(0, 0, axes = FALSE, type = "n", xlab = "", ylab = "")
                    }
                    
                    # layout 5
                    if (isTRUE(rsc)) {
                        par(mar = c(margins[1], 0, 0, 0.5))
                        r.sideColors = row.sideColors[ybottom:ytop, n.rsc:1]
                        sideColors = as.vector(r.sideColors)
                        img = matrix(1:(n.rsc * nc.zoom), nrow = n.rsc,
                                     byrow = TRUE)
                        
                        image(1:n.rsc, 1:nc.zoom, img, col = sideColors,
                              axes = FALSE, xlab = "", ylab = "")
                        abline(v = 1:(n.rsc - 1) + 0.5, lwd = 2,
                               col = ifelse(par("bg") == "transparent", "white",
                                            par("bg")))
                    }
                    
                    # layout 6
                    par(mar = c(margins[1], 0, 0, margins[2]))
                    image(1:nr.zoom, 1:nc.zoom, mat.zoom, col = color,
                          breaks = breaks, axes = FALSE, xlab = "", ylab = "")
                    
                    axis(1, 1:nr.zoom, labels = rlab.zoom, las = 2, 
                         line = -0.5, tick = 0, cex.axis = r.cex)
                    
                    if (!is.null(xlab)) 
                        mtext(xlab, side = 1, line = margins[1] - 1.25)
                    
                    axis(4, 1:nc.zoom, labels = clab.zoom, las = 2, 
                         line = -0.5, tick = 0, cex.axis = c.cex)
                    
                    if (!is.null(ylab)) 
                        mtext(ylab, side = 4, line = margins[2] - 1.25)
                }
                
                par(plot.par)
            }
        }  
        par(op)
    }

################################################################################
# created: 2011
# last modified: 2011
#
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
# SPDX-License-Identifier: AGPL-3.0-only AND (GPL-2.0-only OR GPL-3.0-only)
#
# Code from other packages keeps its own licence. It is listed here and
# marked in the code with SPDX-SnippetBegin and SPDX-SnippetEnd:
#   - ica.def() and ica.par(): adapted from ica.R.def() and ica.R.par()
#     of the fastICA package (version 1.1-13) by J L Marchini, C Heaton
#     and B D Ripley, GPL-2 | GPL-3.
###############################################################################

# SPDX-SnippetBegin
# SPDX-SnippetCopyrightText: J L Marchini, C Heaton, B D Ripley (fastICA)
# SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
ica.def <-
    function (X, ncomp, tol, fun, alpha, max.iter, verbose, w.init)
    {
        n <- nrow(X)
        p <- ncol(X)
        W <- matrix(0, ncomp, ncomp)
        for (i in 1:ncomp) {
            w <- matrix(w.init[i,], ncomp, 1)
            if (i > 1) {
                t <- w
                t[1:length(t)] <- 0
                for (u in 1:(i - 1)) {
                    k <- sum(w * W[u, ])
                    t <- t + k * W[u, ]
                }
                w <- w - t
            }
            w <- w/sqrt(sum(w^2))
            lim <- rep(1000, max.iter)
            it <- 1
            if (fun == "logcosh") {
                while (lim[it] > tol && it < max.iter) {
                    wx <- t(w) %*% X
                    gwx <- tanh(alpha * wx)
                    gwx <- matrix(gwx, ncomp, p, byrow = TRUE)
                    xgwx <- X * gwx
                    v1 <- apply(xgwx, 1, FUN = mean)
                    g.wx <- alpha * (1 - (tanh(alpha * wx))^2)
                    v2 <- mean(g.wx) * w
                    w1 <- v1 - v2
                    w1 <- matrix(w1, ncomp, 1)
                    it <- it + 1
                    if (i > 1) {
                        t <- w1
                        t[1:length(t)] <- 0
                        for (u in 1:(i - 1)) {
                            k <- sum(w1 * W[u, ])
                            t <- t + k * W[u, ]
                        }
                        w1 <- w1 - t
                    }
                    w1 <- w1/sqrt(sum(w1^2))
                    lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
                    w <- matrix(w1, ncomp, 1)
                }
            }
            if (fun == "exp") {
                while (lim[it] > tol && it < max.iter) {
                    wx <- t(w) %*% X
                    gwx <- wx * exp(-(wx^2)/2)
                    gwx <- matrix(gwx, ncomp, p, byrow = TRUE)
                    xgwx <- X * gwx
                    v1 <- apply(xgwx, 1, FUN = mean)
                    g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
                    v2 <- mean(g.wx) * w
                    w1 <- v1 - v2
                    w1 <- matrix(w1, ncomp, 1)
                    it <- it + 1
                    if (i > 1) {
                        t <- w1
                        t[1:length(t)] <- 0
                        for (u in 1:(i - 1)) {
                            k <- sum(w1 * W[u, ])
                            t <- t + k * W[u, ]
                        }
                        w1 <- w1 - t
                    }
                    w1 <- w1/sqrt(sum(w1^2))
                    lim[it] <- Mod(Mod(sum((w1 * w))) - 1)
                    if (verbose)
                        message("Iteration ", it - 1, " tol = ", format(lim[it]))
                    w <- matrix(w1, ncomp, 1)
                }
            }
            W[i, ] <- w
        }
        return(W)
    }

ica.par <- function (X, ncomp, tol, fun, alpha, max.iter, verbose, w.init)
{
    Diag <- function(d) if(length(d) > 1L) diag(d) else as.matrix(d)
    n <- nrow(X)
    p <- ncol(X)
    W <- w.init+matrix(c(1:9),ncomp,ncomp)
    sW <- La.svd(W)
    W <- sW$u %*% Diag(1/sW$d) %*% t(sW$u) %*% W
    W1 <- W
    lim <- rep(1000, max.iter)
    it <- 1
    if (fun == "logcosh") {
        while (lim[it] > tol && it < max.iter) {
            wx <- W %*% X
            gwx <- tanh(alpha * wx)
            v1 <- gwx %*% t(X)/p
            g.wx <- alpha * (1 - (gwx)^2)
            v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
            W1 <- v1 - v2
            sW1 <- La.svd(W1)
            W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
            lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
            W <- W1
            if (verbose)
                message("Iteration ", it, " tol = ", format(lim[it + 1]))
            it <- it + 1
        }
    }
    if (fun == "exp") {
        while (lim[it] > tol && it < max.iter) {
            wx <- W %*% X
            gwx <- wx * exp(-(wx^2)/2)
            v1 <- gwx %*% t(X)/p
            g.wx <- (1 - wx^2) * exp(-(wx^2)/2)
            v2 <- Diag(apply(g.wx, 1, FUN = mean)) %*% W
            W1 <- v1 - v2
            sW1 <- La.svd(W1)
            W1 <- sW1$u %*% Diag(1/sW1$d) %*% t(sW1$u) %*% W1
            lim[it + 1] <- max(Mod(Mod(diag(W1 %*% t(W))) - 1))
            W <- W1
            if (verbose)
                message("Iteration ", it, " tol = ", format(lim[it + 1]))
            it <- it + 1
        }
    }
    return(W)
}
# SPDX-SnippetEnd

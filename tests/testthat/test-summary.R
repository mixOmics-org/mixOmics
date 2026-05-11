context("summary")

test_that("(summary:basic): pls", {
    data(linnerud)
    object <- pls(linnerud$exercise, linnerud$physiological, ncomp = 2)

    out <- summary(object)

    expect_s3_class(out, "summary")
    expect_equal(out$method, "pls")
    expect_equal(out$ncomp, object$ncomp)
    expect_equal(out$mode, object$mode)
    expect_equal(out$keep.var$X, colnames(object$X))
    expect_equal(out$keep.var$Y, colnames(object$Y))

    expect_equal(dim(out$Cm.X$own), c(ncol(object$X), object$ncomp))
    expect_equal(dim(out$Cm.X$opp), c(ncol(object$X), object$ncomp))
    expect_equal(dim(out$Cm.Y$own), c(ncol(object$Y), object$ncomp))
    expect_equal(dim(out$Cm.Y$opp), c(ncol(object$Y), object$ncomp))
    expect_equal(colnames(out$Cm.X$own), paste("comp", seq_len(object$ncomp)))

    expect_equal(dim(out$Rd.X$own), c(object$ncomp, 2))
    expect_equal(dim(out$Rd.Y$own), c(object$ncomp, 2))
    expect_equal(colnames(out$Rd.X$own), c("Proportion", "Cumulative"))
    expect_equal(colnames(out$Rd.Y$own), c("Proportion", "Cumulative"))

    expect_equal(out$VIP, vip(object))
    expect_match(paste(capture.output(print(out)), collapse = "\n"),
                 "PLS mode:")
})

test_that("(summary:parameter): spls, 'what' and 'keep.var'", {
    data(linnerud)
    object <- spls(linnerud$exercise, linnerud$physiological, ncomp = 2,
                   keepX = c(2, 2), keepY = c(2, 2))

    out <- summary(object, what = c("redundancy", "VIP"), keep.var = TRUE)

    selected.X <- apply(object$loadings$X != 0, 1, any)
    selected.Y <- apply(object$loadings$Y != 0, 1, any)

    expect_s3_class(out, "summary")
    expect_equal(out$method, "spls")
    expect_equal(out$keepX, object$keepX)
    expect_equal(out$keepY, object$keepY)
    expect_equal(out$keep.var$X, colnames(object$X)[selected.X])
    expect_equal(out$keep.var$Y, colnames(object$Y)[selected.Y])

    expect_null(out$Cm.X)
    expect_equal(dim(out$Rd.X$own), c(object$ncomp, 2))
    expect_equal(dim(out$VIP), c(sum(selected.X), object$ncomp))
    expect_equal(out$VIP, vip(object)[selected.X, ])
    expect_match(paste(capture.output(print(out)), collapse = "\n"),
                 "sPLS mode:")
})

test_that("(summary:basic): rcc", {
    data(linnerud)
    object <- rcc(linnerud$exercise, linnerud$physiological, ncomp = 2,
                  lambda1 = 0.1, lambda2 = 0.1)

    out <- summary(object)

    expect_s3_class(out, "summary")
    expect_equal(out$method, "rcc")
    expect_equal(out$ncomp, object$ncomp)
    expect_equal(out$cor, setNames(object$cor[seq_len(object$ncomp)],
                                   paste(seq_len(object$ncomp), "th", sep = "")))
    expect_equal(out$can.cor, out$cor)
    expect_null(out$cutoff)
    expect_equal(out$keep$X, colnames(object$X))
    expect_equal(out$keep$Y, colnames(object$Y))

    expect_equal(dim(out$Cm.X$own), c(ncol(object$X), object$ncomp))
    expect_equal(dim(out$Rd.X$own), c(object$ncomp, 2))
    expect_equal(rownames(out$Rd.X$own), paste("comp", seq_len(object$ncomp)))
    expect_match(paste(capture.output(print(out)), collapse = "\n"),
                 "Canonical correlations")

})

test_that("(summary:error): rcc cutoff selects no variables", {
    data(linnerud)
    object <- rcc(linnerud$exercise, linnerud$physiological, ncomp = 2,
                  lambda1 = 0.1, lambda2 = 0.1)

    expect_error(summary(object, cutoff = 2), "No variable was selected")
})

test_that("(summary:basic): pca", {
    data(linnerud)
    object <- pca(linnerud$exercise, ncomp = 2, scale = TRUE)

    out <- summary(object)

    expect_s3_class(out, "summary.prcomp")
    expect_equal(dim(out$importance), c(3, object$ncomp))
    expect_equal(rownames(out$importance),
                 c("Standard deviation",
                   "Proportion of Variance",
                   "Cumulative Proportion"))
    expect_equal(colnames(out$importance), colnames(object$rotation))
    expect_equal(out$importance["Standard deviation", ], object$sdev)
    expect_equal(out$importance["Proportion of Variance", ],
                 round(object$prop_expl_var$X, 5))
    expect_equal(out$importance["Cumulative Proportion", ],
                 round(cumsum(object$prop_expl_var$X), 5))
})

context("vip")

test_that("(vip:basic): pls, multivariate Y", {
    data(linnerud)
    res <- pls(linnerud$exercise, linnerud$physiological, ncomp = 2)

    v <- vip(res)

    expect_true(is.matrix(v))
    expect_equal(dim(v), c(ncol(linnerud$exercise), 2))
    expect_equal(rownames(v), colnames(linnerud$exercise))
    expect_equal(colnames(v), c("comp1", "comp2"))
    # VIP invariant: column-wise sum of squared VIPs equals p
    expect_equal(colSums(v^2),
                 setNames(rep(ncol(linnerud$exercise), 2), c("comp1", "comp2")))
})

test_that("(vip:data): pls, linnerud reference values", {
    data(linnerud)
    res <- pls(linnerud$exercise, linnerud$physiological, ncomp = 2)

    v <- vip(res)

    # Regression baseline — guards against silent changes to the VIP formula
    expected <- matrix(
        c(1.062280, 1.293793, 0.444592,
          0.904423, 1.139641, 0.939807),
        nrow = 3, ncol = 2,
        dimnames = list(c("Chins", "Situps", "Jumps"),
                        c("comp1", "comp2"))
    )
    expect_equal(v, expected, tolerance = 1e-5)
})

test_that("(vip:edge.case): pls, univariate Y (q == 1)", {
    data(linnerud)
    res <- pls(linnerud$exercise, linnerud$physiological[, 1, drop = FALSE],
               ncomp = 2)

    v <- vip(res)

    expect_equal(dim(v), c(ncol(linnerud$exercise), 2))
    expect_equal(colSums(v^2),
                 setNames(rep(ncol(linnerud$exercise), 2), c("comp1", "comp2")))
})

test_that("(vip:basic): spls", {
    data(liver.toxicity)
    res <- spls(liver.toxicity$gene, liver.toxicity$clinic,
                ncomp = 2, keepX = c(50, 50), keepY = c(10, 10))

    v <- vip(res)

    expect_equal(dim(v), c(ncol(liver.toxicity$gene), 2))
    expect_equal(rownames(v), colnames(liver.toxicity$gene))
    expect_equal(colnames(v), c("comp1", "comp2"))
    # sparse loadings are still unit-norm, so the invariant must hold
    expect_equal(colSums(v^2),
                 setNames(rep(ncol(liver.toxicity$gene), 2),
                          c("comp1", "comp2")))
})

test_that("(vip:basic): plsda", {
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- as.factor(breast.tumors$sample$treatment)
    res <- plsda(X, Y, ncomp = 2)

    v <- vip(res)

    expect_equal(dim(v), c(ncol(X), 2))
    expect_equal(colSums(v^2), setNames(rep(ncol(X), 2), c("comp1", "comp2")))
})

test_that("(vip:basic): splsda", {
    data(breast.tumors)
    X <- breast.tumors$gene.exp
    Y <- as.factor(breast.tumors$sample$treatment)
    res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))

    v <- vip(res)

    expect_equal(dim(v), c(ncol(X), 2))
    expect_equal(rownames(v), colnames(X))
    expect_equal(colnames(v), c("comp1", "comp2"))
    expect_equal(colSums(v^2), setNames(rep(ncol(X), 2), c("comp1", "comp2")))
})

test_that("(vip:edge.case): ncomp == 1", {
    data(linnerud)
    res <- pls(linnerud$exercise, linnerud$physiological, ncomp = 1)

    v <- vip(res)

    expect_equal(dim(v), c(ncol(linnerud$exercise), 1))
    expect_equal(colnames(v), "comp1")
    # comp1 is always sqrt(p) * |W[, 1]| regardless of Y dimensionality
    expect_equal(unname(v[, 1]),
                 unname(sqrt(ncol(linnerud$exercise)) * abs(res$loadings$X[, 1])))
})

test_that("(vip:error): unsupported object classes", {
    data(linnerud)
    res.pca <- pca(linnerud$exercise, ncomp = 2)
    expect_error(vip(res.pca),
                 "'vip' is only implemented for the following objects")
    expect_error(vip(list(a = 1)),
                 "'vip' is only implemented for the following objects")
})

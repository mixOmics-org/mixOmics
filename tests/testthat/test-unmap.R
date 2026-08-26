context("unmap and map")

test_that("(unmap:basic): factor input with default args", {
    data(nutrimouse)
    Y <- unmap(nutrimouse$diet)

    expect_true(is.matrix(Y))
    expect_equal(nrow(Y), length(nutrimouse$diet))
    expect_equal(ncol(Y), nlevels(nutrimouse$diet))
    expect_equal(unname(rowSums(Y)), rep(1, nrow(Y)))
    expect_equal(attr(Y, "levels"), levels(nutrimouse$diet))
})

test_that("(unmap:basic): character vector input", {
    Y <- unmap(c("a", "b", "a", "c"))

    expect_equal(dim(Y), c(4, 3))
    expect_equal(unname(rowSums(Y)), rep(1, 4))
    expect_null(attr(Y, "levels"))
})

test_that("(unmap:parameter): groups", {
    cls <- c("a", "b", "a")
    Y <- unmap(cls, groups = c("a", "b"))

    expect_equal(dim(Y), c(3, 2))
    expect_equal(unname(rowSums(Y)), rep(1, 3))
})

test_that("(unmap:error): groups incompatible with classification", {
    expect_error(
        unmap(c("a", "b", "c"), groups = c("a", "b")),
        "groups incompatible with classification"
    )
})

test_that("(unmap:parameter): noise places noise class last", {
    Y <- unmap(c("a", "b", "noise", "a"), noise = "noise")

    expect_equal(nrow(Y), 4)
    expect_equal(colnames(Y), c("a", "b", "noise"))
    expect_equal(unname(Y[, "noise"]), c(0, 0, 1, 0))
})

test_that("(unmap:edge.case): noise reorder preserves encoding", {
    Y <- unmap(c("a", "0", "b"), groups = c("0", "a", "b"), noise = "0")

    expect_equal(colnames(Y), c("a", "b", "noise"))
    expect_equal(unname(Y[, "a"]), c(1, 0, 0))
    expect_equal(unname(Y[, "b"]), c(0, 0, 1))
    expect_equal(unname(Y[, "noise"]), c(0, 1, 0))
})

test_that("(unmap:error): noise incompatible with groups", {
    expect_error(
        unmap(c("a", "b", "c"), noise = "missing_class"),
        "noise incompatible with classification"
    )
})

test_that("(map:basic): one-hot matrix", {
    Y <- diag(3)
    expect_equal(map(Y), 1:3)
})

test_that("(map:edge.case): ties resolved to first column", {
    Y <- rbind(c(1, 1, 0), c(0, 1, 1))
    expect_equal(map(Y), c(1, 2))
})

test_that("(unmap:basic): round-trip via map preserves codes", {
    data(nutrimouse)
    expect_equal(map(unmap(nutrimouse$diet)), as.integer(nutrimouse$diet))
})

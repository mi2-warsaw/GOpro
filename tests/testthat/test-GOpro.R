context("rtcga data")

test_that("rtcga object is valid", {
    expect_true(length(exrtcga) > 2)
    obj <- rownames(exrtcga)
    expect_identical(obj[[1]], obj[[2]])
    expect_identical(obj[[1]], obj[[3]])
})


context("aovTopTest")

groups <- assay(experiments(exrtcga))

test_that("when top is missing or 0 expect warning", {
    expect_warning(aovTopTest(groups, 0))
})

test_that("empty list causes error", {
    expect_error(aovTopTest(list(), 10))
})

context("groupByTukey")

test_that("empty list causes error", {
    expect_error(groupByTukey(list()))
})

context("fisher_test")

test_that("fisher.test/stats and fisher_test/GOpro p-values are equal", {
    f_stats <- fisher.test(matrix(c(3, 1, 2, 7), nrow = 2), alternative = 'greater')[[1]]
    f_GOpro <- fisher_test(3, 1, 2, 7)
    expect_equal(f_stats, f_GOpro)
})

context("unbundleCluster")

test_that("the number of groups = 2*(number of rows)-1", {
    kk <- data.frame(a = c(1, 5, 1), b = c(1, 2, 1))
    rownames(kk) <- c('r1', 'r2', 'r3')
    expect_equal(length(unbundleCluster(hclust(dist(kk)))[[1]]), 2 * nrow(kk) - 1)
})

context("ontologyOver")

test_that("over-represented ontology is correct", {
    expect_identical(ontologyOver(data.frame(GROUP = 1, ONTOLOGY = 'BP')), 'BP')
})

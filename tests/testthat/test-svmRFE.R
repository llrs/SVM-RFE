test_that("multiplication works", {
    set.seed(12345)
    library(e1071)
    data(input)
    res <- svmRFE(input, k = 10, halve.above = 100)
    expect_equal(head(res, 23),
                 c(337,319, 173, 175, 503, 287, 304, 492, 1, 47, 286, 214, 191,
                   402, 313, 177, 187, 407, 255, 468, 256, 259, 391 ))
})

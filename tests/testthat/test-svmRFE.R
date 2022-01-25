test_that("multiplication works", {
    set.seed(12345)
    data(input)
    res <- svmRFE(input[, -1], response = input$DX2, k = 10, halve.above = 100)
    expect_equal(head(res, 23),
                 c(337L, 173L, 313L, 1L, 94L, 503L, 47L, 255L, 147L, 50L, 175L,
                   187L, 205L, 217L, 256L, 497L, 186L, 502L, 319L, 262L, 169L, 287L,
                   481L))
})

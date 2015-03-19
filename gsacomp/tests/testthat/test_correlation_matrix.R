test_that("The block correlation matrix works correctly", {
    M <- matrix(0.11, nrow = 3, ncol = 3)
    diag(M) <- 1
    expect_that(createCorrelationMatrix(size = 3, rho = 0),
                equals(diag(3)))
    expect_that(createCorrelationMatrix(size = 3, rho = 0.11),
                equals(M))
    expect_that(createCorrelationMatrix(size = 3, rho = 2),
                throws_error("rho must be <= 1"))
    expect_that(createBlockCorrelationMatrix(sizes = 3, rho = 0.11),
                is_equivalent_to(M))
})

## test_that("runCamera is reproducible", {
##     set.seed(123)
##     out <- runCamera()
##     pv <- c(0.371034610354304,
##             0.380796198965292,
##             0.672732881527667,
##             0.973447702597027, 
##             0.973447702597027)
##     expect_that(out$FDR, equals(pv))
## })

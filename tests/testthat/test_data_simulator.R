library(testthat)

# test for data simulator
test_that("generate data with binary dependent variable", {
    parameter <- list(
        N = 100,
        P = 500,
        p_ref = 30,
        tau = 0, 
        sigma = 1,
        rho = 0.5,
        link = "logit"
    )
    }
)
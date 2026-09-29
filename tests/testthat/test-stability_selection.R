# Stability selection needs q^2 < p * PFER (cutoff below 1), so the number of
# predictors and q are chosen accordingly.
P <- 20
q <- 5
PFER <- 2

one_data <- function(...) sim_data(n = 1, N = 80, P = P, p_ref = 4, ...)[[1]]

# glmnet.adalasso ------------------------------------------------------------

adalasso_input <- function(seed = 1, N = 100) {
    set.seed(seed)
    x <- matrix(rnorm(N * 8), nrow = N, dimnames = list(NULL, paste0("X", 1:8)))
    y <- rbinom(N, 1, plogis(3 * x[, 1]))
    list(x = x, y = y)
}

test_that("glmnet.adalasso returns selected variables and a selection path", {
    d <- adalasso_input()
    out <- glmnet.adalasso(d$x, d$y, q = 3, l2_lambda = 10, type = "conservative",
                           family = "binomial", standardize = FALSE)

    expect_named(out, c("selected", "path"))
    expect_type(out$selected, "logical")
    expect_named(out$selected, paste0("X", 1:8))
    expect_equal(nrow(out$path), 8)
    expect_true(is.matrix(out$path))
})

test_that("glmnet.adalasso limits the number of selected variables", {
    d <- adalasso_input(N = 200)
    cons <- glmnet.adalasso(d$x, d$y, q = 3, l2_lambda = 10, type = "conservative",
                            family = "binomial", standardize = FALSE)
    # glmnet >= 5.0 deprecates the dfmax argument used by this type
    anti <- suppressWarnings(
        glmnet.adalasso(d$x, d$y, q = 3, l2_lambda = 10, type = "anticonservative",
                        family = "binomial", standardize = FALSE)
    )

    expect_lte(sum(cons$selected), 3)
    expect_lte(sum(anti$selected), 2)
})

test_that("glmnet.adalasso selects a strong true predictor", {
    d <- adalasso_input(N = 300)
    out <- glmnet.adalasso(d$x, d$y, q = 3, l2_lambda = 1, type = "conservative",
                           family = "binomial", standardize = FALSE)

    expect_true(out$selected[["X1"]])
})

test_that("glmnet.adalasso supports unweighted fits and rejects an unknown type", {
    d <- adalasso_input()

    expect_no_error(glmnet.adalasso(d$x, d$y, q = 3, l2_lambda = 10, type = "conservative",
                                    family = "binomial", weighted = FALSE, standardize = FALSE))
    expect_error(glmnet.adalasso(d$x, d$y, q = 3, l2_lambda = 10, type = "other",
                                 family = "binomial", standardize = FALSE))
})

# cpss_adaptive_lasso --------------------------------------------------------

test_that("cpss_adaptive_lasso returns a stabsel object", {
    res <- cpss_adaptive_lasso(one_data(), "y", features = NULL, q = q, PFER = PFER, seed = 1)

    expect_s3_class(res, "stabsel")
    expect_equal(res$q, q)
    expect_lte(res$PFER, PFER)
    expect_equal(nrow(res$phat), P)
})

test_that("cpss_adaptive_lasso is reproducible for a fixed seed", {
    a <- cpss_adaptive_lasso(one_data(), "y", NULL, q, PFER, seed = 1)
    b <- cpss_adaptive_lasso(one_data(), "y", NULL, q, PFER, seed = 1)

    expect_identical(a$phat, b$phat)
})

test_that("cpss_adaptive_lasso respects features and binary_class", {
    d <- one_data()
    sub <- cpss_adaptive_lasso(d, "y", paste0("X", 1:10), q = 3, PFER = PFER, seed = 1)
    expect_equal(nrow(sub$phat), 10)

    renamed <- data.table::setnames(data.table::copy(d), "y", "outcome")
    res <- cpss_adaptive_lasso(renamed, "outcome", NULL, q, PFER, seed = 1)
    expect_equal(nrow(res$phat), P)
})

# cpss_glmboost --------------------------------------------------------------

test_that("cpss_glmboost returns a stabsel object", {
    res <- cpss_glmboost(one_data(), "y", features = NULL, iter = 50, q = q, PFER = PFER, seed = 1)

    expect_s3_class(res, "stabsel")
    expect_equal(res$q, q)
    expect_equal(nrow(res$phat), P + 1)  # one row per base-learner, incl. the intercept
})

test_that("cpss_glmboost is reproducible for a fixed seed", {
    a <- cpss_glmboost(one_data(), "y", NULL, iter = 50, q = q, PFER = PFER, seed = 1)
    b <- cpss_glmboost(one_data(), "y", NULL, iter = 50, q = q, PFER = PFER, seed = 1)

    expect_identical(a$phat, b$phat)
})

# cpss_gamboost --------------------------------------------------------------

gam_one <- function() gam_fixture(1, N = 100)[[1]]

test_that("cpss_gamboost returns a stabsel object and honours iter", {
    short <- cpss_gamboost(gam_one(), "y", features = NULL, iter = 60, q = 3, PFER = PFER, seed = 1)
    long  <- cpss_gamboost(gam_one(), "y", features = NULL, iter = 150, q = 3, PFER = PFER, seed = 1)

    expect_s3_class(short, "stabsel")
    expect_equal(short$q, 3)
    expect_lt(ncol(short$phat), ncol(long$phat))
})

test_that("cpss_gamboost is reproducible for a fixed seed", {
    a <- cpss_gamboost(gam_one(), "y", NULL, iter = 80, q = 3, PFER = PFER, seed = 1)
    b <- cpss_gamboost(gam_one(), "y", NULL, iter = 80, q = 3, PFER = PFER, seed = 1)

    expect_identical(a$phat, b$phat)
})

# run_cpss_* -----------------------------------------------------------------

test_that("run_cpss_adaptive_lasso runs on every data set with q = 10 and PFER = 2", {
    dat <- sim_data(n = 2, N = 60, P = 100, p_ref = 4)
    res <- run_cpss_adaptive_lasso(dat, "y", seed = 1)

    expect_length(res, 2)
    for (r in res) {
        expect_s3_class(r, "stabsel")
        expect_equal(r$q, 10)
        expect_lte(r$PFER, 2)
    }
})

test_that("run_cpss_glmboost runs on every data set with q = 10 and PFER = 2", {
    dat <- sim_data(n = 2, N = 60, P = 100, p_ref = 4)
    res <- run_cpss_glmboost(dat, "y", iter = 200, seed = 1)

    expect_length(res, 2)
    for (r in res) {
        expect_s3_class(r, "stabsel")
        expect_equal(r$q, 10)
        expect_lte(r$PFER, 2)
    }
})

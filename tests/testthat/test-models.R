# Fit small fixtures once and share them between tests.
cv <- cv_fixture(n = 1, N = 80, P = 10, p_ref = 3)[[1]]
train <- get_train_cv_data(cv)
test <- get_test_cv_data(cv)
ridge <- withr::with_seed(10, fit_ridge(train, binary_class = "y"))
adalasso <- withr::with_seed(10, fit_adaptive_lasso(train, binary_class = "y"))
glmb <- withr::with_seed(10, fit_glmboost(train, binary_class = "y", cores = 2))

predictors <- paste0("X", 1:10)

# fit_ridge ------------------------------------------------------------------

test_that("fit_ridge returns a cv.glmnet model", {
    expect_s3_class(ridge, "cv.glmnet")
    expect_equal(ridge$glmnet.fit$dim[1], 10)
})

test_that("fit_ridge does not shrink any coefficient to exactly zero", {
    coefs <- as.vector(coef(ridge, s = "lambda.min"))[-1]

    expect_true(all(coefs != 0))
})

test_that("fit_ridge respects features and binary_class", {
    sub <- withr::with_seed(
        10, fit_ridge(train, binary_class = "y", features = c("X1", "X2"))
    )
    expect_equal(sub$glmnet.fit$dim[1], 2)

    renamed <- data.table::copy(train)
    data.table::setnames(renamed, "y", "outcome")
    renamed_fit <- withr::with_seed(
        10, fit_ridge(renamed, binary_class = "outcome")
    )
    expect_equal(renamed_fit$glmnet.fit$dim[1], 10)
})

# fit_adaptive_lasso ---------------------------------------------------------

test_that("fit_adaptive_lasso returns a finite cv.glmnet model", {
    expect_s3_class(adalasso, "cv.glmnet")
    expect_equal(adalasso$glmnet.fit$dim[1], 10)
    expect_true(all(is.finite(as.vector(coef(adalasso, s = "lambda.min")))))
})

test_that("fit_adaptive_lasso respects features and binary_class", {
    sub <- withr::with_seed(
        10, fit_adaptive_lasso(train, binary_class = "y", features = c("X1", "X2", "X3"))
    )
    expect_equal(sub$glmnet.fit$dim[1], 3)

    renamed <- data.table::copy(train)
    data.table::setnames(renamed, "y", "outcome")
    renamed_fit <- withr::with_seed(
        10, fit_adaptive_lasso(renamed, binary_class = "outcome")
    )
    expect_equal(renamed_fit$glmnet.fit$dim[1], 10)
})

# fit_glmboost ---------------------------------------------------------------

test_that("fit_glmboost returns a model with a valid mstop", {
    expect_s3_class(glmb, "glmboost")
    expect_gte(mboost::mstop(glmb), 1)
    expect_lte(mboost::mstop(glmb), 1000)
})

test_that("fit_glmboost is reproducible for a fixed RNG seed", {
    again <- withr::with_seed(10, fit_glmboost(train, binary_class = "y", cores = 2))

    expect_equal(mboost::mstop(again), mboost::mstop(glmb))
    expect_equal(coef(again), coef(glmb))
})

test_that("fit_glmboost restricts the model to requested features", {
    sub <- withr::with_seed(
        10, fit_glmboost(train, binary_class = "y", features = c("X1", "X2"), cores = 2)
    )
    learners <- setdiff(colnames(mboost::extract(sub, "design")), "(Intercept)")

    expect_setequal(learners, c("X1", "X2"))
})

# make_gamboost_formula and fit_gamboost -------------------------------------

test_that("make_gamboost_formula builds linear, smooth and intercept learners", {
    dat <- gam_fixture(1)[[1]]
    f <- make_gamboost_formula(dat, "y", features = c("age", "sex"))

    expect_equal(as.character(f[[2]]), "y")
    expect_setequal(
        attr(terms(f), "term.labels"),
        c("bols(age, intercept = FALSE)",
          "bols(sex, intercept = FALSE)",
          "bbs(age, knots = 8, degree = 4, df = 1, center = TRUE)",
          "bols(Intercept, intercept = FALSE)")
    )
})

test_that("make_gamboost_formula uses all other columns when features is NULL", {
    dat <- gam_fixture(1)[[1]]
    labels <- attr(terms(make_gamboost_formula(dat, "y")), "term.labels")

    expect_true("bols(x1, intercept = FALSE)" %in% labels)
    expect_true(any(grepl("^bbs\\(x1", labels)))
    expect_false(any(grepl("^bbs\\(sex", labels)))
})

test_that("fit_gamboost fits a model with the requested base-learners", {
    dat <- gam_fixture(1)[[1]]
    fit <- withr::with_seed(
        10, fit_gamboost(dat, binary_class = "y", features = c("age", "sex"), cores = 2)
    )

    expect_s3_class(fit, "gamboost")
    expect_gte(mboost::mstop(fit), 1)
    bnames <- names(fit$baselearner)
    expect_true(any(grepl("^bbs\\(age", bnames)))
    expect_true(any(grepl("^bols\\(sex", bnames)))
    expect_false(any(grepl("^bbs\\(sex", bnames)))
    expect_true(any(grepl("^bols\\(Intercept", bnames)))
})

# extract_adalasso_selected_features -----------------------------------------

test_that("extract_adalasso_selected_features returns selected predictors", {
    sel <- extract_adalasso_selected_features(adalasso)

    expect_s3_class(sel, "data.table")
    expect_named(sel, "variable")
    expect_true(all(sel$variable %in% predictors))
})

test_that("extract_adalasso_selected_features gives no rows for a model that selects nothing", {
    x <- as.matrix(train[, !"y"])
    y <- as.numeric(train$y) - 1
    empty <- withr::with_seed(
        10,
        glmnet::cv.glmnet(x, y, family = "binomial", lambda = c(100, 50),
                          foldid = rep(1:5, length.out = nrow(x)))
    )

    sel <- extract_adalasso_selected_features(empty)
    expect_equal(nrow(sel), 0)
})

# extract_glmboost_varimp ----------------------------------------------------

test_that("extract_glmboost_varimp returns sorted variable importance", {
    vi <- extract_glmboost_varimp(glmb)

    expect_s3_class(vi, "data.table")
    expect_true(all(c("Variable", "blearner", "reduction", "selfreq") %in% names(vi)))
    expect_true(all(diff(vi$reduction) <= 0))
})

test_that("summarise_glmboost_varimp aggregates positive risk reductions", {
    vi <- data.table::data.table(
        Variable = c("x1", "x1", "x2"),
        blearner = c("linear", "linear", "linear"),
        tar_batch = c(1L, 1L, 1L),
        tar_rep = c(1L, 1L, 1L),
        tar_seed = c(1L, 1L, 1L),
        reduction = c(2, 4, -1),
        selfreq = c(0.5, 1, 0)
    )

    out <- summarise_glmboost_varimp(vi)
    expect_s3_class(out, "data.table")
    expect_equal(nrow(out), 1)
    expect_equal(out$mrd, 3)
    expect_equal(out$msf, 0.75)
    expect_equal(out$selfreq, 2L)
})

# test-set predictions -------------------------------------------------------

expect_valid_prediction <- function(pred, newdata) {
    expect_s3_class(pred, "data.table")
    expect_equal(nrow(pred), nrow(newdata))
    expect_true(all(c("y", "response", "predict") %in% names(pred)))
    expect_s3_class(pred$y, "factor")
    expect_true(all(pred$response >= 0 & pred$response <= 1))
    expect_equal(levels(pred$predict), c("0", "1"))
    expect_equal(as.character(pred$predict), ifelse(pred$response > 0.5, "1", "0"))
}

test_that("get_ridge_test_prediction returns response and class for every test row", {
    expect_valid_prediction(get_ridge_test_prediction(ridge, test), test)
})

test_that("get_adalasso_test_prediction returns response and class for every test row", {
    expect_valid_prediction(get_adalasso_test_prediction(adalasso, test), test)
})

test_that("get_glmboost_test_prediction returns response and class for every test row", {
    expect_valid_prediction(get_glmboost_test_prediction(glmb, test), test)
})

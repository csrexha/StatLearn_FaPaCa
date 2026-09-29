# Models are fitted once and shared by the tests below.
set.seed(10)
cv    <- cv_fixture(n = 2, N = 80, P = 10, p_ref = 3)
train <- get_train_cv_data(cv)
test  <- get_test_cv_data(cv)
ridge <- fit_ridge(train, binary_class = "y")
adalasso <- fit_adaptive_lasso(train, binary_class = "y")
glmb  <- fit_glmboost(train, binary_class = "y", seed = 1)

predictors <- paste0("X", 1:10)

# fit_ridge ------------------------------------------------------------------

test_that("fit_ridge returns one cv.glmnet model per data set", {
    expect_length(ridge, 2)
    for (m in ridge) {
        expect_s3_class(m, "cv.glmnet")
        expect_equal(m$glmnet.fit$dim[1], 10)
    }
})

test_that("fit_ridge does not shrink any coefficient to exactly zero", {
    coefs <- as.vector(coef(ridge[[1]], s = "lambda.min"))[-1]

    expect_true(all(coefs != 0))
})

test_that("fit_ridge respects features and binary_class", {
    sub <- fit_ridge(train, binary_class = "y", features = c("X1", "X2"))
    expect_equal(sub[[1]]$glmnet.fit$dim[1], 2)

    renamed <- lapply(train, function(d) data.table::setnames(data.table::copy(d), "y", "outcome"))
    expect_equal(fit_ridge(renamed, binary_class = "outcome")[[1]]$glmnet.fit$dim[1], 10)
})

# fit_adaptive_lasso ---------------------------------------------------------

test_that("fit_adaptive_lasso returns one cv.glmnet model per data set", {
    expect_length(adalasso, 2)
    for (m in adalasso) {
        expect_s3_class(m, "cv.glmnet")
        expect_equal(m$glmnet.fit$dim[1], 10)
        expect_true(all(is.finite(as.vector(coef(m, s = "lambda.min")))))
    }
})

test_that("fit_adaptive_lasso respects features and binary_class", {
    sub <- fit_adaptive_lasso(train, binary_class = "y", features = c("X1", "X2", "X3"))
    expect_equal(sub[[1]]$glmnet.fit$dim[1], 3)

    renamed <- lapply(train, function(d) data.table::setnames(data.table::copy(d), "y", "outcome"))
    expect_equal(fit_adaptive_lasso(renamed, binary_class = "outcome")[[1]]$glmnet.fit$dim[1], 10)
})

# fit_glmboost ---------------------------------------------------------------

test_that("fit_glmboost returns one glmboost model per data set with a valid mstop", {
    expect_length(glmb, 2)
    for (m in glmb) {
        expect_s3_class(m, "glmboost")
        expect_gte(mboost::mstop(m), 1)
        expect_lte(mboost::mstop(m), 1000)
    }
})

test_that("fit_glmboost is reproducible for a fixed seed", {
    again <- fit_glmboost(train, binary_class = "y", seed = 1)

    expect_equal(mboost::mstop(again[[1]]), mboost::mstop(glmb[[1]]))
    expect_equal(coef(again[[1]]), coef(glmb[[1]]))
})

test_that("fit_glmboost restricts the model to the requested features", {
    sub <- fit_glmboost(train, binary_class = "y", features = c("X1", "X2"), seed = 1)
    learners <- setdiff(colnames(mboost::extract(sub[[1]], "design")), "(Intercept)")

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

test_that("fit_gamboost returns one gamboost model per data set", {
    fits <- fit_gamboost(gam_fixture(2), binary_class = "y", features = c("age", "sex"), seed = 1)

    expect_length(fits, 2)
    for (m in fits) {
        expect_s3_class(m, "gamboost")
        expect_gte(mboost::mstop(m), 1)
        bnames <- names(m$baselearner)
        expect_true(any(grepl("^bbs\\(age", bnames)))
        expect_true(any(grepl("^bols\\(sex", bnames)))
        expect_false(any(grepl("^bbs\\(sex", bnames)))
        expect_true(any(grepl("^bols\\(Intercept", bnames)))
    }
})

# extract_adalasso_selected_features -----------------------------------------

test_that("extract_adalasso_selected_features returns Simulation and variable", {
    sel <- extract_adalasso_selected_features(adalasso)

    expect_s3_class(sel, "data.table")
    expect_named(sel, c("Simulation", "variable"))
    expect_true(all(sel$variable %in% predictors))
    expect_true(all(sel$Simulation %in% 1:2))
})

test_that("extract_adalasso_selected_features gives no rows for a model that selects nothing", {
    d <- train[[1]]
    x <- as.matrix(d[, !"y"])
    y <- as.numeric(d$y) - 1
    empty <- glmnet::cv.glmnet(x, y, family = "binomial",
                               lambda = c(100, 50),  # far above lambda.max: no coefficient is non-zero
                               foldid = rep(1:5, length.out = nrow(x)))

    sel <- extract_adalasso_selected_features(list(empty))
    expect_equal(nrow(sel), 0)
})

# extract_glmboost_varimp ----------------------------------------------------

test_that("extract_glmboost_varimp summarises the variable importance per model", {
    vi <- extract_glmboost_varimp(glmb)

    expect_s3_class(vi, "data.table")
    expect_true(all(c("Variable", "blearner", "Simulation", "mrd", "msf", "selfreq") %in% names(vi)))
    expect_true(all(vi$Simulation %in% 1:2))
    expect_true(all(vi$mrd > 0))
    expect_false(is.unsorted(rev(vi$selfreq)))
})

# test-set predictions -------------------------------------------------------

expect_valid_prediction <- function(pred, newdata) {
    expect_s3_class(pred, "data.table")
    expect_equal(nrow(pred), sum(vapply(newdata, nrow, integer(1))))
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
    pred <- get_glmboost_test_prediction(glmb, test)

    expect_valid_prediction(pred, test)
    expect_true("Simulation" %in% names(pred))
    expect_equal(as.vector(table(pred$Simulation)),
                 vapply(test, nrow, integer(1)))
})

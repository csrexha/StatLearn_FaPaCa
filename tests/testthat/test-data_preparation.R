# binary_data_generator ------------------------------------------------------

test_that("binary_data_generator returns y and x of the right shape", {
    set.seed(1)
    out <- binary_data_generator(N = 50, P = 8, p_ref = 3, tau = 0, sigma = 1, rho = 0.5)

    expect_named(out, c("y", "x"))
    expect_length(out$y, 50)
    expect_true(all(out$y %in% 0:1))
    expect_equal(dim(out$x), c(50, 8))
})

test_that("binary_data_generator guarantees at least 4 cases of the minor class", {
    set.seed(1)
    out <- binary_data_generator(N = 30, P = 5, p_ref = 2, tau = -4, sigma = 1, rho = 0.5)

    expect_gte(sum(out$y), 4)
})

test_that("binary_data_generator is reproducible", {
    args <- list(N = 40, P = 6, p_ref = 2, tau = 0, sigma = 1, rho = 0.5)
    set.seed(2); a <- do.call(binary_data_generator, args)
    set.seed(2); b <- do.call(binary_data_generator, args)

    expect_identical(a, b)
})

test_that("binary_data_generator correlates only the relevant predictors", {
    set.seed(3)
    out <- binary_data_generator(N = 3000, P = 6, p_ref = 3, tau = 0, sigma = 1, rho = 0.8)
    cm <- cor(out$x)

    relevant <- mean(cm[1:3, 1:3][lower.tri(diag(3))])
    noise    <- mean(abs(cm[4:6, 4:6][lower.tri(diag(3))]))

    expect_gt(relevant, 0.7)
    expect_lt(noise, 0.1)
})

test_that("binary_data_generator accepts another link and rejects P < p_ref", {
    set.seed(4)
    expect_no_error(
        binary_data_generator(N = 40, P = 6, p_ref = 2, tau = 0, sigma = 1, rho = 0.5,
                              link = "probit")
    )
    expect_error(
        binary_data_generator(N = 20, P = 2, p_ref = 5, tau = 0, sigma = 1, rho = 0.5)
    )
})

# get_simulated_data ---------------------------------------------------------

test_that("get_simulated_data returns n data.tables with factor y and X1..XP", {
    dat <- sim_data(n = 3, N = 40, P = 8, p_ref = 2)

    expect_length(dat, 3)
    for (d in dat) {
        expect_s3_class(d, "data.table")
        expect_named(d, c("y", paste0("X", 1:8)))
        expect_s3_class(d$y, "factor")
        expect_equal(nrow(d), 40)
    }
})

test_that("get_simulated_data is reproducible and depends on the seed", {
    expect_identical(sim_data(seed = 5), sim_data(seed = 5))
    expect_false(identical(sim_data(seed = 5), sim_data(seed = 6)))
})

test_that("get_simulated_data expands vectors of N and P into one data set per value", {
    by_n <- get_simulated_data(N = rep(c(20, 30), each = 2), P = 6, p_ref = 2,
                               tau = -1, sigma = 2, rho = 0.5, seed = 1)
    expect_equal(vapply(by_n, nrow, integer(1)), c(20L, 20L, 30L, 30L))

    by_p <- get_simulated_data(N = 30, P = c(6, 10), p_ref = 2,
                               tau = -1, sigma = 2, rho = 0.5, seed = 1)
    expect_equal(vapply(by_p, ncol, integer(1)), c(7L, 11L))
})

test_that("get_simulated_data expands vectors of p_ref, tau and rho", {
    by_pref <- get_simulated_data(N = 30, P = 8, p_ref = c(2, 4), tau = -1, sigma = 2,
                                  rho = 0.5, seed = 1)
    expect_length(by_pref, 2)

    by_tau <- get_simulated_data(N = 30, P = 8, p_ref = 2, tau = c(0, -1, -2), sigma = 2,
                                 rho = 0.5, seed = 1)
    expect_length(by_tau, 3)

    by_rho <- get_simulated_data(N = 30, P = 8, p_ref = 2, tau = -1, sigma = 2,
                                 rho = c(0.1, 0.5, 0.9), seed = 1)
    expect_length(by_rho, 3)

    by_rho_n <- get_simulated_data(n = 2, N = 30, P = 8, p_ref = 2, tau = -1, sigma = 2,
                                   rho = c(0.1, 0.5, 0.9), seed = 1)
    expect_length(by_rho_n, 6)
})

test_that("get_simulated_data draws a default rho when rho is NULL", {
    dat <- get_simulated_data(n = 2, N = 50, P = 8, p_ref = 2, tau = -1, sigma = 2,
                              rho = NULL, seed = 1)
    expect_length(dat, 2)
})

test_that("get_simulated_data requires a seed", {
    expect_error(get_simulated_data(n = 1))
})

# get_stratified_cv_data -----------------------------------------------------

test_that("get_stratified_cv_data splits every data set into train and test", {
    dat <- sim_data(n = 2, N = 60, P = 12)
    cv  <- get_stratified_cv_data(dat, seed = 1)

    expect_length(cv, 2)
    for (x in cv) {
        expect_named(x, c("train", "test"))
        expect_equal(nrow(x$train) + nrow(x$test), 60)
        expect_gte(nrow(x$train), 40)
        expect_lte(nrow(x$train), 44)
        expect_s3_class(x$train$y, "factor")
        expect_equal(nlevels(droplevels(x$train$y)), 2)
        expect_equal(nlevels(droplevels(x$test$y)), 2)
    }
})

test_that("get_stratified_cv_data standardises with the training statistics only", {
    cv <- get_stratified_cv_data(sim_data(n = 1, N = 60, P = 12), seed = 1)[[1]]
    train <- as.matrix(cv$train[, !"y"])
    test  <- as.matrix(cv$test[, !"y"])

    expect_equal(unname(colMeans(train)), rep(0, 12), tolerance = 1e-8)
    expect_equal(unname(apply(train, 2, sd)), rep(1, 12), tolerance = 1e-8)
    expect_gt(max(abs(colMeans(test))), 1e-3)
})

test_that("get_stratified_cv_data works for any number of predictors", {
    expect_no_error(get_stratified_cv_data(sim_data(n = 1, P = 7, p_ref = 2), seed = 1))
    expect_no_error(get_stratified_cv_data(sim_data(n = 1, P = 30, p_ref = 5), seed = 1))
})

test_that("get_stratified_cv_data is reproducible and leaves the input unchanged", {
    dat <- sim_data()
    dat_copy <- lapply(dat, data.table::copy)

    expect_identical(get_stratified_cv_data(dat, seed = 1),
                     get_stratified_cv_data(dat, seed = 1))
    expect_false(identical(get_stratified_cv_data(dat, seed = 1),
                           get_stratified_cv_data(dat, seed = 2)))
    expect_identical(dat, dat_copy)
})

# get_train_cv_data / get_test_cv_data ---------------------------------------

test_that("get_train_cv_data and get_test_cv_data extract the right element", {
    cv <- cv_fixture(n = 3)

    expect_identical(get_train_cv_data(cv), lapply(cv, function(x) x$train))
    expect_identical(get_test_cv_data(cv), lapply(cv, function(x) x$test))
    expect_length(get_train_cv_data(cv), 3)
})

# read_raw_data --------------------------------------------------------------

test_that("read_raw_data reads a csv and converts character columns to factors", {
    f <- withr::local_tempfile(fileext = ".csv")
    data.table::fwrite(data.frame(id = 1:4, grp = c("a", "b", "a", "c"),
                                  val = c(1.5, 2, NA, 4)), f)
    out <- read_raw_data(f)

    expect_s3_class(out, "data.table")
    expect_s3_class(out$grp, "factor")
    expect_equal(levels(out$grp), c("a", "b", "c"))
    expect_true(is.numeric(out$val))
    expect_true(is.integer(out$id))
})

test_that("read_raw_data works without character columns and fails on a missing file", {
    f <- withr::local_tempfile(fileext = ".csv")
    data.table::fwrite(data.frame(a = 1:3, b = c(0.1, 0.2, 0.3)), f)

    expect_no_error(read_raw_data(f))
    expect_error(read_raw_data(withr::local_tempfile(fileext = ".csv")))
})

# impute_data ----------------------------------------------------------------

impute_fixture <- function() {
    set.seed(5)
    dt <- data.table(a = rnorm(40), b = rnorm(40),
                     f = factor(sample(c("x", "y"), 40, replace = TRUE)))
    dt[c(3, 7), a := NA]
    dt[5, f := NA]
    dt
}

test_that("impute_data fills all missing values and keeps observed ones", {
    dt  <- impute_fixture()
    res <- impute_data(dt, features = c("a", "b", "f"), maxiter = 2, ntree = 10)

    expect_false(anyNA(res))
    expect_equal(dim(res), c(40, 3))
    expect_equal(res$b, dt$b)
    expect_equal(res$a[-c(3, 7)], dt$a[-c(3, 7)])
    expect_equal(levels(res$f), levels(dt$f))
})

test_that("impute_data returns only the selected features", {
    res <- impute_data(impute_fixture(), features = c("a", "b"), maxiter = 2, ntree = 10)

    expect_named(res, c("a", "b"))
})

test_that("impute_data fails for an unknown feature", {
    expect_error(impute_data(impute_fixture(), features = "zzz", maxiter = 2, ntree = 10))
})

# export_imputed_data / read_imputed_data ------------------------------------

test_that("export_imputed_data returns the path and the data survive a round trip", {
    f  <- withr::local_tempfile(fileext = ".csv")
    dt <- data.table(a = c(1.5, 2.5, 3.5), b = c("u", "v", "w"))

    expect_identical(export_imputed_data(dt, f), f)
    expect_true(file.exists(f))

    back <- read_imputed_data(f)
    expect_s3_class(back, "data.table")
    expect_equal(back, dt)
})

# get_fapaca_cv_data ---------------------------------------------------------

test_that("get_fapaca_cv_data returns 4 x 10 train and test folds", {
    cv <- get_fapaca_cv_data(fapaca_like_data(), features_rescale = c("age", "prot1"), seed = 1)

    expect_named(cv, c("train", "test"))
    expect_length(cv$train, 40)
    expect_length(cv$test, 40)
})

test_that("get_fapaca_cv_data folds are disjoint and cover all rows", {
    dt <- fapaca_like_data()
    cv <- get_fapaca_cv_data(dt, features_rescale = c("age", "prot1"), seed = 1)

    for (i in seq_along(cv$train)) {
        expect_length(intersect(cv$train[[i]]$ID, cv$test[[i]]$ID), 0)
        expect_setequal(c(cv$train[[i]]$ID, cv$test[[i]]$ID), dt$ID)
    }
})

test_that("get_fapaca_cv_data centres the features on the training fold only", {
    dt <- fapaca_like_data()
    cv <- get_fapaca_cv_data(dt, features_rescale = c("age", "prot1"), seed = 1)

    for (i in seq_along(cv$train)) {
        expect_equal(mean(cv$train[[i]]$age), 0, tolerance = 1e-8)
        expect_equal(mean(cv$train[[i]]$prot1), 0, tolerance = 1e-8)
        # centring keeps differences between observations
        expect_equal(diff(cv$train[[i]]$age),
                     diff(dt[match(cv$train[[i]]$ID, ID), age]),
                     tolerance = 1e-8)
    }
    # columns that are not rescaled are untouched
    expect_equal(cv$train[[1]]$prot2, dt[match(cv$train[[1]]$ID, ID), prot2])
    expect_equal(cv$test[[1]]$sex, dt[match(cv$test[[1]]$ID, ID), sex])
})

test_that("get_fapaca_cv_data is stratified by Status", {
    dt <- fapaca_like_data()
    cv <- get_fapaca_cv_data(dt, features_rescale = "age", seed = 1)
    overall <- mean(dt$Status == "L")

    for (x in cv$train) {
        expect_equal(mean(x$Status == "L"), overall, tolerance = 0.05)
    }
})

test_that("get_fapaca_cv_data is reproducible and does not modify its input", {
    dt <- fapaca_like_data()
    dt_copy <- data.table::copy(dt)

    a <- get_fapaca_cv_data(dt, features_rescale = "age", seed = 1)
    b <- get_fapaca_cv_data(dt, features_rescale = "age", seed = 1)
    c <- get_fapaca_cv_data(dt, features_rescale = "age", seed = 2)

    expect_identical(a, b)
    expect_false(identical(a, c))
    expect_identical(dt, dt_copy)
})

# get_fapaca_train_data ------------------------------------------------------

test_that("get_fapaca_train_data keeps the requested features in every fold", {
    cv  <- get_fapaca_cv_data(fapaca_like_data(), features_rescale = "age", seed = 1)
    out <- get_fapaca_train_data(cv, features = c("age", "sex"))

    expect_length(out, 40)
    for (x in out) expect_named(x, c("age", "sex"))
})

test_that("get_fapaca_train_data validates its input", {
    cv <- get_fapaca_cv_data(fapaca_like_data(), features_rescale = "age", seed = 1)

    expect_error(get_fapaca_train_data(list(train = "a"), "age"), "list of data.tables")
    expect_error(get_fapaca_train_data(list(), "age"), "list of data.tables")
    expect_error(get_fapaca_train_data(cv, features = "not_a_column"))
})

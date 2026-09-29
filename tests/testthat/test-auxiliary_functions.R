predict_fixture <- function(n = 60, seed = 6) {
    set.seed(seed)
    d <- data.table(ID = seq_len(n), Status = NA, FamilyID = rep(1:6, length.out = n),
                    Clinical_Findings = rep(c("a", "b"), length.out = n),
                    x1 = rnorm(n), x2 = rnorm(n))
    d[, Status := factor(rbinom(n, 1, plogis(1.5 * x1)), levels = 0:1)]
    fit <- mboost::glmboost(Status ~ x1 + x2, data = d, family = mboost::Binomial(),
                            control = mboost::boost_control(mstop = 50))
    list(fit = fit, data = d)
}

test_that("predict_result returns identifiers, offset, base learners, eta and response", {
    fx <- predict_fixture()
    out <- predict_result(fx$fit, fx$data, fx$data, bl = c("x1", "x2"))

    expect_s3_class(out, "data.table")
    expect_equal(nrow(out), nrow(fx$data))
    expect_named(out, c("ID", "Status", "FamilyID", "Clinical_Findings",
                        "offset", "x1", "x2", "eta", "response"))
})

test_that("predict_result computes eta as offset plus the base learner contributions", {
    fx <- predict_fixture()
    out <- predict_result(fx$fit, fx$data, fx$data, bl = c("x1", "x2"))

    expect_equal(out$eta, out$offset + out$x1 + out$x2)
    expect_equal(order(out$response), order(out$eta))
    expect_true(all(out$response > 0 & out$response < 1))
})

test_that("predict_result: response is plogis(2 * eta) when all base-learners are requested", {
    fx <- predict_fixture()
    # mboost's Binomial() works on the half-log-odds scale; the intercept is a base-learner too
    out <- predict_result(fx$fit, fx$data, fx$data, bl = c("(Intercept)", "x1", "x2"))

    expect_equal(out$response, plogis(2 * out$eta))
})

test_that("predict_result keeps only the requested cols and accepts several links", {
    fx <- predict_fixture()

    out <- predict_result(fx$fit, fx$data, fx$data, bl = "x1", cols = "ID")
    expect_named(out, c("ID", "offset", "x1", "eta", "response"))

    expect_no_error(predict_result(fx$fit, fx$data, fx$data, bl = "x1", link = "probit"))
    expect_no_error(predict_result(fx$fit, fx$data, fx$data, bl = "x1", link = plogis))
})

# summarySE ------------------------------------------------------------------

summary_fixture <- function() {
    data.frame(g = rep(c("a", "b"), each = 4),
               v = c(1, 2, 3, 4, 10, 12, 14, 16))
}

test_that("summarySE summarises a variable by group", {
    out <- summarySE(summary_fixture(), measurevar = "v", groupvars = "g")

    expect_equal(nrow(out), 2)
    expect_named(out, c("g", "N", "v", "sd", "se", "ci"))
    expect_equal(out$N, c(4, 4))
    expect_equal(out$v, c(2.5, 13))
    expect_equal(out$sd, c(sd(1:4), sd(c(10, 12, 14, 16))))
    expect_equal(out$se, out$sd / sqrt(out$N))
    expect_equal(out$ci, out$se * qt(0.975, out$N - 1))
})

test_that("summarySE handles missing values and the confidence level", {
    d <- summary_fixture()
    d$v[1] <- NA

    expect_true(is.na(summarySE(d, "v", "g")$v[1]))

    out <- summarySE(d, "v", "g", na.rm = TRUE)
    expect_equal(out$N, c(3, 4))
    expect_equal(out$v[1], mean(2:4))

    narrow <- summarySE(summary_fixture(), "v", "g", conf.interval = 0.8)
    wide   <- summarySE(summary_fixture(), "v", "g", conf.interval = 0.99)
    expect_true(all(narrow$ci < wide$ci))
})

test_that("summarySE without grouping variables summarises all rows", {
    out <- summarySE(summary_fixture(), measurevar = "v")

    expect_equal(nrow(out), 1)
    expect_equal(out$N, 8)
    expect_equal(out$v, mean(summary_fixture()$v))
})

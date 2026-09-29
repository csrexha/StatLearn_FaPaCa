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

test_that("predict_result keeps only the requested cols and accepts several links", {
    fx <- predict_fixture()

    out <- predict_result(fx$fit, fx$data, fx$data, bl = "x1", cols = "ID")
    expect_named(out, c("ID", "offset", "x1", "eta", "response"))

    expect_no_error(predict_result(fx$fit, fx$data, fx$data, bl = "x1", link = "probit"))
    expect_no_error(predict_result(fx$fit, fx$data, fx$data, bl = "x1", link = plogis))
})

#' Perform stability selection with adaptive lasso
#'
#' @description
#' The function runs complementary pairs stability selection using adaptive lasso model as
#' feature selection model.
#'
#' @param data  data
#' @param q     number of (unique) selected variables that are selected on each subsample.
#' @param PFER  upper bound for the per-family error rate.
#'
#' @returns stabsel object
#' @export

cpss_adaptive_lasso <- function(data, q, PFER){

    # convert data to matrix
    xmat <- as.matrix(data[, -1])
    y <- data$y

    # stratified subsampling
    stabs_rsmp <- subsample(rep(1, nrow(xmat)), B = 50, strata = as.factor(y))

    simdata_stabsel <- stabsel(
        x = xmat,
        y = y,
        fitfun = glmnet.adalasso,
        args.fitfun = list(
            type = "conservative",
            family = "binomial",
            standardize = FALSE,
            l2_lambda = 10,
            weighted = TRUE,
            gamma = 2
        ),
        sampling.type = "SS",
        assumption = "unimod",
        B = 50,
        folds  = stabs_rsmp,
        q = q,
        PFER = PFER
    )

    return(simdata_stabsel)
}

#' Perform stability selection with glmboost
#'
#' @description
#' The function runs complementary pairs stability selection using glmboost model as
#' feature selection model.
#'
#' @param data  data
#' @param q     number of (unique) selected variables that are selected on each subsample.
#' @param PFER  upper bound for the per-family error rate.
#'
#' @returns stabsel object
#' @export

cpss_glmboost <- function(data, q = 10, PFER = 2){

    dimnames(data$x) <- list(NULL, paste0("X", 1:ncol(data$x)))

    # fit the logistic boosting model
    mboost_fit <- glmboost(
        x = cbind(Intercept = 1, data$x),
        y = as.factor(data$y),
        family = Binomial(link = "logit"),
        control = boost_control(mstop = 500, nu = 0.1)
    )

    # stratified subsampling
    stabs_rsmp <- subsample(
        model.weights(mboost_fit),
        B = 50,
        strata = mboost_fit$response
    )

    mboost_stabsel <- stabsel(
        mboost_fit,
        q = q,
        PFER = PFER,
        sampling.type = "SS",
        assumption = "unimod",
        folds = stabs_rsmp,
        grid = 0:500
    )

    return(mboost_stabsel)
}


#' Run stability selection with adaptive lasso on a list of data sets
#'
#' @description
#' Run stability selection using adaptive lasso model on a list of data sets with default
#' parameters q=10 and PFER=2
#'
#' @param data  A list of data sets
#'
#' @returns A list of stabsel objects
#' @export

run_cpss_adaptive_lasso <- function(data) {
    data |>
        map(function(x) cpss_adaptive_lasso(x, 10, 2))
}

#' Run stability selection with glmboost on a list of data sets
#'
#' @description
#' Run stability selection using glmboost model on a list of data sets with default
#' parameters q=10 and PFER=2
#'
#' @param data  A list of data sets
#'
#' @returns A list of stabsel objects
#' @export

run_cpss_glmboost <- function(data) {
    data |>
        map(function(x) cpss_glmboost(x, 10, 2))
}

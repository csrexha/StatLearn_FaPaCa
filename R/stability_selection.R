#' Adaptive lasso for stability selection
#'
#' @param x dependent variable matrix
#' @param y independent variable vector
#' @param q Coefficient-count cap: limit the maximum number of variables in the model at any lambda.
#' @param l2_lambda lambda for the ridge regression (to determine the penalty factor)
#' @param type conservative or anticonservative
#' @param family family for the glm
#' @param weighted the observations should be weighted or not. Default TRUE
#' @param gamma gamma for penalty factor. The higher the gamma, the more aggressive the penalty factor. Default 2.
#'
#' @returns list of selected variables and the selection paths
#' @export
glmnet.adalasso <- function(x, y, q, l2_lambda, type = c("conservative", "anticonservative"),
                            family, weighted = TRUE, gamma = 2, ...) {
    require(glmnet)

    # calculate the weights
    if(weighted){
        w0 <- 0.5 * length(y)/sum(y == 0)
        w1 <- 0.5 * length(y)/sum(y == 1)
        w01 <- ifelse(y == 0, w0, w1)
    } else{
        w01 <- rep(1, length(y))
    }

    # determine the penalty factor
    ridgefit <- glmnet(x, y, alpha = 0, weight = w01, intercept = TRUE, family = family, ...)
    weight <- as.vector(1/abs(coef(ridgefit, s = l2_lambda))**gamma)[-1]
    pf <- ncol(x) * weight/sum(weight)

    # fit model
    type <- match.arg(type)
    if (type == "conservative")
        fit <- suppressWarnings(glmnet(x, y, pmax = q, alpha = 1,
                                       weight = w01, family = family, penalty.factor = pf, ...))
    if (type == "anticonservative")
        fit <- glmnet(x, y, dfmax = q - 1, alpha = 1,
                      weight = w01, family = family, penalty.factor = pf, ...)

    # which coefficients are non-zero?
    selected <- predict(fit, type = "nonzero")
    selected <- selected[[length(selected)]]
    ret <- logical(ncol(x))
    ret[selected] <- TRUE
    names(ret) <- colnames(x)

    # compute selection paths
    cf <- fit$beta
    sequence <- as.matrix(cf != 0)

    # return both
    return(list(selected = ret, path = sequence))
}

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

    # fit the logistic boosting model
    mboost_fit <- glmboost(
        y ~ .,
        data = data,
        family = Binomial(type = "adaboost", link = "logit"),
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

#' Adaptive lasso fitting function for stability selection
#'
#' @description
#' Fitting function for [stabs::stabsel()]. A ridge regression at the fixed penalty
#' `l2_lambda` gives adaptive weights `1 / |coefficient|^gamma`, rescaled to sum to
#' `ncol(x)`. An adaptive lasso path with these penalty factors is then fitted and the
#' variables selected at the end of the path are returned.
#'
#' @param x Numeric matrix of predictors, with column names.
#' @param y Numeric binary outcome vector (0 and 1).
#' @param q Maximum number of selected variables. For `type = "conservative"` at most `q`
#'   coefficients may be non-zero (`pmax = q`), for `"anticonservative"` the degrees of
#'   freedom are limited to `q - 1` (`dfmax = q - 1`).
#' @param l2_lambda Penalty of the ridge regression used to compute the adaptive weights.
#' @param type `"conservative"` (default) or `"anticonservative"`.
#' @param family Family of the generalised linear model, e.g. `"binomial"`.
#' @param weighted Whether observations are weighted so that both classes carry the same
#'   total weight. Default `TRUE`.
#' @param gamma Exponent of the adaptive weights. The larger the value, the more
#'   aggressive the penalty factors. Default 2.
#' @param ... Further arguments passed to both [glmnet::glmnet()] calls, e.g.
#'   `standardize`.
#'
#' @returns A list with `selected`, a named logical vector of length `ncol(x)`, and
#'   `path`, a logical matrix (variables in rows, lambda values in columns) marking the
#'   non-zero coefficients along the path.
#' @seealso [stabs::stabsel()]
#' @export
glmnet.adalasso <- function(x, y, q, l2_lambda, type = c("conservative", "anticonservative"),
                            family, weighted = TRUE, gamma = 2, ...) {
    # calculate the weights
    if(weighted){
        w0 <- 0.5 * length(y)/sum(y == 0)
        w1 <- 0.5 * length(y)/sum(y == 1)
        w01 <- ifelse(y == 0, w0, w1)
    } else{
        w01 <- rep(1, length(y))
    }

    # determine the penalty factor
    ridgefit <- glmnet(x, y, alpha = 0, weights = w01, intercept = TRUE, family = family, ...)
    weight <- as.vector(1/abs(coef(ridgefit, s = l2_lambda))**gamma)[-1]
    pf <- ncol(x) * weight/sum(weight)

    # fit model
    type <- match.arg(type)
    # pmax / dfmax are passed via `control` (glmnet >= 5.0 deprecates them as arguments)
    if (type == "conservative")
        fit <- suppressWarnings(glmnet(x, y, alpha = 1, weights = w01, family = family,
                                       penalty.factor = pf, control = list(pmax = q), ...))
    if (type == "anticonservative")
        fit <- glmnet(x, y, alpha = 1, weights = w01, family = family,
                      penalty.factor = pf, control = list(dfmax = q - 1), ...)

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

#' Stability selection with adaptive lasso
#'
#' @description
#' Runs complementary pairs stability selection ([stabs::stabsel()], `sampling.type = "SS"`)
#' with the adaptive lasso ([glmnet.adalasso()]) as the selection method.
#'
#' @details
#' Uses 50 subsamples stratified by the outcome, the unimodality assumption, a ridge
#' penalty of 10 for the adaptive weights, `gamma = 2` and class weights. The data are
#' not standardised inside the fit (`standardize = FALSE`).
#'
#' @param data A data.table with the outcome column `binary_class` (two-level factor with
#'   levels 0 and 1) and the predictors.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. `NULL` uses all columns
#'   except `binary_class`.
#' @param q Number of (unique) variables selected on each subsample. Together with
#'   `PFER` and the number of variables `p` it must satisfy `q^2 < p * PFER`, otherwise
#'   the selection threshold exceeds 1 and [stabs::stabsel()] fails.
#' @param PFER Upper bound for the per-family error rate.
#'
#' @returns A [stabs::stabsel()] object.
#' @seealso [stabs::stabsel()]
#' @export
cpss_adaptive_lasso <- function(data, binary_class, features, q, PFER) {

    # convert data to matrix
    if(!is.null(features)) {
        xmat <- as.matrix(data[, ..features])
    } else {
        xmat <- as.matrix(data[, !..binary_class])
    }

    y <- as.numeric(data[[binary_class]]) - 1

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

#' Stability selection with glmboost
#'
#' @description
#' Runs complementary pairs stability selection ([stabs::stabsel()], `sampling.type = "SS"`)
#' with [mboost::glmboost()] as the selection method.
#'
#' @details
#' Uses 50 subsamples stratified by the outcome, the unimodality assumption, a logistic
#' `Binomial(type = "adaboost")` model with step length `nu = 0.1`, the iteration grid
#' `0:iter` and 2 cores. If `iter` is too small to select `q` base-learners in some
#' subsamples, [stabs::stabsel()] warns; increase `iter` then.
#'
#' @param data A data.table with the outcome column `binary_class` (two-level factor) and
#'   the predictors.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. `NULL` uses all columns
#'   except `binary_class`.
#' @param iter Number of boosting iterations.
#' @param q Number of (unique) variables selected on each subsample. Together with
#'   `PFER` and the number of variables `p` it must satisfy `q^2 < p * PFER`, otherwise
#'   the selection threshold exceeds 1 and [stabs::stabsel()] fails.
#' @param PFER Upper bound for the per-family error rate.
#'
#' @returns A [stabs::stabsel()] object. The intercept counts as a base-learner.
#' @seealso [stabs::stabsel()]
#' @export
cpss_glmboost <- function(data, binary_class, features, iter, q, PFER) {

    if(!is.null(features)) {
        model <- as.formula(paste0(binary_class, "~",
                                   paste(features, collapse = " + ")))
    } else {
        model <- formula(paste0(binary_class, "~ ."))
    }

    # fit the logistic boosting model
    mboost_fit <- glmboost(
        model,
        data = data,
        family = Binomial(type = "adaboost", link = "logit"),
        control = boost_control(mstop = iter, nu = 0.1)
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
        grid = 0:iter,
        mc.cores = 2
    )

    return(mboost_stabsel)
}


#' Stability selection with gamboost
#'
#' @description
#' Runs complementary pairs stability selection ([stabs::stabsel()], `sampling.type = "SS"`)
#' with [mboost::gamboost()] as the selection method.
#'
#' @details
#' The base-learners are the same as in [fit_gamboost()] (linear for every feature, centred
#' smooth for numeric features, intercept). Uses 50 subsamples stratified by the outcome,
#' the unimodality assumption, step length `nu = 0.1`, the iteration grid `0:iter` and 2
#' cores. If `iter` is too small to select `q` base-learners in some subsamples,
#' [stabs::stabsel()] warns; increase `iter` then.
#'
#' @param data A data.table with the outcome column `binary_class` (two-level factor) and
#'   the predictors. It must not contain a column named `Intercept`.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. `NULL` uses all columns
#'   except `binary_class`.
#' @param iter Number of boosting iterations.
#' @param q Number of (unique) variables selected on each subsample. Together with
#'   `PFER` and the number of variables `p` it must satisfy `q^2 < p * PFER`, otherwise
#'   the selection threshold exceeds 1 and [stabs::stabsel()] fails.
#' @param PFER Upper bound for the per-family error rate.
#' @param seed Random seed for the subsampling.
#'
#' @returns A [stabs::stabsel()] object.
#' @seealso [stabs::stabsel()]
#' @export
cpss_gamboost <- function(data, binary_class, features, iter, q, PFER) {

    model <- make_gamboost_formula(data, binary_class, features)

    # fit the logistic boosting model
    fit <- mboost::gamboost(
        formula = model,
        data = cbind(data, Intercept = 1),
        family = mboost::Binomial(link = "logit"),
        control = mboost::boost_control(mstop = iter, nu = 0.1)
        )

    # stratified subsampling
    stabs_rsmp <- subsample(
        model.weights(fit),
        B = 50,
        strata = fit$response
    )

    mboost_stabsel <- stabsel(
        fit,
        q = q,
        PFER = PFER,
        sampling.type = "SS",
        assumption = "unimod",
        folds = stabs_rsmp,
        grid = 0:iter,
        mc.cores = 2
    )

    return(mboost_stabsel)
}

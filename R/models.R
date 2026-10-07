#' Fit ridge regression models
#'
#' @description
#' Fits a ridge-penalised logistic regression ([glmnet::cv.glmnet()] with `alpha = 0`)
#' to one data set.
#'
#' @details
#' The penalty is chosen by 5-fold cross-validation of the binomial deviance over 500
#' lambda values. Observations are weighted so that both classes carry the same total
#' weight. The folds are drawn at random from the current RNG state (there is no `seed`
#' argument).
#'
#' @param data A data.table containing the two-level factor `binary_class` (levels 0 and
#'   1) and the predictors.
#' @param binary_class Name of the binary outcome column in `data`.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#'
#' @returns A [glmnet::cv.glmnet()] object.
#' @export
fit_ridge <- function(data, binary_class, features = NULL) {

    # convert data to matrix
    if(!is.null(features)) {
        xmat <- as.matrix(data[, ..features])
        } else {
            xmat <- as.matrix(data[, !..binary_class])
            }
            
    y <- as.numeric(data[[binary_class]]) - 1

    # assign weights
    w0 <- 0.5*length(y)/sum(y == 0)
    w1 <- 0.5*length(y)/sum(y == 1)
    w01 <- ifelse(y == 0, w0, w1)

    foldid <- createFolds(as.factor(y), k = 5, list = FALSE)
    
    # fit ridge regression
    glmnet::cv.glmnet(x = xmat, 
                      y = y,
                      alpha = 0,
                      nlambda = 500,
                      weights = w01,
                      family = "binomial",
                      type.measure = "deviance",
                      foldid = foldid)
}

#' Fit adaptive lasso models
#'
#' @description
#' Fits an adaptive lasso logistic regression to a data set, in two steps:
#' a ridge regression ([glmnet::cv.glmnet()] with `alpha = 0`) gives the coefficients at
#' `lambda.min`, from which the penalty factors `1 / |coefficient|^2` are computed; then
#' a lasso ([glmnet::cv.glmnet()] with `alpha = 1`) is fitted with these penalty factors.
#'
#' @details
#' Both steps use 500 lambda values, class weights that give both classes the same total
#' weight, and the same 5 cross-validation folds. The folds are drawn at random from the
#' current RNG state (there is no `seed` argument).
#'
#' @param data A data.table containing the two-level factor `binary_class` (levels 0 and
#'   1) and the predictors.
#' @param binary_class Name of the binary outcome column in `data`.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#'
#' @returns A [glmnet::cv.glmnet()] object.
#' @export
fit_adaptive_lasso <- function(data, binary_class, features = NULL) {

    # convert data to matrix
    if (!is.null(features)) {
        xmat <- as.matrix(data[, ..features])
    } else {
        xmat <- as.matrix(data[, !..binary_class])
    }

    y <- as.numeric(data[[binary_class]]) - 1

    # assign weights
    w0 <- 0.5 * length(y) / sum(y == 0)
    w1 <- 0.5 * length(y) / sum(y == 1)
    w01 <- ifelse(y == 0, w0, w1)

    foldid <- createFolds(as.factor(y), k = 5, list = FALSE)

    # estimate penalty factor (adaptive weights) with ridge regression
    l2fit <- glmnet::cv.glmnet(x = xmat,
                               y = y,
                               alpha = 0,
                               nlambda = 500,
                               weights = w01,
                               family = "binomial",
                               type.measure = "deviance",
                               foldid = foldid)

    # obtain penalty factor
    weight <- as.vector(1 / abs(coef(l2fit, s = "lambda.min"))^2)[-1]

    # fit adaptive lasso
    glmnet::cv.glmnet(x = xmat,
                      y = y,
                      alpha = 1,
                      penalty.factor = weight,
                      weights = w01,
                      nlambda = 500,
                      family = "binomial",
                      foldid = foldid)
}


#' Fit glmboost models
#'
#' @description
#' Fits a logistic [mboost::glmboost()] model (`Binomial(type = "adaboost")`, step length
#' `nu = 0.1`) and chooses the number of boosting iterations.
#'
#' @details
#' The number of iterations starts at 50 and is increased in steps of 50 (up to 1000)
#' until it exceeds 1.2 times the AIC-optimal number. The final number of iterations is
#' then chosen by bootstrap resampling stratified by the outcome ([mboost::cvrisk()] with
#' `cores` cores).
#'
#' @param data A data.table containing the two-level factor `binary_class` (levels 0 and
#'   1) and the predictors.
#' @param binary_class Name of the binary outcome column in `data`.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#' @param cores Number of cores used by [mboost::cvrisk()]. Default 4.
#'
#' @returns An [mboost::glmboost()] model.
#' @export
fit_glmboost <- function(data, binary_class, features = NULL, cores = 4) {

    if (!is.null(features)) {
        model <- as.formula(paste0(
            binary_class, " ~ ", paste(features, collapse = " + ")
        ))
    } else {
        model <- formula(paste0(binary_class, " ~ ."))
    }

    # Set the initial number of iterations.
    iter <- 50

    # Fit the logistic boosting model.
    fit <- mboost::glmboost(
        model,
        data = data,
        family = mboost::Binomial(type = "adaboost", link = "logit"),
        control = mboost::boost_control(mstop = iter, nu = 0.1)
    )

    # Adjust the number of iterations using AIC (logit link only).
    aic <- AIC(fit, method = "classical")

    while (iter <= mboost::mstop(aic) * 1.2 && iter < 1000) {
        iter <- iter + 50
        mboost::mstop(fit) <- iter
        aic <- AIC(fit, method = "classical")
    }

    # Set a resampling scheme.
    rsmp <- mboost::cv(
        model.weights(fit),
        type = "bootstrap",
        strata = fit$response
    )

    # Search for the optimal iteration using resampling.
    fit_cvrisk <- mboost::cvrisk(
        fit,
        folds = rsmp,
        mc.cores = cores
    )

    # Set the optimal number of iterations.
    mboost::mstop(fit) <- mboost::mstop(fit_cvrisk)
    
    return(fit)
}

#' Build the formula of a gamboost model
#'
#' @description
#' Every feature gets a linear base-learner (\code{bols}). Numeric features with more
#' than two distinct values also get a centred smooth base-learner (\code{bbs}), so
#' that the smooth term only captures the non-linear part. Binary and factor
#' features get the linear base-learner only. An intercept base-learner
#' (\code{bols(Intercept)}) is always added: the data passed to
#' \code{\link[mboost]{gamboost}} must contain a column \code{Intercept = 1}.
#'
#' @param data A data frame
#' @param binary_class The name of the binary class variable
#' @param features A vector of feature names. Default is NULL, which means all
#'   columns except \code{binary_class}.
#'
#' @returns A formula
#' @noRd
make_gamboost_formula <- function(data, binary_class, features = NULL) {

    if (!binary_class %in% names(data)) {
        stop(sprintf("Column '%s' not found in data.", binary_class))
    }

    if ("Intercept" %in% names(data)) {
        stop("The data must not contain a column named 'Intercept'; it is used internally by fit_gamboost().")
    }

    if (is.null(features)) {
        features <- setdiff(names(data), binary_class)
    }

    is_smooth <- vapply(
        features,
        function(f) is.numeric(data[[f]]) && length(unique(data[[f]])) > 2,
        logical(1)
    )
    smooth <- features[is_smooth]

    blrns <- c(
        paste0("bols(", features, ", intercept = FALSE)"),
        if (length(smooth) > 0) {
            paste0("bbs(", smooth, ", knots = 8, degree = 4, df = 1, center = TRUE)")
        },
        "bols(Intercept, intercept = FALSE)"
    )

    as.formula(paste0(binary_class, " ~ ", paste(blrns, collapse = " + ")))
}

#' Fit gamboost models
#'
#' @description
#' Fits a logistic [mboost::gamboost()] model (`Binomial(link = "logit")`, step length
#' `nu = 0.1`) to a data set and chooses the number of boosting iterations.
#' Every feature gets a linear base-learner and numeric features also a centred smooth
#' base-learner; see the details.
#'
#' @details
#' The base-learners are built by an internal helper: every feature gets a linear
#' `bols(x, intercept = FALSE)` learner, numeric features with more than two distinct
#' values also get a centred smooth `bbs(x, knots = 8, degree = 4, df = 1, center = TRUE)`
#' learner, and an intercept learner `bols(Intercept, intercept = FALSE)` is added. The
#' column `Intercept` is added to the data internally, so the data must not contain a
#' column with that name. Because the linear learners have no intercept, the features
#' should be centred.
#'
#' The number of iterations starts at 50 and is increased in steps of 50 (up to 1000)
#' until it exceeds 1.2 times the AIC-optimal number. The final number of iterations is
#' then chosen by bootstrap resampling stratified by the outcome ([mboost::cvrisk()] with
#' `cores` cores).
#'
#' @param data A data frame or data.table containing the two-level factor `binary_class`
#'   (levels 0 and 1) and the predictors.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#' @param cores Number of cores used by [mboost::cvrisk()]. Default 1.
#'
#' @returns An [mboost::gamboost()] model.
#' @export
fit_gamboost <- function(data, binary_class, features = NULL, cores = 1) {

    model <- make_gamboost_formula(data, binary_class, features)

    # setting initial number of iteration
    iter <- 50

    # fit the logistic boosting model
    fit <- mboost::gamboost(formula = model,
                            data = cbind(data, Intercept = 1),
                            family = mboost::Binomial(link = "logit"),
                            control = mboost::boost_control(mstop = iter, nu = 0.1))

    # adjust the number of iteration using AIC (logit link only)
    aic <- stats::AIC(fit, method = "classical")

    while(iter <= mboost::mstop(aic) * 1.2 && iter < 1000){
        iter <- iter + 50
        mboost::mstop(fit) <- iter
        aic <- stats::AIC(fit, method = "classical")
    }

    # set a resampling scheme.
    rsmp <- mboost::cv(model.weights(fit),
                        type = "bootstrap",
                        strata = fit$response)

    # using resampling to search for the optimal iteration.
    fit_cvrisk <- mboost::cvrisk(fit,
                                    folds = rsmp,
                                    mc.cores = cores)

    # set the mstop to optimal mstop
    mboost::mstop(fit) <- mboost::mstop(fit_cvrisk)

    return(fit)
}

#' Extract the selected features of adaptive lasso models
#'
#' @description
#' Extracts the predictors with a non-zero coefficient at `lambda.min`.
#'
#' @param model A [glmnet::cv.glmnet()] model, e.g. one returned by
#'   [fit_adaptive_lasso()].
#'
#' @returns A data.table with a `variable` column containing the names of selected
#'   predictors. The intercept is excluded; if no predictors are selected, the table
#'   has zero rows.
#' @export
extract_adalasso_selected_features <- function(model) {

    coefs <- as.matrix(coef(model, s = "lambda.min"))
    
    selected <- data.table(variable = setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)"))
    
    return(selected)
}

#' Extract the variable importance of glmboost models
#'
#' @description
#' Extracts variable importance ([mboost::varimp()]) for a fitted model and sorts the
#' results by decreasing risk reduction.
#'
#' @param model A [mboost::glmboost()] model, e.g. one returned by [fit_glmboost()].
#'
#' @returns A data.table containing the columns returned by [mboost::varimp()], sorted
#'   by decreasing `reduction`. The `variable` column is renamed to `Variable`, and
#'   `Variable` and `blearner` are converted to character.
#' @export
extract_glmboost_varimp <- function(model) {

    # extract variable importance of each fold 
    glmboost_varimp <- as.data.table(varimp(model))[order(reduction, decreasing = TRUE)]
    
    setnames(glmboost_varimp, "variable", "Variable")
    
    glmboost_varimp[, `:=` (blearner  = as.character(blearner),
                         Variable  = as.character(Variable))]

    return(glmboost_varimp)
    }

#' Summarise glmboost variable importance
#'
#' @description
#' Summarises positive risk reductions by variable, base-learner and simulation
#' identifiers.
#'
#' @param varimp A data.table containing glmboost variable importance values,
#'   including `reduction`, `selfreq`, `Variable`, `blearner`, `tar_batch`,
#'   `tar_rep` and `tar_seed` columns.
#'
#' @returns A data.table grouped by `Variable`, `blearner`, `tar_batch`,
#'   `tar_rep` and `tar_seed`, containing mean risk reduction (`mrd`), mean
#'   selection frequency (`msf`) and the number of rows summarised (`selfreq`),
#'   sorted by decreasing `selfreq`.
#' @export
summarise_glmboost_varimp <- function(varimp) {

    # summarise the results of all folds
    varimp[reduction > 0, 
        .(mrd = mean(reduction), msf = mean(selfreq), selfreq = .N),
        by = .(Variable, blearner, tar_batch, tar_rep, tar_seed)][order(-selfreq)]
}

#' Predict the test sets with adaptive lasso models
#'
#' @description
#' Predicts a test set with a fitted adaptive lasso model at `lambda.min`.
#'
#' @param model A [glmnet::cv.glmnet()] model, e.g. from [fit_adaptive_lasso()].
#' @param newdata A data.table of test observations. The outcome column `y` must be
#'   first, followed by the predictors used in the model.
#'
#' @returns A data.table with one row per test observation and columns `y` (observed
#'   outcome, converted to a factor), `response` (predicted probability) and `predict`
#'   (predicted class, a factor with levels 0 and 1).
#' @export
get_adalasso_test_prediction <- function(model, newdata) {

    # obtain test results
    data.table(y = as.factor(newdata$y),
               response = as.vector(predict(model, 
                                            s = model$lambda.min,
                                            newx = as.matrix(newdata[, -1]),
                                            type = "response")),
               predict = factor(predict(model,
                                        s = model$lambda.min,
                                        newx = as.matrix(newdata[, -1]),
                                        type = "class"),
                                        levels = c(0, 1)))
}

#' Predict the test sets with glmboost models
#'
#' @description
#' Predicts a test set with a fitted glmboost model.
#'
#' @param model A [mboost::glmboost()] model, e.g. from [fit_glmboost()].
#' @param newdata A data.table of test observations containing the outcome column `y`
#'   and the predictors used in the model.
#'
#' @returns A data.table with one row per test observation and columns `y` (observed
#'   outcome), `response` (predicted probability) and `predict` (predicted class, a
#'   factor with levels 0 and 1, classified using a 0.5 threshold).
#' @export
get_glmboost_test_prediction <- function(model, newdata) {

    # obtain test result
    simdata_glmboost_testpred <- data.table(
        y = newdata$y,
        response = as.vector(predict(model, newdata = newdata, type = "response"))
        )

    # classify the test prediction
    simdata_glmboost_testpred[, predict := factor(ifelse(response > 0.5, 1, 0), levels = c(0, 1))]

    return(simdata_glmboost_testpred)
}

#' Predict the test sets with ridge regression models
#'
#' @description
#' Predicts a test set with a fitted ridge regression model at `lambda.min`.
#'
#' @param model A [glmnet::cv.glmnet()] model, e.g. from [fit_ridge()].
#' @param newdata A data.table of test observations. The outcome column `y` must be
#'   first, followed by the predictors used in the model.
#'
#' @returns A data.table with one row per test observation and columns `y` (observed
#'   outcome, converted to a factor), `response` (predicted probability) and `predict`
#'   (predicted class, a factor with levels 0 and 1).
#' @export
get_ridge_test_prediction <- function(model, newdata) {

    # obtain test results
    data.table(y = as.factor(newdata$y),
               response = as.vector(predict(model, 
                                            s = model$lambda.min,
                                            newx = as.matrix(newdata[, -1]),
                                            type = "response")),
               predict = factor(predict(model,
                                        s = model$lambda.min,
                                        newx = as.matrix(newdata[, -1]),
                                        type = "class"),
                                        levels = c(0, 1)))
}

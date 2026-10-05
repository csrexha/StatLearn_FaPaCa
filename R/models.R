#' Fit ridge regression models
#'
#' @description
#' Fits a ridge-penalised logistic regression ([glmnet::cv.glmnet()] with `alpha = 0`) to
#' each data set in a list.
#'
#' @details
#' The penalty is chosen by 5-fold cross-validation of the binomial deviance over 500
#' lambda values. Observations are weighted so that both classes carry the same total
#' weight. The folds are drawn at random from the current RNG state (there is no `seed`
#' argument).
#'
#' @param data A list of data.tables, one per data set. Each contains the two-level factor
#'   `binary_class` (levels 0 and 1) and the predictors.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#'
#' @returns A list of [glmnet::cv.glmnet()] objects, one per data set.
#' @export
fit_ridge <- function(data, binary_class, features = NULL) {

    ## perform ridge regression / L2-regularisation ####
    data |>
        map (
            function(data) {

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
                glmnet::cv.glmnet(x = xmat, y = y,
                                  alpha = 0,
                                  nlambda = 500,
                                  weights = w01,
                                  family = "binomial",
                                  type.measure = "deviance",
                                  foldid = foldid)
            }
        )
}

#' Fit adaptive lasso models
#'
#' @description
#' Fits an adaptive lasso logistic regression to each data set in a list, in two steps:
#' a ridge regression ([glmnet::cv.glmnet()] with `alpha = 0`) gives the coefficients at
#' `lambda.min`, from which the penalty factors `1 / |coefficient|^2` are computed; then
#' a lasso ([glmnet::cv.glmnet()] with `alpha = 1`) is fitted with these penalty factors.
#'
#' @details
#' Both steps use 500 lambda values, class weights that give both classes the same total
#' weight, and the same 5 cross-validation folds. The folds are drawn at random from the
#' current RNG state (there is no `seed` argument).
#'
#' @param data A list of data.tables, one per data set. Each contains the two-level factor
#'   `binary_class` (levels 0 and 1) and the predictors.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#'
#' @returns A list of [glmnet::cv.glmnet()] objects, one per data set.
#' @export
fit_adaptive_lasso <- function(data, binary_class, features = NULL) {
    data |>
        map(
            function(data) {

                # convert data to matrix
                if(!is.null(features)) {
                    xmat <- as.matrix(data[, ..features])
                } else {
                    xmat <- as.matrix(data[, !..binary_class])
                }

                y <- as.numeric(data[[binary_class]]) - 1

                # assign weights
                w0 <- 0.5 * length(y)/sum(y == 0)
                w1 <- 0.5 * length(y)/sum(y == 1)
                w01 <- ifelse(y == 0, w0, w1)

                foldid <- createFolds(as.factor(y), k = 5, list = FALSE)

                # estimate penalty factor (adpative weights) with ridge regression
                l2fit <- cv.glmnet(x = xmat, y = y,
                                   alpha = 0,
                                   nlambda = 500,
                                   weights = w01,
                                   family = "binomial",
                                   type.measure = "deviance",
                                   foldid = foldid)

                # obtain penalty factor
                weight <- as.vector(1/abs(coef(l2fit, s = "lambda.min"))^2)[-1]

                # fit adalasso
                cv.glmnet(x = xmat, y = y,
                          alpha = 1,
                          penalty.factor = weight,
                          weights = w01,
                          nlambda = 500,
                          family = "binomial",
                          foldid = foldid)
            }
        )
}


#' Fit glmboost models
#'
#' @description
#' Fits a logistic [mboost::glmboost()] model (`Binomial(type = "adaboost")`, step length
#' `nu = 0.1`) to each data set in a list and chooses the number of boosting iterations.
#'
#' @details
#' The number of iterations starts at 50 and is increased in steps of 50 (up to 1000)
#' until it exceeds 1.2 times the AIC-optimal number. The final number of iterations is
#' then chosen by bootstrap resampling stratified by the outcome ([mboost::cvrisk()] with
#' `cores` cores).
#'
#' @param data A list of data.tables, one per data set. Each contains the two-level factor
#'   `binary_class` (levels 0 and 1) and the predictors.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#' @param cores Number of cores used by [mboost::cvrisk()]. Default 4.
#'
#' @returns A list of [mboost::glmboost()] models, one per data set.
#' @export
fit_glmboost <- function(data, binary_class, features = NULL, cores = 4) {

    fit_list <- data |>
        map(
            function(data) {

                if(!is.null(features)) {
                    model <- as.formula(paste0(binary_class, "~",
                                              paste(features, collapse = " + ")))
                } else {
                    model <- formula(paste0(binary_class, "~ ."))
                }

                # setting initial number of iteration
                iter <- 50

                # fit the logistic boosting model
                fit <- glmboost(model,
                                data = data,
                                family = Binomial(type = "adaboost", link = "logit"),
                                control = boost_control(mstop = iter, nu = 0.1))

                # adjust the number of iteration using AIC (logit link only)
                aic <- AIC(fit, method = "classical")

                while(iter <= mstop(aic) * 1.2 && iter < 1000){
                    iter <- iter + 50
                    mstop(fit) <- iter
                    aic <- AIC(fit, method = "classical")
                }

                # set a resampling scheme.
                rsmp <- cv(model.weights(fit),
                           type = "bootstrap",
                           strata = fit$response)

                # using resampling to search for the optimal iteration.
                fit_cvrisk <- cvrisk(fit,
                                     folds = rsmp,
                                     mc.cores = cores)

                # obtain the optimal model according to mstop
                mstop(fit) <- mstop(fit_cvrisk)

                # return fitted model
                return(fit)
            }
        )
    return(fit_list)
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
#' `nu = 0.1`) to each data set in a list and chooses the number of boosting iterations.
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
#' @param data A list of data.tables, one per data set. Each contains the two-level factor
#'   `binary_class` (levels 0 and 1) and the predictors.
#' @param binary_class Name of the binary outcome column.
#' @param features Character vector of the predictors to use. Default `NULL` uses all
#'   columns except `binary_class`.
#' @param cores Number of cores used by [mboost::cvrisk()]. Default 10.
#'
#' @returns A list of [mboost::gamboost()] models, one per data set.
#' @export
fit_gamboost <- function(data, binary_class, features = NULL, cores = 10) {

    fit_list <- data |>
        purrr::map(
            function(data) {

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
        )

    return(fit_list)
}

#' Extract the selected features of adaptive lasso models
#'
#' @description
#' Extracts, for each fitted model, the predictors with a non-zero coefficient at
#' `lambda.min`.
#'
#' @param model A list of [glmnet::cv.glmnet()] models, e.g. from [fit_adaptive_lasso()].
#'
#' @returns A data.table with the columns `Simulation` (position of the model in the list)
#'   and `variable` (name of a selected predictor). A model that selects nothing
#'   contributes no rows.
#' @export
extract_adalasso_selected_features <- function(model) {

    model |>
        imap(
            function(x, idx){
                coefs <- as.matrix(coef(x, s = "lambda.min"))
                selected <- setdiff(rownames(coefs)[coefs[, 1] != 0], "(Intercept)")
                return(data.table(Simulation = rep(idx, length(selected)),
                                  variable = selected))
            }
        ) |>
        rbindlist()
}



#' Extract the variable importance of glmboost models
#'
#' @description
#' Extracts the variable importance ([mboost::varimp()]) of each fitted model and keeps the
#' variables with a positive risk reduction.
#'
#' @details
#' The importance is summarised by variable, base-learner and model. glmboost has one
#' base-learner per variable, so each group contains a single row and `selfreq` is 1.
#'
#' @param model A list of [mboost::glmboost()] models, e.g. from [fit_glmboost()].
#'
#' @returns A data.table sorted by decreasing `selfreq`, with the columns `Variable`,
#'   `blearner`, `Simulation` (position of the model in the list), `mrd` (mean risk
#'   reduction), `msf` (mean relative selection frequency) and `selfreq` (number of rows
#'   summarised).
#' @export
extract_glmboost_varimp <- function(model) {

    ## extract variable importance of each fold ####
    glmboost_varimp <- imap(
        model,
        function(x, idx){
            varimp_logit <- as.data.table(varimp(x))[order(reduction, decreasing = TRUE)]
            setnames(varimp_logit, "variable", "Variable")
            varimp_logit[, `:=` (blearner  = as.character(blearner),
                                 Variable  = as.character(Variable),
                                 Simulation = idx)]
            return(varimp_logit)
        }) |>
        rbindlist()

    ## summarise the results of all folds ####
    glmboost_varimp_summary <- glmboost_varimp[
        reduction > 0,
        .(mrd = mean(reduction), msf = mean(selfreq), selfreq = .N),
        by = .(Variable, blearner, Simulation)][order(-selfreq)]

    return(glmboost_varimp_summary)
}

#' Predict the test sets with adaptive lasso models
#'
#' @description
#' Predicts each test set with the corresponding fitted adaptive lasso model at
#' `lambda.min`.
#'
#' @param model A list of [glmnet::cv.glmnet()] models, e.g. from [fit_adaptive_lasso()].
#' @param newdata A list of data.tables, one per model in the same order. The outcome
#'   column `y` must be the first column, followed by the predictors used in the model.
#'
#' @returns A data.table with one row per test observation, stacked over the models, with
#'   the columns `y` (observed outcome, factor), `response` (predicted probability) and
#'   `predict` (predicted class, factor with levels 0 and 1, threshold 0.5).
#' @export
get_adalasso_test_prediction <- function(model, newdata) {

    # obtain test results
    map2(model, newdata,
         function(x, y){
             data.table(y = as.factor(y$y),
                        response = as.vector(predict(x, s = x$lambda.min,
                                                     newx = as.matrix(y[, -1]),
                                                     type = "response")),
                        predict = factor(predict(x,
                                                 s = x$lambda.min,
                                                 newx = as.matrix(y[, -1]),
                                                 type = "class"),
                                         levels = c(0, 1)))
         }) |>
        rbindlist()
}

#' Predict the test sets with glmboost models
#'
#' @description
#' Predicts each test set with the corresponding fitted glmboost model.
#'
#' @param model A list of [mboost::glmboost()] models, e.g. from [fit_glmboost()].
#' @param newdata A list of data.tables, one per model in the same order, containing the
#'   outcome column `y` and the predictors used in the model.
#'
#' @returns A data.table with one row per test observation, stacked over the models, with
#'   the columns `Simulation` (position of the model in the list), `y` (observed outcome),
#'   `response` (predicted probability) and `predict` (predicted class, factor with
#'   levels 0 and 1, threshold 0.5).
#' @export
get_glmboost_test_prediction <- function(model, newdata) {

    # obtain test result
    simdata_glmboost_testpred <- pmap(
        list(model, newdata, seq_along(model)),
        function(x, y, z) data.table(
            Simulation = z,
            y = y$y,
            response = as.vector(predict(x, newdata = y, type = "response"))
        )
    )  |>
        rbindlist()

    # classify the test prediction
    simdata_glmboost_testpred[, predict := factor(ifelse(response > .5, 1, 0), levels = c(0, 1))]

    return(simdata_glmboost_testpred)
}

#' Predict the test sets with ridge regression models
#'
#' @description
#' Predicts each test set with the corresponding fitted ridge regression model at
#' `lambda.min`.
#'
#' @param model A list of [glmnet::cv.glmnet()] models, e.g. from [fit_ridge()].
#' @param newdata A list of data.tables, one per model in the same order. The outcome
#'   column `y` must be the first column, followed by the predictors used in the model.
#'
#' @returns A data.table with one row per test observation, stacked over the models, with
#'   the columns `y` (observed outcome, factor), `response` (predicted probability) and
#'   `predict` (predicted class, factor with levels 0 and 1, threshold 0.5).
#' @export
get_ridge_test_prediction <- function(model, newdata) {

    # obtain test results
    map2(model, newdata,
         function(x, y){
             data.table(y = as.factor(y$y),
                        response = as.vector(predict(x, s = x$lambda.min,
                                                     newx = as.matrix(y[, -1]),
                                                     type = "response")),
                        predict = factor(predict(x,
                                                 s = x$lambda.min,
                                                 newx = as.matrix(y[, -1]),
                                                 type = "class"),
                                         levels = c(0, 1)))
         }
    ) |>
        rbindlist()
}

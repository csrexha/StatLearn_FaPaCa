#' Fit ridge regression model
#'
#' @description
#' This function fits a ridge regression model to each dataset in the provided list of data.
#'
#' @param data A list of data frames
#'
#' @returns A list of fitted ridge regression models
#' @export

fit_ridge <- function(data) {

    ## perform ridge regression / L2-regularisation ####
    data |>
        map (
            function(data) {

                # convert data to matrix
                xmat <- as.matrix(data[, -1])
                y <- data$y

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

#' Fit adaptive lasso model
#'
#' @description
#' The function fits an adaptive lasso model to each dataset in the provided
#' list of data. It first estimates the penalty factors using ridge regression
#' and then fits the adaptive lasso model using these weights.
#' The function returns a list of fitted models for each dataset.
#'
#' @param data A list of data frames
#'
#' @returns A list of fitted adaptive lasso models
#' @export

fit_adaptive_lasso <- function(data) {
    data |>
        map(
            function(data) {
                # convert data to matrix
                xmat <- as.matrix(data[, -1])
                y <- data$y

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


#' Fit glmboost model
#'
#' @description
#' This function fits a glmboost model to each dataset in the provided list
#' of simulation data. It first sets an initial number of iterations,
#' fits the logistic boosting model, and then adjusts the number of iterations using AIC.
#' The function also uses resampling to search for the optimal iteration.
#' Finally, it returns a list of fitted models for each dataset.
#'
#' @param data A list of data frames
#'
#' @returns A list of fitted glmboost models
#' @export

fit_glmboost <- function(data) {

    data |>
        map(
            function(data) {
                # setting initial number of iteration
                iter <- 50

                # fit the logistic boosting model
                fit <- glmboost(y ~ .,
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
                                     mc.cores = 4)

                # obtain the optimal model according to mstop
                mstop(fit) <- mstop(fit_cvrisk)

                # return fitted model
                return(fit)
            }
        )
}

#' Extract selected features from adaptive lasso model
#'
#' @description
#' This function extracts the selected features from each fitted adaptive lasso model
#' in the provided list of models. It retrieves the coefficients at the optimal lambda value
#' and returns a data table containing the simulation index and the names of
#' the selected variables.
#'
#' @param model A list of fitted adaptive lasso models
#'
#' @returns A data table containing the selected features for each model
#' @export

extract_adalasso_selected_features <- function(model) {

    model |>
        imap(
            function(x, idx){
                coef <- coef(x, s = "lambda.min")
                return(data.table(Simulation = idx,
                                  variable = rownames(coef)[c(summary(coef)$i)[-1]]))
            }
        ) |>
        rbindlist()
}



#' Extract variable importance from glmboost model
#'
#' @description
#' This function extracts the variable importance from each fitted glmboost model.
#'
#' @param model A list of fitted glmboost models
#'
#' @returns A data table containing the variable importance for each model
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

#' Estimate test predictions from fitted adpative lasso model
#'
#' @description
#' This function estimates the test predictions from each fitted adaptive lasso model.
#'
#' @param model A list of fitted adaptive lasso models
#' @param newdata A list of new data tables for prediction
#'
#' @returns A data table containing the test predictions for each model
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

#' Estimate test predictions from fitted glmboost model
#'
#' @description
#' This function estimates the test predictions from each fitted glmboost model.
#'
#' @param model A list of fitted glmboost models
#' @param newdata A list of new data tables for prediction
#'
#' @returns A data table containing the test predictions for each model
#' @export

get_glmboost_test_prediction <- function(model, newdata) {

    # obtain test result
    simdata_glmboost_testpred <- pmap(
        list(model, newdata, 1:100),
        function(x, y, z) data.table(
            Simultaion = z,
            y = y$y,
            response = as.vector(predict(x, newdata = y, type = "response"))
        )
    )  |>
        rbindlist()

    # classify the test prediction
    simdata_glmboost_testpred[, predict := factor(ifelse(response > .5, 1, 0))]

    return(simdata_glmboost_testpred)
}

#' Estimate test predictions from fitted ridge regression model
#'
#' @description
#' This function estimates the test predictions from each fitted ridge regression model.
#'
#' @param model A list of fitted ridge regression models
#' @param newdata A list of new data tables for prediction
#'
#' @returns A data table containing the test predictions for each model
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

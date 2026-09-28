# Load packages required to define the pipeline:
library(targets)
library(crew)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Set target options:
tar_option_set(
  packages = c("data.table", "purrr", "magrittr", "caret", "glmnet", "mboost", "stabs"),
  format = "qs",
  controller = crew_controller_local(workers = 4)
)



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

extract_adalasso_selected_features <- function(fit) {
    fit |>
        imap(
            function(x, idx){
                coef <- coef(x, s = "lambda.min")
                return(data.table(Simulation = idx,
                                  variable = rownames(coef)[c(summary(coef)$i)[-1]]))
                }
            ) |>
        rbindlist()
}

get_adalasso_test_prediction <- function(fit, newdata) {
    map2(fit, newdata,
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
                                control = boost_control(mstop = iter, nu = .1))

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

extract_glmboost_varimp <- function(fit) {

    ## extract variable importance of each fold ####
    glmboost_varimp <- imap(
        fit,
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

get_glmboost_test_prediction <- function(fit, newdata) {
    ## obtain test result ####
    simdata_glmboost_testpred <- pmap(
        list(fit, newdata, 1:100),
        function(x, y, z) data.table(
            Simultaion = z,
            y = y$y,
            response = as.vector(predict(x, newdata = y, type = "response"))
            )
        )  |>
        rbindlist()

    ## classify the test prediction ####
    simdata_glmboost_testpred[, predict := factor(ifelse(response > .5, 1, 0))]

    return(simdata_glmboost_testpred)
}

fit_ridge <- function(data) {

    ## perform ridge regression / L2-regularisation ####
    ridge_fit <- data |>
        map (
         function(data) {

             # convert data to matrix
             xmat <- as.matrix(data[, -1])
             y <- data$y

             # assign weights
             w0 <- .5*length(y)/sum(y == 0)
             w1 <- .5*length(y)/sum(y == 1)
             w01 <- ifelse(y == 0, w0, w1)

             # fit ridge regression
             glmnet::cv.glmnet(x = xmat, y = y,
                               alpha = 0,
                               nlambda = 500,
                               weights = w01,
                               family = "binomial",
                               type.measure = "deviance",
                               nfolds = 10)
         }
        )
}

get_ridge_test_prediction <- function(fit, newdata) {
    ## obtain test results ####
    map2(fit, newdata,
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

# stability selection function with adaptive lasso
cpss_adaptive_lasso <- function(data, q, PFER){

    # stratified subsampling
    stabs_rsmp <- stabs::subsample(rep(1, nrow(data$x)), B = 50, strata = as.factor(data$y))

    simdata_stabsel <- stabs::stabsel(
        x = data$x,
        y = data$y,
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

generate_simulated_data <- function(
        N = 50,
        P = 500,
        p_ref = 10,
        tau = -9,
        sigma = 10,
        rho = runif(10, 0.05, 0.95),
        seed) {

    set.seed(seed)

    CJ(
        N = N,
        P = P,
        p_ref = p_ref,
        tau = tau,
        sigma = sigma,
        rho = rho
    ) |>
        pmap(binary_data_generator)
}

run_cpss_adaptive_lasso <- function(data) {
    data |>
        map(function(data) cpss_adaptive_lasso(data, 10, 2))
}

# function for stability selection using glmboost
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

run_cpss_glmboost <- function(data) {
    data |>
        map(function(data) cpss_glmboost(data, 10, 2))
}

# pipelines
list(
    # simulation experiment 1 --------------------------------------------------------------

    # simulate data
    tar_target(simulated_data, get_simulated_data(n = 100, seed = 119752361)),
    tar_target(scaled_cv_data, get_stratified_cv_data(simulated_data, seed = 873)),
    tar_target(train_data, get_train_cv_data(scaled_cv_data)),
    tar_target(test_data, get_test_cv_data(scaled_cv_data)),

    # fit ridge, adaptive lasso, glmboost regression to training data
    tar_target(ridge, fit_ridge(train_data)),
    tar_target(adaptive_lasso, fit_adaptive_lasso(train_data)),
    tar_target(glmboost, fit_glmboost(train_data)),

    # estimate the test prediction
    tar_target(test_pred_ridge,
               get_ridge_test_prediction(glmboost, test_data)),
    tar_target(test_pred_adaptive_lasso,
               get_adalasso_test_prediction(adaptive_lasso, test_data)),
    tar_target(test_pred_glmboost,
               get_glmboost_test_prediction(glmboost, test_data)),

    # extract the selected features in adpative lasso and glmboost
    tar_target(selected_features_adaptive_lasso,
               extract_adalasso_selected_features(adaptive_lasso)),
    tar_target(glmboost_varimp, extract_glmboost_varimp(glmboost)),

    # simulation experiment 2 --------------------------------------------------------------

    # simulate data with different scenarios (different N, P, p_ref, tau and rho)
    tar_target(simulated_data_N,
               generate_simulated_data(N = rep(c(30, 50, 100), each = 100),
                                       seed = 39374)),
    tar_target(simulated_data_P,
               generate_simulated_data(P = rep(c(100, 500, 5000), each = 100),
                                       seed = 44234374)),
    tar_target(simulated_data_pref,
               generate_simulated_data(p_ref = rep(c(2, 5, 10, 20), each = 100),
                                       seed = 123474)),
    tar_target(simulated_data_tau,
               generate_simulated_data(tau = rep(c(0, -2, -4, -6), each = 100),
                                       seed = 964896)),
    tar_target(simulated_data_rho,
               generate_simulated_data(rho = rep(c(0.1, 0.3, 0.5, 0.7), each = 100),
                                       seed = 6418241)),

    # stability selection using adaptive lasso
    tar_target(stabsel_adaptive_lasso_N, run_cpss_adaptive_lasso(simulated_data_N)),
    tar_target(stabsel_adaptive_lasso_P, run_cpss_adaptive_lasso(simulated_data_P)),
    tar_target(stabsel_adaptive_lasso_pref, run_cpss_adaptive_lasso(simulated_data_pref)),
    tar_target(stabsel_adaptive_lasso_tau, run_cpss_adaptive_lasso(simulated_data_tau)),
    tar_target(stabsel_adaptive_lasso_rho, run_cpss_adaptive_lasso(simulated_data_rho)),

    # stability selection using glmboost
    tar_target(stabsel_glmboost_N, run_cpss_glmboost(simulated_data_N)),
    tar_target(stabsel_glmboost_P, run_cpss_glmboost(simulated_data_P)),
    tar_target(stabsel_glmboost_pref, run_cpss_glmboost(simulated_data_pref)),
    tar_target(stabsel_glmboost_tau, run_cpss_glmboost(simulated_data_tau)),
    tar_target(stabsel_glmboost_rho, run_cpss_glmboost(simulated_data_rho))
)

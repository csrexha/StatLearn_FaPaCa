# Load packages required to define the pipeline:
library(targets)
library(crew)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Set target options:
tar_option_set(
  packages = c("data.table", "purrr", "magrittr", "caret", "glmnet", "mboost"),
  format = "qs",
  controller = crew_controller_local(workers = 2)
)

get_simulated_data <- function(n, seed) {

    set.seed(seed)

    map(1:n, function(x){
        simdat <- binary_data_generator(N = 50, P = 500, p_ref = 10,
                                        tau = -9, sigma = 10, rho = runif(10, .05, .95))
        simdat <- simdat %$% data.table(as.factor(y), x)
        setnames(simdat, c("y", paste0("X", 1:500)))

        return(simdat)
    })
}

get_stratified_cv_data <- function(data, seed) {

    set.seed(seed)

    cv_ind <- map(data, function(x) createDataPartition(x$y, p = .7, list = FALSE))

    ## standardise the training and test data ####
    vars <- paste0("X", 1:500)

    simulated_data_rescaled <- map2(
        cv_ind,
        data,
        function(x, y){
            dat <- copy(y)[, c(vars):=lapply(.SD, function(a) (a - mean(a[x]))/sd(a[x])), .SDcols=c(vars)]
            return(list(train = dat[x], test = dat[-x]))
            }
        )

    return(simulated_data_rescaled)
}

get_train_cv_data <- function(data) {
    map(data, function(x) x$train)
}

get_test_cv_data <- function(data) {
    map(data, function(x) x$test)
}

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

# Replace the target list below with your own:
list(
  tar_target(simulated_data, get_simulated_data(n = 100, seed = 119752361)),
  tar_target(scaled_cv_data, get_stratified_cv_data(simulated_data, seed = 873)),
  tar_target(train_data, get_train_cv_data(scaled_cv_data)),
  tar_target(test_data, get_test_cv_data(scaled_cv_data)),
  tar_target(adaptive_lasso, fit_adaptive_lasso(train_data)),
  tar_target(test_pred_adaptive_lasso,
             get_adalasso_test_prediction(adaptive_lasso, test_data)),
  tar_target(selected_features_adaptive_lasso,
             extract_adalasso_selected_features(adaptive_lasso))
)

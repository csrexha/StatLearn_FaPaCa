# Load packages required to define the pipeline:
library(targets)
library(crew)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Set target options:
tar_option_set(
  packages = c("data.table", "purrr", "magrittr", "caret", "glmnet", "doParallel", "foreach"),
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

# Replace the target list below with your own:
list(
  tar_target(simulated_data, get_simulated_data(n = 100, seed = 119752361)),
  tar_target(scaled_cv_data, get_stratified_cv_data(simulated_data, seed = 873)),
  tar_target(train_data, get_train_cv_data(scaled_cv_data)),
  tar_target(test_data, get_test_cv_data(scaled_cv_data))
)

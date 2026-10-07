#' Simulate a binary data set
#'
#' @description
#' Simulates a binary outcome and a matrix of predictors. The outcome depends on a latent
#' standard normal factor `f` through `P(y = 1) = linkinv(tau + sigma * f)`. The first
#' `p_ref` predictors (the relevant ones) are noisy copies of `f` with pairwise
#' correlation `rho`; the remaining `P - p_ref` predictors are independent standard
#' normal noise.
#'
#' @details
#' The outcome is redrawn until it contains at least four cases with `y = 1`, so that the
#' minor class can be split and resampled.
#'
#' @param N Number of samples.
#' @param P Total number of predictors.
#' @param p_ref Number of relevant predictors (the first `p_ref` columns of `x`). Must
#'   not exceed `P`.
#' @param tau Intercept of the linear predictor on the link scale. Negative values make
#'   `y = 1` rare.
#' @param sigma Scale (standard deviation) of the latent factor in the linear predictor.
#'   Larger values make the outcome more predictable from the relevant predictors.
#' @param rho Correlation between the relevant predictors, between 0 and 1. A vector is
#'   recycled over the `N * p_ref` values of the relevant block.
#' @param link Name of a link function accepted by [stats::make.link()]. Default
#'   `"logit"`.
#'
#' @returns A list with `y`, an integer vector of length `N` with values 0 and 1, and `x`,
#'   a numeric `N x P` matrix.
#' @export
binary_data_generator <- function(N, P, p_ref, tau, sigma, rho, link = "logit"){

    linkinv <- make.link(link)$linkinv

    f <- rnorm(N, 0, 1)
    y <- rbinom(N, 1, linkinv(tau + sigma*f))

    # regenerate if total number of minor group smaller than 4
    while(sum(y) < 4){
        y <- rbinom(N, 1, linkinv(tau + sigma*f))
    }

    x <- cbind(sqrt(1 - rho) * matrix(rnorm(N*p_ref, 0, 1), nrow = N) +
                   sqrt(rho) * matrix(rep(f, p_ref), nrow = N),
               matrix(rnorm(N * (P - p_ref), 0, 1), nrow = N))

    return(list(y = y, x = x))
}

#' Generate a data set for the simulation study
#'
#' @description
#' Generates one data set with [binary_data_generator()], converting the outcome to
#' a factor and the predictor matrix to named columns.
#'
#' @param N,P,p_ref,tau,sigma,rho Arguments passed to
#'   [binary_data_generator()].
#'
#' @returns A data.table with factor column `y` and numeric predictor columns
#'   `X1`, ..., `XP`.
#' @export
get_simulated_data <- function(N, P, p_ref, tau, sigma, rho) {

    simulated_data <- binary_data_generator(N, P, p_ref, tau, sigma, rho) |>
        as.data.table()

    simulated_data[, y := as.factor(y)]
    setnames(simulated_data, old = paste0("x.V", 1:P), new = paste0("X", 1:P))

    return(simulated_data)
}

#' Split a simulated data set into standardised training and test sets
#'
#' @description
#' Splits a single data set 70/30 into a training and a test set, stratified by
#' `y` using [caret::createDataPartition()]. All predictors (every column except
#' `y`) are then standardised using the mean and standard deviation of the
#' training set, and the same scaling is applied to the test set so that no
#' information from the test data enters the centering or scaling.
#'
#' @param data A data.table returned by [get_simulated_data()], with a factor
#'   column `y` and numeric predictor columns.
#'
#' @returns A list with two elements: `train` and `test`, each a data.table.
#' @export
get_stratified_cv_data <- function(data) {

    cv_ind <- createDataPartition(data$y, p = 0.7, list = FALSE)

    # standardise the training and test data
    vars <- setdiff(names(data), "y")
    data_rescaled <- copy(data)[, c(vars) := lapply(.SD,
        function(a) (a - mean(a[cv_ind])) / sd(a[cv_ind])), .SDcols = c(vars)]

    return(list(train = data_rescaled[cv_ind], test = data_rescaled[-cv_ind]))
}

#' Extract the training set
#'
#' @description
#' Returns the training component from the output of [get_stratified_cv_data()].
#'
#' @param data A list returned by [get_stratified_cv_data()] with `train` and
#'   `test` elements.
#'
#' @returns The training data.table.
#' @export
get_train_cv_data <- function(data) {
    return(data$train)
}

#' Extract the test set
#'
#' @description
#' Returns the test component from the output of [get_stratified_cv_data()].
#'
#' @param data A list returned by [get_stratified_cv_data()] with `train` and
#'   `test` elements.
#'
#' @returns The test data.table.
#' @export
get_test_cv_data <- function(data) {
    return(data$test)
}

#' Read the raw data file
#'
#' @description
#' Reads a CSV file with [data.table::fread()] and converts character columns to factors,
#' which is required for imputation with [missForest::missForest()].
#'
#' @param file_path Path to the raw CSV file.
#'
#' @returns A data.table with the raw data, character columns converted to factors.
#' @export
read_raw_data <- function(file_path) {
    dt <- data.table::fread(file_path)

    # Ensure character columns are converted to factors for missForest
    char_cols <- names(dt)[sapply(dt, is.character)]
    if (length(char_cols) > 0) {
        dt[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]
    }

    return(dt)
}


#' Impute missing values with missForest
#'
#' @description
#' Imputes missing values non-parametrically with random forests
#' ([missForest::missForest()]). Continuous and categorical variables are handled
#' automatically.
#'
#' @param data A data.table with missing values (`NA`). Categorical variables must be
#'   factors (see [read_raw_data()]).
#' @param features Character vector of the columns to impute. Only these columns are used
#'   for the imputation and returned.
#' @param maxiter Maximum number of iterations if the stopping criterion is not met.
#'   Default 10.
#' @param ntree Number of trees in each forest. Default 100.
#'
#' @returns A data.frame with the `features` columns and no missing values.
#' @export
impute_data <- function(data, features, maxiter = 10, ntree = 100) {

    # Convert the selected features to a data frame for missForest
    xmis <- as.data.frame(data[, ..features])

    # missForest returns a list; $ximp contains the imputed data frame
    imputed_result <- missForest::missForest(
        xmis = xmis,
        maxiter = maxiter,
        ntree = ntree
    )
    return(imputed_result$ximp)
}


#' Export imputed data to CSV
#'
#' @description
#' Writes the imputed data to a CSV file with [data.table::fwrite()] and returns the path,
#' so that the file can be tracked by a `targets` pipeline.
#'
#' @param data A data.frame or data.table with the imputed data.
#' @param output_path Path of the CSV file to write.
#'
#' @returns `output_path`, as a character string.
#' @export
export_imputed_data <- function(data, output_path) {
    data.table::fwrite(data, output_path)
    return(output_path) # Return file path for target tracking
}

#' Read the imputed data file
#'
#' @description
#' Reads the CSV file written by [export_imputed_data()] with [data.table::fread()].
#'
#' @param file_path Path to the imputed CSV file.
#'
#' @returns A data.table with the imputed data.
#' @export
read_imputed_data <- function(file_path) {
    data.table::fread(file_path)
}

#' Create repeated cross-validation splits with centred features
#'
#' @description
#' Creates 10 repeats of stratified 4-fold cross-validation ([caret::createMultiFolds()],
#' stratified by `Status`) and centres the selected features with the mean of the training
#' part of each fold. The features are centred only, not scaled.
#'
#' @details
#' The training mean is subtracted from both the training and the test part of a fold, so
#' no information from the test part enters the centring. The result has 40 folds named
#' `Fold1.Rep01`, ..., `Fold4.Rep10`.
#'
#' @param data A data.table with a factor column `Status` and the features to be centred.
#' @param features_rescale Character vector of the columns to centre.
#'
#' @returns A list with two elements:
#' \describe{
#'   \item{train}{A list of 40 training data.tables, one per fold.}
#'   \item{test}{A list of 40 test data.tables, one per fold, in the same order.}
#' }
#' @export
get_fapaca_cv_data <- function(data, features_rescale) {

    # Generate 10 repeats of 4-fold stratified cross-validation indices based on 'Status'
    cv_ind <- caret::createMultiFolds(data$Status, 4, 10)

    # Rescale specified features using fold-specific training means
    data_cv_rescaled <- purrr::map(
        cv_ind,
        function(x) {
            # Extract unique training row indices for the current fold
            ind <- unique(x)

            # Create a deep copy of data and mean-center selected features using training fold mean
            data_rescaled <- copy(data)[, c(features_rescale) := lapply(.SD, function(y) (y - mean(y[ind]))),
                                        .SDcols = c(features_rescale)]

            # Split rescaled dataset into train and test folds
            return(list(train = data_rescaled[ind], test = data_rescaled[-ind]))
        }
    )

    # Extract all training and testing splits into separate nested lists
    data_train <- purrr::map(data_cv_rescaled, function(x) x$train)
    data_test <- purrr::map(data_cv_rescaled, function(x) x$test)

    # Return structured train/test cross-validation dataset lists
    return(list(train = data_train, test = data_test))
}


#' Select features from the cross-validation training sets
#'
#' @description
#' Keeps only the given feature columns in every training fold.
#'
#' @param data A list with an element `train`, a list of data.tables, as returned by
#'   [get_fapaca_cv_data()].
#' @param features Character vector of the columns to keep.
#'
#' @returns A list of data.tables, one per training fold, with only the columns `features`.
#' @export
get_fapaca_train_data <- function(data, features) {

    # Validate that input list contains a valid 'train' element of type list
    stopifnot("Input data must be a list of data.tables" = is.list(data$train))

    # Subset specified features from each data.table in the training folds list
    purrr::map(data$train, function(x) x[, ..features])
}

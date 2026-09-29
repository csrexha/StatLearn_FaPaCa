# FaPaCa study pipeline: imputation of the proteomics data, cross-validation splits and
# glmboost / gamboost models.
# Run with source("make.R") (needs data/raw_data.csv) or with
# targets::tar_make(script = "analyses/fapaca_study.R",
#                   store = "outputs/fapaca_study")

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

# Specify the features to impute, e.g. age, sex, bmi, proteins etc.
features_impute <- c("age", "sex")

# Specify the features to be rescaled, e.g. age, bmi, proteins etc.
features_rescale <- c("age")

# Specify the features for modelling, e.g. age, sex, bmi, proteins etc.
features_model <- c("age", "sex")

# Pipeline
list(
    # Data preparation ---------------------------------------------------------------------

    # 1.1. Track changes of the raw data file (place the raw data into the data/ folder)
    tar_target(
        name = raw_file,
        command = "data/raw_data.csv",
        format = "file"
    ),

    # 1.2. Read the raw data
    tar_target(
        name = raw_data,
        command = read_raw_data(raw_file)
    ),

    # 1.3. Impute missing values using missForest
    tar_target(
        name = imputed_data,
        command = impute_data(raw_data,
                              features = features_impute,
                              maxiter = 10,
                              ntree = 100)
    ),

    # 1.4. Save the imputed data to disk and track the output file
    tar_target(
        name = output_file,
        command = export_imputed_data(imputed_data, "data/imputed_data.csv"),
        format = "file"
    ),

    # 1.5. Read the exported imputed data back into the pipeline
    tar_target(
        name = reloaded_imputed_data,
        command = read_imputed_data(output_file)
    ),

    # 1.6. Prepare the data for cross-validation
    tar_target(
        name = fapaca_cv_data,
        command = get_fapaca_cv_data(reloaded_imputed_data,
                                     features_rescale = features_rescale)
    ),

    # 1.7. Prepare the training data for model fitting
    tar_target(
        name = fapaca_train_data,
        command = get_fapaca_train_data(fapaca_cv_data, features = features_model)
    ),


    # Model fitting ------------------------------------------------------------------------

    # 2.1. Fit glmboost models to the training data
    tar_target(
        name = fapaca_glmboost,
        command = fit_glmboost(fapaca_train_data,
                               binary_class = "Status",
                               features = features_model,
                               seed = 1006613948)
    ),

    # 2.2. Fit gamboost models to the training data
    tar_target(
        name = fapaca_gamboost,
        command = fit_gamboost(fapaca_train_data,
                               binary_class = "Status",
                               features = features_model,
                               seed = 1006613948)
    )
)

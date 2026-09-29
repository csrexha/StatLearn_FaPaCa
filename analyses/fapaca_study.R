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

# specify the features to impute, e.g. age, sex, bmi, proteins etc.
features_impute <- c("age", "sex")

# specify the features for modelling, e.g. age, sex, bmi, proteins etc.
features_model <- c("age", "sex")


list(
    # 1. Track raw input file changes
    tar_target(
        name = raw_file,
        command = "data/raw_data.csv", # place the raw data into the data/ folder and specify the path here
        format = "file"
    ),

    # 2. Read raw dataset
    tar_target(
        name = raw_data,
        command = read_raw_data(raw_file)
    ),

    # 3. Impute missing values using missForest
    tar_target(
        name = imputed_data,
        command = impute_data(raw_data,
                              features = features_impute,
                              maxiter = 10,
                              ntree = 100)
    ),

    # 4. Save imputed data to disk and track output file
    tar_target(
        name = output_file,
        command = export_imputed_data(imputed_data, "data/imputed_data.csv"),
        format = "file"
    ),

    # 5. Read the exported imputed data file back into the pipeline
    tar_target(
        name = reloaded_imputed_data,
        command = read_imputed_data(output_file)
    ),

    # 6. Prepare data for cross-validation
    tar_target(
        name = fapaca_cv_data,
        command = get_fapaca_cv_data(reloaded_imputed_data,
                                     features_rescale = features_rescale)
    ),

    # 7. Prepare training data for model fitting
    tar_target(
        name = fapaca_train_data,
        command = get_fapaca_train_data(fapaca_cv_data,
                                        features = features_model)
    ),

    # 8. Fit glmboost model to training data
    tar_target(
        name = fapaca_glmboost,
        command = fit_glmboost(fapaca_train_data,
                               binary_class = "Status",
                               features = features_model,
                               seed = 1006613948)
    ),

    # 9. Fit gamboost model to training data
    tar_target(
        name = fapaca_gamboost,
        command = fit_gamboost(fapaca_train_data,
                               binary_class = "Status",
                               features = features_model,
                               seed = 1006613948)
    )






)

# Simulation study pipeline: Experiment 1 (ridge, adaptive lasso and glmboost on one
# simulated data set) and Experiment 2 (stability selection with varying parameters).
# Run with source("make.R") or with
# targets::tar_make(script = "analyses/simulation_study_pipeline.R",
#                   store = "outputs/simulation_study")

# Load packages required to define the pipeline:
library(targets)
library(crew)
library(tarchetypes)
library(data.table)

# Run the R scripts in the R/ folder with your custom functions:
tar_source()

# Set target options:
tar_option_set(
    packages = c("data.table", "purrr", "magrittr", "caret", "glmnet", "mboost", "stabs"),
    format = "qs",
    controller = crew_controller_local(workers = 4),
    seed = 119752361
)

# Name of the binary outcome column in the simulated data
binary_class <- "y"

scenario_grid <- make_simulation_senario(
    default = data.table(N = 50, P = 500, p_ref = 10, tau = -9, sigma = 10,
                         rho=list(runif(10, 0.05, 0.95))),
    values = list(
        N = c(30L, 50L, 100L),
        P = c(100L, 500L, 5000L),
        p_ref = c(2L, 5L, 10L, 20L),
        tau = c(0L, -2L, -4L, -6L),
        rho = c(0.1, 0.3, 0.5, 0.7)
    )
)

# Pipeline
list(
    # Simulation experiment 1 --------------------------------------------------------------

    # 1.1. Generate simulation data and rescale the data for cross-validation
    tar_rep(
        name = e1_simulated_data,
        command = list(data = get_simulated_data(
            N = 50,
            P = 10,
            p_ref = 10,
            tau = -9,
            sigma = 10,
            rho=runif(10, 0.05, 0.95)
            )
            ),
        batches = 2,
        reps = 2,
        iteration = "list"
    ),
    
    tar_rep2(
        name = e1_scaled_cv_data,
        command = list(
            data = get_stratified_cv_data(e1_simulated_data$data)
            ),
        e1_simulated_data,
        iteration = "list"
    ),

    # 1.2. Split the scaled data into training and test sets for cross-validation
    tar_rep2(
        name = e1_train_test_data,
        command = list(
            train = get_train_cv_data(e1_scaled_cv_data$data),
            test = get_test_cv_data(e1_scaled_cv_data$data)
            ),
        e1_scaled_cv_data,
        iteration = "list"
    ),

    # 1.3. Fit ridge, adaptive lasso and glmboost models to the training data
    tar_rep2(
        name = e1_ridge,
        command = list(model = fit_ridge(
            e1_train_test_data$train,
            binary_class = binary_class
        )),
        e1_train_test_data,
        iteration = "list"
    ),
    tar_rep2(
        name = e1_adaptive_lasso,
        command = list(model = fit_adaptive_lasso(
            e1_train_test_data$train,
            binary_class = binary_class
        )),
        e1_train_test_data,
        iteration = "list"
    ),
    tar_rep2(
        name = e1_glmboost,
        command = list(model = fit_glmboost(
            e1_train_test_data$train,
            binary_class = binary_class
        )),
        e1_train_test_data,
        iteration = "list"
    ),

    # 1.4. Get the test set predictions for ridge, adaptive lasso and glmboost
    tar_rep2(
        name = e1_test_pred_ridge,
        command = get_ridge_test_prediction(
            e1_ridge$model,
            e1_train_test_data$test
            ),
        e1_ridge,
        e1_train_test_data
    ),

    tar_rep2(
        name = e1_test_pred_adaptive_lasso,
        command = get_adalasso_test_prediction(
            e1_adaptive_lasso$model,
            e1_train_test_data$test
            ),
        e1_adaptive_lasso,
        e1_train_test_data
    ),

    tar_rep2(
        name = e1_test_pred_glmboost,
        command = get_glmboost_test_prediction(
            e1_glmboost$model,
            e1_train_test_data$test
            ),
        e1_glmboost,
        e1_train_test_data
    ),

    # 1.5. Extract the selected features (adaptive lasso) and the variable importance (glmboost)
    tar_rep2(
        name = e1_selected_features_adaptive_lasso,
        command = extract_adalasso_selected_features(e1_adaptive_lasso$model),
        e1_adaptive_lasso
    ),
    tar_rep2(
        name = e1_glmboost_varimp,
        command = extract_glmboost_varimp(e1_glmboost$model) ,
        e1_glmboost
    ),

    tar_target(
        name = e1_glmboost_varimp_summary,
        command = summarise_glmboost_varimp(e1_glmboost_varimp)
    ),


    # Simulation experiment 2 --------------------------------------------------------------

    tar_map(
        values = scenario_grid,
        names = tidyselect::any_of("scenario"),

        # 2.1. Generate simulation data with varying parameters N, P, p_ref, tau and rho
        tar_rep(
            name = e2_simulated_data,
            command = list( data = get_simulated_data(
                N = N,
                P = P,
                p_ref = p_ref,
                tau = tau,
                sigma = sigma,
                rho = rho
            )),
            batches = 5,
            reps = 20,
            iteration = "list"
        ),

        # 2.2. Stability selection using adaptive lasso and glmboost
        tar_rep2(
            name = e2_stabsel_adaptive_lasso,
            command = cpss_adaptive_lasso(
                e2_simulated_data$data,
                binary_class = binary_class,
                features = NULL,
                q = 10,
                PFER = 2
            ),
            e2_simulated_data,
            iteration = "list"
        ),

        tar_rep2(
            name = e2_stabsel_glmboost,
            command = cpss_glmboost(
                e2_simulated_data$data,
                binary_class = binary_class,
                features = NULL,
                iter = 500,
                q = 10,
                PFER = 2
            ),
            e2_simulated_data,
            iteration = "list"
        )
    )
)

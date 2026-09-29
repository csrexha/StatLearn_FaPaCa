# Simulation study pipeline: Experiment 1 (ridge, adaptive lasso and glmboost on one
# simulated data set) and Experiment 2 (stability selection with varying parameters).
# Run with source("make.R") or with
# targets::tar_make(script = "analyses/simulation_study_pipeline.R",
#                   store = "outputs/simulation_study")

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

# Name of the binary outcome column in the simulated data
binary_class <- "y"

# Pipeline
list(
    # Simulation experiment 1 --------------------------------------------------------------

    # 1.1. Generate simulation data and rescale the data for cross-validation
    tar_target(
        name = e1_simulated_data,
        command = get_simulated_data(n = 100, seed = 119752361)
    ),
    tar_target(
        name = e1_scaled_cv_data,
        command = get_stratified_cv_data(e1_simulated_data, seed = 873)
    ),

    # 1.2. Split the scaled data into training and test sets for cross-validation
    tar_target(
        name = e1_train_data,
        command = get_train_cv_data(e1_scaled_cv_data)
    ),
    tar_target(
        name = e1_test_data,
        command = get_test_cv_data(e1_scaled_cv_data)
    ),

    # 1.3. Fit ridge, adaptive lasso and glmboost models to the training data
    tar_target(
        name = e1_ridge,
        command = fit_ridge(e1_train_data, binary_class = binary_class)
    ),
    tar_target(
        name = e1_adaptive_lasso,
        command = fit_adaptive_lasso(e1_train_data, binary_class = binary_class)
    ),
    tar_target(
        name = e1_glmboost,
        command = fit_glmboost(e1_train_data, binary_class = binary_class)
    ),

    # 1.4. Get the test set predictions for ridge, adaptive lasso and glmboost
    tar_target(
        name = e1_test_pred_ridge,
        command = get_ridge_test_prediction(e1_ridge, e1_test_data)
    ),
    tar_target(
        name = e1_test_pred_adaptive_lasso,
        command = get_adalasso_test_prediction(e1_adaptive_lasso, e1_test_data)
    ),
    tar_target(
        name = e1_test_pred_glmboost,
        command = get_glmboost_test_prediction(e1_glmboost, e1_test_data)
    ),

    # 1.5. Extract the selected features (adaptive lasso) and the variable importance (glmboost)
    tar_target(
        name = e1_selected_features_adaptive_lasso,
        command = extract_adalasso_selected_features(e1_adaptive_lasso)
    ),
    tar_target(
        name = e1_glmboost_varimp,
        command = extract_glmboost_varimp(e1_glmboost)
    ),


    # Simulation experiment 2 --------------------------------------------------------------

    # 2.1. Generate simulation data with varying parameters N, P, p_ref, tau and rho
    tar_target(
        name = e2_simulated_data_N,
        command = get_simulated_data(N = rep(c(30, 50, 100), each = 100),
                                     seed = 39374)
    ),
    tar_target(
        name = e2_simulated_data_P,
        command = get_simulated_data(P = rep(c(100, 500, 5000), each = 100),
                                     seed = 44234374)
    ),
    tar_target(
        name = e2_simulated_data_pref,
        command = get_simulated_data(p_ref = rep(c(2, 5, 10, 20), each = 100),
                                     seed = 123474)
    ),
    tar_target(
        name = e2_simulated_data_tau,
        command = get_simulated_data(tau = rep(c(0, -2, -4, -6), each = 100),
                                     seed = 964896)
    ),
    tar_target(
        name = e2_simulated_data_rho,
        command = get_simulated_data(rho = rep(c(0.1, 0.3, 0.5, 0.7), each = 100),
                                     seed = 6418241)
    ),

    # 2.2. Stability selection using adaptive lasso
    tar_target(
        name = e2_stabsel_adaptive_lasso_N,
        command = run_cpss_adaptive_lasso(e2_simulated_data_N,
                                          binary_class = binary_class,
                                          features = NULL,
                                          q = 10,
                                          PFER = 2,
                                          seed = 87435)
    ),
    tar_target(
        name = e2_stabsel_adaptive_lasso_P,
        command = run_cpss_adaptive_lasso(e2_simulated_data_P,
                                          binary_class = binary_class,
                                          features = NULL,
                                          q = 10,
                                          PFER = 2,
                                          seed = 87435)
    ),
    tar_target(
        name = e2_stabsel_adaptive_lasso_pref,
        command = run_cpss_adaptive_lasso(e2_simulated_data_pref,
                                          binary_class = binary_class,
                                          features = NULL,
                                          q = 10,
                                          PFER = 2,
                                          seed = 87435)
    ),
    tar_target(
        name = e2_stabsel_adaptive_lasso_tau,
        command = run_cpss_adaptive_lasso(e2_simulated_data_tau,
                                          binary_class = binary_class,
                                          features = NULL,
                                          q = 10,
                                          PFER = 2,
                                          seed = 87435)
    ),
    tar_target(
        name = e2_stabsel_adaptive_lasso_rho,
        command = run_cpss_adaptive_lasso(e2_simulated_data_rho,
                                          binary_class = binary_class,
                                          features = NULL,
                                          q = 10,
                                          PFER = 2,
                                          seed = 87435)
    ),

    # 2.3. Stability selection using glmboost
    tar_target(
        name = e2_stabsel_glmboost_N,
        command = run_cpss_glmboost(e2_simulated_data_N,
                                    binary_class = binary_class,
                                    features = NULL,
                                    iter = 500,
                                    q = 10,
                                    PFER = 2,
                                    seed = 29374)
    ),
    tar_target(
        name = e2_stabsel_glmboost_P,
        command = run_cpss_glmboost(e2_simulated_data_P,
                                    binary_class = binary_class,
                                    features = NULL,
                                    iter = 500,
                                    q = 10,
                                    PFER = 2,
                                    seed = 29374)
    ),
    tar_target(
        name = e2_stabsel_glmboost_pref,
        command = run_cpss_glmboost(e2_simulated_data_pref,
                                    binary_class = binary_class,
                                    features = NULL,
                                    iter = 500,
                                    q = 10,
                                    PFER = 2,
                                    seed = 29374)
    ),
    tar_target(
        name = e2_stabsel_glmboost_tau,
        command = run_cpss_glmboost(e2_simulated_data_tau,
                                    binary_class = binary_class,
                                    features = NULL,
                                    iter = 500,
                                    q = 10,
                                    PFER = 2,
                                    seed = 29374)
    ),
    tar_target(
        name = e2_stabsel_glmboost_rho,
        command = run_cpss_glmboost(e2_simulated_data_rho,
                                    binary_class = binary_class,
                                    features = NULL,
                                    iter = 500,
                                    q = 10,
                                    PFER = 2,
                                    seed = 29374)
    )
)

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

# pipelines
list(
    # simulation experiment 1 --------------------------------------------------------------

    # simulate data
    tar_target(e1_simulated_data, get_simulated_data(n = 100, seed = 119752361)),
    tar_target(e1_scaled_cv_data, get_stratified_cv_data(e1_simulated_data, seed = 873)),
    tar_target(e1_train_data, get_train_cv_data(e1_scaled_cv_data)),
    tar_target(e1_test_data, get_test_cv_data(e1_scaled_cv_data)),

    # fit ridge, adaptive lasso, glmboost regression to training data
    tar_target(e1_ridge, fit_ridge(e1_train_data)),
    tar_target(e1_adaptive_lasso, fit_adaptive_lasso(e1_train_data)),
    tar_target(e1_glmboost, fit_glmboost(e1_train_data)),

    # estimate the test prediction
    tar_target(e1_test_pred_ridge,
               get_ridge_test_prediction(e1_ridge, e1_test_data)),
    tar_target(e1_test_pred_adaptive_lasso,
               get_adalasso_test_prediction(e1_adaptive_lasso, e1_test_data)),
    tar_target(e1_test_pred_glmboost,
               get_glmboost_test_prediction(e1_glmboost, e1_test_data)),

    # extract the selected features in adpative lasso and glmboost
    tar_target(e1_selected_features_adaptive_lasso,
               extract_adalasso_selected_features(e1_adaptive_lasso)),
    tar_target(e1_glmboost_varimp, extract_glmboost_varimp(e1_glmboost)),

    # simulation experiment 2 --------------------------------------------------------------

    # simulate data with different scenarios (different N, P, p_ref, tau and rho)
    tar_target(e2_simulated_data_N,
               get_simulated_data(N = rep(c(30, 50, 100), each = 100),
                                  seed = 39374)),
    tar_target(e2_simulated_data_P,
               get_simulated_data(P = rep(c(100, 500, 5000), each = 100),
                                  seed = 44234374)),
    tar_target(e2_simulated_data_pref,
               get_simulated_data(p_ref = rep(c(2, 5, 10, 20), each = 100),
                                  seed = 123474)),
    tar_target(e2_simulated_data_tau,
               get_simulated_data(tau = rep(c(0, -2, -4, -6), each = 100),
                                  seed = 964896)),
    tar_target(e2_simulated_data_rho,
               get_simulated_data(rho = rep(c(0.1, 0.3, 0.5, 0.7), each = 100),
                                  seed = 6418241)),

    # stability selection using adaptive lasso
    tar_target(e2_stabsel_adaptive_lasso_N, run_cpss_adaptive_lasso(e2_simulated_data_N)),
    tar_target(e2_stabsel_adaptive_lasso_P, run_cpss_adaptive_lasso(e2_simulated_data_P)),
    tar_target(e2_stabsel_adaptive_lasso_pref, run_cpss_adaptive_lasso(e2_simulated_data_pref)),
    tar_target(e2_stabsel_adaptive_lasso_tau, run_cpss_adaptive_lasso(e2_simulated_data_tau)),
    tar_target(e2_stabsel_adaptive_lasso_rho, run_cpss_adaptive_lasso(e2_simulated_data_rho)),

    # stability selection using glmboost
    tar_target(e2_stabsel_glmboost_N, run_cpss_glmboost(e2_simulated_data_N)),
    tar_target(e2_stabsel_glmboost_P, run_cpss_glmboost(e2_simulated_data_P)),
    tar_target(e2_stabsel_glmboost_pref, run_cpss_glmboost(e2_simulated_data_pref)),
    tar_target(e2_stabsel_glmboost_tau, run_cpss_glmboost(e2_simulated_data_tau)),
    tar_target(e2_stabsel_glmboost_rho, run_cpss_glmboost(e2_simulated_data_rho))
)

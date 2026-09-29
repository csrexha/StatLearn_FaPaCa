# analyses

The two `targets` pipelines of this project. Both are run by
[`make.R`](../make.R), or individually with `targets::tar_make()`.

| Script | Purpose | Store |
| :----- | :------ | :---- |
| [`simulation_study_pipeline.R`](simulation_study_pipeline.R) | Simulation study. Experiment 1: ridge, adaptive lasso and glmboost on one simulated data set, with test-set predictions, selected features and variable importance. Experiment 2: stability selection (adaptive lasso, glmboost) while varying `N`, `P`, `p_ref`, `tau` and `rho`. | `outputs/simulation_study` |
| [`fapaca_study.R`](fapaca_study.R) | FaPaCa proteomics study: read, impute (`missForest`), prepare the cross-validation data and fit glmboost and gamboost. Requires `data/raw_data.csv`. | `outputs/fapaca_study` |

The functions used by the pipelines are in [`R/`](../R).

Note: `fapaca_study.R` still contains placeholder feature sets and an
undefined `features_rescale`; set them before running it on real data.

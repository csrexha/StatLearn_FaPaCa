#' StatLearn_FaPaCa: A Research Compendium
#'
#' @description
#' Runs the simulation study and the FaPaCa study of "Proteomics biomarker
#' discovery for individualized prevention of familial pancreatic cancer using
#' statistical learning" as two targets pipelines.
#'
#' @author Chung Shing Rex Ha \email{hachungshingrex@gmail.com}
#'
#' @date 2026/09/24



## Install Dependencies (listed in rproject.toml) ----

if (!nzchar(Sys.which("rv"))) {
    stop("rv is not installed. Install it from https://github.com/A2-ai/rv")
}
system2("rv", c("sync", "-c", here::here("rproject.toml")))


## Load Project Addins (R Functions and Packages) ----

devtools::load_all(here::here())


## Global Variables ----

run_simulation <- TRUE
run_fapaca     <- file.exists(here::here("data", "raw_data.csv"))  # data are not distributed


## Run Project ----

if (run_simulation) {
    targets::tar_make(script = here::here("analyses", "simulation_study_pipeline.R"),
                      store  = here::here("outputs", "simulation_study"))
}

if (run_fapaca) {
    targets::tar_make(script = here::here("analyses", "fapaca_study.R"),
                      store  = here::here("outputs", "fapaca_study"))
} else {
    message("data/raw_data.csv not found: skipping the FaPaCa study")
}

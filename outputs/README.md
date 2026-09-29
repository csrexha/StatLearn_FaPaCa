# outputs

Results produced by the `targets` pipelines in [`analyses/`](../analyses).
This folder is not tracked by git (see `.gitignore`) because the results are
large and can be rebuilt.

- `simulation_study/`: `targets` store of `analyses/simulation_study_pipeline.R`
  (simulation Experiments 1 and 2)
- `fapaca_study/`: `targets` store of `analyses/fapaca_study.R`
  (requires `data/raw_data.csv`)

To recreate the results, run `source("make.R")` from the project root.
Experiment 2 of the simulation study (stability selection) takes a long time.
Use `targets::tar_read(<name>, store = "outputs/simulation_study")` to load a
result.

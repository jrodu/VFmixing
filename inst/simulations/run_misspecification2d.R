# Misspecification 2-D simulation study
#
# Runs runSimulationStudy_Misspecification2D() once per seed and saves each replicate to its
# own file, so a long run can be stopped and restarted without redoing finished
# seeds. Then combines every replicate's Results table into one data frame.
#
# Find this script with:
#   system.file("simulations", "run_misspecification2d.R", package = "VFmixing")
#
# Each replicate covers 3 error levels x 3 trajectory complexities x 4 method
# pairs and takes roughly a minute, so 1000 seeds is on the order of 15 hours
# on one core. Seeds are independent, so they can be split across machines or
# cores (for example, run seeds 1:250 on one and 251:500 on another).

library(VFmixing)

seeds = 1:1000
out_dir = "misspecification2d_results"   # one .rds file per seed is written here
keep_fits = FALSE                       # TRUE keeps every fitted path and trajectory (~25 MB per seed)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for(seed in seeds){

  file = file.path(out_dir, sprintf("misspecification2d_seed%04d.rds", seed))

  if(file.exists(file)){
    next
  }

  run = runSimulationStudy_Misspecification2D(seed = seed, keep_fits = keep_fits)
  saveRDS(run, file)

  message(sprintf("seed %d done (%s)", seed, format(Sys.time(), "%H:%M:%S")))

}

# Combine the results of every finished seed
files = list.files(out_dir, pattern = "^misspecification2d_seed\\d+\\.rds$", full.names = TRUE)
results = do.call(rbind, lapply(files, function(f){readRDS(f)$Results}))

saveRDS(results, file.path(out_dir, "misspecification2d_results_all.rds"))

# Example summary: mean of each metric by error level, complexity and method
metric_cols = setdiff(names(results)[vapply(results, is.numeric, logical(1))], c("seed", "error_sd", "n_true_regions"))
summary_table = aggregate(results[, metric_cols], by = results[, c("error_level", "complexity", "method")], FUN = mean, na.rm = TRUE)
print(summary_table)

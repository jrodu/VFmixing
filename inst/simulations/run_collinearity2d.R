# Collinearity2D simulation study
#
# Runs runSimulationStudy_Collinearity2D() once per seed and saves each replicate to its
# own file, so a long run can be stopped and restarted without redoing finished
# seeds. Then combines every replicate's Results table into one data frame.
#
# Find this script with:
#   system.file("simulations", "run_collinearity2d.R", package = "VFmixing")
#
# Each seed draws one angle theta and fits the low error, medium complexity
# setting with 2x BIC + bootstrap (B = 20).
# Each seed takes roughly 1 minute on one core. Seeds are independent,
# so they can be split across machines or cores (for example, seeds 1:50 on one
# and 51:100 on another, or one block of seeds per task of a job array).

library(VFmixing)

seeds = 1:100
out_dir = "collinearity2d_results"   # one .rds file per seed is written here
keep_fits = FALSE        # TRUE also keeps every fitted path, trajectory and variance piece

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for(seed in seeds){

  file = file.path(out_dir, sprintf("collinearity2d_seed%04d.rds", seed))

  if(file.exists(file)){
    next
  }

  run = runSimulationStudy_Collinearity2D(seed = seed, keep_fits = keep_fits)
  saveRDS(run, file)

  message(sprintf("seed %d done (%s)", seed, format(Sys.time(), "%H:%M:%S")))

}

# Combine the results of every finished seed
files = list.files(out_dir, pattern = "^collinearity2d_seed[0-9]+[.]rds$", full.names = TRUE)
results = do.call(rbind, lapply(files, function(f){readRDS(f)$Results}))

saveRDS(results, file.path(out_dir, "collinearity2d_results_all.rds"))

# Example summary: mean of each metric by error level, complexity and variant (to study
# theta, model the metrics in results against results$theta instead)
metric_cols = setdiff(names(results)[vapply(results, is.numeric, logical(1))], c("seed", "error_sd", "n_true_regions", "theta"))
summary_table = aggregate(results[, metric_cols], by = results[, c("error_level", "complexity", "variant")], FUN = mean, na.rm = TRUE)
print(summary_table)

# Monte Carlo standard error of the mean coverage of the first trajectory, by the same groups
mc_se = aggregate(results["traj_coverage_Trajectory1"], by = results[, c("error_level", "complexity", "variant")],
                  FUN = function(x){stats::sd(x, na.rm = TRUE)/sqrt(sum(!is.na(x)))})
print(mc_se)

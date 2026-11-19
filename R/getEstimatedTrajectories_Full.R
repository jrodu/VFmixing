#' Full Method: Estimated Path and Trajectories with Optional Bootstrap
#'
#' @description
#' Runs the whole method on the data: a path is estimated with a covariance
#' structure (FullFours by default), then the trajectories are estimated from it
#' with [getEstimatedTrajectories_OneShot()]. Both node searches use the "2x BIC"
#' undersmoothing rule by default. Optionally, a bootstrap adds the uncertainty
#' from the node selection to OneShot's covariance through the law of total
#' variance.
#'
#' @details
#' **Undersmoothing (`undersmooth`).** Each node search returns the candidate
#' with twice as many regions as the one its criterion selects. This reduces the
#' bias of the path's velocity and of the trajectory spline, so the trajectory
#' intervals are centred closer to the truth. See [getEstimatedPaths_FullFours()].
#'
#' **Refining the trajectory (`traj_refine`).** With low positional error and
#' complex trajectories, most of the remaining bias comes from trajectory nodes
#' that sit away from where the true trajectories change. The path can only
#' support so many regions, but the trajectory can be split further: each
#' selected trajectory region is split into `traj_refine` equal regions and refit
#' (see the `refine` argument of [getEstimatedTrajectories_OneShot()]). This is a
#' heuristic that trades wider intervals for less bias. In the bootstrap, every
#' replicate's trajectory is refined the same way.
#'
#' **Bootstrap (`bootstrap`).** OneShot's covariance treats the selected nodes as
#' fixed. The nodes are chosen from the data, so they vary from data set to data
#' set, and that variability is missing from it. By the law of total variance,
#' \deqn{Var(\hat\beta) = E[Var(\hat\beta \mid nodes)] + Var(E[\hat\beta \mid nodes]).}
#' The first term is OneShot's covariance (the within-node part). The second
#' (the selection part) is estimated with a residual bootstrap:
#' 1. The positions of the fitted path at the data times are taken as the truth,
#'    with the centred residuals as the noise. The residuals are scaled by
#'    \eqn{\sqrt{n/(n - p)}}, with \eqn{p} the path's free parameters per
#'    dimension, so the resampled errors have the same variance as the path's
#'    own estimate, \eqn{RSS/(n - p)}: least-squares residuals are smaller than
#'    the errors.
#' 2. For each of `B` replicates, the residuals are resampled with replacement
#'    and added to the fitted positions, and the whole method (both node
#'    searches) is rerun, giving that replicate's path and trajectory nodes.
#' 3. The replicate's nodes are then fit, without any noise, to the fitted
#'    positions. This is the conditional mean \eqn{m_b = E[\hat\beta \mid nodes]}
#'    for those nodes, a C1 cubic spline on the replicate's own trajectory nodes.
#' 4. The selection part at each time is the sample variance (divisor
#'    \eqn{B - 1}) of the conditional means there,
#'    \eqn{V_{sel}(t) = \frac{1}{B - 1}\sum_b (m_b(t) - \bar m(t))^2}.
#'    The conditional means are kept as they are, with no projection onto the
#'    original nodes.
#'
#' The pointwise intervals are then
#' \eqn{\hat\beta(t) \pm z \sqrt{s_w^2(t) + V_{sel}(t)}}, with \eqn{s_w^2(t)} the
#' within-node variance from OneShot's covariance. The returned trajectory keeps
#' OneShot's covariance as `Cov` and carries the conditional means as
#' `SelectionMeans`; [visualizeEstTraj()], [coverageRateTraj()] and
#' [avgCIWidthTraj()] add the selection part whenever `SelectionMeans` is
#' present. Set it to `NULL` for intervals from OneShot's covariance alone.
#'
#' **Parallel replicates (`cores`).** The replicates are independent, so with
#' `cores` greater than 1 they run in parallel. Every replicate's resampling is
#' drawn before any of them run, so the result is the same for any number of
#' cores, and a seed gives the same answer on any machine.
#'
#' The replicates run in separate R processes (a socket cluster), on any
#' system, with nothing to set up: each starts with one math-library thread (so
#' the processes don't compete for the cores) and with this session's package
#' libraries. They load VFmixing from the installed package (so install the
#' version you are using; a `devtools::load_all()` copy isn't seen by them),
#' and the global variables and functions that `baseVectorFields` uses are
#' copied to them. On a cluster scheduler, request as many cores as you pass.
#'
#' Compiled code can't be copied to other processes: if `baseVectorFields`
#' uses a function compiled with Rcpp (`Rcpp::sourceCpp()` or
#' `Rcpp::cppFunction()`), rebuild it in each of them with `worker_setup`, for
#' example `worker_setup = function(){Rcpp::sourceCpp("interp.cpp")}` (it runs
#' after the global objects are copied, so what it defines replaces them).
#' Functions from installed packages are fine, compiled or not.
#'
#' On Linux and macOS, `options(VFmixing.parallel = "fork")` uses forked
#' processes ([parallel::mclapply()]) instead: copies of this session, which
#' start faster and need no `worker_setup`, but can hang silently depending on
#' how R and its math library were built (for example after threaded OpenBLAS
#' has run, or inside RStudio). A forked process can need as much memory as the
#' whole session, so request enough for about `cores` copies; if they can't be
#' started for lack of memory, fewer are used (halving down to one), with a
#' warning, and if a worker dies, its replicates are rerun in the session.
#'
#' With `quiet = FALSE`, a timestamped line is printed as each replicate
#' finishes. In a terminal or batch log these appear as they happen; in RStudio
#' and other consoles, the lines from separate R processes aren't shown, but
#' the start and end lines are.

#'
#' Because the conditional means are fit without noise, the selection part
#' settles with few replicates: 20 is usually plenty. This part still conditions
#' on the base vector fields evaluated along the estimated path, so the error
#' from the path position in those evaluations is not included.
#'
#' @param Data a numeric n x (d+1) data frame or matrix. The first column is time with the following each dimension's position.
#' @param baseVectorFields a baseVectorField function.
#' @param path_nsplits a positive integer, or NULL for a quarter of the number of data points.
#' The maximum number of splits for the path search.
#' @param traj_nsplits a positive integer. The maximum number of splits for the trajectory search.
#' @param path_method a string. The path method: "FullFours", "FullRegression" or "Projection".
#' @param V_Smooth a boolean. TRUE for C2 smoothness of the path and trajectory splines, FALSE for C1.
#' @param Model_Selection_Type a string. The criterion for both node searches, "BIC" or "AIC".
#' @param patience NULL or a positive integer. Passed to both node searches.
#' @param undersmooth a boolean. If TRUE (the default), both node searches use twice the selected number of regions. See Details.
#' @param traj_refine a positive whole number. Each selected trajectory region is split into this many equal
#' regions and refit. The default, 1, keeps the selected regions. See Details.
#' @param bootstrap a boolean. If TRUE, add the node-selection variance to the trajectory covariance. See Details. Default FALSE.
#' @param B a positive integer. The number of bootstrap replicates.
#' @param quad_nodes a positive integer. Passed to [getEstimatedTrajectories_OneShot()].
#' @param t_grid_size a positive number. The time grid for the pointwise variances in `Variance$Pointwise`.
#' @param cores a positive whole number. The number of cores to run the bootstrap replicates on. See Details.
#' @param worker_setup NULL or a function with no arguments, run once in each separate R process the replicates
#' use (on Windows and inside RStudio) before they start. Use it to rebuild what can't be copied to those
#' processes, such as functions compiled with Rcpp. See Details.
#' @param quiet a boolean. If TRUE, the progress output of the estimation functions is suppressed.
#'
#' @returns a list with:
#' * `EstPath`: the EstimatedPath object.
#' * `EstTraj`: the EstimatedTrajectory object, with OneShot's covariance as `Cov`. If `bootstrap`, it also has
#'   `SelectionMeans`: the conditional means, one per successful replicate (each a list with `Traj` and `TimeSplits`).
#' * `Variance`: if `bootstrap`, a list with `Within` (OneShot's covariance per region) and `Pointwise`, a data
#'   frame with the time, trajectory, and the within, selection and total variances of the trajectory value with
#'   the selection share of the total. NULL otherwise.
#' * `Bootstrap`: if `bootstrap`, a list with `B`, `Successful` (the number of replicates that fit), `Nodes`
#'   (each replicate's path and trajectory nodes) and `ResidualScale` (the residuals' scale factor). NULL otherwise.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' fit = getEstimatedTrajectories_Full(ExData, baseVectorFields, traj_nsplits = 50,
#'                                     bootstrap = TRUE, B = 10, quiet = TRUE)
#' head(fit$Variance$Pointwise)
#' # Coverage of the true trajectories with the total variance
#' coverageRateTraj(fit$EstTraj, ExTraj$Traj, ExTraj$TimeSplits)
#' }
#'
getEstimatedTrajectories_Full = function(Data, baseVectorFields, path_nsplits = NULL, traj_nsplits = 100,
                                         path_method = "FullFours", V_Smooth = F, Model_Selection_Type = "BIC",
                                         patience = 10, undersmooth = T, traj_refine = 1, bootstrap = F, B = 20,
                                         quad_nodes = 20, t_grid_size = 0.01, cores = 1, worker_setup = NULL, quiet = F){

  checkCores(cores)

  if(!is.null(worker_setup) && !is.function(worker_setup)){
    stop("worker_setup must be NULL or a function with no arguments.", call. = FALSE)
  }

  if(!(path_method %in% c("FullFours", "FullRegression", "Projection"))){
    stop("path_method must be \"FullFours\", \"FullRegression\" or \"Projection\": OneShot needs a path with a covariance structure.",
         call. = FALSE)
  }

  if(!(Model_Selection_Type %in% c("BIC", "AIC"))){
    stop("Model_Selection_Type must be \"BIC\" or \"AIC\".", call. = FALSE)
  }

  checkPatience(patience)
  checkRefine(traj_refine)

  if(bootstrap && (length(B) != 1 || !is.numeric(B) || B < 2 || B != round(B))){
    stop("B must be a whole number of at least 2.", call. = FALSE)
  }

  Data = data.frame(Data)
  n_dim = ncol(Data) - 1
  names(Data) = c("t", stringr::str_c("X", 1:n_dim))

  if(is.null(path_nsplits)){
    path_nsplits = floor(nrow(Data)/4)
  }

  runQuietly = function(expr, silent){
    if(silent){
      res = NULL
      utils::capture.output(res <- expr)
      res
    } else{
      expr
    }
  }

  fitPath = function(d){
    switch(path_method,
           FullFours = getEstimatedPaths_FullFours(d, nsplits = path_nsplits, V_Smooth = V_Smooth, Model_Selection_Type = Model_Selection_Type,
                                                   patience = patience, undersmooth = undersmooth),
           FullRegression = getEstimatedPaths_FullRegression(d, nsplits = path_nsplits, V_Smooth = V_Smooth, Model_Selection_Type = Model_Selection_Type,
                                                             patience = patience, undersmooth = undersmooth),
           Projection = getEstimatedPaths_Projection(d, nsplits = path_nsplits, V_Smooth = V_Smooth, Model_Selection_Type = Model_Selection_Type,
                                                     patience = patience, undersmooth = undersmooth))
  }

  fitTraj = function(EstPath, return_cov = TRUE){
    getEstimatedTrajectories_OneShot(EstPath, nsplits = traj_nsplits, baseVectorFields = baseVectorFields, V_Smooth = V_Smooth,
                                     Model_Selection_Type = Model_Selection_Type, quad_nodes = quad_nodes,
                                     patience = patience, undersmooth = undersmooth, refine = traj_refine, return_cov = return_cov)
  }

  # ---- The fit on the data ---------------------------------------------------

  EstPath = runQuietly(fitPath(Data), quiet)
  EstTraj = runQuietly(fitTraj(EstPath), quiet)

  if(!bootstrap){
    return(list(EstPath = EstPath, EstTraj = EstTraj, Variance = NULL, Bootstrap = NULL))
  }

  # ---- The selection part of the variance ------------------------------------

  fitted = getEstPathPosition(Data$t, EstPath)
  Data_fit = Data
  Data_fit[, -1] = fitted
  resid = as.matrix(Data[, -1]) - fitted
  resid = sweep(resid, 2, colMeans(resid))

  # Least-squares residuals are smaller than the errors: rescale them so the
  # resampled errors have variance RSS/(n - p), with p the path's free
  # parameters per dimension, matching the path's own variance estimate.
  n_obs = nrow(resid)
  p_path = splineFreeParams(length(EstPath$Path), inferPathVSmooth(EstPath, default = V_Smooth))
  residual_scale = if(n_obs > p_path) sqrt(n_obs/(n_obs - p_path)) else 1
  resid = resid*residual_scale

  # Every replicate's resampling is drawn here, before any replicate runs, so
  # the results are the same however many cores run them.
  draws = lapply(seq_len(B), function(b){sample(nrow(resid), replace = TRUE)})

  runReplicate = function(b){

    Data_b = Data_fit
    Data_b[, -1] = fitted + resid[draws[[b]], , drop = FALSE]

    tryCatch(runQuietly({
      EstPath_b = fitPath(Data_b)
      EstTraj_b = fitTraj(EstPath_b, return_cov = FALSE)   # only its nodes are used

      # The conditional mean for this replicate's nodes: the same nodes fit to the noise-free fitted positions
      EstPath_m = runEstPathStep_FullRegression(Data_fit, EstPath_b$TimeSplits, V_Smooth = V_Smooth)
      cond_traj = oneShotFixedNodes(EstPath_m, EstTraj_b$TimeSplits, baseVectorFields, V_Smooth, quad_nodes, with_cov = FALSE)

      # Kept as it is, on the replicate's own nodes: no projection
      list(Mean = list(Traj = cond_traj$Traj, TimeSplits = cond_traj$TimeSplits),
           Nodes = list(Path = EstPath_b$TimeSplits, Traj = EstTraj_b$TimeSplits))
    }, TRUE), error = function(e){list(Mean = NULL, Nodes = NULL, Error = conditionMessage(e))})   # this replicate failed to fit

  }

  if(!quiet){
    message(sprintf("[%s] Running %d bootstrap replicates%s...", format(Sys.time(), "%H:%M:%S"), B,
                    if(cores > 1) stringr::str_c(" on ", min(cores, B), " cores (", parallelType(), " workers)") else ""))
  }
  boot_start = Sys.time()

  replicates = runInParallel(seq_len(B), runReplicate, cores, baseVectorFields, progress = !quiet, label = "Bootstrap replicate",
                             worker_setup = worker_setup)

  if(!quiet){
    message(sprintf("[%s] Bootstrap done in %.1f min", format(Sys.time(), "%H:%M:%S"),
                    as.numeric(difftime(Sys.time(), boot_start, units = "mins"))))
  }

  # A replicate that failed to fit has no conditional mean (and its error message)
  errors = unlist(lapply(replicates, function(r){
    if(is.list(r)) r$Error else if(inherits(r, "try-error")) as.character(r) else NULL
  }))
  replicates = lapply(replicates, function(r){if(is.list(r) && !is.null(r$Mean)) r else NULL})
  means = lapply(replicates, function(r){r$Mean})
  nodes = lapply(replicates, function(r){r$Nodes})

  ok = !vapply(means, is.null, logical(1))
  n_ok = sum(ok)

  if(n_ok < 2){
    socket = cores > 1 && parallelType() == "socket"
    compiled = length(errors) > 0 && any(grepl("symbol address|external pointer|NativeSymbol", errors))
    stop("Fewer than two bootstrap replicates could be fit.",
         if(length(errors) > 0) stringr::str_c(" The first error was: ", errors[1]),
         if(socket && compiled)
           stringr::str_c(" The parallel workers are separate R processes, and compiled code that baseVectorFields uses ",
                          "(such as an Rcpp function) can't be copied to them. Rebuild it in each worker with worker_setup, ",
                          "for example worker_setup = function(){Rcpp::sourceCpp(\"your_file.cpp\")}, or, on Linux or macOS, use forked ",
                          "workers (options(VFmixing.parallel = \"fork\")).")
         else if(socket)
           " The parallel workers use the installed VFmixing, so check it is the version you are using.",
         call. = FALSE)
  }

  # OneShot's covariance stays the within part; the conditional means give the
  # selection part pointwise, and the interval functions add the two.
  EstTraj$SelectionMeans = means[ok]

  list(EstPath = EstPath,
       EstTraj = EstTraj,
       Variance = list(Within = EstTraj$Cov, Pointwise = pointwiseTrajVariance(EstTraj, t_grid_size)),
       Bootstrap = list(B = B, Successful = n_ok, Nodes = nodes[ok], ResidualScale = residual_scale))

}

#' Refine the Trajectories of a Full-Method Fit
#'
#' @description
#' Splits every trajectory region of a [getEstimatedTrajectories_Full()] fit
#' into `refine` equal regions and refits, reusing the fit's node selections
#' instead of searching again. With a bootstrap, every replicate's conditional
#' mean is refit on its own refined trajectory nodes too, so the result keeps
#' the node-selection variance.
#'
#' @details
#' Refinement happens after node selection, so every node search (on the data
#' and in each bootstrap replicate) would select the same nodes with or without
#' it. This gives the same result as rerunning [getEstimatedTrajectories_Full()]
#' with `traj_refine` multiplied by `refine` on the same bootstrap replicates,
#' at the cost of fixed-node refits only. Comparing a fit with its refined
#' version is also a paired comparison: both use the same replicates.
#'
#' @param fit a list returned by [getEstimatedTrajectories_Full()].
#' @param Data the data the fit was made from.
#' @param baseVectorFields the baseVectorField function the fit was made with.
#' @param refine a positive whole number. Each trajectory region is split into this many equal regions.
#' @param V_Smooth,quad_nodes,t_grid_size the values the fit was made with (the defaults match those of
#' [getEstimatedTrajectories_Full()]).
#'
#' @returns a list like that of [getEstimatedTrajectories_Full()], with the refined trajectories.
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' fit = getEstimatedTrajectories_Full(ExData, baseVectorFields, traj_nsplits = 50,
#'                                     bootstrap = TRUE, B = 10, quiet = TRUE)
#' refined = refineEstimatedTrajectories_Full(fit, ExData, baseVectorFields, refine = 2)
#' length(fit$EstTraj$Traj)
#' length(refined$EstTraj$Traj)
#' }
#'
refineEstimatedTrajectories_Full = function(fit, Data, baseVectorFields, refine = 2, V_Smooth = F, quad_nodes = 20,
                                            t_grid_size = 0.01){

  checkRefine(refine)

  if(refine == 1){
    return(fit)
  }

  Data = data.frame(Data)
  names(Data) = c("t", stringr::str_c("X", seq_len(ncol(Data) - 1)))

  EstTraj = oneShotFixedNodes(fit$EstPath, refineTimeSplits(fit$EstTraj$TimeSplits, refine), baseVectorFields, V_Smooth, quad_nodes)

  if(is.null(fit$Bootstrap)){
    return(list(EstPath = fit$EstPath, EstTraj = EstTraj, Variance = NULL, Bootstrap = NULL))
  }

  # Each replicate's conditional mean on its own refined trajectory nodes, with
  # its path nodes fit to the noise-free fitted positions as before
  Data_fit = Data
  Data_fit[, -1] = getEstPathPosition(Data$t, fit$EstPath)

  nodes = lapply(fit$Bootstrap$Nodes, function(n){list(Path = n$Path, Traj = refineTimeSplits(n$Traj, refine))})

  EstTraj$SelectionMeans = lapply(nodes, function(n){
    EstPath_m = runEstPathStep_FullRegression(Data_fit, n$Path, V_Smooth = V_Smooth)
    cond_traj = oneShotFixedNodes(EstPath_m, n$Traj, baseVectorFields, V_Smooth, quad_nodes, with_cov = FALSE)
    list(Traj = cond_traj$Traj, TimeSplits = cond_traj$TimeSplits)
  })

  Bootstrap = fit$Bootstrap
  Bootstrap$Nodes = nodes

  list(EstPath = fit$EstPath,
       EstTraj = EstTraj,
       Variance = list(Within = EstTraj$Cov, Pointwise = pointwiseTrajVariance(EstTraj, t_grid_size)),
       Bootstrap = Bootstrap)

}

# Validate a cores argument: a positive whole number.
checkCores = function(cores){

  if(length(cores) != 1 || !is.numeric(cores) || cores < 1 || cores != round(cores)){
    stop("cores must be a positive whole number.", call. = FALSE)
  }

  invisible(cores)
}

# lapply(X, FUN), on up to `cores` cores. Forks on Linux and macOS
# (parallel::mclapply). On Windows, which can't fork, it starts a socket
# cluster: the workers load VFmixing from the same libraries, and any global
# variables and functions that baseVectorFields uses are copied to them.
runInParallel = function(X, FUN, cores, baseVectorFields = NULL, progress = FALSE, label = "Job", worker_setup = NULL){

  cores = min(cores, length(X))
  job = parallelJob(X, FUN, progress, label, limit_threads = cores > 1)

  if(cores <= 1){
    return(lapply(seq_along(X), job))
  }

  if(parallelType() == "fork"){
    # Keep this session's math library to one thread while forking: forked
    # children of a session that has used threaded OpenBLAS can hang
    restore = limitSessionThreads()
    on.exit(restore(), add = TRUE)
    return(forkWithFallback(seq_along(X), job, cores))
  }

  # Separate R processes, started with what they need from the first
  # instruction: one math-library thread each (a new process's OpenBLAS would
  # otherwise start a thread per core on the whole machine, and the workers
  # would fight over the cores), and this session's package libraries (so they
  # load the same VFmixing). New processes inherit these environment
  # variables; this session's own values are restored once they are running.
  # Their output (the progress lines) goes to this session's standard output,
  # which a terminal or batch log shows.
  restore_env = setEnvTemporarily(c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
                                    VECLIB_MAXIMUM_THREADS = "1",
                                    R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep)))
  cl = tryCatch(parallel::makeCluster(cores, outfile = ""), finally = restore_env())
  on.exit(parallel::stopCluster(cl), add = TRUE)

  # The same library paths, set in each worker too (with a function that
  # doesn't need VFmixing loaded to run)
  set_paths = function(paths){.libPaths(paths); invisible(NULL)}
  environment(set_paths) = baseenv()
  parallel::clusterCall(cl, set_paths, .libPaths())

  # The global objects the field function and the setup function use
  globals = unique(c(if(!is.null(baseVectorFields)) fieldGlobals(baseVectorFields),
                     if(!is.null(worker_setup)) fieldGlobals(worker_setup)))
  if(length(globals) > 0){
    parallel::clusterExport(cl, globals, envir = globalenv())
  }

  # Anything the workers must build themselves, such as compiled (Rcpp)
  # functions, whose pointers into compiled code don't survive being copied
  if(!is.null(worker_setup)){
    parallel::clusterCall(cl, worker_setup)
  }

  parallel::parLapplyLB(cl, seq_along(X), job)

}

# How to run in parallel: "socket" (separate R processes, the default) or
# "fork". Forked processes are copies of this session, which is faster to start
# but can hang silently depending on how R and its math library were built (for
# example after threaded OpenBLAS has run, or inside RStudio), so they are only
# used on request: options(VFmixing.parallel = "fork"). Never on Windows, which
# can't fork.
parallelType = function(){

  if(.Platform$OS.type != "windows" && identical(getOption("VFmixing.parallel"), "fork")){
    return("fork")
  }
  "socket"

}

# Set environment variables, returning a function that restores the old values.
setEnvTemporarily = function(values){

  old = Sys.getenv(names(values), unset = NA, names = TRUE)
  do.call(Sys.setenv, as.list(values))

  function(){
    for(name in names(old)){
      if(is.na(old[[name]])) Sys.unsetenv(name) else do.call(Sys.setenv, stats::setNames(list(old[[name]]), name))
    }
    invisible(NULL)
  }

}

# The function each parallel job runs, built here so it carries only what it
# needs: it limits the math library to one thread, runs FUN on its element and
# reports when it is done.
parallelJob = function(X, FUN, progress, label, limit_threads){

  n = length(X)
  force(FUN); force(progress); force(label); force(limit_threads)

  function(i){
    if(limit_threads){
      limitWorkerThreads()
    }
    t0 = Sys.time()
    out = FUN(X[[i]])
    if(progress){
      message(sprintf("[%s] %s %d/%d done (%.1f s)", format(Sys.time(), "%H:%M:%S"), label, i, n,
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    out
  }

}

# One math-library thread in a parallel worker, so the workers don't compete
# for the cores (needs the RhpcBLASctl package for BLAS libraries that are
# already running; the environment variables cover anything started later).
limitWorkerThreads = function(){

  Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
  if(requireNamespace("RhpcBLASctl", quietly = TRUE)){
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
  }
  invisible(NULL)

}

# One math-library thread in this session, returning a function that restores
# the previous number.
limitSessionThreads = function(){

  if(!requireNamespace("RhpcBLASctl", quietly = TRUE)){
    return(function(){invisible(NULL)})
  }
  old = RhpcBLASctl::blas_get_num_procs()
  RhpcBLASctl::blas_set_num_threads(1)
  function(){RhpcBLASctl::blas_set_num_threads(old); invisible(NULL)}

}

# parallel::mclapply(X, FUN) on `cores` forked processes, made robust to memory
# limits. If the processes can't be started (a fork can need as much memory as
# the whole session, so many cores can exceed a job's memory), it retries on
# half as many cores, down to running in this session. If a worker dies while
# running (for example, killed for using too much memory), its elements are
# rerun here. FUN must not return NULL itself: NULL marks a dead worker.
forkWithFallback = function(X, FUN, cores){

  while(cores > 1){

    res = tryCatch(parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE), error = function(e){e})

    if(!inherits(res, "error")){
      lost = which(vapply(res, function(r){is.null(r) || inherits(r, "try-error")}, logical(1)))
      if(length(lost) > 0){
        warning(length(lost), " of ", length(X), " parallel jobs didn't finish (a worker may have run out of memory); ",
                "running them here instead.", call. = FALSE)
        res[lost] = lapply(X[lost], FUN)
      }
      return(res)
    }

    if(!grepl("fork|memory", conditionMessage(res), ignore.case = TRUE)){
      stop(res)
    }

    new_cores = max(1, floor(cores/2))
    warning("Couldn't start ", cores, " processes (", conditionMessage(res), "); trying ", new_cores, ".", call. = FALSE)
    cores = new_cores

  }

  lapply(X, FUN)

}

# The global variables and functions a vector field function defined at the
# top level uses, followed through the global functions it calls.
fieldGlobals = function(f){

  if(!is.function(f) || !identical(environment(f), globalenv()) || !requireNamespace("codetools", quietly = TRUE)){
    return(character(0))
  }

  found = character(0)
  queue = list(f)

  while(length(queue) > 0){
    g = queue[[1]]
    queue = queue[-1]
    names_g = codetools::findGlobals(g)
    names_g = names_g[vapply(names_g, exists, logical(1), envir = globalenv(), inherits = FALSE)]
    new = setdiff(names_g, found)
    found = c(found, new)
    for(nm in new){
      obj = get(nm, envir = globalenv())
      if(is.function(obj) && identical(environment(obj), globalenv())){
        queue[[length(queue) + 1]] = obj
      }
    }
  }

  found

}

# OneShot trajectory with fixed nodes (no node search) on a given path. With
# with_cov = FALSE, only the coefficients and errors (no covariance).
oneShotFixedNodes = function(EstPath, TimeSplits, baseVectorFields, V_Smooth = F, quad_nodes = 20, with_cov = TRUE){

  n_dim = nrow(EstPath$Path[[1]])
  n_models = ncol(baseVectorFields(min(EstPath$TimeSplits), rep(0, n_dim)))
  gl = gaussLegendre(quad_nodes)

  Omega_p = constructUnboundedCubicSplineTransitionMatrix_Traj(EstPath$TimeSplits, n_dim, inferPathVSmooth(EstPath, default = V_Smooth))

  caches = lapply(seq_len(length(TimeSplits) - 1), function(r){
    oneShotRegionCache(TimeSplits[r], TimeSplits[r + 1], EstPath, baseVectorFields, n_dim, n_models, gl)
  })

  runEstTrajStep_OneShot(EstPath, TimeSplits, caches, Omega_p, n_dim, n_models, V_Smooth, with_cov = with_cov)

}

# The trajectory values (one column per model) at times t.
evaluateTrajAt = function(t, EstTraj){

  ts = EstTraj$TimeSplits
  idx = pmin(pmax(findInterval(t, ts, left.open = TRUE), 1), length(ts) - 1)
  u = t - ts[idx]
  n_models = nrow(EstTraj$Traj[[1]])

  out = matrix(NA_real_, length(t), n_models)
  for(i in seq_along(t)){
    out[i, ] = EstTraj$Traj[[idx[i]]] %*% c(u[i]^3, u[i]^2, u[i], 1)
  }
  out

}

# The within-node variance of each trajectory's value at times t, from
# OneShot's covariance (one row per time, one column per trajectory).
trajWithinVariance = function(EstTraj, t){

  TimeSplits = EstTraj$TimeSplits
  n_models = nrow(EstTraj$Traj[[1]])
  idx = pmin(pmax(findInterval(t, TimeSplits, left.open = TRUE), 1), length(TimeSplits) - 1)
  u = t - TimeSplits[idx]
  B = cbind(u^3, u^2, u, 1)

  variance = matrix(0, length(t), n_models)
  for(r in unique(idx)){
    sel = idx == r
    Bs = B[sel, , drop = FALSE]
    for(m in seq_len(n_models)){
      block = EstTraj$Cov[[r]][1:4 + (m - 1)*4, 1:4 + (m - 1)*4]
      variance[sel, m] = rowSums((Bs %*% block) * Bs)
    }
  }
  pmax(variance, 0)

}

# The node-selection variance of each trajectory's value at times t: the sample
# variance (divisor B - 1) of the bootstrap conditional means at each time.
# Zero if the trajectory has no bootstrap (one row per time, one column per trajectory).
trajSelectionVariance = function(EstTraj, t){

  n_models = nrow(EstTraj$Traj[[1]])

  if(is.null(EstTraj$SelectionMeans)){
    return(matrix(0, length(t), n_models))
  }

  values = vapply(EstTraj$SelectionMeans, function(m){evaluateTrajAt(t, m)}, matrix(0, length(t), n_models))
  values = array(values, c(length(t), n_models, length(EstTraj$SelectionMeans)))
  apply(values, c(1, 2), stats::var)

}

# Within, selection and total variances of each trajectory's value on a time grid.
pointwiseTrajVariance = function(EstTraj, t_grid_size){

  TimeSplits = EstTraj$TimeSplits
  t_grid = seq(min(TimeSplits), max(TimeSplits), by = t_grid_size)
  within = trajWithinVariance(EstTraj, t_grid)
  selection = trajSelectionVariance(EstTraj, t_grid)

  rows = lapply(seq_len(ncol(within)), function(m){
    total = within[, m] + selection[, m]
    data.frame(t = t_grid, trajectory = m, var_within = within[, m], var_selection = selection[, m],
               var_total = total, selection_share = ifelse(total > 0, selection[, m]/total, 0))
  })

  do.call(rbind, rows)

}

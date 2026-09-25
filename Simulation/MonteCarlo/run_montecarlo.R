# Monte Carlo comparison of Diff-in-Disc bandwidth selectors.
#
# Usage (from the project root):  Rscript Simulation/MonteCarlo/run_montecarlo.R [R] [cores]
# Writes Simulation/MonteCarlo/results/mc_results.rds, read by MonteCarlo.Rmd.
suppressMessages(library(rdrobust))
source(here::here("Bandwidth/SGD/BW_Functions.R"))
source(here::here("Simulation/MonteCarlo/dgp.R"))
source(here::here("Simulation/MonteCarlo/estimators.R"))

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1) as.integer(args[1]) else 1000
cores <- if (length(args) >= 2) as.integer(args[2]) else max(1, parallel::detectCores() - 2)
designs <- c("sharp", "smooth")
fixed_h <- c(50, 100, 150, 200, 250, 300, 400, 500, 600, 800, 1000)
R_check <- min(10, R)  # replications where DiRD_BW itself is run, to check it finds the same minimum

# Replication r uses seed 10000 + r in both designs (common random numbers)
one_rep <- function(r){
  out <- lapply(designs, function(dsg){
    res <- estimate_all(simulate_sim_data(seed = 10000 + r, design = dsg), fixed_h)
    cbind(rep = r, design = dsg, res)
  })
  do.call(rbind, out)
}

cat(sprintf("Running %d replications x %d designs on %d cores...\n", R, length(designs), cores))
t_main <- system.time(reps <- parallel::mclapply(seq_len(R), one_rep, mc.cores = cores))
failed <- vapply(reps, inherits, logical(1), what = "try-error")
if (any(failed)) {
  cat(sprintf("%d replications failed, first error:\n", sum(failed)))
  print(reps[[which(failed)[1]]])
}
results <- do.call(rbind, reps[!failed])
cat(sprintf("Done in %.1f min\n", t_main[["elapsed"]] / 60))

# DiRD_BW (the SGD itself) vs the grid + golden-section minimizer of the same objective
check_rep <- function(job){
  data <- simulate_sim_data(seed = 10000 + job$rep, design = job$design)
  d <- data.frame(ID = data$PropertyID, y = data$Price_m2, x = data$Distance_to_Border,
                  time_var = data$Year)
  J <- function(h) penalized_MSE(d, h, c = 0, t0 = 3, N = nrow(d), kernel = job$kernel)
  set.seed(job$rep)
  sgd <- DiRD_BW(y = d$y, x = d$x, c = 0, t0 = 3, time_var = d$time_var, ID = d$ID,
                 kernel = job$kernel)
  h_grid <- bw_sgd_criterion(data, job$kernel)
  data.frame(job, h_sgd = sgd$optimal_h, J_sgd = sgd$optimal_mse,
             h_grid = h_grid, J_grid = J(h_grid))
}
jobs <- expand.grid(rep = seq_len(R_check), design = designs,
                    kernel = c("triangular", "uniform"), stringsAsFactors = FALSE)
t_check <- system.time(
  sgd_check <- do.call(rbind, parallel::mclapply(split(jobs, seq_len(nrow(jobs))), check_rep,
                                                 mc.cores = cores)))
cat(sprintf("SGD check done in %.1f min\n", t_check[["elapsed"]] / 60))

saveRDS(list(results = results, sgd_check = sgd_check, R = R, fixed_h = fixed_h, tau = 350,
             created = Sys.time(), session = sessionInfo()),
        here::here("Simulation/MonteCarlo/results/mc_results.rds"))

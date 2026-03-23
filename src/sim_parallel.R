source("src/functions.R")
library(parallel)

start <- Sys.time() # track simulation time taken

# simulation parameters
n_values <- c(100, 500, 1000)

rho_values  <- c(0, 0.3, 0.6, 0.9)
sigma2 <- 0.5

B  <- 10 # has to be 500
l_values <- c(1, 2, 5)
g_values <- c(3, 5, 10)
h_grid <- seq(0.01, 0.2, by = 0.001)
p <- 1

# grid of sim combinations instead of two for loops
scenarios <- expand.grid(n = n_values, rho = rho_values)

# cluster
n_cores <- detectCores() - 1 
cl      <- makeCluster(n_cores)

clusterEvalQ(cl, library(MASS))

# export everything we need to each core
clusterExport(cl, varlist = c("simulate_bandwidth", "bandwidth_cv","bandwidth_mcv",
                              "bandwidth_pcv","bandwidth_cdpi", "local_poly_est",
                              "generate_ar1_errors", "r_true", "l_values", 
                              "g_values", "h_grid", "p", "sigma2", "B", "scenarios"))

clusterSetRNGStream(cl, 123) 

# run in parallel over scenarios
results_list <- parLapply(cl, seq_len(nrow(scenarios)), function(i) {
  n   <- scenarios$n[i]
  rho <- scenarios$rho[i]
  cat("Running: n =", n, "| rho =", rho, "\n")
  simulate_bandwidth(n, rho, sigma2, B, l_values, g_values, h_grid, p)
})

stopCluster(cl)

results <- do.call(rbind, lapply(results_list, function(x) x$summary))

end <- Sys.time()
print(end - start)

saveRDS(results_list, file = "data/sim_results_list_20260323.rds")
saveRDS(results,      file = "data/sim_results_20260323.rds")

h_mats <- lapply(results_list, function(x) x$h_mat)
saveRDS(h_mats, file = "data/sim_h_mats_20260323.rds")
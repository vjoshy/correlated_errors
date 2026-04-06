source("src/functions.R")
library(parallel)
library(pbapply)

start <- Sys.time()

# Simulation Parameters
n_values <- c(1000) 
rho_values  <- c(0, 0.3, 0.6, 0.9)
sigma2 <- 0.5

B  <- 500
l_values  <- c(1, 2, 5)
g_values  <- c(3, 5, 10)
h_grid <- seq(0.01, 0.2, by = 0.001)
p <- 1


n_cores <- detectCores() - 1
cl      <- makeCluster(n_cores)


clusterEvalQ(cl, library(MASS))

clusterExport(cl, varlist = c("bandwidth_cv", "bandwidth_mcv", "bandwidth_pcv", 
                              "bandwidth_cdpi", "local_poly_est", 
                              "generate_ar1_errors", "r_true"))

clusterSetRNGStream(cl, 123) 


results_list <- vector("list", length(n_values) * length(rho_values))
idx  <- 1


for (n in n_values) {
  for (rho in rho_values) {
    
    check_point <- Sys.time()
    
    cat("\n========================================\n")
    cat("Starting Scenario: n =", n, "| rho =", rho, "\n")
    cat("========================================\n")
    
    results_list[[idx]] <- simulate_bandwidth(n, rho, sigma2, B,
                                              l_values, g_values,
                                              h_grid, p, cl = cl)
    idx <- idx + 1
    
    saveRDS(results_list, file = "data/sim_results_list_partial.rds")
    
    cat("Time taken: ", Sys.time() - check_point)
    
  }
}


stopCluster(cl)


results <- do.call(rbind, lapply(results_list, function(x) x$summary))

end <- Sys.time()
print(end - start)


saveRDS(results, file = "data/sim_results_20260323.rds")

h_mats <- lapply(results_list, function(x) x$h_mat)
saveRDS(h_mats, file = "data/sim_h_mats_20260323.rds")
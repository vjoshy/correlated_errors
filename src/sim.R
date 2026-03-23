source("src/functions.R")

start <- Sys.time()

# simulation parameters
#n_values  <- c(10)
n_values <- c(10, 10, 10)

rho_values  <- c(0, 0.3, 0.6, 0.9)
sigma2 <- 0.5

B  <- 10 # has to be 500
l_values  <- c(1, 2, 5)
g_values  <- c(3, 5, 10)
h_grid <- seq(0.01, 0.2, by = 0.001)
p <- 1

set.seed(123)
results_list <- vector("list", length(n_values) * length(rho_values))
idx          <- 1

for (n in n_values) {
  for (rho in rho_values) {
    cat("Running: n =", n, "| rho =", rho, "\n")
    results_list[[idx]] <- simulate_bandwidth(n, rho, sigma2, B,
                                              l_values, g_values,
                                              h_grid, p)
    idx <- idx + 1
    
    saveRDS(results_list, file = "data/sim_results_list_partial.rds")
  }
}

results <- do.call(rbind, results_list)

end <- Sys.time()

print(end - start)

saveRDS(results, file = "data/sim_results.rds")

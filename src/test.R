source("src/functions.R")

set.seed(123)

n <- 1000
sigma2 <- 0.5
x <- ((1:n)/n) 

errors_iid <- rnorm(n, 0, sigma2) 
y_true <- r_true(x)

p <- 1

errors_ar1 <- generate_ar1_errors(n, 0.9, sigma2)
h_grid <- seq(0.0001, 0.1, by = 0.001)

# CV bandwidth
#y_baseline  <- y_true + errors_iid
y_baseline  <- y_true + errors_ar1

h_baseline  <- bandwidth_cv(x, y_baseline, p = p,h_grid = h_grid)$minimum

h_mcv  <- bandwidth_mcv(x, y_baseline, p = p, l = 1,h_grid = h_grid)$minimum
h_pcv  <- bandwidth_pcv(x, y_baseline, p = p, g = 5, h_grid = h_grid)$minimum

h_cdpi <- bandwidth_cdpi(x,y_baseline)$minimum

y_baseline_hat <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                        p = p, h = h_baseline)$beta[1])
y_mcv_hat <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                        p = p, h = h_mcv)$beta[1])
y_pcv_hat <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                        p = p, h = h_pcv)$beta[1])

y_cdpi_hat <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                p = p, h = h_cdpi)$beta[1])

plot(x, y_baseline, col = "grey", pch = 19, cex = 0.3)
lines(x, y_true, col = "red", lwd = 2)
lines(x, y_baseline_hat, col = "blue", lwd = 2)
lines(x, y_mcv_hat, col = "black", lwd = 2)
lines(x, y_pcv_hat, col = "green", lwd = 2)
lines(x, y_cdpi_hat, col = "pink", lwd = 2)
legend("topright", legend = c("True", "NW fit"), col = c("red", "blue"), lty = 1)


# plot ----------------------------------------------------------------------------------------


source("src/functions.R")

set.seed(123)

n <- 1000
sigma2 <- 0.5
x <- ((1:n)/n) 

errors_iid <- rnorm(n, 0, sigma2) 
y_true <- r_true(x)

p <- 1

cor <- c(0, 0.3, 0.6, 0.9)
h_grid <- seq(0.00001, 0.1, by = 0.00001)

y_list <- list()
h_list <- list()

for(e in 1:length(cor)){
  
  errors_ar1 <- generate_ar1_errors(n,cor[e], sigma2)
  
  if (e == 1){
    h_grid_temp <- h_grid[h_grid < 0.07 & h_grid > 0.05]
  } else {
    h_grid_temp <- h_grid[h_grid < 0.01]
  }
  
  y_list[[e]]  <- y_true + errors_ar1
  h_list[[e]]  <- bandwidth_cv(x, y_list[[e]], p = p, h_grid = h_grid_temp)$minimum
  
}





library(ggplot2)
library(gridExtra)

# Generate fits for each rho
plot_list <- list()
rho_labels <- c("ρ = 0", "ρ = 0.3", "ρ = 0.6", "ρ = 0.9")

for(e in 1:length(cor)){
  
  h_selected <- h_list[[e]]
  y_obs <- y_list[[e]]
  
  # Get fitted values at each x using selected bandwidth
  y_hat <- sapply(x, function(x0){
    fit <- local_poly_est(x0, x, y_obs, p, h_selected, dnorm)
    fit$beta[1]
  })
  
  df <- data.frame(
    x      = x,
    y_obs  = y_obs,
    y_true = y_true,
    y_hat  = y_hat
  )
  
  plot_list[[e]] <- ggplot(df, aes(x = x)) +
    geom_point(aes(y = y_obs), color = "grey70", size = 0.3, alpha = 0.5) +
    geom_line(aes(y = y_true, color = "True function"), linewidth = 0.8) +
    geom_line(aes(y = y_hat, color = "LOOCV fit"), linewidth = 0.8) +
    scale_color_manual(
      values = c("True function" = "black", "LOOCV fit" = "red")
    ) +
    labs(
      title = rho_labels[e],
      subtitle = paste0("h = ", round(h_selected, 5)),
      x = "x", y = "y",
      color = NULL
    ) +
    theme_bw() +
    theme(
      legend.position  = "bottom",
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 9)
    )
}

# Combine into 2x2 grid
grid.arrange(
  plot_list[[1]], plot_list[[2]],
  plot_list[[3]], plot_list[[4]],
  ncol = 2)

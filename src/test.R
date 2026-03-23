source("src/functions.R")

set.seed(123)

n <- 1000
sigma2 <- 0.5
x <- ((1:n)/n) 

errors_iid <- rnorm(n, 0, sigma2) 
y_true <- r_true(x)

p <- 1

errors_ar1 <- generate_ar1_errors(n, 0.3, sigma2)

# CV bandwidth
#y_baseline  <- y_true + errors_iid
y_baseline  <- y_true + errors_ar1

h_baseline  <- bandwidth_cv_fast(x, y_baseline, p = p,
                         h_grid = seq(0.01, 0.05, by = 0.001))$minimum

h_mcv  <- bandwidth_mcv(x, y_baseline, p = p, l = 2,
                                 h_grid = seq(0.01, 0.05, by = 0.001))$minimum

h_pcv  <- bandwidth_pcv(x, y_baseline, p = p, g = 5,
                        h_grid = seq(0.01, 0.05, by = 0.001))$minimum

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
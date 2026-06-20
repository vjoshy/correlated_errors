source("src/functions.R")
source("src/mixed_functions.R")

#' Bandwidth Leave one out - Cross-validation (From A2)
#'
#' @param x vector of inputs
#' @param y vector of response
#' @param p degree of regression
#' @param h_grid vector of bandwidth sizes
#' @param kernel kernel function, default = standard Gaussian
#'
#' @returns list of optimal bandwidth size and a vector of MSE corresponding to each bandwidth in h_grid
#' @export
#'
bandwidth_cv <- function(x, y, p, h_grid, kernel = dnorm) {
  n        <- length(x)
  mse_grid <- numeric(length(h_grid))
  for (i in seq_along(h_grid)) {
    avg_mse <- numeric(n)
    for (j in 1:n) {
      y_hat      <- local_poly_est(x[j], x[-j], y[-j], p, h_grid[i], kernel)
      avg_mse[j] <- (y_hat$beta[1] - y[j])^2
    }
    mse_grid[i] <- mean(avg_mse)
  }
  list(minimum = h_grid[which.min(mse_grid)], mse = mse_grid)
}


bandwidth_mcv <- function(x, y, p, l, h_grid, kernel = dnorm) {
  n        <- length(x)
  mse_grid <- numeric(length(h_grid))
  for (i in seq_along(h_grid)) {
    cv_vals <- numeric(n)
    for (j in 1:n) {
      keep   <- which(abs((1:n) - j) > l)  # indices where |i - j| > l
      fit    <- local_poly_est(x[j], x[keep], y[keep], p, h_grid[i], kernel)
      cv_vals[j] <- (fit$beta[1] - y[j])^2
    }
    mse_grid[i] <- mean(cv_vals)
  }
  list(minimum = h_grid[which.min(mse_grid)], mse = mse_grid)
}


bandwidth_pcv <- function(x, y, p, g, h_grid, kernel = dnorm) {
  n        <- length(x)
  mse_mat  <- matrix(0, nrow = g, ncol = length(h_grid))
  
  for (k in 1:g) {
    idx          <- seq(k, n, by = g)
    mse_mat[k, ] <- bandwidth_cv(x[idx], y[idx], p, h_grid, kernel)$mse
  }
  
  avg_mse  <- colMeans(mse_mat)
  h_star   <- h_grid[which.min(avg_mse)]
  h_pcv    <- g^(-1/5) * h_star
  
  list(minimum = h_pcv, h_star = h_star, mse = avg_mse)
}


bandwidth_cdpi <- function(x, y, kernel = dnorm) {
  n <- length(x)
  
  # local linear pilot 
  lambda <- n^(-1/5)
  y_pilot <- sapply(x, function(x0)
    local_poly_est(x0, x, y, p = 1, lambda, kernel)$beta[1])
  
  # c(eps) estmimate
  resid <- y - y_pilot
  sigma2 <- mean(resid^2)
  rho <- sum(resid[-n] * resid[-1]) / sum(resid^2)
  c_eps <- sigma2 * (1 + rho) / (1 - rho)
  
  # theta_2,2 (cubic fit )
  g <- n^(-1/7)
  f2_hat  <- sapply(x, function(x0)
    2 * local_poly_est(x0, x, y, p = 3, g, kernel)$beta[3])
  theta22 <- mean(f2_hat^2)
  
  # plug in 
  C01K  <- (1 / (2 * sqrt(pi)))^(1/5) # constant from Ruppert (1995)
  h_dep  <- C01K * (c_eps / (n * theta22))^(1/5) 
  
  list(minimum = h_dep, c_eps = c_eps, theta22 = theta22, rho = rho)
}


bandwidth_batch_cv <- function(x, y, batch, p, h_grid, kernel = dnorm) {
  n <- length(x)
  mse_grid <- numeric(length(h_grid))
  
  for (i in seq_along(h_grid)) {
    avg_mse <- numeric(n)
    
    for (j in 1:n) {
      keep <- which(batch != batch[j])
      fit <- local_poly_est(x[j], x[keep], y[keep], p, h_grid[i], kernel)
      avg_mse[j] <- (fit$beta[1] - y[j])^2
    }
    mse_grid[i] <- mean(avg_mse)
  }
  
  list(minimum = h_grid[which.min(mse_grid)], mse = mse_grid)
}

simulate_bandwidth <- function(n, rho, sigma2, B, 
                               l_values, g_values,
                               h_grid, p = 1, cl = NULL) {
  
  x      <- (1:n) / n
  y_true <- r_true(x)
  
  method_names <- c("LOOCV",
                    paste0("MCV_l", l_values),
                    paste0("PCV_g", g_values),
                    "CDPI")
  
  clusterExport(cl, varlist = c("n", "rho", "sigma2", "x", "y_true", 
                                "l_values", "g_values", "h_grid", "p", 
                                "method_names"), 
                envir = environment())
  
  sim_results <- pblapply(1:B, function(b) {
    
    errors <- if (rho == 0) rnorm(n, 0, sqrt(sigma2)) else generate_ar1_errors(n, rho, sigma2)
    y <- y_true + errors
    
    h_loocv <- bandwidth_cv(x, y, p, h_grid)$minimum
    h_mcv   <- sapply(l_values, function(l) bandwidth_mcv(x, y, p, l, h_grid)$minimum)
    h_pcv   <- sapply(g_values, function(g) bandwidth_pcv(x, y, p, g, h_grid)$minimum)
    h_cdpi  <- bandwidth_cdpi(x, y)$minimum
    
    all_h <- c(h_loocv, h_mcv, h_pcv, h_cdpi)
    
    mse_vals <- numeric(length(all_h))
    for (m in seq_along(all_h)) {
      y_hat       <- sapply(x, function(x0) local_poly_est(x0, x, y, p, all_h[m])$beta[1])
      mse_vals[m] <- mean((y_hat - y_true)^2)
    }
    
    list(h = all_h, mse = mse_vals)
    
  }, cl = cl) 
  
  h_mat   <- do.call(rbind, lapply(sim_results, `[[`, "h"))
  mse_mat <- do.call(rbind, lapply(sim_results, `[[`, "mse"))
  
  colnames(h_mat)   <- method_names
  colnames(mse_mat) <- method_names
  
  return(
    list(
      h_mat = h_mat,
      summary = data.frame(
        n         = n,
        rho       = rho,
        method    = method_names,
        mean_h    = colMeans(h_mat),
        sd_h      = apply(h_mat, 2, sd),
        mean_mse  = colMeans(mse_mat),
        sd_mse    = apply(mse_mat, 2, sd)
      )
    )
  )
}

simulate_np_mixed_bandwidth <- function(k, nk,
                                        sigma,
                                        tau_intercept,
                                        B,
                                        l_values,
                                        g_values,
                                        h_grid,
                                        tau_slope = 0,
                                        rho_within_group = 0,
                                        effect_type = "random_intercept",
                                        design = "balanced",
                                        t_range = c(0, 1),
                                        eta_fun = r_true,
                                        p = 1,
                                        cl = NULL){
  
  method_names <- c("LOOCV",
                    "BatchCV",
                    paste0("MCV_l", l_values),
                    paste0("PCV_g", g_values),
                    "CDPI")
  if (!is.null(cl)) {
    
    clusterExport(
      cl,
      varlist = c("k", "nk", "sigma", "tau_intercept", "tau_slope",
                  "rho_within_group", "effect_type", "design", "t_range",
                  "eta_fun", "B", "l_values", "g_values", "h_grid",
                  "p", "method_names",
                  "r_true", "generate_ar1_errors", "generate_np_mixed_data",
                  "local_poly_est", "bandwidth_cv", "bandwidth_batch_cv",
                  "bandwidth_mcv", "bandwidth_pcv", "bandwidth_cdpi"),
      envir = environment()
    )
  }
  sim_results <- pblapply(1:B, function(b) {
    
    # Generate one nonparametric mixed-effects dataset
    dat <- generate_np_mixed_data(k, nk, sigma, tau_intercept, tau_slope,
                                  rho_within_group, effect_type, design, 
                                  t_range, eta_fun)
    
    x <- dat$t
    y <- dat$y
    batch <- dat$batch
    # Evaluation grid for eta(t)
    x_eval <- seq(t_range[1], t_range[2], length.out = 100)
    y_true <- eta_fun(x_eval)
    
    h_loocv <- bandwidth_cv(x, y, p, h_grid)$minimum
    h_mcv   <- sapply(l_values, function(l) bandwidth_mcv(x, y, p, l, h_grid)$minimum)
    h_pcv   <- sapply(g_values, function(g) bandwidth_pcv(x, y, p, g, h_grid)$minimum)
    h_cdpi  <- bandwidth_cdpi(x, y)$minimum
    h_batchcv <- bandwidth_batch_cv(x, y, batch, p, h_grid)$minimum
    
    all_h <- c(h_loocv, h_batchcv, h_mcv, h_pcv, h_cdpi)
    
    # Estimate eta(t) using each bandwidth
    mse_vals <- numeric(length(all_h))
    
    for (m in seq_along(all_h)) {
      y_hat <- sapply(x_eval, function(x0) local_poly_est(x0, x, y, p, all_h[m])$beta[1])
      mse_vals[m] <- mean((y_hat - y_true)^2)
    }
    list(h = all_h, mse = mse_vals)
    
  }, cl = cl)
  
  h_mat <- do.call(rbind, lapply(sim_results, `[[`, "h"))
  mse_mat <- do.call(rbind, lapply(sim_results, `[[`, "mse"))
  
  colnames(h_mat) <- method_names
  colnames(mse_mat) <- method_names
  
  return(
    list(
      h_mat = h_mat,
      mse_mat = mse_mat,
      summary = data.frame(
        k = k,
        nk = nk,
        n = k * nk,
        sigma = sigma,
        tau_intercept = tau_intercept,
        tau_slope = tau_slope,
        rho_within_group = rho_within_group,
        effect_type = effect_type,
        method = method_names,
        mean_h = colMeans(h_mat),
        sd_h = apply(h_mat, 2, sd),
        mean_mse = colMeans(mse_mat),
        sd_mse = apply(mse_mat, 2, sd)
      )
    )
  )
  
}


# # simulation
# simulate_bandwidth <- function(n, rho, sigma2, B, 
#                                l_values, g_values,
#                                h_grid, p = 1) {
#   x     <- (1:n) / n
#   y_true <- r_true(x)
#   
#   # storage: methods are LOOCV, MCV(l) for each l, PCV(g) for each g, CDPI
#   method_names <- c("LOOCV",
#                     paste0("MCV_l", l_values),
#                     paste0("PCV_g", g_values),
#                     "CDPI")
#   
#   mse_mat <- matrix(NA, nrow = B, ncol = length(method_names))
#   colnames(mse_mat) <- method_names
#   
#   h_mat <- matrix(NA, nrow = B, ncol = length(method_names))
#   colnames(h_mat) <- method_names
#   
#   for (b in 1:B) {
#     if (b %% 50 == 0) cat("n =", n, "| rho =", rho, "| rep", b, "/", B, "\n")
#     
#     # generate data
#     errors <- if (rho == 0) rnorm(n, 0, sqrt(sigma2)) else
#       generate_ar1_errors(n, rho, sigma2)
#     y <- y_true + errors
#     
#     #  LOOCV
#     h_loocv <- bandwidth_cv(x, y, p, h_grid)$minimum
#     
#     #  MCV 
#     h_mcv <- sapply(l_values, function(l)
#       bandwidth_mcv(x, y, p, l, h_grid)$minimum)
#     
#     #  PCV 
#     h_pcv <- sapply(g_values, function(g)
#       bandwidth_pcv(x, y, p, g, h_grid)$minimum)
#     
#     #  CDPI 
#     h_cdpi <- bandwidth_cdpi(x, y)$minimum
#     
#     #  compute fits and MSE 
#     all_h <- c(h_loocv, h_mcv, h_pcv, h_cdpi)
#     h_mat[b, ]   <- all_h
#     
#     for (m in seq_along(all_h)) {
#       y_hat        <- sapply(x, function(x0)
#         local_poly_est(x0, x, y, p, all_h[m])$beta[1])
#       mse_mat[b, m] <- mean((y_hat - y_true)^2)
#     }
#   }
#   
#   # summarize
#   
#   return(
#     list(
#       h_mat = h_mat,
#       summary = data.frame(
#         n      = n,
#         rho    = rho,
#         method = method_names,
#         mean_h   = colMeans(h_mat),
#         sd_h     = apply(h_mat, 2, sd),
#         mean_mse  = colMeans(mse_mat),
#         sd_mse    = apply(mse_mat, 2, sd)
#       )
#     )
#   )
#   
# }



r_true <- function(x) {
  x + 2 * exp(-16 * x^2)
}


generate_ar1_errors <- function(n, rho, sigma2) {
  errors <- numeric(n)
  errors[1] <- rnorm(1, 0, (sigma2 / (1 - rho^2)))  
  for (i in 2:n) {
    errors[i] <- rho * errors[i-1] + rnorm(1, 0, sigma2 * (1 - rho^2))
  }
  errors
}


#' Local polynomial regression estimation (From A2)
#'
#' @param x0 constant fixed
#' @param x vector of input values
#' @param y vector of responses
#' @param p degree of regression, when p = 1 it is simply nadaraya watson
#' @param h bandwidth size
#' @param kernel kernel function, default = standard Gaussian
#'
#' @return list of regression coefficients and smoother matrix
#' @export
#'
local_poly_est <- function(x0, x, y, p, h, kernel = dnorm) {
  stopifnot(length(x0) == 1)
  
  dx <- x - x0
  x_mat <- matrix(1, nrow = length(x), ncol = p + 1)
  
  if (p > 0) {
    for (j in 1:p) x_mat[, j + 1] <- (dx^j) / factorial(j)
  }
  w  <- kernel(dx / h) / h
  WX <- x_mat * w
  XtWX <- crossprod(x_mat, WX)

  inv_XtWX <- tryCatch(
    solve(XtWX), 
    error = function(e) MASS::ginv(XtWX)
  )
  wt_x <- inv_XtWX %*% t(WX)
  
  list(beta = as.numeric(wt_x %*% y), smooth = wt_x[1, ])
}

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
      keep       <- which(abs((1:n) - j) > l)  # indices where |i - j| > l
      fit        <- local_poly_est(x[j], x[keep], y[keep], p, h_grid[i], kernel)
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


simulate_bandwidth <- function(n, rho, sigma2, B, 
                               l_values, g_values,
                               h_grid, p = 1, cl = NULL) {
  
  x      <- (1:n) / n
  y_true <- r_true(x)
  
  method_names <- c("LOOCV",
                    paste0("MCV_l", l_values),
                    paste0("PCV_g", g_values),
                    "CDPI")
  
  # Export the scenario-specific variables to the cluster workers
  clusterExport(cl, varlist = c("n", "rho", "sigma2", "x", "y_true", 
                                "l_values", "g_values", "h_grid", "p", 
                                "method_names"), 
                envir = environment())
  
  # Run the B iterations in parallel with a progress bar!
  sim_results <- pblapply(1:B, function(b) {
    
    # Generate data
    errors <- if (rho == 0) rnorm(n, 0, sqrt(sigma2)) else generate_ar1_errors(n, rho, sigma2)
    y <- y_true + errors
    
    # Bandwidth Selection
    h_loocv <- bandwidth_cv(x, y, p, h_grid)$minimum
    h_mcv   <- sapply(l_values, function(l) bandwidth_mcv(x, y, p, l, h_grid)$minimum)
    h_pcv   <- sapply(g_values, function(g) bandwidth_pcv(x, y, p, g, h_grid)$minimum)
    h_cdpi  <- bandwidth_cdpi(x, y)$minimum
    
    all_h <- c(h_loocv, h_mcv, h_pcv, h_cdpi)
    
    # Compute Fits and MSE 
    mse_vals <- numeric(length(all_h))
    for (m in seq_along(all_h)) {
      y_hat       <- sapply(x, function(x0) local_poly_est(x0, x, y, p, all_h[m])$beta[1])
      mse_vals[m] <- mean((y_hat - y_true)^2)
    }
    
    # Return the row results for this specific iteration
    list(h = all_h, mse = mse_vals)
    
  }, cl = cl) # Execute on the cluster
  
  # Reconstruct the h_mat and mse_mat from the parallelized list results
  h_mat   <- do.call(rbind, lapply(sim_results, `[[`, "h"))
  mse_mat <- do.call(rbind, lapply(sim_results, `[[`, "mse"))
  
  colnames(h_mat)   <- method_names
  colnames(mse_mat) <- method_names
  
  # Summarize and return
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
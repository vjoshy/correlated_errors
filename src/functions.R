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

fit_local_poly_curve <- function(x, y, p, h, grid, kernel = dnorm) {
  
  y_hat <- numeric(length(grid))
  
  for (i in seq_along(grid)) {
    fit <- local_poly_est(grid[i], x, y, p, h, kernel)
    y_hat[i] <- fit$beta[1]
  }
  
  data.frame(x = grid, y_hat = y_hat)
}


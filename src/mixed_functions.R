source("src/functions.R")

#' Generate data for one simulation under a linear mixed model
#' 
#' Under the linear mixed model, covariates x are generated from a standard normal distribution. Errors e are generated from a normal distribution with mean 0 and variance \eqn{\sigma^2}. Batch effect is generated from a normal distribution with mean 0 and variance \eqn{\tau^2} We then generate the response variables y using the formula y = beta0 + beta1 * x + batch effect + e, where beta0 and beta1 are the fixed intercept and slope, respectively.
#' 
#' @param b0 beta0 value for the linear model
#' @param b1 beta1 value for the linear model
#' @param sigma standard deviation of the error terms
#' @param tau standard deviation of the batch effect
#' @param k total number of batches (default is nobs/10)
#' @param nk number of observations per group
#' 
#' @return A dataframe containing the following columns: y, x, batch, e

one_sim_lmm <- function(b0=0, b1, sigma=1, tau, k, nk) {
  
  nobs <- k*nk  # k clusters each with nk obs
  
  # generate data
  x <- rnorm(nobs, 0, 1) # univariate continuous covariate
  re_batch <- rnorm(k, 0, tau) # random intercepts for each cluster
  batch <- as.factor(sort(rep(1:k, length.out = nobs)))
  re <- rep(re_batch, each = nk)
  e <- rnorm(nobs, 0, sigma) # errors
  y <- b0 + b1 * x + re + e
  
  results <- data.frame(y=y,
                        x=x, 
                        batch=batch, 
                        e=e, 
                        re=re)
}

generate_np_mixed_data <- function(k, nk, sigma = 1,
                                   tau_intercept = 1,
                                   tau_slope = 0,
                                   rho_within_group = 0,
                                   effect_type = c("random_intercept",
                                                   "random_slope"),
                                   design = c("balanced", "random"),
                                   t_range = c(0, 1),
                                   eta_fun = r_true,
                                   slope_reference = NULL){
  
  effect_type <- match.arg(effect_type)
  design <- match.arg(design)
  nobs <- k * nk
  
  batch <- factor(rep(seq_len(k), each = nk))
  
  if (is.null(slope_reference)) {
    slope_reference <- mean(t_range)
  }
  if (design == "balanced") {
    # Every batch has the same t values
    t_single <- seq(t_range[1], t_range[2], length.out = nk)
    t <- rep(t_single, times = k)
    
  } else {
    # Each batch has its own randomly generated t values
    t <- unlist(
      replicate(
        k,
        sort(runif(nk, min = t_range[1], max = t_range[2])),
        simplify = FALSE
      )
    )
  }
  
  # Evaluate the true population curve eta(t_ij)
  eta <- eta_fun(t)
  
  # Batch level random effects 
  random_int_batch <- rnorm(k, mean = 0, sd = tau_intercept)
  random_slope_batch <- rnorm(k, mean = 0, sd = tau_slope)
  
  random_int <- rep(random_int_batch, each = nk)
  random_slope <- rep(random_slope_batch, each = nk)
  
  # Batchs specific random effect curve v_i(t_ij)
  if (effect_type == "random_intercept") {
    # v_i(t_ij) = random_int_i
    random_effect <- random_int
    
  } else if (effect_type == "random_slope") {
    # v_i(t_ij) = random_int_i + random_slope_i * (t_ij - t_ref)
    t_for_slope <- t - slope_reference
    
    random_effect <- random_int + random_slope * t_for_slope
  }
  
  # Generate e_ij
  error <- numeric(nobs)

  for (i in seq_len(k)) {
    idx <- which(batch == i)
    ni <- length(idx)
    if (rho_within_group == 0) {
      error[idx] <- rnorm(ni, mean = 0, sd = sigma)
    } else {
      error[idx] <- generate_ar1_errors(
        n = ni,
        rho = rho_within_group,
        sigma2 = sigma^2
      )
    }
  }
  # Observed response
  y <- eta + random_effect + error
  
  data.frame(
    y = y,
    t = t,
    batch = batch,
    eta = eta,
    random_effect = random_effect,
    error = error,
    random_int = random_int,
    random_slope = random_slope
  )
}
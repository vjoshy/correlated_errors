source("src/functions.R")
source("src/mixed_functions.R")
source("src/bandwidth_selectors.R")

set.seed(123)

# Common simulation settings
B <- 20
k <- 20
nk <- 15
sigma <- 0.3

l_values <- c(1, 2, 3)
g_values <- c(2, 3, 4)
h_grid <- seq(0.02, 0.30, by = 0.01)
p <- 1

res_intercept_balanced <- simulate_np_mixed_bandwidth(
  k = k,
  nk = nk,
  sigma = sigma,
  tau_intercept = 0.7,
  B = B,
  l_values = l_values,
  g_values = g_values,
  h_grid = h_grid,
  tau_slope = 0,
  rho_within_group = 0,
  effect_type = "random_intercept",
  design = "balanced",
  p = p
)
res_intercept_balanced$summary

res_intercept_random <- simulate_np_mixed_bandwidth(
  k = k,
  nk = nk,
  sigma = sigma,
  tau_intercept = 0.7,
  B = B,
  l_values = l_values,
  g_values = g_values,
  h_grid = h_grid,
  tau_slope = 0,
  rho_within_group = 0,
  effect_type = "random_intercept",
  design = "random",
  p = p
)

res_intercept_random$summary

res_slope_balanced <- simulate_np_mixed_bandwidth(
  k = k,
  nk = nk,
  sigma = sigma,
  tau_intercept = 0.4,
  B = B,
  l_values = l_values,
  g_values = g_values,
  h_grid = h_grid,
  tau_slope = 0.8,
  rho_within_group = 0,
  effect_type = "random_slope",
  design = "balanced",
  p = p
)

res_slope_balanced$summary

res_slope_random <- simulate_np_mixed_bandwidth(
  k = k,
  nk = nk,
  sigma = sigma,
  tau_intercept = 0.4,
  B = B,
  l_values = l_values,
  g_values = g_values,
  h_grid = h_grid,
  tau_slope = 0.8,
  rho_within_group = 0,
  effect_type = "random_slope",
  design = "random",
  p = p
)

res_slope_random$summary

res_intercept_ar1_balanced <- simulate_np_mixed_bandwidth(
  k = k,
  nk = nk,
  sigma = sigma,
  tau_intercept = 0.7,
  B = B,
  l_values = l_values,
  g_values = g_values,
  h_grid = h_grid,
  tau_slope = 0,
  rho_within_group = 0.6,
  effect_type = "random_intercept",
  design = "balanced",
  p = p
)

res_intercept_ar1_balanced$summary


set.seed(123)

dat_check <- generate_np_mixed_data(
  k = 5,
  nk = 15,
  sigma = 0.3,
  tau_intercept = 0.4,
  tau_slope = 0.8,
  rho_within_group = 0,
  effect_type = "random_slope",
  design = "balanced"
)

head(dat_check)

truth <- data.frame(
  t = sort(unique(dat_check$t)),
  eta = r_true(sort(unique(dat_check$t)))
)

ggplot(dat_check, aes(x = t, y = y, group = batch, color = batch)) +
  geom_line(alpha = 0.7) +
  geom_point(size = 1.5) +
  geom_line(
    data = truth,
    aes(x = t, y = eta),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 1.2
  ) +
  theme_bw() +
  labs(
    title = "Balanced random-slope mixed-effects data",
    subtitle = "Black curve is the population mean eta(t); colored curves are batches",
    x = "t",
    y = "y"
  )

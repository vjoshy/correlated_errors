library(tidyverse)

summary_files <- c(
  "data/sim_results_20260323_100_500.rds",
  "data/sim_results_20260323_1000.rds"
)

hmat_files <- c(
  "data/sim_h_mats_20260323_100_500.rds",
  "data/sim_h_mats_20260323_1000.rds"
)

summary_results <- summary_files %>%
  map(readRDS) %>%
  bind_rows()

h_mats <- hmat_files %>%
  map(readRDS)

method_levels <- c("LOOCV", "MCV_l1", "MCV_l2", "MCV_l5",
                   "PCV_g3", "PCV_g5", "PCV_g10", "CDPI")

summary_results$method <- factor(summary_results$method, levels = method_levels)

# sort for cleaner plotting
summary_results <- summary_results %>%
  arrange(n, rho, method)

dir.create("plots", showWarnings = FALSE)

# 1. MCV sensitivity to l
mcv_plot_df <- summary_results %>%
  filter(str_detect(method, "^MCV_l")) %>%
  mutate(l = as.numeric(str_extract(as.character(method), "\\d+")))

p_mcv <- ggplot(mcv_plot_df,
                aes(x = l, y = mean_mse, group = factor(rho), color = factor(rho))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ n) +
  labs(
    title = "MCV sensitivity to l",
    x = "Block size l",
    y = "Mean integrated MSE",
    color = expression(rho)
  ) +
  theme_bw()

ggsave(
  filename = "plots/mcv_sensitivity_l.png",
  plot = p_mcv,
  width = 8,
  height = 5,
  dpi = 300
)

# 2. PCV sensitivity to g
pcv_plot_df <- summary_results %>%
  filter(str_detect(method, "^PCV_g")) %>%
  mutate(g = as.numeric(str_extract(as.character(method), "\\d+")))

p_pcv <- ggplot(pcv_plot_df,
                aes(x = g, y = mean_mse, group = factor(rho), color = factor(rho))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~ n) +
  labs(
    title = "PCV sensitivity to g",
    x = "Partition parameter g",
    y = "Mean integrated MSE",
    color = expression(rho)
  ) +
  theme_bw()

ggsave(
  filename = "plots/pcv_sensitivity_g.png",
  plot = p_pcv,
  width = 8,
  height = 5,
  dpi = 300
)

# 3. Mean MSE across rho
summary_small <- summary_results %>%
  filter(method %in% c("LOOCV", "MCV_l5", "PCV_g5", "CDPI"))

p_mse <- ggplot(summary_small,
                aes(x = rho, y = mean_mse, color = method, group = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  facet_wrap(~ n) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9))+
  labs(
    title = "Mean integrated MSE across dependence levels",
    x = expression(rho),
    y = "Mean integrated MSE",
    color = "Method"
  ) +
  theme_bw()

ggsave(
  filename = "plots/mean_mse_across_rho.png",
  plot = p_mse,
  width = 8,
  height = 5,
  dpi = 300
)

# 4. Mean h across rho
p_h <- ggplot(summary_small,
              aes(x = rho, y = mean_h, color = method, group = method)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  facet_wrap(~ n) +
  scale_x_continuous(breaks = c(0, 0.3, 0.6, 0.9))+
  labs(
    title = "Average selected bandwidth across dependence levels",
    x = expression(rho),
    y = "Mean selected bandwidth",
    color = "Method"
  ) +
  theme_bw()

ggsave(
  filename = "plots/mean_h_across_rho.png",
  plot = p_h,
  width = 8,
  height = 5,
  dpi = 300
)

bandwidth_table <- summary_small %>%
  filter(n == 1000, rho %in% c(0.3, 0.6, 0.9)) %>%
  select(rho, method, mean_h) %>%
  pivot_wider(names_from = method, values_from = mean_h)

bandwidth_table

fit_curve <- function(x, y, h, p = 1) {
  sapply(x, function(x0) {
    local_poly_est(x0, x, y, p = p, h = h)$beta[1]
  })
}

n <- 1000
sigma2 <- 0.5
p <- 1
x <- (1:n) / n
y_true <- r_true(x)

make_fit_plot <- function(rho_val) {
  set.seed(123)
  bw_row <- bandwidth_table %>% filter(rho == rho_val)
  
  errors_ar1 <- generate_ar1_errors(1000, rho_val, 0.5)
  y_baseline  <- y_true + errors_ar1
  
  h_loocv  <- bandwidth_cv(x, y_baseline, p = p,
                              h_grid = seq(0.01, 0.05, by = 0.001))$minimum
  
  h_mcv  <- bandwidth_mcv(x, y_baseline, p = p, l = 5,
                          h_grid = seq(0.01, 0.05, by = 0.001))$minimum
  
  h_pcv  <- bandwidth_pcv(x, y_baseline, p = p, g = 5,
                          h_grid = seq(0.01, 0.05, by = 0.001))$minimum
  
  h_cdpi <- bandwidth_cdpi(x,y_baseline)$minimum
  
  fit_loocv <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                          p = p, h = h_loocv)$beta[1])
  fit_mcv <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                     p = p, h = h_mcv)$beta[1])
  fit_pcv <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                     p = p, h = h_pcv)$beta[1])
  
  fit_cdpi <- sapply(x, function(x0) local_poly_est(x0, x, y_baseline,
                                                      p = p, h = h_cdpi)$beta[1])
  
  # put into one data frame
  plot_df <- tibble(
    x = x,
    y = y_baseline,
    true = y_true,
    loocv = fit_loocv,
    mcv = fit_mcv,
    pcv = fit_pcv,
    cdpi = fit_cdpi
  )
  
  ggplot(plot_df, aes(x = x)) +
    geom_point(aes(y = y_baseline), color = "grey70", alpha = 0.5, size = 0.7) +
    geom_line(aes(y = true, color = "True"), linewidth = 1.2) +
    geom_line(aes(y = loocv, color = "LOOCV"),  alpha=0.7,linewidth = 1) +
    geom_line(aes(y = mcv, color = "MCV_l5"), linewidth = 1) +
    geom_line(aes(y = pcv, color = "PCV_g5"),  alpha=0.7,linewidth = 1) +
    geom_line(aes(y = cdpi, color = "CDPI"), linewidth = 1) +
    labs(
      title = paste("Fit comparison for rho =", rho_val),
      x = "x",
      y = "y",
      color = NULL
    ) +
    scale_color_manual(values = c(
      "True" = "red",
      "LOOCV" = "purple",
      "MCV_l5" = "darkgreen",
      "PCV_g5" = "orange",
      "CDPI" = "skyblue"
    ))+
    theme_bw()
}

p_r03 <- make_fit_plot(0.3)
p_r06 <- make_fit_plot(0.6)
p_r09 <- make_fit_plot(0.9)

p_r03
p_r06
p_r09

ggsave("plots/fit_compare_rho03.png", p_r03, width = 8, height = 5, dpi = 300)
ggsave("plots/fit_compare_rho06.png", p_r06, width = 8, height = 5, dpi = 300)
ggsave("plots/fit_compare_rho09.png", p_r09, width = 8, height = 5, dpi = 300)



make_mcv_compare_plot <- function(rho, seed = 123) {
  set.seed(seed)
  
  errors <- generate_ar1_errors(n, rho, sigma2)
  y <- y_true + errors
  
  # bandwidths
  h_loocv <- bandwidth_cv(x, y, p = p, h_grid = h_grid)$minimum
  h_l1    <- bandwidth_mcv(x, y, p = p, l = 1, h_grid = h_grid)$minimum
  h_l2    <- bandwidth_mcv(x, y, p = p, l = 2, h_grid = h_grid)$minimum
  h_l5    <- bandwidth_mcv(x, y, p = p, l = 5, h_grid = h_grid)$minimum
  
  # fits
  y_loocv <- sapply(x, function(x0) {
    local_poly_est(x0, x, y, p = p, h = h_loocv)$beta[1]
  })
  
  y_l1 <- sapply(x, function(x0) {
    local_poly_est(x0, x, y, p = p, h = h_l1)$beta[1]
  })
  
  y_l2 <- sapply(x, function(x0) {
    local_poly_est(x0, x, y, p = p, h = h_l2)$beta[1]
  })
  
  y_l5 <- sapply(x, function(x0) {
    local_poly_est(x0, x, y, p = p, h = h_l5)$beta[1]
  })
  
  df_lines <- tibble(
    x = x,
    `True regression` = y_true,
    `LOOCV` = y_loocv,
    `MCV l=1` = y_l1,
    `MCV l=2` = y_l2,
    `MCV l=5` = y_l5
  ) %>%
    pivot_longer(-x, names_to = "method", values_to = "fit")
  
  list(
    points = x,
    lines = df_lines,
    h = c(
      LOOCV = h_loocv,
      `MCV l=1` = h_l1,
      `MCV l=2` = h_l2,
      `MCV l=5` = h_l5
    )
  )
}

plot_mcv_compare <- function(obj, rho_value) {
  ggplot() +
    geom_point(
      data = obj$points,
      aes(x = x, y = y),
      color = "grey70", alpha = 0.6, size = 0.7
    ) +
    geom_line(
      data = obj$lines,
      aes(x = x, y = fit, color = method),
      linewidth = 1
    ) +
    scale_color_manual(
      values = c(
        "True regression" = "red",
        "LOOCV" = "blue",
        "MCV l=1" = "black",
        "MCV l=2" = "darkgreen",
        "MCV l=5" = "purple"
      )
    ) +
    labs(
      title = paste0("LOOCV vs MCV for rho = ", rho_value),
      subtitle = paste0(
        "h_LOOCV = ", round(obj$h["LOOCV"], 3),
        ", h_l1 = ", round(obj$h["MCV l=1"], 3),
        ", h_l2 = ", round(obj$h["MCV l=2"], 3),
        ", h_l5 = ", round(obj$h["MCV l=5"], 3)
      ),
      x = "x",
      y = "y",
      color = NULL,
      linetype = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )
}

obj_r03 <- make_mcv_compare_plot(rho = 0.3, seed = 123)
obj_r06 <- make_mcv_compare_plot(rho = 0.6, seed = 123)
obj_r09 <- make_mcv_compare_plot(rho = 0.9, seed = 123)

p_r03_mcv <- plot_mcv_compare(obj_r03, 0.3)
p_r06_mcv <- plot_mcv_compare(obj_r06, 0.6)
p_r09_mcv <- plot_mcv_compare(obj_r09, 0.9)

p_r03_mcv
p_r06_mcv
p_r09_mcv
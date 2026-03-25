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


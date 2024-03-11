library(tidyverse)
load(file = "data/sim_pipe.rda")

sim_data <- sim_sine_6d

set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(n_jellies = c(20, 50, 100),
                      max_tries = c(50, 100),
                      d = c(6, 8, 10, 12),
                      sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

tour_level_best_basis <- sim_data |>
  group_by(id) |>
  filter(index_val == max(index_val)) |>
  filter(row_number() == 1)


# columns
# setup:      number of experiment (1-4),
# tour:       number of tour run in each setup (1-50)
# J:          number of jellyfish,
# P_J:        proportion of jellyfish found
# B:          max tries set,
# B_star:     max.tries used (doens't look right here)
# D:          the coefficient in c_t > .5 (the 3 in best_jelly - 3 * runif(1) ...)
# I_max:      best index across all the jellyfish in each trial run,
# basis:      projection basis
tour_level <- sim_data |>
  group_by(id) |>
  reframe(I_max = max(index_val),
          P_J = sum(abs(I_max - index_val) <= 0.05)/n()) |>
  left_join(tour_level_best_basis)


# columns
# P_J_hat:    proportion of trials found
# I_max_max:  best index of each setup
setup_level <- tour_level |>
  group_by(n_jellies, max_tries, d) |>
  reframe(P_J_hat = sum(abs(1 - I_max) <= 0.05)/n(), # TODO: use a better rule
          I_max_max = max(I_max))

setup_level |>
  ggplot(aes( x= as.factor(n_jellies), y = as.factor(max_tries))) +
  geom_tile(aes(fill = P_J_hat)) +
  facet_wrap(vars(d)) +
  colorspace::scale_fill_continuous_sequential(palette = "Sunset")

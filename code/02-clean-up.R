library(tidyverse)
load(file = "data/sim_pipe.rda")

tour_level_best_basis <- sim_pipe |>
  group_by(id) |>
  filter(index_val == max(index_val)) |>
  filter(row_number() == 1)

setup_level_best_basis <- sim_pipe |>
  group_by(n_jellies, max_tries, d) |>
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
tour_level_pipe <- sim_pipe |>
  group_by(id) |>
  reframe(I_max = max(index_val),
          P_J = sum(abs(I_max - index_val) <= 0.05)/n()) |>
  left_join(tour_level_best_basis)


# columns
# P_J_hat:    proportion of trials found
# I_max_max:  best index of each setup
setup_level_pipe <- tour_level_pipe |>
  group_by(n_jellies, max_tries, d) |>
  reframe(
    I_max_max = max(I_max),
    P_J_hat = sum(abs(I_max_max - I_max) <= 0.05)/n() # TODO: use a better rule
    )

setup_level_pipe |>
  ggplot(aes( x= as.factor(n_jellies), y = P_J_hat,
              group = d, color = as.factor(d))) +
  geom_line() +
  facet_wrap(vars(max_tries)) +
  #facet_wrap(vars(d)) +
  scale_color_brewer(palette = "Dark2")

#############################################################################
load(here::here("data/sim_sine_6d_dcor2d.rda"))
load(here::here("data/sim_sine_6d.rda")) # loess
load(here::here("data/sim_sine_68d_TICMIC.rda"))
load(here::here("data/sim_sine_6d_spline.rda"))

sim_data_sine <- bind_rows(sim_sine_6d_dcor2d,
                      sim_sine_6d |> mutate(idx_f = "loess2d"),
                      sim_sine_68d_TICMIC,
                      sim_sine_6d_spline |> mutate(idx_f = "spline"))


tour_level_best_basis <- sim_data_sine |>
  group_by(idx_f, id) |>
  filter(index_val == max(index_val)) |>
  filter(row_number() == 1)

setup_level_best_basis <- sim_data_sine |>
  group_by(idx_f, n_jellies, max_tries, d) |>
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
tour_level_sine <- sim_data_sine |>
  group_by(id) |>
  reframe(I_max = max(index_val),
          P_J = sum(abs(I_max - index_val) <= 0.05)/n()) |>
  left_join(tour_level_best_basis)


# columns
# P_J_hat:    proportion of trials found
# I_max_max:  best index of each setup
setup_level_sine <- tour_level_sine |>
  group_by(idx_f, n_jellies, max_tries, d) |>
  reframe(
    I_max_max = max(I_max),
    P_J_hat = sum(abs(I_max_max - I_max) <= 0.05)/n() # TODO: use a better rule
  )


sim_summary <- setup_level_pipe |> mutate(idx_f = "holes") |> bind_rows(setup_level_sine)
save(sim_summary, file = here::here("data/sim_summary.rda"))



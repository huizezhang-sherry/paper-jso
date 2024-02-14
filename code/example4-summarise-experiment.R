library(tidyverse)
load(here::here("data/pipe_raw.rda")) # loaded as b
raw <- b |> unnest(result)
raw_clean <- raw |> rename(setup = dim, tour = sim, id = loop, dim = d) |>
  mutate(setup = as.numeric(setup), tour = as.numeric(tour), id = as.numeric(id))

jelly_level <- raw_clean |> select(-c(info: alpha))

tour_level_best_basis <- jelly_level |>
  group_by(setup, tour) |>
  filter(index_val == max(index_val)) |> select(setup, tour, basis)

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
tour_level <- jelly_level |>
  group_by(setup, tour) |>
  reframe(J = 1000, B = 300, D = 3, B_star = n()/J, I_max = max(index_val),
          # TODO: use a better rule to define jellyfish found
          P_J = sum(abs(I_max - index_val) <= 0.05)/J) |>
  left_join(tour_level_best_basis) |>
  select(setup, tour, J, P_J, B, B_star, D, I_max, basis)

# columns
# P_J_hat:    proportion of trials found
# I_max_max:  best index of each setup
setup_level <- tour_level |>
  group_by(setup) |>
  reframe(P_J_hat = sum(abs(1 - I_max) <= 0.05)/n(), # TODO: use a better rule
          I_max_max = max(I_max))
save(tour_level, file = here::here("data/tour_level.rda"))
save(setup_level, file = here::here("data/setup_level.rda"))

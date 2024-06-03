library(tidyverse)
load(file = "data/sim_pipe.rda")

# all the pipe simulation are conducted in one go,
# so id is different for each n_jellies, max_tries, d combination
pipe_run_df <- sim_pipe |> get_best(group = id)
pipe_setup_df <- sim_pipe |>
  mutate(id2 = paste0(n_jellies, max_tries, d)) |>
  get_best(group = id2)

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
pipe_run_summ <- sim_pipe |>
  mutate(idx_f = "holes") |>
  group_by(id) |>
  reframe(I_max = max(index_val),
          P_J = sum(abs(I_max - index_val) <= 0.05)/n(),
          time = mean(time))

# columns
# P_J_hat:    proportion of trials found
# I_max_max:  best index of each setup
pipe_setup_summ <- pipe_run_df |>
  mutate(idx_f = "holes") |>
  group_by(n_jellies, max_tries, d, idx_f) |>
  reframe(
    I_max_max = max(index_val),
    P_J_hat = sum(abs(I_max_max - index_val) <= 0.05)/n(),
    time = mean(time)
  )



pipe_setup |>
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
load(here::here("data/sim_sine_6d_stringy.rda"))

sim_data_sine <- bind_rows(sim_sine_6d_dcor2d,
                      sim_sine_6d |> mutate(idx_f = "loess2d"),
                      sim_sine_68d_TICMIC,
                      sim_sine_6d_spline |> mutate(idx_f = "spline"),
                      sim_sine_6d_stringy |> mutate(idx_f = "stringy"),
                      )

sine_run_df <- sim_data_sine |>
  mutate(id2 = paste0(idx_f, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2)

sine_setup_df <- sim_data_sine |>
  mutate(id2 = paste0(idx_f, n_jellies, max_tries, d)) |>
  get_best(group = id2)


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
sine_run_summ <- sim_data_sine |>
  group_by(idx_f, n_jellies, max_tries, d, sim) |>
  reframe(I_max = max(index_val),
          P_J = sum(abs(I_max - index_val) <= 0.05)/n(),
          time = mean(time))


# columns
# P_J_hat:    proportion of trials found
# I_max_max:  best index of each setup
sine_setup_summ <- sine_run_df |>
  group_by(idx_f, n_jellies, max_tries, d) |>
  reframe(
    I_max_max = max(index_val),
    P_J_hat = sum(abs(I_max_max - index_val) <= 0.05)/n(),
    time = mean(time)
  )


sim_summary <- pipe_setup_summ |>  bind_rows(sine_setup_summ) |>
  rename(index = idx_f) |>
  mutate(index = ifelse(index == "spline", "splines2d", index))

############################################################################
############################################################################
load(here::here("data", "squintability.rda"))
load(here::here("data", "smoothness.rda"))
sim_df <- sim_summary |>
  left_join(smoothness |> select(n, index, smoothness) |> rename(d = n)) |>
  left_join(squintability |> select(index, n, squint) |> rename(d = n, squintability = squint)) |>
  select(index, d, I_max_max, P_J_hat, n_jellies, max_tries, smoothness, squintability, time)
save(sim_df, file = here::here("data", "sim_df.rda"))



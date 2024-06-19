library(tidyverse)

# calculate success rate for pipe-finding problem
load(file = "data-raw/sim_pipe.rda")
# all the pipe simulation are conducted in one go,
# so id is different for each n_jellies, max_tries, d combination
pipe_run_df <- sim_pipe |> get_best(group = id)
pipe_setup_df <- sim_pipe |>
  mutate(id2 = paste0(n_jellies, max_tries, d)) |>
  get_best(group = id2)

# columns
# I_max:      best index across all the simulations under each setting (n_jellies, max_tries, d)
# P_J:        success rate, proportion of simulation, out of 50 that find close (< 0.05) index value
# time:       average time across 50 simulations
pipe_setup_summ <- pipe_run_df |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(abs(I_max - index_val) <= 0.05)/n(),
    time = mean(time)
  )

#############################################################################
#############################################################################
load(here::here("data-raw/sim_sine_6d_dcor2d.rda"))
load(here::here("data-raw/sim_sine_6d_loess2d.rda"))
load(here::here("data-raw/sim_sine_68d_TICMIC.rda"))
load(here::here("data-raw/sim_sine_6d_splines2d.rda"))
load(here::here("data-raw/sim_sine_6d_stringy.rda"))

sim_data_sine <- bind_rows(
  sim_sine_6d_dcor2d, sim_sine_6d_loess2d, sim_sine_68d_TICMIC,
  sim_sine_6d_splines2d, sim_sine_6d_stringy)

sine_run_df <- sim_data_sine |>
  mutate(id2 = paste0(index, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2)
# columns
# I_max:      best index across all the simulations under each setting (n_jellies, max_tries, d)
# P_J:        success rate, proportion of simulation, out of 50 that find close (< 0.05) index value
# time:       average time across 50 simulations
sine_setup_summ <- sine_run_df |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(abs(I_max - index_val) <= 0.05)/n(),
    time = mean(time)
  )


#############################################################################
#############################################################################
# combine sine and pipe success rate
sim_summary <- pipe_setup_summ |>  bind_rows(sine_setup_summ)
save(sim_summary, file = here::here("data", "sim_summary.rda"))

############################################################################
############################################################################
load(here::here("data", "squintability.rda"))
load(here::here("data", "smoothness.rda"))
sim_df <- sim_summary |>
  left_join(smoothness |> select(n, index, smoothness) |> rename(d = n)) |>
  left_join(squintability |> select(index, n, squint) |> rename(d = n, squintability = squint))
save(sim_df, file = here::here("data", "sim_df.rda"))



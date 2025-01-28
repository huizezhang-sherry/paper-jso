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
#load(here::here("data-raw/sim_sine_6d_stringy.rda"))
load(here::here("data-raw/sim_sine_6d_skinny.rda"))
load(here::here("data-raw/sim_sine_8d_splines2d.rda"))
load(here::here("data-raw/sim_sine_8d_loess2d.rda"))
load(here::here("data-raw/sim_sine_4d_stringy.rda"))

sim_data_sine <- bind_rows(
  sim_sine_6d_dcor2d, sim_sine_6d_loess2d |> mutate(index = "loess2d"),
  sim_sine_68d_TICMIC,
  sim_sine_6d_splines2d,
  sim_sine_8d_splines2d, sim_sine_8d_loess2d)

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
    P_J = sum(abs(I_max - index_val) <= 0.05)/n()
  )

#############################################################################
#############################################################################
# process stringy and skinny
row1 <- sim_sine_6d_stringy |> mutate(index = "stringy2") |>
  mutate(id2 = paste0(index, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2) |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(index_val > 0.46)/n()
  )

row4 <- sim_sine_4d_stringy |> mutate(index = "stringy2") |>
  mutate(id2 = paste0(index, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2) |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(index_val > 0.595)/n()
  )


row2 <- sim_sine_6d_skinny |> mutate(index = "skinny") |>
  mutate(id2 = paste0(index, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2) |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(index_val > 0.67)/n()
  )

row3 <- sim_sine_6d_splines2d_100 |> mutate(index = "splines2d") |>
  mutate(id2 = paste0(index, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2) |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(abs(I_max - index_val) <= 0.05)/n()
  )

#############################################################################
#############################################################################
# process 4d data
load(here::here("data-raw/sim_pipe_4d.rda"))
load(here::here("data-raw/sim_sine_4d_dcor2d.rda"))
load(here::here("data-raw/sim_sine_4d_loess2d.rda"))
load(here::here("data-raw/sim_sine_4d_MICTIC.rda"))
load(here::here("data-raw/sim_sine_4d_splines2d.rda"))
pipe_setup_summ_4d <- sim_pipe_4d |>
  mutate(id2 = paste0(index, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2) |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(abs(I_max - index_val) <= 0.05)/n(),
    time = mean(time)
  )

sim_data_sine_4d <- bind_rows(
  sim_sine_4d_dcor2d |> mutate(index = "dcor2d_2") , sim_sine_4d_loess2d, sim_sine_4d_MICTIC,
  sim_sine_4d_splines2d)

sine_run_df_4d <- sim_data_sine_4d |>
  mutate(id2 = paste0(index, n_jellies, max_tries, sim, d)) |>
  get_best(group = id2)

sine_setup_summ_4d <- sine_run_df_4d |>
  group_by(index, d, n_jellies, max_tries) |>
  reframe(
    I_max = max(index_val),
    P_J = sum(abs(I_max - index_val) <= 0.05)/n(),
    time = mean(time)
  )

#############################################################################
#############################################################################
# combine sine and pipe success rate
sim_summary <- sine_setup_summ |> filter(index != "stringy") |>
  bind_rows(row1, row2, row3, row4, pipe_setup_summ_4d, sine_setup_summ_4d,
            pipe_setup_summ) |> select(-time) |>
  arrange(index, d)
save(sim_summary, file = here::here("data", "sim_summary.rda"))

############################################################################
############################################################################
load(here::here("data", "squintability.rda"))
load(here::here("data", "smoothness.rda"))
sim_df <- sim_summary |>
  left_join(smoothness |> mutate(smoothness = rank(smoothness)) |> select(n, index, smoothness) |> rename(d = n)) |>
  left_join(squintability |> mutate(squint = rank(squint)) |> select(index, n, squint) |> rename(d = n, squintability = squint))
save(sim_df, file = here::here("data", "sim_df.rda"))



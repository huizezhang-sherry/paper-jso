library(tidyverse)
library(ferrn)
library(tourr)

sim <- function(d = d, n_jellies = n_jellies, max.tries = max.tries,
                optim_seed = seed){
  if (d == 6){
    pipe1000 <- ferrn::pipe1000_6d
  } else if (d == 8){
    pipe1000 <- ferrn::pipe1000_8d
  } else if (d == 10){
    pipe1000 <- ferrn::pipe1000_10d
  } else if (d == 12){
    pipe1000 <- ferrn::pipe1000_12d
  }


  cat("n_jellies: ", n_jellies, "\n")
  cat("max.tries: ", max.tries, "\n")
  cat("d: ", d, "\n")
  res <- map_dfr(optim_seed, function(seed){
    t1 <- Sys.time()
    set.seed(seed)
    res <- animate_xy(pipe1000, guided_tour(
      holes(), d = 2, n_jellies = n_jellies,
      search_f =  search_jellyfish,
      max.tries = max.tries
    ), max_frames = max.tries) |> dplyr::select(-id)
    t2 <- Sys.time()
    res <- res |> mutate(time = t2 - t1)
    return(res)
  })

  return(res)
}

set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(n_jellies = c(20, 50, 100),
                      max_tries = c(50, 100),
                      d = c(6, 8, 10, 12),
                      sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, n_jellies, max_tries, optim_seed = seed))) |>
  ungroup() |>
  unnest(res)

sim_pipe_raw <- sim_res |> select(index, d, n_jellies, max_tries, sim:time)
sim_pipe <- tibble(index = "holes") |> bind_cols(sim_pipe_raw)
save(sim_pipe, file = "data-raw/sim_pipe.rda")

################################################################################
################################################################################
# 4D data is simulated in 10-sim-4d.R
load(here::here("data-raw/sim_pipe_4d.rda"))
sim_pipe_run_best <- bind_rows(sim_pipe_4d, sim_pipe) |>
  group_by(id, n_jellies, max_tries, d) |>
  summarise(I_max = max(index_val))
save(sim_pipe_run_best, file = "data/sim_pipe_run_best.rda")

pipe_jellyfish <- sim_pipe |>
  filter(n_jellies == 100, max_tries == 100) |>
  get_best(group = id) |>
  ungroup() |>
  mutate(optimiser = "jellyfish")
save(pipe_jellyfish, file = "data/pipe_jellyfish.rda")



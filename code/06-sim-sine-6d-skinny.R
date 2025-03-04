library(tidyverse)
library(ferrn)
library(tourr)

################################################################################
# d = 6
skinny <- function(){
  function(mat){
    cassowaryr::sc_skinny(mat[,1], mat[,2])
  }
}

sim <- function(d = d, n_jellies = n_jellies, max.tries = max.tries,
                optim_seed = seed, sim = sim){

  if (d == 4){
    sine1000 <- ferrn::sine1000_4d
  } else if (d == 6){
    sine1000 <- ferrn::sine1000_6d
  }

  cat("sim: ", sim, "\n")
  cat("n_jellies: ", n_jellies, "\n")
  cat("max.tries: ", max.tries, "\n")
  cat("d: ", d, "\n")
  res <- map_dfr(optim_seed, function(seed){
      t1 <- Sys.time()
      set.seed(seed)
      res <- animate_xy(sine1000, guided_tour(
        skinny(), d = 2, n_jellies = n_jellies,
        search_f =  search_jellyfish,
        max.tries = max.tries
      ), max_frames = max.tries) |> dplyr::select(-id)
      t2 <- Sys.time()
      res <- res |> mutate(time = t2 - t1)
      return(res)
  })

  return(res)
}

set.seed(1234)
seed <- sample(1000:10000, size = 50)
sim_setup <- tidyr::crossing(n_jellies = 30, max_tries = 50, d = 6) |>
  tidyr::crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
t1
set.seed(123)
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1


sim_sine_6d_skinny <- bind_rows(sim_sine_6d_skinny_0105, sim_sine_6d_skinny_0610,
                                sim_sine_6d_skinny_1120, sim_sine_6d_skinny_2130,
                                sim_sine_6d_skinny_3140, sim_res3233,
                                sim_sine_6d_skinny_4150, sim_res48) |>
  arrange(sim)
save(sim_sine_6d_skinny, file = here::here("data/sim_sine_6d_skinny.rda"))

##################################################
set.seed(1234)
seed <- sample(1000:10000, size = 50)
sim_setup <- tidyr::crossing(n_jellies = 30,
                             max_tries = 50,
                             d = 4) |>
  tidyr::crossing(sim = 31:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
t1
set.seed(123)
sim_sine_4d_skinny_1130 <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1
save(sim_sine_4d_skinny_3150, file = here::here("data-raw/sim_sine_4d_skinny_3150.rda"))

sim_sine_4d_skinny <- bind_rows(
  sim_sine_4d_skinny_0110, sim_sine_4d_skinny_1130, sim_sine_4d_skinny_3150
  ) |> arrange(sim)
save(sim_sine_4d_skinny, file = here::here("data-raw/sim_sine_4d_skinny.rda"))
save(sim_sine_4d_skinny_0110, file = here::here("data-raw/sim_sine_4d_skinny_0110.rda"))
save(sim_sine_4d_skinny_1130, file = here::here("data-raw/sim_sine_4d_skinny_1130.rda"))

##################################################
# diagnostics
library(ferrn)
sim_sine_4d_skinny |> get_best(sim) |> plot_projection(sine1000_4d)

sim_sine_6d_skinny |> filter(sim == 2) |> plot_projection(sine1000, id = loop, animate_along = tries)
anim_save(filename = "skinny_sim1.gif")



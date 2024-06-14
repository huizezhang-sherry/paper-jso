library(tidyverse)
library(spinebil)
library(tourr)

################################################################################
# d = 6
sim <- function(d = d, n_jellies = n_jellies, max.tries = max.tries,
                data_seed = 123456, optim_seed = seed, sim = sim){
  set.seed(data_seed)
  sine1000 <- sinData(d, 1000) %>% scale() %>% as_tibble()
  colnames(sine1000) <- paste0("V", 1:d)

  cat("sim: ", sim, "\n")
  cat("n_jellies: ", n_jellies, "\n")
  cat("max.tries: ", max.tries, "\n")
  cat("d: ", d, "\n")
  res <- map_dfr(optim_seed, function(seed){
    t1 <- Sys.time()
    set.seed(seed)
    res <- animate_xy(sine1000, guided_tour(
      tourr::stringy(), d = 2, n_jellies = n_jellies,
      search_f =  search_jellyfish,
      max.tries = max.tries
    ), max_frames = max.tries) |> dplyr::select(-id)
    t2 <- Sys.time()
    res <- res |> mutate(time = t2 - t1)
    return(res)
  })

  return(res)
}

# separate into sim 1-30 and 31-50 to run, each takes < 2 hrs
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(n_jellies = 50,
                      max_tries = 50,
                      d = 6) |>
  crossing(sim = 21:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
set.seed(123)
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1

#sim_sine_6d_stringy12 <- sim_res |> select(-alpha)
#sim_sine_6d_stringy310 <- sim_res |> select(-alpha)
#sim_sine_6d_stringy1120 <- sim_res |> select(-alpha)
sim_sine_6d_stringy2150 <- sim_res |> select(-alpha)

sim_sine_6d_stringy <- bind_rows(sim_sine_6d_stringy12, sim_sine_6d_stringy310,
                                 sim_sine_6d_stringy1120, sim_sine_6d_stringy2150)

sim_sine_6d_stringy <- tibble(index = "stringy") |>
  bind_cols(sim_sine_6d_stringy) |>
  select(index, d, n_jellies, max_tries, sim:time)
save(sim_sine_6d_stringy, file = here::here("data-raw/sim_sine_6d_stringy.rda"))



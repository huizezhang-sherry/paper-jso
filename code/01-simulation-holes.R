library(tidyverse)
library(spinebil)
library(tourr)

sim <- function(d = d, n_jellies = n_jellies, max.tries = max.tries,
                data_seed = 123456, optim_seed = seed){
  set.seed(data_seed)
  pipe1000 <- pipeData(d, 1000) %>% scale() %>% as_tibble()
  colnames(pipe1000) <- paste0("V", 1:d)

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

# sim_pipe <- sim_res |> select(-alpha)
# save(sim_pipe, file = "data/sim_pipe.rda")


pipe_jellyfish <- sim_pipe |>
  filter(n_jellies == 100, max_tries == 100) |>
  get_best(group = id) |>
  ungroup() |>
  mutate(optimiser = "jellyfish")
# save(pipe_jellyfish, file = "data/pipe_jellyfish.rda")

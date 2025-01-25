library(tidyverse)
library(spinebil)
library(tourr)

skinny <- function(){
  function(mat){
    cassowaryr::sc_skinny(mat[,1], mat[,2])
  }
}

stringy2 <- function(){
  function(mat){
    x <- cassowaryr::scree(mat[,1], mat[,2])
    b <- cassowaryr:::gen_mst(x$del, x$weights)
    diameter <- length(igraph::get_diameter(b))
    length <-  length(b) - 1
    diameter / length
  }
}


sim_pipe <- function(d = d, n_jellies = n_jellies, max.tries = max.tries,
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

t1 <- Sys.time()
set.seed(123)
seed <- sample(1000:10000, size = 50)
sim_setup <- tidyr::crossing(n_jellies = c(20, 50, 100), max_tries = c(50, 100),
                             d = 4, sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim_pipe(d, n_jellies, max_tries, optim_seed = seed))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1

sim_pipe <- sim_res |> mutate(index = "holes") |> select(index, d, n_jellies, max_tries, sim:time)
save(sim_pipe, file = "data-raw/sim_pipe_4d.rda")

################################################################################
################################################################################
sim <- function(index, d = d, n_jellies = n_jellies, max.tries = max.tries,
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
      get(index)(), d = 2, n_jellies = n_jellies,
      search_f =  search_jellyfish,
      max.tries = max.tries, verbose = TRUE,
    ), max_frames = max.tries) |> dplyr::select(-id)
    t2 <- Sys.time()
    res <- res |> mutate(time = t2 - t1)
    return(res)
  })

  return(res)
}

set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- tidyr::crossing(n_jellies = 20, max_tries = 50, d = 4) |>
  tidyr::crossing(sim = 1:2) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
set.seed(123)
sim_res <- sim_setup |>
  mutate(index = "splines2d") |>
  rowwise() |>
  mutate(res = list(sim(index = index, d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1
sim_sine_4d_splines2d <- sim_res
save(sim_sine_4d_splines2d, file = "data-raw/sim_sine_4d_splines2d.rda")

################################################################################
################################################################################
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- tidyr::crossing(n_jellies = c(10, 20), max_tries = c(20, 30), d = 4) |>
  tidyr::crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
set.seed(123)
sim_res <- sim_setup |>
  mutate(index = "loess2d") |>
  rowwise() |>
  mutate(res = list(sim(index = index, d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1
sim_sine_4d_loess2d <- sim_res
save(sim_sine_4d_loess2d, file = "data-raw/sim_sine_4d_loess2d.rda")

################################################################################
################################################################################
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- tidyr::crossing(index = c("MIC", "TIC"),
                      d = 4,
                      n_jellies = c(50, 100),
                      max_tries = c(30, 50)) |>
  tidyr::crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
set.seed(123)
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(index = index, d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1
sim_sine_4d_MICTIC <- sim_res
save(sim_sine_4d_MICTIC, file = "data-raw/sim_sine_4d_MICTIC.rda")

################################################################################
################################################################################
# 50 mins
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- tidyr::crossing(index = "dcor2d", d = 4, n_jellies = c(20), max_tries = c(30, 50)) |>
  tidyr::crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
set.seed(123)
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(index = index, d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1
sim_sine_4d_dcor2d <- sim_res
save(sim_sine_4d_dcor2d, file = "data-raw/sim_sine_4d_dcor2d.rda")


library(cassowaryr)
library(igraph)
library(tidyverse)
library(spinebil)
library(tourr)

set.seed(123456)
sine1000 <- sinData(6, 1000) %>% scale() %>% as_tibble()
colnames(sine1000) <- paste0("V", 1:6)
################################################################################
# d = 6
stringy2 <- function(){
  function(mat){
    x <- cassowaryr::scree(mat[,1], mat[,2])
    b <- cassowaryr:::gen_mst(x$del, x$weights)
    diameter <- length(igraph::get_diameter(b))
    length <-  length(b) - 1
    diameter / length
  }
}

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
        stringy2(), d = 2, n_jellies = n_jellies,
        search_f =  search_jellyfish,
        max.tries = max.tries, verbose = TRUE
      ), max_frames = max.tries, verbose = TRUE) |> dplyr::select(-id)
      t2 <- Sys.time()
      res <- res |> mutate(time = t2 - t1)
      return(res)
  })

  return(res)
}

set.seed(1234)
seed <- sample(1000:10000, size = 50)
sim_setup <- tidyr::crossing(n_jellies = 50,
                      max_tries = 50,
                      d = 6) |>
  tidyr::crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
t1
set.seed(123)
sim_res_150 <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1

sim_sine_6d_stringy <- sim_res_150 |> mutate(index = "stringy")
save(sim_sine_6d_stringy, file = here::here("data/sim_sine_6d_stringy.rda"))


##################################################
# diagnostics
sim_sine_6d_stringy |> get_best(sim) |> plot_projection(sine1000)

sim_sine_6d_stringy |> filter(sim == 2) |> plot_projection(sine1000, id = loop, animate_along = tries)
anim_save(filename = "stringy_sim1.gif")



library(tidyverse)
library(ferrn)
library(tourr)

################################################################################
# d = 6
# d = 4 data is simulated in 10-sim-4d.R
sim <- function(d = d, idx_f = idx_f, n_jellies = n_jellies, max.tries = max.tries,
                data = ferrn::sine1000_6d, optim_seed = seed, sim = sim){

  cat("sim:", sim, "| n_jellies:", n_jellies, "| max.tries:", max.tries,
  "| idx_f:", idx_f, "| d:", d, "\n")

  idx_f <- eval(sym(idx_f))
  res <- map_dfr(optim_seed, function(seed){
    t1 <- Sys.time()
    set.seed(seed)
    res <- animate_xy(data, guided_tour(
      idx_f(), d = 2, n_jellies = n_jellies,
      search_f =  search_jellyfish,
      max.tries = max.tries
    ), max_frames = max.tries) |> dplyr::select(-id)
    t2 <- Sys.time()
    res <- res |> mutate(time = t2 - t1)
    return(res)
  })

  return(res)
}

dcor2d_2 <- function() {
  function(mat) {
    xy <- na.omit(data.frame(x = mat[, 1], y = mat[, 2]))
    measure <- with(xy, energy::dcor2d(x, y, type = "U"))
    return(measure)
  }
}

#  n_jellies | max_tries
#         20 |       50
#         20 |      100
#         50 |       50
#         50 |      100
#        100 |       50
#        100 |      100
# each 50 repetitions
# take less than 2 hrs
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(idx_f = c("dcor2d_2"), d = 6,
                      n_jellies = c(20, 50),
                      max_tries = c(50, 100)) |>
  filter(!(n_jellies == 50 & max_tries == 100)) |>
  crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, idx_f, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1


sim_sine_6d_dcor2d <- sim_res |> select(-alpha) |> rename(index = idx_f)
save(sim_sine_6d_dcor2d , file = "data-raw/sim_sine_6d_dcor2d.rda")

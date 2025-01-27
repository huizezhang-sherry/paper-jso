library(tidyverse)
library(spinebil)
library(tourr)
library(minerva)

################################################################################
# d = 6
sim <- function(d = d, idx_f = idx_f, n_jellies = n_jellies, max.tries = max.tries,
                data_seed = 123456, optim_seed = seed, sim = sim){
  set.seed(data_seed)
  sine1000 <- sinData(d, 1000) %>% scale() %>% as_tibble()
  colnames(sine1000) <- paste0("V", 1:d)

  cat("sim: ", sim, "\n")
  cat("n_jellies: ", n_jellies, "\n")
  cat("max.tries: ", max.tries, "\n")
  cat("idx_f: ", idx_f, "\n")
  cat("d: ", d, "\n")

  idx_f <- eval(sym(idx_f))
  res <- map_dfr(optim_seed, function(seed){
    t1 <- Sys.time()
    set.seed(seed)
    res <- animate_xy(sine1000, guided_tour(
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

MIC <- function(){
  function(mat){
    minerva::mine(mat[,1], mat[,2], alpha = 0.3,  est = "mic_e")$MIC
  }
}

TIC <- function(){
  function(mat){
    minerva::mine(mat[,1], mat[,2], est = "mic_e", alpha = 0.3)$TIC
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
# take 2 hrs
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(idx_f = c("MIC", "TIC"),
                      d = c(6, 8),
                      n_jellies = c(20, 50),
                      max_tries = c(50, 100)) |>
  crossing(sim = 1) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, idx_f, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1

sim_sine_68d_TICMIC <- sim_res |> select(-alpha) |> rename(index = idx_f)
save(sim_sine_68d_TICMIC, file = "data-raw/sim_sine_68d_TICMIC.rda")


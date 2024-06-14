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
      loess2d(), d = 2, n_jellies = n_jellies,
      search_f =  search_jellyfish,
      max.tries = max.tries
    ), max_frames = max.tries) |> dplyr::select(-id)
    t2 <- Sys.time()
    res <- res |> mutate(time = t2 - t1)
    return(res)
  })

  return(res)
}

loess2d <- function() {
  function(mat) {
    mat <- as.data.frame(mat)
    colnames(mat) <- c("x", "y")
    loess_fit <- loess(y ~ x, data = mat, span = 0.05)
    loess_fit2 <- loess(x ~ y, data = mat, span = 0.05)
    measure <- max(1 - var(residuals(loess_fit), na.rm = T) / var(mat$y, na.rm = T),
                   1 - var(residuals(loess_fit2), na.rm = T) / var(mat$y, na.rm = T)
    )
    return(measure)
  }
}

#  n_jellies | max_tries
#         20 |       50
#         20 |      100
#         50 |       50
# each 50 repetitions
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(n_jellies = c(20, 50),
                      max_tries = c(50, 100),
                      d = 6) |>
  filter(row_number() != 4) |>
  crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)

sim_sine_6d <- sim_res |> select(-alpha)
sim_sine_6d_loess2d <- tibble(index = "loess") |>
  bind_cols(sim_sine_6d) |>
  select(index, d, n_jellies, max_tries, sim:time)
save(sim_sine_6d_loess2d, file = "data-raw/sim_sine_6d_loess2d.rda")



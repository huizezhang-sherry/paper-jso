library(tidyverse)
library(spinebil)
library(tourr)

set.seed(123)
seed <- sample(1000: 10000, size = 50)
historyarray2dt <- function(array, index_f){
  stopifnot(inherits(array, "history_array"))
  data <- attr(array, "data")
  n <- dim(array)[3]
  res <- purrr::map_dfr(1:n, ~tibble(basis = list(unclass(array)[,,.x])))
  projs <-  map(res$basis, ~data %*% .x)
  res <- res |>
    dplyr::mutate(index_val = map_dbl(projs, ~index_f(.x)),
                  id = dplyr::row_number())

  return(res)

}

sim <- function(d, optimiser = c("search_better_random", "search_better", "search_jellyfish"),
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


  if (optimiser == "search_better_random"){
    res <- map_dfr(optim_seed, function(seed){
      set.seed(seed)
      spiral_better <- save_history(
        pipe1000, guided_tour(
          holes(), search_f = search_better_random, max.tries = 1000
        )
      )
      historyarray2dt(spiral_better[,,-1], holes())
    }, .id = "sim") |>
      mutate(d = d)
  } else if (optimiser == "search_better"){
    res <- map_dfr(optim_seed, function(seed){
      set.seed(seed)
      spiral_better <- save_history(
        pipe1000, guided_tour(
          holes(), search_f = search_better, max.tries = 1000
        )
      )
      historyarray2dt(spiral_better[,,-1], holes())
    }, .id = "sim") |>
      mutate(d = d)
  }
  return(res)
}

data <- tibble(x = list(ferrn::pipe1000_4d, ferrn::pipe1000_6d,
                        ferrn::pipe1000_10d, ferrn::pipe1000_12d),
               dim = c(6, 8, 10, 12))
# a <- map_dfr(c(6, 8, 10, 12),
#              ~safe_sim(d = .x, optimiser = "search_better",
#                        optim_seed = seed),
#              .id = "dim"
#              )
# take around an hour to run
c <- map_dfr(c(6, 8, 10, 12),
             ~sim(d = .x, optimiser = "search_better_random",
                       optim_seed = seed),
             .id = "dim"
             )


pipe_better <- a |>
  group_by(sim, dim) |>
  filter(index_val == max(index_val), row_number() == max(row_number())) |>
  left_join(data, c("d" = "dim")) |>
  mutate(proj = list(as_tibble(x[[1]] %*% basis[[1]])), optimiser = "better") |>
  dplyr::select(-x, -id, -d) |>
  ungroup()

pipe_better_random <- c |>
  group_by(sim, dim) |>
  filter(index_val == max(index_val)) |>
  filter(row_number() == max(row_number())) |>
  left_join(data, c("d" = "dim")) |>
  mutate(proj = list(as_tibble(x[[1]] %*% basis[[1]])), optimiser = "better_random") |>
  dplyr::select(-x, -id, -d) |>
  ungroup()
# save(pipe_better, file = "data/pipe_better.rda")
# save(pipe_better_random, file = here::here("data/pipe_better_random.rda"))

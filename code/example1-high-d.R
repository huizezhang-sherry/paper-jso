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

sim <- function(d, optimiser = c("search_better", "search_jellyfish"), data_seed = 123456, optim_seed = seed){
  set.seed(data_seed)
  pipe1000 <- pipeData(d, 1000) %>% scale() %>% as_tibble()
  colnames(pipe1000) <- paste0("V", 1:d)

  if (optimiser == "search_better"){
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
  } else if (optimiser == "search_jellyfish"){

    res <- map_dfr(optim_seed, function(seed){
      set.seed(seed)
        animate_xy(pipe1000, guided_tour(
          holes(), d = 2, n_jellies = 1000,
          search_f =  search_jellyfish,
          max.tries = 300, min.tries = 100
        ), max_frames = 50)
    }, .id = "sim") |>
      mutate(d = d)

  }

  return(res)
}

data <- map_dfr(c(6, 8, 10, 12), ~{
  set.seed(123456)
  tibble(x = list(pipeData(.x, 1000) %>% scale()))
}) |>
  mutate(dim = c(6, 8, 10, 12))
safe_sim <- safely(sim)
#a <- map_dfr(c(6, 8, 10, 12), ~safe_sim(d = .x, optimiser = "search_better", optim_seed = seed), .id = "dim")
# take around an hour to run
b <- map_dfr(
  c(6, 8, 10, 12),
  ~safe_sim(d = .x, optimiser = "search_jellyfish", optim_seed = seed),
  .id = "dim")

#save(b, file = "data/pipe_raw.rda")
# pipe_better <- a |>
#   unnest(result) |>
#   group_by(sim, dim) |>
#   filter(index_val == max(index_val), row_number() == max(row_number())) |>
#   left_join(data, c("d" = "dim")) |>
#   mutate(proj = list(as_tibble(x[[1]] %*% basis[[1]])), optimiser = "better") |>
#   dplyr::select(-x, -id, -d) |>
#   ungroup()
# save(pipe_better, file = "data/pipe_better.rda")
# pipe_jellyfish <- b |>
#   unnest(result) |>
#   group_by(sim, dim) |>
#   filter(index_val == max(index_val)) |>
#   left_join(data, c("d" = "dim")) |>
#   mutate(proj = list(as_tibble(x[[1]] %*% basis[[1]])), optimiser = "jellyfish") |>
#   dplyr::select(dim, sim, basis, index_val, proj, optimiser) |>
#   ungroup()
# save(pipe_jellyfish, file = "data/pipe_jellyfish.rda")
#
#

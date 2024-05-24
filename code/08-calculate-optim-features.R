library(ferrn)
library(tourr)
library(tidyverse)

smoothness_holes <- tibble::tibble(
  n = c(6, 8, 10, 12),
  index = rep("holes", 4),
  data = list(pipe1000, pipe1000_8d, pipe1000_10d, pipe1000_12d),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE),
              matrix(c(rep(0, 12), 1, 0, 0, 1), nrow = 8, byrow = TRUE),
              matrix(c(rep(0, 16), 1, 0, 0, 1), nrow = 10, byrow = TRUE),
              matrix(c(rep(0, 20), 1, 0, 0, 1), nrow = 12, byrow = TRUE))) |>
  dplyr::mutate(smoothness = purrr::pmap(
    list(data, best, n),
    function(data, best, n){calc_smoothness("holes", data = data, n = n, best = best)}))

holes_tidy <- smoothness_holes |> tidyr::unnest_wider(smoothness) |> unnest(measure)

idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "stringy", "splines2d")
smoothness_sine <- tibble::tibble(
  n = 6, index = idx_names, data = list(sine1000),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE))) |>
  rowwise() |>
  dplyr::mutate(smoothness = list(calc_smoothness(index)))


sine_tidy <- smoothness_sine |> tidyr::unnest_wider(smoothness) |> unnest(measure)

smoothness <- bind_rows(sine_tidy, holes_tidy)
save(smoothness, file = here::here("data", "smoothness.rda"))


################################################################################
################################################################################
load(here::here("data", "sim_summary.rda"))
sim_summary |> rename(index = idx_f) |>
  left_join(smoothness |> select(n, index, smoothness) |> rename(d = n))




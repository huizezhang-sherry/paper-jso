library(ferrn)
library(tourr)
library(tidyverse)

sine_8d_tbl <- function(index){
  tibble::tibble(
    n = 8, index = index, data = list(sine1000_8d),
    best = list(matrix(c(rep(0, 12), 1, 0, 0, 1), nrow = 8, byrow = TRUE)))
}

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
  dplyr::bind_rows(sine_8d_tbl("MIC")) |>
  dplyr::bind_rows(sine_8d_tbl("TIC")) |>
  rowwise() |>
  dplyr::mutate(smoothness = list(calc_smoothness(index)))


sine_tidy <- smoothness_sine |> tidyr::unnest_wider(smoothness) |> unnest(measure)

smoothness <- bind_rows(sine_tidy, holes_tidy) |>
  mutate(index = factor(index, levels = c("holes", "MIC", "TIC", "dcor2d_2",
                                          "loess2d", "splines2d", "stringy"))) |>
  arrange(index) |>
  select(index, n, variance:nugget)
save(smoothness, file = here::here("data", "smoothness.rda"))

################################################################################
################################################################################
# squintability
sq_holes_basis_df <- tibble::tibble(
  n = c(6, 8, 10, 12),
  index = rep("holes", 4),
  data = list(pipe1000, pipe1000_8d, pipe1000_10d, pipe1000_12d),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE),
              matrix(c(rep(0, 12), 1, 0, 0, 1), nrow = 8, byrow = TRUE),
              matrix(c(rep(0, 16), 1, 0, 0, 1), nrow = 10, byrow = TRUE),
              matrix(c(rep(0, 20), 1, 0, 0, 1), nrow = 12, byrow = TRUE))) |>
  dplyr::mutate(basis_df = purrr::pmap(
    list(data, best, n), function(data, best, n){
      calc_squintability("holes", data = data, n = n, best = best, return_early = TRUE)}))

idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "splines2d", "stringy")
sq_sine_basis_df <- tibble::tibble(
  n = 6, index = idx_names, data = list(sine1000),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE))) |>
  dplyr::bind_rows(sine_8d_tbl("MIC")) |>
  dplyr::bind_rows(sine_8d_tbl("TIC")) |>
  dplyr::rowwise() |>
  dplyr::mutate(basis_df = list(calc_squintability(
    index, data = data, n = n, best = best, return_early = TRUE)
  ))

sq_basis_df <- sq_holes_basis_df |>
  dplyr::bind_rows(sq_sine_basis_df) |>
  unnest(basis_df) |>
  group_by(n, index) |>
  mutate(
    holes = ifelse(
      index == "holes",
      (holes - min(holes, na.rm = TRUE))/
        (max(holes, na.rm = TRUE) - min(holes, na.rm = TRUE)),
      NA),
    stringy = ifelse(
      index == "stringy",
      (stringy - min(stringy, na.rm = TRUE))/
        (max(stringy, na.rm = TRUE) - min(stringy, na.rm = TRUE)),
      NA))|>
  nest(basis_df = c(id:stringy)) |>
  ungroup()
save(sq_basis_df, file = here::here("data/sq_basis_df.rda"))


start_prms1 <- c(theta1 = 1, theta2 = 1, theta3 = 3, theta4 = 0)
start_prms2 <- c(theta1 = 1, theta2 = 0.01, theta3 = 50, theta4 = 0)
squintability <- sq_basis_df |>
  dplyr::mutate(start = ifelse(index == "stringy", list(start_prms2),
                               list(start_prms1))) |>
  dplyr::rowwise() |>
  dplyr::mutate(res = list(ferrn:::fit_nls(basis_df, idx = index,
                                           nls_params = list(start)))) |>
  tidyr::unnest_wider(res) |>
  mutate(index = factor(index, levels = c("holes", "MIC", "TIC", "dcor2d_2",
                                          "loess2d", "splines2d", "stringy"))) |>
  arrange(index) |>
  select(index, n, theta1: theta4) |>
  mutate(squint = abs(theta1 * theta2 * theta3 / 4))
save(squintability, file = here::here("data", "squintability.rda"))



############################################################################
############################################################################
sq_basis_dist_idx <- sq_basis_df |>
  dplyr::select(basis_df, n) |>
  unnest(basis_df) |>
  pivot_longer(holes:stringy, names_to = "index", values_to = "idx",
               values_drop_na = TRUE) |>
  mutate(dist = ceiling(dist / 0.005) * 0.005,
         index = factor(index, levels = c("holes", "MIC", "TIC", "dcor2d_2",
                                          "loess2d", "splines2d", "stringy"))) |>
  group_by(n, index, dist) |>
  summarise(idx = mean(idx), .groups = "drop")
save(sq_basis_dist_idx, file = here::here("data", "sq_basis_dist_idx.rda"))

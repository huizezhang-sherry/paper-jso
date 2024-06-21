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
  dplyr::mutate(basis_df = purrr::pmap(
    list(data, best),
    function(data, best, n){sample_bases("holes", data = data, best = best)}))

holes_tidy <- smoothness_holes |>
  rowwise() |>
  dplyr::mutate(smooth = calc_smoothness(basis_df)) |>
  unnest(smooth)

idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "stringy", "splines2d")
smoothness_sine <- tibble::tibble(
  n = 6, index = idx_names, data = list(sine1000),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE))) |>
  dplyr::bind_rows(sine_8d_tbl("MIC")) |>
  dplyr::bind_rows(sine_8d_tbl("TIC")) |>
  rowwise() |>
  dplyr::mutate(basis_df = list(sample_bases(index)))

sine_tidy <- smoothness_sine |>
  rowwise() |>
  dplyr::mutate(smooth = calc_smoothness(basis_df)) |>
  unnest(smooth)

smoothness <- bind_rows(sine_tidy, holes_tidy) |>
  mutate(index = factor(index, levels = c("holes", "MIC", "TIC", "dcor2d_2",
                                          "loess2d", "splines2d", "stringy"))) |>
  arrange(index) |>
  select(index, n, variance:nugget)
save(smoothness, file = here::here("data", "smoothness.rda"))

################################################################################
################################################################################
# squintability - sample bases
# 2-3 mins
sq_holes_basis_df <- tibble::tibble(
  n = c(6, 8, 10, 12),
  index = rep("holes", 4),
  data = list(pipe1000, pipe1000_8d, pipe1000_10d, pipe1000_12d),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE),
              matrix(c(rep(0, 12), 1, 0, 0, 1), nrow = 8, byrow = TRUE),
              matrix(c(rep(0, 16), 1, 0, 0, 1), nrow = 10, byrow = TRUE),
              matrix(c(rep(0, 20), 1, 0, 0, 1), nrow = 12, byrow = TRUE))) |>
  dplyr::mutate(basis_df = purrr::pmap(
    list(data, best), function(data, best){
      sample_bases("holes", data = data, n_basis = 50, min_proj_dist = 1.5,
                   step_size = 0.005, best = best)}))
save(sq_holes_basis_df, file = here::here("data/sq_holes_basis_df.rda"))

idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "splines2d", "stringy")
# about 40 mins
sq_sine_basis_df <- tibble::tibble(
  n = 6, index = idx_names, data = list(sine1000),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE))) |>
  dplyr::bind_rows(sine_8d_tbl("MIC")) |>
  dplyr::bind_rows(sine_8d_tbl("TIC")) |>
  dplyr::rowwise() |>
  dplyr::mutate(basis_df = list(sample_bases(
    idx = index, data = data, n_basis = 50, min_proj_dist = 1.5,
    step_size = 0.005, best = best, parallel = TRUE)
  ))
save(sq_sine_basis_df, file = here::here("data/sq_sine_basis_df.rda"))

# calculate squintability
res_holes <- sq_holes_basis_df |>
  rowwise() |>
  mutate(res = calc_squintability(
    basis_df, method = "nls", bin_width = 0.005, scale = TRUE,
    other_params = list(start = c(theta1 = 1, theta2 = 1, theta3 = 3, theta4 = 0)))) |>
  unnest(res) |>
  select(index, n, theta1: squint)


param_tbl <- tibble(other_params =  c(
  rep(list(start = list(theta1 = 1, theta2 = 1, theta3 = 2, theta4 = 0)), 5),
  list(start = list(theta1 = 1, theta2 = 0, theta3 = 100, theta4 = 0)),
  rep(list(start = list(theta1 = 1, theta2 = 1, theta3 = 2, theta4 = 0)), 2)))

res_sine <- sq_sine_basis_df|>
  bind_cols(other_params = param_tbl) |>
  bind_cols(scale = c(FALSE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)) |>
  mutate(res = calc_squintability(
    basis_df, method = "nls", bin_width = 0.005, scale = scale,
    other_params = list(start = other_params))) |>
  unnest(res) |>
  select(index, n, theta1:squint)


squintability <- res_holes |> bind_rows(res_sine) |> select(index, n, theta1:theta4, squint)
save(squintability, file = here::here("data", "squintability.rda"))
############################################################################
############################################################################
# TODO
sq_basis_dist_idx <- bind_rows(sq_holes_basis_df, sq_sine_basis_df) |>
  dplyr::select(basis_df, n, index) |>
  rename(index_name = index) |>
  unnest(basis_df) |>
  mutate(dist = ceiling(dist / 0.005) * 0.005,
         index_name = factor(index_name, levels = c("holes", "MIC", "TIC", "dcor2d_2",
                                          "loess2d", "splines2d", "stringy"))) |>
  group_by(n, index_name, dist) |>
  summarise(index = mean(index, na.rm = TRUE), .groups = "drop")


sq_basis_dist_idx <- sq_basis_dist_idx |>
  filter(index_name %in% c("holes", "TIC", "stringy")) |>
  group_by(n, index_name) |>
  mutate(index = (index - min(index)) / (max(index) - min(index))) |>
  ungroup() |>
  bind_rows(sq_basis_dist_idx |>
              filter(!index_name %in% c("holes", "TIC", "stringy")))
save(sq_basis_dist_idx, file = here::here("data", "sq_basis_dist_idx.rda"))

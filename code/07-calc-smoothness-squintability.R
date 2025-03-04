library(ferrn)
library(tourr)
library(tidyverse)

sine_8d_tbl <- function(index){
  tibble::tibble(
    n = 8, index = index, data = list(as.matrix(sine1000_8d)),
    best = list(matrix(c(rep(0, 12), 1, 0, 0, 1), nrow = 8, byrow = TRUE)))
}

smoothness_holes <- tibble::tibble(
  n = c(4, 6, 8, 10, 12),
  index = rep("holes", 5),
  data = list(pipe1000_4d, pipe1000_6d, pipe1000_8d, pipe1000_10d, pipe1000_12d),
  best = list(matrix(c(rep(0, 4), 1, 0, 0, 1), nrow = 4, byrow = TRUE),
              matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE),
              matrix(c(rep(0, 12), 1, 0, 0, 1), nrow = 8, byrow = TRUE),
              matrix(c(rep(0, 16), 1, 0, 0, 1), nrow = 10, byrow = TRUE),
              matrix(c(rep(0, 20), 1, 0, 0, 1), nrow = 12, byrow = TRUE))) |>
  dplyr::mutate(basis_df = purrr::pmap(
    list(data, best),
    function(data, best, n){sample_bases("holes", data = data, best = best, n_basis = 500)}))

holes_tidy <- smoothness_holes |>
  rowwise() |>
  dplyr::mutate(smooth = calc_smoothness(basis_df)) |>
  unnest(smooth)

##################################################
idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "splines2d", "skinny", "stringy2")
smoothness_sine <- tibble::tibble(
  n = 6, index = idx_names, data = list(sine1000),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE))) |>
  dplyr::bind_rows(sine_8d_tbl("MIC")) |>
  dplyr::bind_rows(sine_8d_tbl("TIC")) |>
  dplyr::bind_rows(sine_8d_tbl("loess2d")) |>
  dplyr::bind_rows(sine_8d_tbl("splines2d")) |>
  rowwise() |>
  dplyr::mutate(basis_df = list(sample_bases(index, n_basis = 500, best = best,
                                             data = as.matrix(data))))

idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "splines2d", "skinny", "stringy2")
smoothness_4d <- tibble::tibble(
  n = 4, index = idx_names, data = list(sine1000_4d),
  best = list(matrix(c(rep(0, 4), 1, 0, 0, 1), nrow = 4, byrow = TRUE))) |>
  rowwise() |>
  dplyr::mutate(basis_df = list(sample_bases(index, n_basis = 500, best = best,
                                             data = as.matrix(data))))

sine_tidy <- bind_rows(smoothness_sine, smoothness_4d) |>
  dplyr::mutate(smooth = calc_smoothness(basis_df)) |>
  unnest(smooth)

smoothness <- bind_rows(sine_tidy, holes_tidy) |>
  mutate(index = factor(index, levels = c("holes", "MIC", "TIC", "dcor2d_2",
                                          "loess2d", "splines2d", "skinny", "stringy2"))) |>
  arrange(index) |>
  select(index, n, variance:nugget)
save(smoothness, file = here::here("data", "smoothness.rda"))
################################################################################
################################################################################
# squintability - sample bases
# 2-3 mins
sq_holes_basis_df <- tibble::tibble(
  n = c(4, 6, 8, 10, 12),
  index = rep("holes", 5),
  data = list(pipe1000_4d, pipe1000_6d, pipe1000_8d, pipe1000_10d, pipe1000_12d),
  best = list(matrix(c(rep(0, 4), 1, 0, 0, 1), nrow = 4, byrow = TRUE),
              matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE),
              matrix(c(rep(0, 12), 1, 0, 0, 1), nrow = 8, byrow = TRUE),
              matrix(c(rep(0, 16), 1, 0, 0, 1), nrow = 10, byrow = TRUE),
              matrix(c(rep(0, 20), 1, 0, 0, 1), nrow = 12, byrow = TRUE))) |>
  dplyr::mutate(basis_df = purrr::pmap(
    list(data, best), function(data, best){
      sample_bases("holes", data = data, n_basis = 50, min_proj_dist = 1.5,
                   step_size = 0.005, best = best)}))
save(sq_holes_basis_df, file = here::here("data/sq_holes_basis_df.rda"))

idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "splines2d", "skinny", "stringy2")
# about 40 mins
sq_sine_68s_basis_df <- tibble::tibble(
  n = 6, index = idx_names, data = list(sine1000),
  best = list(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE))) |>
  dplyr::bind_rows(sine_8d_tbl("MIC")) |>
  dplyr::bind_rows(sine_8d_tbl("TIC")) |>
  dplyr::bind_rows(sine_8d_tbl("loess2d")) |>
  dplyr::bind_rows(sine_8d_tbl("splines2d")) |>
  dplyr::rowwise() |>
  dplyr::mutate(basis_df = list(sample_bases(
    idx = index, data = as.matrix(data), n_basis = 50, min_proj_dist = 1.5,
    step_size = 0.005, best = best, parallel = TRUE)
  ))
# for stringy2, run separately
# stringy$basis_df[[2]] <- stringy$basis_df[[2]] |> filter(!is.na(index))
sq_sine_basis_df2 <- sq_sine_basis_df |> filter(index != "stringy") |> bind_rows(stringy)

# about 5 mins
idx_names <- c("dcor2d_2", "loess2d", "MIC", "TIC", "splines2d", "skinny", "stringy2")
sq_4d_basis_df <- tibble::tibble(
  n = 4, index = idx_names, data = list(sine1000_4d),
  best = list(matrix(c(rep(0, 4), 1, 0, 0, 1), nrow = 4, byrow = TRUE))) |>
  dplyr::rowwise() |>
  dplyr::mutate(basis_df = list(sample_bases(
    idx = index, data = as.matrix(data), n_basis = 50, min_proj_dist = 1.5,
    step_size = 0.005, best = best, parallel = TRUE)
  ))

sq_sine_basis_df <- bind_rows(sq_sine_68s_basis_df, sq_4d_basis_df) |> arrange(index, n)
save(sq_sine_basis_df, file = here::here("data/sq_sine_basis_df.rda"))

# calculate squintability
res_holes <- sq_holes_basis_df |>
  rowwise() |>
  mutate(res = calc_squintability(
    basis_df, method = "nls", bin_width = 0.005, scale = TRUE,
    other_params = list(start = c(theta1 = 1, theta2 = 1, theta3 = 3, theta4 = 0)))) |>
  unnest(res) |>
  select(index, n, theta1:squint)

param_tbl <- tibble(other_params =  c(
  rep(list(start = list(theta1 = 1, theta2 = 1, theta3 = 2, theta4 = 0)), 11),
  rep(list(start = list(theta1 = 1, theta2 = 0, theta3 = 10, theta4 = 0.1)), 2),
  rep(list(start = list(theta1 = 1, theta2 = 1, theta3 = 2, theta4 = 0)), 3),
  rep(list(start = list(theta1 = 1, theta2 = 0, theta3 = 10, theta4 = 0.15)), 2)
  ))

# TIC, stringy2 and skinny need scale
res_sine <- sq_sine_basis_df |>
  bind_cols(other_params = param_tbl) |>
  bind_cols(scale = c(rep(FALSE, 3), rep(TRUE, 3), rep(FALSE, 5), rep(FALSE, 2), rep(FALSE, 3), rep(TRUE, 2))) |>
  mutate(res = calc_squintability(
    basis_df, method = "nls", bin_width = 0.005, scale = scale,
    other_params = list(start = other_params))) |>
  unnest(res) |>
  select(index, n, theta1:squint)

squintability <- bind_rows(res_holes, res_sine) |> select(index, n, theta1:theta4, squint)
save(squintability, file = here::here("data", "squintability.rda"))
############################################################################
############################################################################
# check smoothness against squintability
smoothness |> select(n, index, smoothness) |>
  left_join(squintability |> select(n, index, squint) |> rename(squintability = squint)) |>
  ggplot(aes(x = smoothness, y = squintability, group = index, color = as.factor(index))) +
  geom_point() +
  ggrepel::geom_label_repel(aes(label = index), nudge_x = 0.01, nudge_y = 0.01) +
  theme(aspect.ratio = 1)

############################################################################
############################################################################
# TODO
sq_basis_df_6d_TIC_spline_skinny <- sq_sine_basis_df |>
  filter(index %in% c("TIC", "splines2d", "skinny"),n  == 6) |>
  dplyr::select(basis_df, n, index) |>
  rename(index_name = index) |>
  unnest(basis_df) |>
  filter(id <= 5) |>
  mutate(index_name = ifelse(index_name == "splines2d", "splines", index_name),
         index_name = factor(index_name, levels = c("TIC", "splines",  "skinny"))) |>
  group_by(id, index_name) |>
  mutate(index = (index - min(index)) / (max(index) - min(index))) |>
  ungroup()

save(sq_basis_df_6d_TIC_spline_skinny, file = here::here("data", "sq_basis_df_6d_TIC_spline_skinny.rda"))



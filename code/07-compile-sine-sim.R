library(tidyverse)
library(minerva)
library(progress)
library(furrr)
library(cassowaryr)
load(here::here("data/sim_sine_6d_dcor2d.rda"))
load(here::here("data/sim_sine_6d.rda")) # loess
load(here::here("data/sim_sine_68d_TICMIC.rda"))

################################################################
################################################################
set.seed(123456)
sine1000 <- spinebil::sinData(6, 1000) %>% scale()
colnames(sine1000) <- paste0("V", 1:6)

sim_res_all <- bind_rows(sim_sine_6d_dcor2d,
                         sim_sine_6d |> mutate(idx_f = "loess2d"),
                         sim_sine_68d_TICMIC |> filter(d == 6)) |>
  filter(n_jellies == 50, max_tries == 50)


sim_res_wide <- sim_res_all |>
  select(idx_f:index_val, tries, loop) |>
  mutate(row = row_number()) |>
  pivot_wider(names_from = idx_f, values_from = index_val) |>
  group_split(row)

add_idx_val <- function(data, idx_f){
  pb$tick()
  idx_f_fn <- eval(sym(idx_f))
  data |>
    mutate(!!sym(idx_f) := ifelse(is.na(!!sym(idx_f)),
                                  idx_f_fn()(sine1000 %*% basis[[1]]),
                                  !!sym(idx_f)))
}

# Declaimer: in this script, I save the intermediate objects after each > 30 mins
# run. But to save storage space, only the final data `sim_sine` is saved in the
# data folder

################################################################
################################################################
# Part 1
# for the dcor2d, loess, TIC, MIC data, compute all the four indexes + stringy

# add dcor: over an hour
pb <- progress_bar$new(total = length(sim_res_wide))
impute_dcor2d <- sim_res_wide |> map_dfr(~add_idx_val(.x, "dcor2d_2")) |> rename(loess2d = loess)

# add loess in three batches: over an hour
pb <- progress_bar$new(total = nrow(impute_dcor2d))
impute_loess1 <- impute_dcor2d |> filter(row < 125001) |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))
pb <- progress_bar$new(total = 125000)
impute_loess2 <- impute_dcor2d |> filter(between(row, 250001, 375000)) |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))
pb <- progress_bar$new(total = 125000)
impute_loess3 <- impute_dcor2d |> filter(row > 375000) |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))
exist_loess <- impute_dcor2d |> filter(!is.na(loess2d)) |> pull(loess2d)
impute_loess2d_vec <- c(impute_loess1$loess2d, exist_loess, impute_loess2$loess2d, impute_loess3$loess2d)

# add MIC: ~40 mins
pb <- progress_bar$new(total = nrow(impute_dcor2d))
impute_MIC <- impute_dcor2d |> group_split(row) |> map_dfr(~add_idx_val(.x, "MIC"))
impute_MIC_vec <- impute_MIC$MIC

# add TIC: ~40 mins
pb <- progress_bar$new(total = nrow(impute_dcor2d))
impute_TIC <- impute_dcor2d |> group_split(row) |> map_dfr(~add_idx_val(.x, "TIC"))
impute_TIC_vec <- impute_TIC$TIC

sim_full <- impute_dcor2d |>
  mutate(loess2d = impute_loess2d_vec, MIC = impute_MIC_vec, TIC = impute_TIC_vec)

# add stringy in parallel: about 2 hours
impute_stringy <- sim_full |>
  pull(basis) |>
  future_map(function(x){stringy()(as.matrix(sine1000) %*% x)})
sim_full2 <- sim_full |> mutate(stringy = unlist(impute_stringy))

# add spline in parallel
impute_spline2d <- sim_full2 |>
  pull(basis) |>
  future_map(function(x){tourr::splines2d()(as.matrix(sine1000) %*% x)})
sim_sine_less_spline <- sim_full2 |> mutate(splines2d = unlist(impute_spline2d))

################################################################
################################################################
# Part 2:
# for the spline data, compute all the five index values: dcor2d, loess, MIC,
# TIC, and stringy
load(here::here("data/sim_sine_6d_spline.rda"))
sim_spline <- sim_sine_6d_spline |>
  mutate(row = row_number() + 500000) |>
  rename(splines2d = index_val) |>
  select(d, n_jellies, max_tries, sim, seed, id, basis, tries, loop, row, splines2d) |>
  mutate(dcor2d = NA, loess2d = NA, MIC = NA, TIC = NA, stringy = NA)

# add MIC
pb <- progress_bar$new(total = 125000)
sim_spline1 <- sim_spline |> group_split(row) |> map_dfr(~add_idx_val(.x, "MIC"))
# add TIC
pb <- progress_bar$new(total = 125000)
sim_spline2 <- sim_spline1 |> group_split(row) |> map_dfr(~add_idx_val(.x, "TIC"))
# add loess2d
pb <- progress_bar$new(total = 125000)
sim_spline3 <- sim_spline2 |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))
# add dcor2d
impute_dcor2d <- sim_spline3 |>
  pull(basis) |>
  future_map(function(x){dcor2d_2()(as.matrix(sine1000) %*% x)})
sim_spline3_dcor2d <- sim_spline3 |> mutate(dcor2d = unlist(impute_dcor2d))
# add stringy
impute_stringy <- sim_spline3_dcor2d |>
  pull(basis) |>
  future_map(function(x){stringy()(as.matrix(sine1000) %*% x)})
sim_sine_spline <- sim_spline3_dcor2d |> mutate(stringy = unlist(impute_stringy))

################################################################
################################################################
# combine everything together: sim_sine

best <-matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE)
sim_sine <- sim_sine_less_spline |>
  bind_rows(sim_sine_spline) |>
  rowwise( ) |>
  mutate(proj_dist = tourr::proj_dist(best, basis)) |>
  ungroup() |>
  mutate(orig_idx = rep(c("dcor2d", "loess2d", "MIC", "TIC", "splines2d"), 125000))
save(sim_sine, file = here::here("data/sim_sine.rda"))


################################################################
stringy <- function(){
  function(mat){
    cassowaryr::sc_stringy(mat[,1], mat[,2])
  }
}

MIC <- function(){
  function(mat){
    mine(mat[,1], mat[,2], alpha = 0.3,  est = "mic_e")$MIC
  }
}

TIC <- function(){
  function(mat){
    mine(mat[,1], mat[,2], est = "mic_e", alpha = 0.3)$TIC
  }
}

dcor2d_2 <- function() {
  function(mat) {
    xy <- na.omit(data.frame(x = mat[, 1], y = mat[, 2]))
    measure <- with(xy, energy::dcor2d(x, y, type = "U"))
    return(measure)
  }
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

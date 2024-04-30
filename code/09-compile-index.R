library(tidyverse)
library(minerva)
library(progress)
load(here::here("data/sim_sine_6d_dcor2d.rda"))
load(here::here("data/sim_sine_6d.rda"))
load(here::here("data/sim_sine_68d_TICMIC.rda"))

set.seed(123456)
sine1000 <- spinebil::sinData(6, 1000) %>% scale()
colnames(sine1000) <- paste0("V", 1:6)

sim_res_all <- bind_rows(
  sim_sine_6d_dcor2d,
  sim_sine_6d |> mutate(idx_f = "loess2d"),
  sim_sine_68d_TICMIC |> filter(d == 6)
) |>
  filter(n_jellies == 50, max_tries == 50)


sim_res_wide <- sim_res_all |>
  select(idx_f:index_val, tries, loop) |>
  mutate(row = row_number()) |>
  pivot_wider(names_from = idx_f, values_from = index_val)

sim_res_wide2 <- sim_res_wide |> group_split(row)

add_idx_val <- function(data, idx_f){
  pb$tick()
  idx_f_fn <- eval(sym(idx_f))
  data |>
    mutate(!!sym(idx_f) := ifelse(is.na(!!sym(idx_f)),
                                  idx_f_fn()(sine1000 %*% basis[[1]]),
                                  !!sym(idx_f)))
}

# over an hour
pb <- progress_bar$new(total = length(sim_res_wide2))
impute_dcor2d <- sim_res_wide2 |> map_dfr(~add_idx_val(.x, "dcor2d_2")) |> rename(loess2d = loess)

# over an hour
pb <- progress_bar$new(total = nrow(impute_dcor2d))
impute_loess1 <- impute_dcor2d |> filter(row < 125001) |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))
pb <- progress_bar$new(total = 125000)
impute_loess2 <- impute_dcor2d |> filter(between(row, 250001, 375000)) |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))
pb <- progress_bar$new(total = 125000)
impute_loess3 <- impute_dcor2d |> filter(row > 375000) |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))
exist_loess <- impute_dcor2d |> filter(!is.na(loess2d)) |> pull(loess2d)
impute_loess2d_vec <- c(impute_loess1$loess2d, exist_loess, impute_loess2$loess2d, impute_loess3$loess2d)

# ~40 mins
pb <- progress_bar$new(total = nrow(impute_dcor2d))
impute_MIC <- impute_dcor2d |> group_split(row) |> map_dfr(~add_idx_val(.x, "MIC"))
impute_MIC_vec <- impute_MIC$MIC

# ~40 mins
pb <- progress_bar$new(total = nrow(impute_dcor2d))
impute_TIC <- impute_dcor2d |> group_split(row) |> map_dfr(~add_idx_val(.x, "TIC"))
impute_TIC_vec <- impute_TIC$TIC

sim_full <- impute_dcor2d |>
  mutate(loess2d = impute_loess2d_vec,
         MIC = impute_MIC_vec,
         TIC = impute_TIC_vec)
save(sim_full, file = here::here("data/sim_full.rda"))

library(furrr)
# about 2 hours
impute_stringy <- sim_full |>
  pull(basis) |>
  future_map(function(x){stringy()(as.matrix(sine1000) %*% x)})
impute_stringy_vec <- unlist(impute_stringy)
sim_full2 <- sim_full |> mutate(stringy = impute_stringy_vec)
save(sim_full2, file = here::here("data/sim_full2.rda"))


# t1 <- Sys.time()
# map(impute_dcor2d$basis[1:10], ~stringy()(sine1000 %*% .x))
# t2 <- Sys.time()
# t2 - t1

load(here::here("data/sim_sine_6d_spline.rda"))
sim_spline <- sim_sine_6d_spline |>
  mutate(row = row_number() + 500000) |>
  rename(splines2d = index_val) |>
  select(d, n_jellies, max_tries, sim, seed, id, basis, tries, loop, row, splines2d) |>
  mutate(dcor2d = NA, loess2d = NA, MIC = NA, TIC = NA, stringy = NA)

set.seed(123456)
sine1000 <- spinebil::sinData(6, 1000) %>% scale()
colnames(sine1000) <- paste0("V", 1:6)

pb <- progress_bar$new(total = 125000)
sim_spline1 <- sim_spline |> group_split(row) |> map_dfr(~add_idx_val(.x, "MIC"))
pb <- progress_bar$new(total = 125000)
sim_spline2 <- sim_spline1 |> group_split(row) |> map_dfr(~add_idx_val(.x, "TIC"))
pb <- progress_bar$new(total = 125000)
sim_spline3 <- sim_spline2 |> group_split(row) |> map_dfr(~add_idx_val(.x, "loess2d"))

save(sim_spline3, file = here::here("data/sim_spline3.rda"))

########### not yet
# need parallel
t1 <- Sys.time()
pb <- progress_bar$new(total = 125000)
sim_spline4 <- sim_spline3 |> group_split(row) |> map_dfr(~add_idx_val(.x, "stringy"))
t2 <- Sys.time()
t2 - t1

pb <- progress_bar$new(total = 125000)
sim_spline5 <- sim_spline4|> group_split(row) |> map_dfr(~add_idx_val(.x, "dcor2d_2"))



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

library(tidyverse)
library(ferrn)
library(tourr)


################################################################################
# d = 6
# d = 4 and 8 data is simulated in 10-sim-4d.R
sim <- function(d = d, n_jellies = n_jellies, max.tries = max.tries,
                optim_seed = seed, sim = sim){
  sine1000 <- ferrn::sine1000_6d

  cat("sim: ", sim, "\n")
  cat("n_jellies: ", n_jellies, "\n")
  cat("max.tries: ", max.tries, "\n")
  cat("d: ", d, "\n")
  res <- map_dfr(optim_seed, function(seed){
    t1 <- Sys.time()
    set.seed(seed)
    res <- animate_xy(sine1000, guided_tour(
      tourr::splines2d(), d = 2, n_jellies = n_jellies,
      search_f =  search_jellyfish,
      max.tries = max.tries
    ), max_frames = max.tries) |> dplyr::select(-id)
    t2 <- Sys.time()
    res <- res |> mutate(time = t2 - t1)
    return(res)
  })

  return(res)
}

# separate into sim 1-30 and 31-50 to run, each takes < 2 hrs
# 50/ 50
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(n_jellies = 50,
                      max_tries = 100,
                      d = 6) |>
  crossing(sim = 1:50) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
set.seed(123)
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1

sim_res |> get_best(sim) |> plot_projection(sine1000)
sim_sine_6d_splines2d_100 <- sim_res
save(sim_sine_6d_splines2d_100, file = here::here("data/sim_sine_6d_splines2d_100.rda"))

sim_sine_6d_spline130 <- sim_res |> select(-alpha)
sim_sine_6d_spline3150 <- sim_res |> select(-alpha)

sim_sine_6d_spline <- bind_rows(
  sim_sine_6d_spline130, sim_sine_6d_spline3150 |> mutate(id = id + 30)
  )
sim_sine_6d_splines2d <- tibble(index = "splines2d") |>
  bind_cols(sim_sine_6d_spline) |>
  select(index, d, n_jellies, max_tries, sim:time)
save(sim_sine_6d_splines2d, file = here::here("data-raw/sim_sine_6d_splines2d.rda"))

################################################################################
################################################################################
# additional smaller data
sim_sine_6d_spline_head <- sim_sine_6d_splines2d |>
  filter(n_jellies == 50, max_tries == 50) |>
  head(10)
save(sim_sine_6d_spline_head, file = here::here("data/sim_sine_6d_spline_head.rda"))
sim_sine_6d_spline_best <- sim_sine_6d_splines2d |>
  filter(n_jellies == 50, max_tries == 50) |>
  get_best(group = sim) |> ungroup()
save(sim_sine_6d_spline_best, file = here::here("data/sim_sine_6d_spline_best.rda"))

mmat <- matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE)
sim_sine_6d_spline_projdist <- sim_sine_6d_splines2d |>
  filter(n_jellies == 50, max_tries == 50) |>
  rowwise() |>
  mutate(proj_dist = tourr::proj_dist(basis, mmat)) |>
  select(-d, -id, -info, -method, -time, -basis) |>
  ungroup()
save(sim_sine_6d_spline_projdist, file = here::here("data/sim_sine_6d_spline_projdist.rda"))


################################################################################
################################################################################
set.seed(123456)
sine1000 <- ferrn::sine1000_6d
sum <- sim_sine_6d_spline |>
  group_by(id) |> filter(index_val == max(index_val)) |> filter(row_number() == 1)
b <- sum |> pull(basis)
dt <- map_dfr(b, ~as_tibble(as.matrix(sine1000) %*% .x), .id = "id")
dt |>
  ggplot(aes(x = V2 ,y = V1)) +
  geom_point(size = 0.3) +
  theme(aspect.ratio = 1) +
  facet_wrap(vars(id))

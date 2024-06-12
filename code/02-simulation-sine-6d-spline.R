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
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(n_jellies = 50,
                      max_tries = 50,
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

sim_sine_6d_spline130 <- sim_res |> select(-alpha)
sim_sine_6d_spline3150 <- sim_res |> select(-alpha)

sim_sine_6d_spline <- bind_rows(sim_sine_6d_spline130, sim_sine_6d_spline3150 |> mutate(id = id + 30))
save(sim_sine_6d_spline, file = here::here("data/sim_sine_6d_spline.rda"))


sim_sine_6d_spline_head <- sim_sine_6d_spline |>
  filter(n_jellies == 50, max_tries == 50) |>
  head(10)
save(sim_sine_6d_spline_head, file = here::here("data/sim_sine_6d_spline_head.rda"))
sim_sine_6d_spline_best <- sim_sine_6d_spline |>
  filter(n_jellies == 50, max_tries == 50) |>
  get_best(group = sim) |> ungroup()
save(sim_sine_6d_spline_best, file = here::here("data/sim_sine_6d_spline_best.rda"))

sim_sine_6d_spline_projdist <- sim_sine_6d_spline |>
  filter(n_jellies == 50, max_tries == 50) |>
  rowwise() |>
  mutate(proj_dist = tourr::proj_dist(basis, matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE))) |>
  select(-d, -id, -info, -method, -time, -basis) |>
  ungroup()
save(sim_sine_6d_spline_projdist, file = here::here("data/sim_sine_6d_spline_projdist.rda"))

set.seed(123456)
sine1000 <- spinebil::sinData(6, 1000) %>% scale() %>% as_tibble()
colnames(sine1000) <- paste0("V", 1:6)
sum <- sim_sine_6d_spline |> group_by(id) |> filter(index_val == max(index_val)) |> filter(row_number() == 1)
b <- sum |> pull(basis)
dt <- map_dfr(b, ~as_tibble(as.matrix(sine1000) %*% .x), .id = "id")
dt |>
  ggplot(aes(x = V2 ,y = V1)) +
  geom_point(size = 0.3) +
  theme(aspect.ratio = 1) +
  facet_wrap(vars(id))

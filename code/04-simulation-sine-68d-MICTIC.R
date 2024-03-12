library(tidyverse)
library(spinebil)
library(tourr)
library(minerva)

################################################################################
# d = 6
sim <- function(d = d, idx_f = idx_f, n_jellies = n_jellies, max.tries = max.tries,
                data_seed = 123456, optim_seed = seed, sim = sim){
  set.seed(data_seed)
  sine1000 <- sinData(d, 1000) %>% scale() %>% as_tibble()
  colnames(sine1000) <- paste0("V", 1:d)

  cat("sim: ", sim, "\n")
  cat("n_jellies: ", n_jellies, "\n")
  cat("max.tries: ", max.tries, "\n")
  cat("idx_f: ", idx_f, "\n")
  cat("d: ", d, "\n")

  idx_f <- eval(sym(idx_f))
  res <- map_dfr(optim_seed, function(seed){
    t1 <- Sys.time()
    set.seed(seed)
    res <- animate_xy(sine1000, guided_tour(
      idx_f(), d = 2, n_jellies = n_jellies,
      search_f =  search_jellyfish,
      max.tries = max.tries
    ), max_frames = max.tries) |> dplyr::select(-id)
    t2 <- Sys.time()
    res <- res |> mutate(time = t2 - t1)
    return(res)
  })

  return(res)
}

MIC <- function(){
  function(mat){
    minerva::mine(mat[,1], mat[,2], alpha = 0.3,  est = "mic_e")$MIC
  }
}

TIC <- function(){
  function(mat){
    minerva::mine(mat[,1], mat[,2], est = "mic_e", alpha = 0.3)$TIC
  }
}

#  n_jellies | max_tries
#         20 |       50
#         20 |      100
#         50 |       50
#         50 |      100
#        100 |       50
#        100 |      100
# each 50 repetitions
# take 2 hrs
set.seed(123)
seed <- sample(1000: 10000, size = 50)
sim_setup <- crossing(idx_f = c("MIC", "TIC"),
                      d = c(6, 8),
                      n_jellies = c(20, 50, 100),
                      max_tries = c(50, 100)) |>
  filter(!(n_jellies == 100 & max_tries == 100)) |>
  crossing(sim = 1) |>
  mutate(seed = seed[sim], id = row_number())

t1 <- Sys.time()
sim_res <- sim_setup |>
  rowwise() |>
  mutate(res = list(sim(d, idx_f, n_jellies, max_tries, optim_seed = seed, sim = sim))) |>
  ungroup() |>
  unnest(res)
t2 <- Sys.time()
t2 - t1

# both indexes are similar
# I_max_max shows that most combination (max_tries x n_jellies x d x idx_f)


# sim_sine_68d_TICMIC <- sim_res |> select(-alpha)
# save(sim_sine_68d_TICMIC, file = "data/sim_sine_68d_TICMIC.rda")
#
#
set.seed(123456)
sine1000 <- spinebil::sinData(6, 1000) %>% scale() %>% as_tibble()
colnames(sine1000) <- paste0("V", 1:6)
sum <- sim_res |> group_by(id) |> filter(index_val == max(index_val)) |> filter(row_number() == 1)
#tour_level |> filter(id == 167) |> pull(basis)
# tour_level |>
#   filter(n_jellies == 50, max_tries == 100, idx_f == "MIC", d == 6) |>
#   filter(id == 164) |> pull(basis) -> b
b <- setup_level_best_basis |> filter(id == 717) |> pull(basis)
ids <- setup_level_best_basis |> filter(d == 6) |> pull(id)
res <- map_dfr(ids, ~{
  b <- setup_level_best_basis |> filter(id == .x) |> pull(basis)
  dt <- as_tibble(as.matrix(sine1000) %*% b[[1]])
  return(dt)
}, .id = "id")

labels_df <- setup_level_best_basis |>
  ungroup() |>
  filter(d == 6) |>
  select(idx_f:max_tries, index_val) |>
  mutate(id = as.character(1:10))

library(ggh4x)
res |>
  left_join(labels_df) |>
  ggplot(aes(x = V2 ,y = V1)) +
  geom_point() +
  geom_text(data = labels_df,
            aes(x = -3, y = 3, label = round(index_val, 2)),
            nudge_x = 0.1, nudge_y = 0.1) +
  facet_nested(idx_f ~ n_jellies + max_tries, labeller = label_both )

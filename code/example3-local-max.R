library(tidyverse)
library(spinebil)
library(tourr)

set.seed(123456)
pipe1000 <- pipeData(6, 1000) |> scale() |> as_tibble()
sine1000 <- sinData(6, 1000) |> scale()
dt <-pipe1000 |> bind_cols(sine1000[,5:6])
colnames(dt) <- paste0("V", 1:(6+ 2))


better_res <- animate_xy(dt, guided_tour(
  holes(), d = 2,
  search_f =  search_better, alpha = 0.7,
  max.tries = 500
))

set.seed(1234)
jellyfish_res <- animate_xy(dt, guided_tour(
  holes(), d = 2, n_jellies = 50,
  search_f =  search_jellyfish,
  max.tries = 50
))

proj1 <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 1, 0, 0), ncol = 2)
proj2 <- matrix(c(0, 0, 0, 0, 0, 0, 1, 0,
                  0, 0, 0, 0, 0, 0, 0, 1), ncol = 2)
proj3 <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 0), ncol = 2)
proj4 <- matrix(c(0, 0, 0, 0, 1, 0, 0, 0,
                  0, 0, 0, 0, 0, 1, 0, 1), ncol = 2)
proj5 <- matrix(c(0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 1, 0), ncol = 2)
proj6 <- matrix(c(0, 0, 0, 0, 0, 1, 0, 0,
                  0, 0, 0, 0, 0, 0, 0, 1), ncol = 2)
proj_list <- list(proj1, proj2, proj3, proj4, proj5, proj6)
res <- map_dfr(proj_list, ~ as_tibble(as.matrix(dt) %*% .), .id = "id")
idx <- map(proj_list, ~holes()(as.matrix(dt) %*% .x))

res |>
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  theme(aspect.ratio = 1) +
  facet_wrap(vars(id))


best_jellyfish <- jellyfish_res #|> filter(info == "current_best")
basis_mtx <- ferrn::get_basis_matrix(best_jellyfish) |>
  bind_cols(V17= best_jellyfish$index_val)
km_res <- kmeans(basis_mtx, centers = 4, nstart = 10)
plot(basis_mtx, col = km_res$cluster)
basis_mtx |> mutate(cluster = as.factor(km_res$cluster)) |>
  GGally::ggpairs(columns = 1:16, aes(colour = cluster))

best_b <- jellyfish_res |> filter(index_val == max(index_val)) |> pull(basis)
res <- as_tibble(as.matrix(dt) %*% best_b[[1]])
res |>
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  theme(aspect.ratio = 1)


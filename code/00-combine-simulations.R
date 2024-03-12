#####################################################################
# simulation sine 68d TICMIC
load(file = "data/sim_sine_68d_TICMIC.rda")
sim_data <- sim_sine_68d_TICMIC
setup_level_best_basis <- sim_data |>
  group_by(n_jellies, max_tries, d, idx_f) |>
  filter(index_val == max(index_val)) |>
  filter(row_number() == 1)

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

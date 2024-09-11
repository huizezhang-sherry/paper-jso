library(palmerpenguins)
library(tidyverse)
library(tourr)
library(ferrn)
library(GpGp)
library(fields)
################################################################################
# Using the JSO in a PPGT
# clean up the penguins data
penguins_cl <- palmerpenguins::penguins |>
  filter(!is.na(bill_length_mm)) |>
  rename(bl = bill_length_mm, bd = bill_depth_mm,
         fl = flipper_length_mm, bm = body_mass_g)
penguins_df <- penguins_cl |> select(bl:bm) |> rescale()
# run the projection pursuit with JSO
res <- animate_xy(penguins_df,
                  guided_tour(lda_pp(cl = penguins_cl$species),
                              search_f = search_jellyfish,
                              n_jellies = 20, max.tries = 50))

# extract the bases from the first jellyfish
bases <- res |> filter(loop == 1) |> pull(basis) |> check_dup(0.1)
# visualise its path with a planned tour
animate_xy(data = penguins_df, tour_path = planned_tour(bases),
           col = penguins_cl$species)

################################################################################
# calculate smoothness and squintability
# smoothness
basis_smoothness <- sample_bases(
  idx = "holes", data = sine1000, n_basis = 300,
  best = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1), nrow = 6))
basis_smoothness |> head(3)
calc_smoothness(basis_smoothness)

# squintability
basis_squint <- sample_bases(
  idx = "holes", data = sine1000, n_basis = 100,
  best = matrix(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1), nrow = 6),
  step_size = 0.1, min_proj_dist = 1.5)
calc_squintability(basis_squint, method = "nls", bin_width = 0.01)

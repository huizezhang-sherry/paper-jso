## ----setup, echo = FALSE, message=FALSE, warning=FALSE-----------------------------
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggh4x)
library(broom)
library(kableExtra)
library(tourr)
library(PPtreeViz)
library(ferrn)
library(colorspace)
load(here::here("data/sim_pipe_run_best.rda"))
load(here::here("data/sim_sine_6d_spline_head.rda"))
load(here::here("data/sim_sine_6d_spline_best.rda"))
load(here::here("data/sim_sine_6d_spline_projdist.rda"))
load(here::here("data/pipe_better.rda"))
load(here::here("data/pipe_jellyfish.rda"))
load(here::here("data/smoothness.rda"))
load(here::here("data/squintability.rda"))
load(here::here("data/sq_basis_dist_idx.rda"))
load(here::here("data/sim_df.rda"))


## ----------------------------------------------------------------------------------
#| label: fig-example-functions
#| echo: false
#| fig-width: 9
#| fig-height: 6
#| out-width: 80%
#| fig-cap: "Examples of PP indexes with large (top row) and small (bottom row) squint angles, shown with a Huber plot, and histogram of the projected data corresponding to the optimal projection. A Huber plot shows the PP index values for all 1D data projections in polar coordinates."
#proj <- animate_xy(flea[,1:6], guided_tour(lda_pp(flea$species)))
best_proj <- matrix(c(-0.44157005,  0.80959831,
                      0.14622089, -0.08228193,
                      -0.08454314, -0.23234287,
                      0.70367852,  0.28185021,
                      0.27876825,  0.44819691,
                      0.45123453,  0.05896646), ncol=2, byrow=T)

flea2D <- as.matrix(flea[,1:6]) %*% best_proj
colnames(flea2D) <- c("X1", "X2")
flea2D <- as.data.frame(flea2D)
flea_huber <- prep_huber(flea2D, index = skewness())
sqa1 <- ggplot() +
  geom_huber(data = flea_huber$idx_df, aes(x = x, y = y)) +
  geom_point(data = flea2D, aes(x = X1, y = X2)) +
  geom_abline(slope = flea_huber$slope, intercept = 0) + 
  coord_fixed() +
  theme_huber() +
  ggtitle("(a) Skewness index")

sqa2 <- ggplot(flea_huber$proj_df, aes(x=x)) + 
  geom_histogram(binwidth=0.25, colour="white") +
  xlab("") + ylab("") +
  ggtitle("(b) Optimally projected data") +
  theme_bw() +
  theme(axis.text.y = element_blank())

data(randu)
randu_std <- as.data.frame(apply(randu, 2, function(x) (x-mean(x))/sd(x)))
randu_std$yz <- sqrt(35)/6*randu_std$y-randu_std$z/6
randu_df <- randu_std[c(1,4)]
randu_huber <- prep_huber(randu_df, index = norm_bin(nr = nrow(randu_df)))

sqa3 <- ggplot() +
  geom_huber(data = randu_huber$idx_df, aes(x = x, y = y)) +
  geom_point(data = randu_df, aes(x = x, y = yz)) +
  geom_abline(slope = randu_huber$slope, intercept = 0) +
  theme_huber() +
  coord_fixed() + 
  ggtitle("(c) Binned normality index")

sqa4 <- ggplot(randu_huber$proj_df, aes(x = x)) +
  geom_histogram(breaks = seq(-2.2, 2.4, 0.12)) +
  xlab("") + ylab("") +
  theme_bw() +
  theme(axis.text.y = element_blank()) + 
  ggtitle("(d) Optimally projected data")

sqa1 + sqa2 + sqa3 + sqa4 + plot_layout(ncol=2, widths=c(2,2))



## ----------------------------------------------------------------------------------
#| label: fig-matern-1d
#| out.width: 100%
#| fig.cap: |
#|   Five random simulations from a Gaussian Process defined on $\mathbb{R}$ with zero mean and Matérn-$\nu$ covariance function, with $\nu=1$ (left), $\nu=2$ (middle), and $\nu=4$ (right), showing that higher values of $\nu$ produce smoother curves.
knitr::include_graphics(here::here("figures", "matern_simulation_1d.png"))


## ----------------------------------------------------------------------------------
#| label: fig-matern-2d
#| out.width: 100%
#| fig.cap: |
#|   One random simulation from a Gaussian Process defined on $\mathbb{R}^2$ with zero mean and Matérn-$\nu$ covariance function, with $\nu=1$ (left), $\nu=2$ (middle), and $\nu=4$ (right), showing that higher values of $\nu$ produce smoother surfaces.
knitr::include_graphics(here::here("figures", "matern_simulation_2d.png"))


## ----------------------------------------------------------------------------------
#| label: fig-success-rate
#| fig.cap: "Illustration of success rate calculation: Final projections based on projection pursuit to find the pipe shape in 8D data using the holes index, optimised by CRS, in 50 simulations. The 50 final projections are sorted by their index values. The highest index value found across all simulations is 0.969. Out of the 50 simulations, 43 achieved an index value within 0.05 of the best, resulting in a success rate of 0.86 (43/50)."
# pipe_jellyfish_12d <- pipe_better |> 
#   filter(dim == 2) |> 
#   mutate(dim = as.numeric(dim) * 8, sim = as.numeric(sim)) 
# proj_dt <- map_dfr(pipe_jellyfish_12d |> pull(basis),
#               ~as_tibble(pipe1000_8d %*% .x), .id = "sim") |>
#   mutate(sim = as.numeric(sim)) |>
#   left_join(pipe_jellyfish_12d |> select(sim, index_val), by = "sim")
# 
# idx_val_sorted <- proj_dt |> pull(index_val) |> unique() |> sort()
# proj_dt |>
#   ggplot(aes(x = V1, y = V2)) +
#   geom_point(size = 0.5) +
#   geom_text(aes(label = round(index_val, 3)), x = 0, y = 3.7,
#             size = 8) +
#   xlim(-4, 4) + ylim(-4, 4) +
#   facet_wrap(~fct_reorder(as.factor(sim), -index_val), ncol = 10) +
#   theme_bw() +
#   theme(aspect.ratio = 1, axis.ticks = element_blank(),
#         axis.text = element_blank(), axis.title = element_blank(),
#         panel.grid = element_blank(), strip.background = element_blank(),
#         strip.text = element_blank()) 
knitr::include_graphics(here::here("figures", "success-rate.png"))


## ----------------------------------------------------------------------------------
#| label: fig-smoothness
#| out.width: 100%
#| fig.cap: "Illustration of steps to calculate smoothness. For a given projection pursuit problem defined by the shape to find, data dimension and the index function, 1) sample random bases given the orthonormality contraint, 2) calculate the projection distance and the index value for each random basis, and 3) fit a Gaussian process model of index values against projection distances to obtain the smoothness measure. "
# basis_smoothness <- sample_bases(idx = "splines2d", n_basis = 1000)
# basis_smoothness |>
#    ggplot(aes(x = dist, y = index)) +
#    geom_point() + 
#    labs(x = "Projection distance", y = "Index value") + 
#    theme_bw() +   
#   theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30))
knitr::include_graphics(here::here("figures", "smoothness.png"))


## ----------------------------------------------------------------------------------
#| label: fig-squintability
#| out.width: 100%
#| fig.cap: |
#|   Illustration of steps to calculate squintability. For a given projection pursuit problem defined by the shape to find, data dimension and the index function, 1) sample random bases given the orthonormality and projection distance contraint, 2) interpolate the sampled bases to the optimal basis and calculate the projection distance and the index value for each interpolated basis. 3) bin the index values by projection distances to obtain the average index value for each bin, 4) fit the scaled sigmoid function in equation \eqref{eq-squintability} to the binned index values against projection distances using non-linear least square, 5) calculate the squintability measure using equation (@eq-squintability-parametric) with parameters estimated from the model.
# p1 <- sim_sine_6d_spline_projdist |> 
#   ggplot(aes(x = proj_dist, y = index_val)) +
#   geom_point(size = 0.1, alpha = 0.5) + 
#   labs(x = "Projection distance", y = "Index value") + 
#   theme_bw() + 
#   theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30))
# 
# 
# dist_bin <- ceiling(sim_sine_6d_spline_projdist[["proj_dist"]] / 0.1) * 0.1
# dt <- sim_sine_6d_spline_projdist |>
#   dplyr::bind_cols(dist_bin = dist_bin) |>
#   dplyr::group_by(dist_bin) |>
#   dplyr::summarise(index_val = mean(index_val, na.rm = TRUE))
# 
# p2 <-  dt |> 
#   ggplot(aes(x = dist_bin, y = index_val)) + 
#   geom_line() +
#   geom_point() + 
#   labs(x = "Projection distance", y = "Index value") + 
#   theme_bw() + 
#   theme(axis.text = element_text(size = 20), axis.title = element_text(size = 30))
knitr::include_graphics(here::here("figures", "squintability.png"))


## ----------------------------------------------------------------------------------
#| fig.height = 7.5,
#| fig.width = 10,
#| label = fig-proj,
#| fig.align = "center",
#| fig.cap = "Projections found by the JSO and CRS at each 10th quantile across 50 simulations. The projection pursuit problem is to find the pipe shape using the holes index in the 6, 8, 10, and 12-dimensional spaces. The JSO uses 100 jellyfishes and a maximum number of tries of 100. The CRS uses a maximum of 1000 tries in each step of random sampling step before the algorithm terminates. In the 6-D data space, JSO always finds a clear pipe shape while the CRS also finds the pipe shape but with a wide rim. At higher data dimensions, JSO finds a higher index value and a clearer pipe shape across all the quantiles than the CRS"

aaa <- map_dfr(
  list(pipe1000, pipe1000_8d, pipe1000_10d, pipe1000_12d), 
  ~{
    d_dim = dim(.x)[2]
    pipe_jellyfish |> 
      filter(d == d_dim) |> 
      arrange(index_val) |>
      filter(row_number() %in% c(1, seq(5, 50, 5))) |> 
      mutate(id = factor(
        paste0(c(0, seq(5, 50, 5)) / 50 * 100, "th"),
        levels = paste0(c(0, seq(5, 50, 5)) / 50 * 100, "th"))) |> 
      compute_projection(data = .x, col = c("d", "index_val", "id"))
  }
)

p1 <- aaa |> 
  ggplot() + 
  geom_point(aes(x = V1, y = V2), size = 0.05) + 
  geom_text(data = aaa |> select(-V1, -V2) |> unique(), 
            aes(label = round(index_val,3)), x = 0, y = 3.7, 
            size = 3) + 
  facet_grid(d ~ id) + 
  xlim(-4, 4) + ylim(-4, 4) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

bbb <- pipe_better |> 
  group_by(dim) |> 
  arrange(index_val) |> 
  filter(row_number() %in% c(1, seq(5, 50, 5))) |>  
  mutate(id = factor(
    paste0(c(0, seq(5, 50, 5)) / 50 * 100, "th"),
    levels = paste0(c(0, seq(5, 50, 5)) / 50 * 100, "th")), 
    d = 4 + 2 * as.numeric(dim)) |> 
  unnest(proj)

p2 <- bbb |> 
  ggplot() + 
  geom_point(aes(x = V1, y = V2), size = 0.05) + 
  geom_text(data = bbb |> select(-V1, -V2) |> unique(), 
            aes(label = round(index_val,3)), x = 0, y = 3.7, 
            size = 3) + 
  facet_grid(d ~ id) + 
  xlim(-4, 4) + ylim(-4, 4) +
  theme_bw() +
  theme(aspect.ratio = 1, axis.ticks = element_blank(),
        axis.text = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

(p1 + ggtitle("a. JSO"))/(p2 + ggtitle("b. CRS"))


## ----------------------------------------------------------------------------------
#| label = fig-proportion,
#| fig.width = 8,
#| fig.height = 3.5,
#| fig.align = "center",
#| fig.cap = "The proportion of simulations that reach near-optimal index values in the pipe-finding problem using the holes index. The proportion is calculated based on the number of simulations, out of 50, that achieve an index value within 0.05 of the best-performing simulation. To quantify uncertainty, Bootstrap samples with 500 are generated. The thin lines represent the proportion for each of the 500 bootstrap replicates, while the thicker lines represent the mean of these bootstrap samples, connnected by lines. As the dimensionality increases, the proportion of simulations reaching the optimal index value increases."
pipe_sim_best <- sim_pipe_run_best |> 
  group_by(n_jellies, max_tries, d) |>
  summarise(I_best = max(I_max))
  
pipe_res <- sim_pipe_run_best |> 
  left_join(pipe_sim_best) |> 
  mutate(diff = abs(I_max - I_best) < 0.05) |> 
  group_by(n_jellies, max_tries, d) |> 
  summarise(proportion = sum(diff) / n()) |> 
  mutate(d = as.factor(d), n_jellies = as.factor(n_jellies))

res_df <- sim_pipe_run_best |> 
  left_join(pipe_sim_best |> ungroup() |> mutate(case_id = row_number()),
            by = c("n_jellies", "max_tries", "d")) |> 
  ungroup() |> 
  group_split(case_id)

resampling <- function(dt, size = size, replace = TRUE){
  row_id <- sample(1:50, size = size, replace = replace)
  
  dt |> filter(row_number() %in% row_id) |> 
    mutate(diff = abs(I_max - I_best) < 0.05) |> 
    summarise(proportion = sum(diff) / n()) |> 
    pull(proportion)
  
}

res_interval_df <- tibble(df = res_df) |> 
  crossing(rep = 1:500) |> 
  rowwise() |> 
  mutate(proportion = resampling(df, size = 50)) 

res_quantile <- res_interval_df |> 
  unnest(df) |> 
  group_by(n_jellies, max_tries, case_id, d) |> 
  select(-I_max, -id) |> 
  distinct() |> 
  mutate(d = as.factor(d),
         n_jellies = ifelse(d == 10, n_jellies - 1, n_jellies),
         n_jellies = ifelse(d == 12, n_jellies + 1, n_jellies)) |> 
  rename(`max # of tries` = max_tries) 
  
res_quantile |> 
  ungroup() |> 
  ggplot(aes(x = n_jellies - 0.5, y = proportion, 
             group = d, color = d)) + 
  geom_segment(aes(xend = n_jellies + 0.5), size = 0.05) + 
  geom_segment(aes(x = n_jellies - 1, xend = n_jellies + 1), 
               data = res_quantile |> summarise(proportion = mean(proportion)), 
               size = 2) + 
  geom_line(aes(x = n_jellies), 
               data = res_quantile |> summarise(proportion = mean(proportion)), 
            alpha = 0.3) + 
  scale_color_brewer(palette = "Dark2", name = "dimension") +
  scale_x_continuous(breaks = c(20, 50, 100), labels = c("20", "50", "100")) +
  facet_wrap(vars(`max # of tries`), 
             labeller = "label_both") + 
  theme_bw() +
  theme(panel.grid.minor = element_blank()) + 
  xlab("Number of jellyfish") + 
  ylab("Proportion of success")


## ----------------------------------------------------------------------------------
sim_df2 <- sim_df |>   
  mutate(n_jellies = n_jellies/10,
         max_tries = max_tries/10,
         long_time = ifelse(time > 20, 1, 0),
         P_J = ifelse(index == "stringy", 0, P_J)) |> 
  rename(proportion = P_J)

mod1 <- glm(proportion ~ smoothness + squintability + d + long_time + n_jellies + max_tries, data = sim_df2, family = binomial)


## ----fig-idx-proj-dist, eval = FALSE-----------------------------------------------
#| fig.width = 9,
#| fig.height = 5,
#| fig.cap = "Index values versus projection distance for the 12 pipe/sine-wave finding problem, after the binning procedure for calculating the squintability measure. The index values, averaged at bin width of 0.005, are scaled from 0 to 1 for comparison (holes, TIC, and stringy). The `MIC` and `TIC` index curves are convex while others are concave. The stringy curve shows an instantaneous jump to the optimum when aproaching the best basis."
## sq_basis_dist_idx |>
##   ggplot(aes(x = dist, y = index, group = index_name)) +
##   geom_line() +
##   facet_nested_wrap(~index_name + n, nrow = 3, labeller = label_both) +
##   xlab("projection distance") + ylab("index value") +
##   theme_bw()


## ----------------------------------------------------------------------------------
#| label: tbl-smoothness-squintability
#| tbl-cap: Parameters estimated from the Gaussian process (outputscale $\eta$, lengthscale $\ell$, smoothness $\nu$, and nugget $\sigma$) and scaled logistic function ($\theta_1$ to $\theta_4$ and $\varsigma$) for the pipe-finding and sine-wave finding problems.  The columns $\nu$ and $\varsigma$ represent the smoothness and squintability measures respectively.
dt <- tibble(shape = c(rep("pipe", 4), rep("sine", 8)), smoothness) |>
  left_join(squintability) |> 
  rename(d = n, smooth = smoothness) |> 
  mutate(index = ifelse(index == "dcor2d_2", "dcor", index),
         index = ifelse(index == "loess2d", "loess", index),
         index = ifelse(index == "splines2d", "splines", index),
         variance = sqrt(variance)) 
dt |> 
  knitr::kable(
    digits = 2, format = "latex", 
    col.names = c(
      "shape", "index", "d",
      "$\\eta$", "$\\ell$", "$\\nu$", "$\\sigma$",
      "$\\theta_1$", "$\\theta_2$", "$\\theta_3$", "$\\theta_4$", "$\\varsigma$"), linesep = "",
    booktabs = TRUE, escape = FALSE, row.names = TRUE, align = "c") |> 
  column_spec(1, border_left = TRUE) |> 
  column_spec(c(4, 8, 13), border_right = TRUE) |> 
  kable_styling(font_size = 10) 


## ----------------------------------------------------------------------------------
#| label: tbl-mod-output
#| tbl-cap: "Jellyfish success rate relative to index properties and jellyfish hyper-parameters. This is the summary from a logistic regression fit to smoothness, squintability, dimension, running time, number of jellyfish and maximum number of tries. Interestingly, squintability and dimension strongly affect jellyfish optimisation success. The number of jellies marginally affects success, but index smoothness, running longer and increasing the number of tries do not."
#| fig-pos: t
library(gtsummary)
set_gtsummary_theme(theme_gtsummary_journal("jama"))
mod1 |> tbl_regression(exponentiate=TRUE, 
                       pvalue_fun = label_style_pvalue(digits = 3),
                       label = list(smoothness = "Smoothness", 
                                    squintability = "Squintability",
                                    d = "Dimension", 
                                    long_time = "Long runtime", 
                                    n_jellies = "Number of jellyfish", 
                                    max_tries = "Maximum number of tries")) |>
  bold_p(t = 0.1)


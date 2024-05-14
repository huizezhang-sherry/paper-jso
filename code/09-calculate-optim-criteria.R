library(tourr)
library(minerva)
library(cassowaryr)
library(energy)
library(GpGp)

calc_smoothness <- function(idx){

  idx <- sym(idx)
  set.seed(123)
  seed <- sample(1: 10000, size = 300)
  res <- map_dfr(1:300, ~{set.seed(seed[.x]); tibble(basis = list(tourr::basis_random(n = 6, d = 2)))}) |>
    rowwise() |>
    mutate(proj_dist = tourr::proj_dist(matrix(c(rep(0, 8), 1, 0, 0, 1), nrow = 6, byrow = TRUE), basis),
           index_val = get(idx)()(as.matrix(sine1000) %*% basis))

  fit <- fit_model(res$index_val, res$proj_dist, #start_parms = c(0.001, 0.5, 2, 2),
                   as.matrix(rep(1,nrow(res))), "matern_isotropic")
  # check convergence
  res <- as_tibble_row(fit$covparms, .name_repair = "unique")
  colnames(res) <- c("variance", "range", "smoothness", "nugget", "convergence")
  res |> mutate(convergence = fit$conv, idx = as.character(idx))

}

idx_names <- c("dcor2d", "loess2d", "MIC", "TIC", "stringy", "splines2d")
res <- map_dfr(idx_names, calc_smoothness)

################################################################
stringy <- function(){
  function(mat){
    cassowaryr::sc_stringy(mat[,1], mat[,2])
  }
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


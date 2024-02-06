library(tourr)
library(MaxSkew)
library(MultiSkew)
library(tidyverse)
x <- matrix(c(
  9160, 1759, 2691,
    12, 7218, 7655,
  4624, 4735, 1887,
  4243, 1527, 2875,
  4609, 3411,  911,
  7702, 6074, 5762,
  3225, 1917, 6834,
  7847, 7384, 5466,
  4714, 2428, 4257,
   358, 9174, 6444
), ncol = 3, byrow = TRUE)
vec <- c(0.816, 0.446, 0.367)
vec <- get_init(x) # hack into the code .MaxSkewThree to get .linear
res1 <- MaxSkew(x, iterations = 3, components = 1, plot = FALSE)
res1 <- x %*% vec
skewness()(res1)

i <- 1
j <- 1
for (i in 1:3){
  for (j in 1:3){
    y<-x[,j]#j-th column of the data matrix
    vec[j]<-c(0) #removes the j-th variable from the linear combination
    MaxSkew:::.MaxSkewBiv(x%*%vec,y)
    s<-.skewdesideredBIV #skewness of the new linear combination
    vec[j]<-.linearBIV[2]/.linearBIV[1]
    res1 <- x %*% vec
  }
  print(skewness()(res1))

}

get_init <- function(data) {
  n <- nrow(data)
  d <- ncol(data)

  provaconcentrationmatrix <- solve(cov(data) * (n - 1) / n)
  provaSS <- provaconcentrationmatrix
  provaaut <- eigen(provaSS)
  sort(provaaut$values, decreasing = FALSE)
  provaVV <- matrix(c(0), nrow = d)
  provaDD <- diag(sort(provaaut$values, decreasing = FALSE))

  for (i in ncol(provaaut$vectors):1) {
    provaVV <- cbind(provaVV, provaaut$vectors[, i])
  }

  provaVV <- provaVV[, 2:ncol(provaVV)]
  provaDD.sqrt <- solve(sqrt(provaDD))
  provaAA <- provaVV %*% provaDD.sqrt %*% t(provaVV)
  provaQ <- solve(provaAA)

  x.mean <- colMeans(data)#vector mean
  m <- sweep(data, 2, x.mean)#centered data
  Z <- m %*% provaQ

  A <- matrix(c(0), nrow = d * d, ncol = d)#initialization of the matrix which sums the tensor products

  for (i in 1:n) {
    A <- A + kronecker(kronecker(Z[i, ], t(Z[i, ])), Z[i, ])#update of the matrix summing the tensor products
  }

  #INITIALIZATION
  s <- svd(A / n)#singular value decomposition of the third standardized cumulant
  c <- matrix(-s$v[, 1])#first right singular vector of the third multivariate cumulant
  v <- provaQ %*% c #starting value of "linear"
  return(as.matrix(v))
}

skewness <- function(){
  function(mat){
    n <- length(mat)
    s3<-sqrt(var(mat)*(n-1)/n)^3 #we compute the skewness of projection...
    mx<-mean(mat)
    sk<-sum((mat-mx)^3)/s3
    M<-sk/n#
    return(abs(M))
  }
}

res <- animate_dist(
  x, guided_tour(
    skewness(), d = 1, n_jellies = 50,
    search_f =  search_jellyfish,
    max.tries = 100, min.tries = 100
  )
)

data("OLYMPIC_DECATHLON_2016")
dt <- OLYMPIC_DECATHLON_2016[,-c(1:3)] |> as.matrix()
dt2 <- OLYMPIC_DECATHLON_2016[-c(11, 23),-c(1:3)]
prcomp_res <- prcomp(dt, scale = TRUE)
res0 <- prcomp_res$x[,1:2] |>
  bind_cols(OLYMPIC_DECATHLON_2016[,c(1:3)]) |>
  as_tibble()
res0 |>
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point() +
  #geom_label(aes(label = ATHLETE)) +
  theme(aspect.ratio = 1)


vec <- prcomp_res$rotation[,1]
i <- 1
j <- 1
for (i in 1:20){
  for (j in 1:ncol(dt)){
    y<-dt[,j]#j-th column of the data matrix
    vec[j]<-c(0) #removes the j-th variable from the linear combination
    MaxSkew:::.MaxSkewBiv(dt %*%vec,y)
    s<-.skewdesideredBIV #skewness of the new linear combination
    vec[j]<-.linearBIV[2]/.linearBIV[1]
    res1 <- dt %*% vec
  }
  print(skewness()(res1))

}
mx<-colMeans(dt) # vector mean
n <- nrow(dt)
.projection <- res1
mp<-mean(.projection) #mean of projection
ssquarep<-c(var(.projection)*(n-1)/n) #variance of projection
spx<-cov(.projection,dt)*(n-1)/n #covariance (projection and data)
pen<-spx/ssquarep #slope
intercept<-mx-pen*mp #intercept
teo<-matrix(1,n,1)%*%intercept+.projection%*%pen #matrix of predicted values
res<-dt-teo #matrix of  the residuals
covres<-cov(res)*(n-1)/n #covariance of the residuals
eigen(covres)#spectral decomposition  of covres
o <- order(eigen(covres)$values, decreasing=FALSE)#we reorder
eigen(covres)$values[o]#we reorder
V<-eigen(covres)$vectors[,o]# spectral decomposition of covres,eigenvector ordered in ascending order
h<-ncol(dt)-1+1
proiezione<-res%*%V[,2:h]#data projection orthogonal to skewed components
dt<-proiezione

vec <- prcomp(dt, scale = TRUE)$rotation[,1]
i <- 1
j <- 1
for (i in 1:20){
  for (j in 1:ncol(dt)){
    y<-dt[,j]#j-th column of the data matrix
    vec[j]<-c(0) #removes the j-th variable from the linear combination
    MaxSkew:::.MaxSkewBiv(dt %*%vec,y)
    s<-.skewdesideredBIV #skewness of the new linear combination
    vec[j]<-.linearBIV[2]/.linearBIV[1]
    res1 <- dt %*% vec
  }
  print(skewness()(res1))

}

res1 <- MaxSkew(as.matrix(dt), iterations = 10, components = 2, plot = FALSE)
athlete <- gsub(".*\\s+", "", OLYMPIC_DECATHLON_2016[,1:3]$ATHLETE)
athlete[12] <- "HELCELET"
athlete[19] <- "DISTELBERGER"
athlete_vec <- tolower(athlete)
skewness()(res1[,1])
skewness()(res1[,2])
as_tibble(res1) |>
  mutate(athlete = athlete_vec) |>
  ggplot(aes(x = V1, y = V2)) +
  geom_point() +
  geom_text(aes(label = athlete),  hjust = 1.2, vjust = 0.4) +
  theme(aspect.ratio = 1) +
  xlim(c(-11, -4)) +
  ylim(c(-2, 4))

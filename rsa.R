create_rsa <- function(x, initial_size = 100, quantile = 0.5, hn = function(n) n^(-0.5), a = 0.25) { ## Recursive Stochastic Approximator
  
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
  }
  
  rsa <- list(
    data = x,
    initial_size = initial_size,
    quantile = quantile,
    hn = hn,
    a = 0.25,
    initialized = FALSE
  )
  
  class(rsa) <- "rsa"
  
  if (nrow(rsa$data) >= rsa$initial_size) {rsa <- initialize.rsa(rsa)}
  
  return(rsa)
  
}


update_xi <- function(xi, dn1, n, data, quantile) {
  return(xi - (dn1 / n)*((data <= xi) - quantile))
} 

dn <- function(fnxi, d0, n, a) {
  return(min(fnxi^(-1), d0*n^a))
}

initialize_fnxi <- function(xi, data, hn, invscale) {
  return((abs(xi - data) <= invscale*hn(1)) / (2*invscale*hn(1)))
}

update_fnxi <- function(fnxi, n, xi, data, hn, invscale) {
  return((1 / n)*((n-1)*fnxi + (abs(xi - data) <= invscale*hn(n)) / (2*invscale*hn(n))))
}

density_estimate <- function(xi, data) {
  tmp <- density(data)
  return((tmp$y[which.min(abs(xi - tmp$x))])^(-1))
}

initialize.rsa <- function(rsa) {
  
  ## Set initial values
  rsa$n <- 0
  sx <- apply(rsa$data, 2, sort)
  rsa$xi0 <- sx[floor(nrow(sx) * rsa$quantile)]
  rsa$d0 <- apply(rsa$data, 2, IQR)^(-1)
  
  # rsa$invscale <- sapply(1:length(rsa$xi), function(x) density_estimate(rsa$xi[x], rsa$data[,x]))
  rsa$invscale <- rsa$d0
  # print(rsa$d0)
  # print(rsa$invscale)
  
  tmp <- rsa$data[1,]
  rsa$data <- rsa$data[-1,,drop=FALSE]
  rsa$n <- rsa$n + 1
  rsa$fnxi <- initialize_fnxi(rsa$xi0, tmp, rsa$hn, rsa$invscale)
  rsa$xi <- update_xi(rsa$xi0, rsa$d0, rsa$n, tmp, rsa$quantile)
  rsa$dn <- dn(rsa$fnxi, rsa$d0, rsa$n, rsa$a)

  # Update with rest of data
  iter <- 0
  while (nrow(rsa$data) > 0 & iter < 1e4) {
    tmp <- rsa$data[1,]
    rsa$data <- rsa$data[-1,,drop=FALSE]
    rsa$n <- rsa$n + 1
    rsa$fnxi <- update_fnxi(rsa$fnxi, rsa$n, rsa$xi, tmp, rsa$hn, rsa$invscale)
    rsa$xi <- update_xi(rsa$xi, rsa$dn, rsa$n, tmp, rsa$quantile)
    rsa$dn <- dn(rsa$fnxi, rsa$d0, rsa$n, rsa$a)
    iter <- iter + 1
  }
  rsa$initialized <- TRUE
  return(rsa)
  
}

update.rsa <- function(rsa, x) { ## Recursive Stochastic Approximator
  
  rsa$data <- rbind(rsa$data, as.matrix(x, ncol = ncol(rsa$data)))
  
  if (!rsa$initialized) {
    if (nrow(rsa$data) >= rsa$initial_size) {rsa <- initialize.rsa(rsa)}
    return(rsa)
  }
  
  if (rsa$initialized) {
    iter <- 0
    while (nrow(rsa$data) > 0 & iter < 1e4) {
      tmp <- rsa$data[1,]
      rsa$data <- rsa$data[-1,,drop=FALSE]
      rsa$n <- rsa$n + 1
      rsa$fnxi <- update_fnxi(rsa$fnxi, rsa$n, rsa$xi, tmp, rsa$hn, rsa$invscale)
      rsa$xi <- update_xi(rsa$xi, rsa$dn, rsa$n, tmp, rsa$quantile)
      rsa$dn <- dn(rsa$fnxi, rsa$d0, rsa$n, rsa$a)
      iter <- iter + 1
    }
    return(rsa)
  }
  
}



## Test
test1 <- rnorm(1000, mean = 2, sd = 2)
quantile(test1, 0.5)

initial <- test1[1:50]
update1 <- test1[51:101]
update2 <- test1[102:200]
final <- test1[201:1000]

rsa <- create_rsa(initial, quantile = 0.5)
rsa <- update(rsa, update1)
rsa <- update(rsa, update2)
rsa <- update(rsa, final)
rsa$xi


## Simulation
set.seed(1234)
sss <- c(100, 200, 500)
qntl <- 0.75
niter <- 1000
initial_size <- 10
res <- matrix(nrow = 3, ncol = 3)

for (j in 1:length(sss)) {
  
  sq <- numeric(niter)
  sa <- numeric(niter)
  ss <- sss[j]
  
  for (i in 1:niter) {
    tmp <- rnorm(ss)
    sq[i] <- quantile(tmp, qntl)
    rsa <- create_rsa(tmp[1:initial_size], initial_size = initial_size, quantile = qntl)
    rsa <- update(rsa, tmp[(initial_size + 1):ss])
    sa[i] <- rsa$xi
  }
  
  res[j,1] <- ss*mean((sa - qnorm(qntl))^2)
  res[j,2] <- ss*mean((sq - qnorm(qntl))^2)
  res[j,3] <- ss*mean((sa - sq)^2)

}
res


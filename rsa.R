create_rsa <- function(x, initial_size = 100, quantiles = c(0.5), hn = function(n) n^(-0.5), a = 0.25) { ## Recursive Stochastic Approximator
  
  if (is.vector(x)) {
    x <- as.matrix(x, ncol = 1)
  }
  
  rsa <- list(
    data = x,
    initial_size = initial_size,
    quantiles = quantiles,
    hn = hn,
    a = 0.25,
    initialized = FALSE
  )
  
  class(rsa) <- "rsa"
  
  if (nrow(rsa$data) > rsa$initial_size) {rsa <- initialize.rsa(rsa)}
  
  return(rsa)
  
}


update_xi <- function(xi, dn, n, data, quantiles) {
  return(xi - (dn / (n + 1))*((data <= xi) - quantiles))
} 

dn <- function(fnxi, d0, n, a) {
  return(min(fnxi^(-1), d0*n^a))
}

initialize.rsa_fnxi <- function(xi, data, hn) {
  return((abs(xi - data) <= hn(1)) / (2*hn(1)))
}

update_fnxi <- function(fnxi, n, xi, data, hn) {
  return((1 / (n + 1))*(n*fnxi + (abs(xi - data) <= hn(n)) / (2*hn(n))))
}

initialize.rsa <- function(rsa) {
  
  ## Set initial values
  rsa$n <- 0
  sx <- apply(rsa$data, 2, sort)
  rsa$xi <- sx[floor(rsa$initial_size * rsa$quantiles)]
  rsa$d0 <- apply(rsa$data, 2, IQR)^(-1)
  
  tmp <- rsa$data[1,]
  rsa$data <- rsa$data[-1,,drop=FALSE]
  rsa$n <- rsa$n + 1
  rsa$fnxi <- initialize.rsa_fnxi(rsa$xi, tmp, rsa$hn)
  rsa$dn <- dn(rsa$fnxi, rsa$d0, rsa$n, rsa$a)
  rsa$xi <- update_xi(rsa$xi, rsa$dn, rsa$n, tmp, rsa$quantiles)

  # Update with rest of data
  iter <- 0
  while (nrow(rsa$data) > 0 & iter < 1e6) {
    tmp <- rsa$data[1,]
    rsa$data <- rsa$data[-1,,drop=FALSE]
    rsa$n <- rsa$n + 1
    rsa$fnxi <- update_fnxi(rsa$fnxi, rsa$n, rsa$xi, tmp, rsa$hn)
    rsa$dn <- dn(rsa$fnxi, rsa$d0, rsa$n, rsa$a)
    rsa$xi <- update_xi(rsa$xi, rsa$dn, rsa$n, tmp, rsa$quantiles)
    iter <- iter + 1
  }
  rsa$initialized <- TRUE
  return(rsa)
  
}

update.rsa <- function(rsa, x) { ## Recursive Stochastic Approximator
  
  rsa$data <- rbind(rsa$data, as.matrix(x, ncol = ncol(rsa$data)))
  
  if (!rsa$initialized) {
    if (nrow(rsa$data) > rsa$initial_size) {rsa <- initialize.rsa(rsa)}
    return(rsa)
  }
  
  if (rsa$initialized) {
    iter <- 0
    while (nrow(rsa$data) > 0 & iter < 1e6) {
      tmp <- rsa$data[1,]
      rsa$data <- rsa$data[-1,,drop=FALSE]
      rsa$n <- rsa$n + 1
      rsa$fnxi <- update_fnxi(rsa$fnxi, rsa$n, rsa$xi, tmp, rsa$hn)
      rsa$dn <- dn(rsa$fnxi, rsa$d0, rsa$n, rsa$a)
      rsa$xi <- update_xi(rsa$xi, rsa$dn, rsa$n, tmp, rsa$quantiles)
      iter <- iter + 1
    }
    return(rsa)
  }
  
}



## Test
set.seed(1234)
test1 <- rnorm(1000, mean = 2)
median(test1)

intial <- test1[1:50]
update1 <- test1[51:101]
update2 <- test1[102:200]
final <- test1[201:1000]

rsa <- create_rsa(intial)
rsa <- update(rsa, update1)
rsa <- update(rsa, update2)
rsa <- update(rsa, final)
rsa$xi


require(matrixStats)
rbridgemixing <- function(n, phi, k2 = 100){
  mat = matrix(rgeom(n*k2, prob = 1-phi^2) + 1, k2, n)
  expscalemat = (2/phi^2)*matrixStats::colCumsums(mat)^(-2)
  colSums(matrix(rexp(n*k2),k2, n)*expscalemat)
}
rmbridge <- function(n, phi, R){
  p = ncol(R)
  lambda = rbridgemixing(n, phi)
  x = matrix(rnorm(n*p), n, p)%*%chol(R)
  return(sqrt(lambda)*x)
}

dmbridge <- function(x, phi, R, log = F, kmax = 200){
  x = as.numeric(x)
  Ckphi = (1:kmax)-0.5 + (-1)^(1:kmax)*(phi-0.5)
  d = ncol(R)
  Rinv = solve(R)
  Rinvx = Rinv%*%x # vector
  xRx = sum(x*Rinvx)
  summand = (-1)^(2:(kmax+1))*Ckphi/(pi^2*Ckphi^2/phi^2 + xRx)^((d+1)/2)
  logsum = log(sum(summand))
  logdensity = logsum + lgamma((d+1)/2) - 2*log(phi) - (d-1)/2*log(pi) - 0.5*determinant(R, logarithm = T)$modulus
  logdensity = as.numeric(logdensity)
  if(log){
    return(logdensity)
  }else{
    return(exp(logdensity))
  }
}



find_var_bridge <- function(phi){
  pi^2*(phi^(-2)-1)/3
}

find_phi_bridge <- function(var){
  pi/sqrt(3*var+pi^2)
}


# Below are copied from "bridgedist" R package by Bruce Swihart
# source: https://github.com/swihart/bridgedist/blob/master/R/Bridge.R

dbridge <- function(x, phi = 1/2, log = FALSE){
  maxlength <- max(length(x), length(phi))
  x <- rep(x    , length=maxlength)
  phi <- rep(phi, length=maxlength)
  d=1/(2*pi) * sin(phi*pi) / (cosh(phi*x) + cos(phi*pi))
  if(log[1]) d = log(d)
  d
}

pbridge <- function(q, phi = 1/2, lower.tail = TRUE, log.p = FALSE){
  maxlength <- max(length(q), length(phi))
  q <- rep(q    , length=maxlength)
  phi <- rep(phi, length=maxlength)
  p=1 - 1/(pi*phi) * (pi/2 - atan( (exp(phi*q) + cos(phi*pi)) / sin(phi*pi) ))
  if(!lower.tail[1]) p = 1-p
  if(log.p[1]) p = log(p)
  p
}

qbridge <- function(p, phi = 1/2, lower.tail = TRUE, log.p = FALSE){
  maxlength <- max(length(p), length(phi))
  p <- rep(p    , length=maxlength)
  phi <- rep(phi, length=maxlength)
  if(log.p[1]) p = exp(p)
  if(!lower.tail[1]) p = 1-p
  q=1/phi * log( sin(phi*pi*p) / sin(phi*pi*(1-p)) )
  q
}

rbridge <- function(n, phi = 1/2){
  if(length(n) == 1){
    r <- qbridge(stats::runif(n), phi[1])
  }else{
    r <- qbridge(stats::runif(length(n)), rep(phi,length=length(n)))
  }
  r
}

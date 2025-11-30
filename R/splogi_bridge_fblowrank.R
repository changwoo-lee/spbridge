#' Spatial logistic model with bridge process random effect with low-rank structure (fully Bayes)
#'
#' \deqn{
#' \operatorname{logit}\left[ \Pr(y_{ij} = 1 \mid X_{ij}, u(s_i)) \right]
#'   = X_{ij}^\top \beta + u(s_i)
#' }
#' where u(s) ~ Bridge process with Matern correlation kernel with low-rank structure, i=1,...n corresponds to n spatial locations and
#' j=1,... N_i correspond to n_i responses at location i,
#' resulting data of size N = N_1 + ... N_n.
#'
#' Priors are specified by "priors" argument, which is a list of hyperparameters.
#' Specifically, we set zero centered normal or t prior for beta,
#' uniform prior for Matern range parameter rho (see fields::Matern),
#' and for phi, we use prior that induces half-Cauchy prior on the standard deviation of u.
#'
#' priors is a named list with the following possible elements:
#' \describe{
#'   \item{beta_intercept_scale}{scale of intercept parameter (default 10)}
#'   \item{beta_scale}{scale of other beta parameters (default 2.5)}
#'   \item{beta_df}{degrees of freedom for t prior on beta (default Inf, normal prior)}
#'   \item{logpriorphi}{function for log-prior on phi (default prior that induces half-Cauchy prior on the standard deviation of u.)}
#'   \item{rho_lb}{lower bound for range parameter rho (default min distance between coords)}
#'   \item{rho_ub}{upper bound for range parameter rho (default max distance between coords)}
#' }
#'
#' @param y N x 1 binary vector
#' @param X N x p fixed-effect design matrix, including intercept
#' @param id N x 1 vector of spatial location id. When N=n, it is point-referenced data
#' @param coords n x 2 matrix of spatial coordinates
#' @param coords_knot q x 2 matrix of knot coordinates
#' @param priors list of prior hyperparameters, see details
#' @param smoothness postive numeric, Matern smoothness parameter
#' @param nburn number of burn-in interation
#' @param nsave number of posterior samples
#' @param nthin thin-in rate
#' @param nparticle number of particles in particle marginal Metropolis-Hastings
#' @param verbose logical, whether to print progress
#'
#' @return list, including posterior samples of beta, phi, rho, and u saved in "post_save".
#' Population-averaged coefficient, which is beta*phi, is saved separately as "betam_save"
#' @export
#'
#' @examples
#' \donttest{
#' library(spbridge)
#' data(gambia)
#' N = length(gambia$pos) # 2035
#' y = gambia$pos # binary response, N = 2035 by 1 vector
#' # define id based on spatial coords unique values
#' id = as.numeric(factor(paste(gambia$x, gambia$y)))
#' n = length(unique(id)) # 65
#' coords = unique(cbind(gambia$x, gambia$y)/1000) # n by 2 matrix, in km
#
#' # standardized covariates following Gelman et al (2008)
#' intercept = rep(1,N) # intercept
#' age = scale(gambia$age/365, scale = 2*sd(gambia$age/365)) # in years
#' netuse = gambia$netuse - mean(gambia$netuse)
#' treated = gambia$treated - mean(gambia$treated)
#' green = scale(gambia$green, scale = 2*sd(gambia$green))
#' green2 = scale(gambia$green^2, scale = 2*sd(gambia$green^2))
#' healthctr = gambia$phc - mean(gambia$phc)
#'
#' X = cbind(intercept, age, netuse, treated, green, green2, healthctr)
#' colnames(X) = c("(Intercept)", "age", "netuse", "treated", "green", "green2", "healthctr")
#'
#' centers = c(attr(age, "scaled:center"), mean(gambia$netuse), mean(gambia$treated),
#'             attr(green, "scaled:center"), attr(green2, "scaled:center"), mean(gambia$phc))
#' scales = c(attr(age, "scaled:scale"), 1, 1, attr(green, "scaled:scale"), attr(green2, "scaled:scale"), 1)
#' fit_bridge = spbridge::splogi_bridge_fblowrank(y = y,
#'                                         X = X,
#'                                         id = id,
#'                                         priors = list(beta_intercept_scale = 10,
#'                                                       beta_scale = 2.5, beta_df = Inf,
#'                                                       rho_lb = 0.01, rho_ub = 100),
#'                                         coords = coords,
#'                                         smoothness = 0.5, nburn = 1000, nsave = 10000, nthin = 1)
#' }
#'
#'
splogi_bridge_fblowrank <- function(y, X, id,
                                     coords, coords_knot,
                                     priors = list(beta_intercept_scale = 10,
                                                   beta_scale = 2.5, beta_df = Inf,
                                                   logpriorphi = NULL,
                                                   rho_lb = NULL, rho_ub = NULL),
                                     smoothness = 1.5,
                                     nburn = 100, nsave = 1000, nthin = 1, nparticle = 20, verbose=TRUE){

  t_start = Sys.time()
  #############################################
  n = length(unique(id)) # this is number of spatial locations, not total obs
  if(nrow(coords)!=n) stop("nrow(coords) must be equal to number of unique elements of id")
  y = as.numeric(y)
  X = as.matrix(X)
  N = nrow(X)
  p = ncol(X)
  if(length(y)!=N ) stop("length(y) must be equal to nrow(X)")
  if(length(id)!=N) stop("length(id) must be equal to nrow(X)")
  # sanity check for range and smoothness which should be positive reals
  if(smoothness <= 0) stop("smoothness must be positive reals")
  # sanity check for z which should be a vector of binary
  if(!is.vector(y) || !all(y %in% c(0,1))) stop("y must be a vector of binary")
  if(N>n){
    Z = lme4::getME(lme4::lmer(y~(1|id)), "Z")
  }else{
    Z = Matrix(diag(n))
  }


  distmat = fields::rdist(coords)

  # set hyperparameters
  if(is.null(priors$beta_df)){
    beta_df = Inf; priors$beta_df = Inf
  }else{
    beta_df = priors$beta_df
  }
  if(is.null(priors$beta_intercept_scale)){
    beta_intercept_scale = 10; priors$beta_intercept_scale = 10
  }else{
    beta_intercept_scale = priors$beta_intercept_scale
  }
  if(is.null(priors$beta_scale)){
    beta_scale = 2.5; priors$beta_scale = 2.5
  }else{
    beta_scale = priors$beta_scale
  }
  if(is.null(priors$logpriorphi)){
    logpriorphi = function(x) log(2*sqrt(3)/((pi^2 - (pi^2-3)*x^2)*sqrt(1-x^2))) # corresponding to half-cauchy on stdev
    priors$logpriorphi = logpriorphi
  }else{
    # check whether priors$phi is a function
    if(!is.function(priors$logpriorphi)) stop("priors$logpriorphi must be a function")
    logpriorphi = priors$logpriorphi
  }
  if(is.null(priors$rho_lb)){ #uniform prior, lower bound
    rho_lb = min(distmat[distmat>0]); priors$rho_lb = rho_lb
  }else{
    rho_lb = priors$rho_lb
  }
  if(is.null(priors$rho_ub)){
    rho_ub = max(distmat); priors$rho_ub = rho_ub
  }else{
    rho_ub = priors$rho_ub
  }

  # Initialize
  phi = 0.5
  rho = (rho_lb + rho_ub)/2
  beta_s = c(beta_intercept_scale, rep(beta_scale,p-1))


  q = nrow(coords_knot)
  #distmat_nn = fields::rdist(coords)
  distmat_qq = fields::rdist(coords_knot)
  distmat_nq = fields::rdist(coords, coords_knot)

  Rqq = fields::Matern(distmat_qq, range = rho, smoothness = smoothness)
  Rnq = fields::Matern(distmat_nq, range = rho, smoothness = smoothness)
  #dnn = 1-diag(Rnq%*%solve(Rqq)%*%t(Rnq))
  dnn = 1-rowSums(t(solve(Rqq, t(Rnq))) * Rnq)
  Dnn = Matrix::Diagonal(n, dnn)
  Dnninv = solve(Dnn)
  #if(any(dnn > 0.5)) warning("too small initial rho?")
  #R = Sig_pp + diag(dnn) # modified predictive process
  #if(any(diag(R)!=1)) stop("diag(R)!=1")
  #Rinv = solve(R)




  #R = fields::Matern(distmat, range = rho, smoothness = smoothness)
  #Rinv = solve(R)
  lambda_particles = rbridgemixing(nparticle, phi = phi) # each row is a set of particles
  lambda = lambda_particles[1]
  beta = rep(0,p)
  gamma = rep(2.5,p) # scale mixture
  linpred = X%*%beta
  omega = rep(1,N)
  #omega_grp = Rfast::group(omega, id, method = "sum")
  omega_grp = as.numeric(rowsum(as.matrix(omega), group = id))
  # Saving objects
  nmcmc = nburn + nsave*nthin
  rho_save = array(0, dim = c(nsave))
  beta_save = array(0, dim = c(nsave, p))
  phi_save = array(0, dim = c(nsave))
  lambda_save = array(0, dim = c(nsave))
  u_save = array(0, dim =c(nsave, n))
  loglik_save = array(0, dim = c(nsave, N))
  acc_save = numeric(nmcmc)

  #t1  = t2 = t3 = t4 = 0
  # adaptive MH tuning (Harrio)
  MH_eps = 0.001
  MH_s_d = (2.38)^2/2 # denominator 2 corresponds to dimension (phi, rho)
  C0 = MH_s_d*diag(2)
  start_adapt = 100 # adapt after 100 iterations

  # Run MCMC
  # pre-calculate
  Xtym0.5 = crossprod(X, (y - 0.5))
  Ztym0.5 = crossprod(Z, (y - 0.5))
  # initialize
  ZtOmegaX = (Matrix::t(Z*omega)%*%X)
  XtOmegaX = (Matrix::t(X*omega)%*%X)


  t_end = Sys.time()
  t_premcmc = difftime(t_end, t_start, units = "secs")
  if(verbose){
    pb <- txtProgressBar(style=3)
  }
  isave = 1
  ##### Start of MCMC #####
  t_start = Sys.time()
  for(imcmc in 1:(nburn + nsave*nthin)){
    if(verbose){
      setTxtProgressBar(pb, imcmc/(nburn + nsave*nthin))
    }

    ##### Step 1: update beta #####
    #t1_start = Sys.time()

    d1 = 1/(lambda*omega_grp + 1/dnn)
    D1 = Matrix::Diagonal(n, d1)
    d2 = 1/dnn - d1/dnn^2
    D2 = Matrix::Diagonal(n, d2)
    #all.equal(t(Rnq)%*%D2%*%Rnq, crossprod(sqrt(d2)*Rnq))

    nnmat_inv = lambda*D1 + lambda*D1%*%Dnninv%*%Rnq%*%solve(Rqq + crossprod(Rnq*sqrt(d2)),t(Rnq))%*%Dnninv%*%D1
    #nnmat_inv = solve(1/lambda*Rinv + diag(omega_grp))

    XtSigma_invX = XtOmegaX - Matrix::t(ZtOmegaX)%*%nnmat_inv%*%ZtOmegaX
    XtSigma_invY = Xtym0.5 - Matrix::t(ZtOmegaX)%*%nnmat_inv%*%Ztym0.5

    if(!is.infinite(beta_df)){ # normal prior
      Q_beta = XtSigma_invX + diag(1/gamma, p)
    }else{ # t prior
      Q_beta = XtSigma_invX + diag(1/beta_s^2, p)
    }
    b_beta = XtSigma_invY # assuming prior mean is zero
    beta = as.numeric(spam::rmvnorm.canonical(1, b_beta, Q_beta))
    #t1 = t1 + as.numeric(difftime(Sys.time(), t1_start, units = "secs"))

    # update Xbeta
    Xbeta = X%*%beta



    # update gamma
    if(!is.infinite(beta_df)){
      gamma = 1/rgamma(p, shape = beta_df/2 + 1/2, rate = beta_s^2*beta_df/2 + beta^2/2)
    }


    #t2_start = Sys.time()
    ##### Step 2-1: update lambda and phi, and range parameter
    # transform phi\in (0,1) to real based on logistic transform
    # transform rho \in rho_lb, rho_ub to real based on logistic transform
    # on this transformed space, run random walk with bivariate normal proposal

    phi_trans = glogit(phi, xmin = 0, xmax = 1)
    rho_trans = glogit(rho, xmin = rho_lb, xmax = rho_ub)

    if(imcmc < start_adapt){
      proposal = c(phi_trans, rho_trans) + spam::rmvnorm(1, rep(0,2), C0)
    }else{
      proposal = c(phi_trans, rho_trans) + spam::rmvnorm(1, rep(0,2), Ct)
    }
    phi_trans_star = proposal[1]; phi_star = inv_glogit(phi_trans_star, 0, 1)
    rho_trans_star = proposal[2]; rho_star = inv_glogit(rho_trans_star, rho_lb, rho_ub)

    lambda_particles_star = rbridgemixing(nparticle, phi = phi_star)

    #D3 = Matrix::Diagonal(n, d3)


    y0.5_Xbomega = (y - 0.5) - Xbeta*omega
    #y0.5_Xbomega_grp = Rfast::group(y0.5_Xbomega, id, method = "sum") # t(Z)%*%y0.5_Xbomega
    y0.5_Xbomega_grp = as.numeric(rowsum(as.matrix(y0.5_Xbomega), group = id))
    linpred_proxy = y0.5_Xbomega_grp/omega_grp # location in L(lambda)


    # collapsed likelihood evaluations
    logweights = rep(0, nparticle)
    for(iparticle in 1:nparticle){
      d3inv = (lambda_particles[iparticle]*dnn + 1/omega_grp)
      #mvnfast::dmvn(linpred_proxy, mu = rep(0,n), sigma = diag(1/omega_grp)+lambda_particles[iparticle]*R, log = T)
      logweights[iparticle] = dmvn_lowrankstr(linpred_proxy, mu = rep(0,n),
                                              X = Rnq, Dinv = Rqq/lambda_particles[iparticle],
                                              Dlogdet = NULL, Ediag = d3inv, log = T)
    }


    Rqq_star = fields::Matern(distmat_qq, range = rho_star, smoothness = smoothness)
    Rnq_star = fields::Matern(distmat_nq, range = rho_star, smoothness = smoothness)
    dnn_star = 1-rowSums(t(solve(Rqq_star, t(Rnq_star))) * Rnq_star )
    logweights_star = rep(0, nparticle)
    for(iparticle in 1:nparticle){
      d3inv = (lambda_particles_star[iparticle]*dnn_star + 1/omega_grp)
      logweights_star[iparticle] = dmvn_lowrankstr(linpred_proxy, mu = rep(0,n),
                                              X = Rnq_star, Dinv = Rqq_star/lambda_particles_star[iparticle],
                                              Dlogdet = NULL, Ediag = d3inv, log = T)
    }

    acc_ratio = log(phi_star) + log(1-phi_star) - log(phi) - log(1-phi) +
      (log(rho_star-rho_lb) + log(rho_ub - rho_star)) - (log(rho-rho_lb) + log(rho_ub - rho)) + # log(rho_ub - rho_lb) terms or cancelled out # https://www.wolframalpha.com/input?i=d%2Fdx+log%28%28x-a%29%2F%28b-a%29%2F%281-%28x-a%29%2F%28b-a%29%29%29
      logpriorphi(phi_star) - logpriorphi(phi) + # uniform prior on rho
      matrixStats::logSumExp(logweights_star) - matrixStats::logSumExp(logweights)

    if(log(runif(1)) < acc_ratio){
      phi = phi_star; phi_trans = glogit(phi, xmin = 0, xmax = 1)
      rho = rho_star; rho_trans = glogit(rho, xmin = rho_lb, xmax = rho_ub)
      Rqq = Rqq_star
      Rnq = Rnq_star
      dnn = dnn_star
      Dnn = Matrix::Diagonal(n, dnn)
      Dnninv = solve(Dnn)

      lambda = lambda_particles[sample(1:nparticle, size = 1, prob = exp(logweights - max(logweights)))]
      acc_save[imcmc] = 1
      lambda_particles = lambda_particles_star
      logweights = logweights_star
    }

    # mut and Ct are recursively updated
    if(imcmc == 1){
      mut = c(phi_trans, rho_trans)
      Ct = MH_s_d*MH_eps*diag(2)
    }else{
      tmpmu = (mut*(imcmc-1)+c(phi_trans, rho_trans))/imcmc
      # eq (3) of Haario et al. 2001
      Ct = (imcmc-1)*Ct/imcmc+MH_s_d/imcmc*(imcmc*tcrossprod(mut)-
                                              (imcmc+1)*tcrossprod(tmpmu) +
                                              tcrossprod(c(phi_trans, rho_trans))+
                                              MH_eps*diag(2))
      mut = tmpmu
    }


    #t2 = t2 + as.numeric(difftime(Sys.time(), t2_start, units = "secs"))

    ##### Step 2-2: update random effects u
    #t3_start = Sys.time()
    d1 = 1/(lambda*omega_grp + 1/dnn)
    D1 = Matrix::Diagonal(n, d1)
    d2 = 1/dnn - d1/dnn^2
    nnmat_inv = lambda*D1 + lambda*D1%*%Dnninv%*%Rnq%*%solve(Rqq + crossprod(Rnq*sqrt(d2)),t(Rnq))%*%Dnninv%*%D1


    utilde = as.numeric(spam::rmvnorm.prec(1, rep(0, q), Rqq + crossprod(sqrt(d2)*Rnq)) )
    u = as.numeric(nnmat_inv%*%y0.5_Xbomega_grp + sqrt(lambda)*D1%*%Dnninv%*%Rnq%*%utilde + rnorm(n, 0, sd = sqrt(lambda*d1)))
    #t3 = t3 + as.numeric(difftime(Sys.time(), t3_start, units = "secs"))


    #### Step 4: update omega
    #t4_start = Sys.time()
    linpred = as.numeric(Xbeta + Z%*%u)
    omega = BayesLogit::rpg.devroye(N, h = 1, z = linpred)
    #t4 = t4 + as.numeric(difftime(Sys.time(), t4_start, units = "secs"))

    # update
    #omega_grp = Rfast::group(omega, id, method = "sum")
    omega_grp = as.numeric(rowsum(as.matrix(omega), group = id))
    ZtOmegaX = (Matrix::t(Z*omega)%*%X)
    XtOmegaX = (Matrix::t(X*omega)%*%X)

    # loglikelihood
    loglik = dbinom(y, size = 1, prob = 1/(1+exp(-linpred)), log = T)

    # save
    if((imcmc > nburn)&&((imcmc-nburn)%%nthin==0)){
      beta_save[isave,] = beta
      lambda_save[isave] = lambda
      phi_save[isave] = phi
      u_save[isave,] = u
      rho_save[isave] = rho
      loglik_save[isave,] = loglik
      isave = isave + 1
    }
  }
  t_end = Sys.time()
  t_mcmc = difftime(t_end, t_start, units = "secs")


  # waic
  out = list()
  if(is.null(colnames(X))){
    colnames(beta_save) = colnames(model.matrix(y~ -1+X))
  }else{
    colnames(beta_save) = colnames(X)
  }
  phi_save = matrix(phi_save); colnames(phi_save) = "phi"
  lambda_save = matrix(lambda_save); colnames(lambda_save) = "lambda"
  rho_save = matrix(rho_save); colnames(rho_save) = "rho"

  out$post_save = coda::mcmc(cbind(beta_save, phi_save, lambda_save, rho_save))
  betam_save = beta_save*as.numeric(phi_save)
  out$betam_save = coda::mcmc(betam_save)
  colnames(u_save) = unique(id)
  out$u_save = coda::mcmc(u_save)

  out$smoothness = smoothness
  out$acc_save = acc_save
  out$Ct = Ct
  out$loglik_save = loglik_save
  out$coords = coords
  out$nsave = nsave
  out$coords_knot = coords_knot
  out$q = q

  out$priors = priors

  out$t_mcmc = t_mcmc
  out$t_premcmc = t_premcmc

  return(out)
}


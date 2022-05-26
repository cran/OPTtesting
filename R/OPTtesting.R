#' @import rootSolve
#' @import quantreg
#' @import RSpectra
#' @import stats
NULL
#> NULL






Parameter_estimation_0 <- function(Z,B,p,proportion){
  Z_p = data.frame(Z = Z,B)
  number_of_Z = as.integer(p*proportion)
  order_Z = order(abs(Z_p$Z))[1:number_of_Z]
  Z_set = Z_p[order_Z,]
  B_1 = as.matrix(Z_set[,-1])
  Z_1 = Z_set$Z
  W_hat <- quantreg::rq(Z_1 ~ B_1, 0.5)$coef
  BW_hat = B %*% c(W_hat[-1])
  Hi_hat = Z - BW_hat
  return(list(Hi_hat = Hi_hat,BW_hat = BW_hat))}




Terms_for_pdf_2.1.1<- function(Sigma,eig,tol){
  # later we will consider how to choose dim_W
  p = ncol(Sigma)
  ev = eig
  dim_W = length(which(ev$values > tol))
  d_w = length(which(ev$values > 0))
  nc = min(500, d_w)
  B = ev$vectors[,1:dim_W] %*% diag(sqrt(ev$values[1:dim_W]),nrow = dim_W,ncol = dim_W)
  Br = ev$vectors[,(dim_W + 1):min(p,nc)] %*% diag(sqrt(ev$values[(dim_W + 1):nc]),nrow = (nc - dim_W ),ncol = (nc - dim_W))
  a_i_n2 = 1 - (B^2) %*% c(rep(1,dim_W))
  a_i_n2[which(a_i_n2 < 0 )] = 0
  table = list(B = B, a_i_n2 = a_i_n2, Br= Br, nc = nc)
  return(table)
}




Parameter_estimation_fund_1_c <- function(Hi_hat,a_i_n2,p,tau_sqr_1,tau_sqr_2,nu_0){
  Hi_hat = Hi_hat - nu_0
  a = sum(Hi_hat)/p
  b1 = sum(Hi_hat^2 - tau_sqr_1 - a_i_n2)/p
  c1 = sum(Hi_hat^3 - 3*(tau_sqr_1 + a_i_n2)*Hi_hat)/p
  d1 = sum(Hi_hat^4 - 3*(tau_sqr_1 + a_i_n2)^2 - 6*(tau_sqr_1 + a_i_n2)*(Hi_hat^2 - tau_sqr_1 - a_i_n2))/p
  b2 = sum(Hi_hat^2 - tau_sqr_2 - a_i_n2)/p
  c2 = sum(Hi_hat^3 - 3*(tau_sqr_2 + a_i_n2)*Hi_hat)/p
  d2 = sum(Hi_hat^4 - 3*(tau_sqr_2 + a_i_n2)^2 - 6*(tau_sqr_2 + a_i_n2)*(Hi_hat^2 - tau_sqr_2 - a_i_n2))/p
  return(list(a = a, b1 = b1, b2 = b2, c1 = c1, c2 = c2, d1 = d1, d2 = d2))
}




Parameter_estimation_fund_c <- function(coeff, tau_sqr_1,tau_sqr_2,l,u,nu_0){
  delta = tau_sqr_2 - tau_sqr_1
  f <- function(pi){
    pi_1 = - tau_sqr_2*pi/delta + (1 - (coeff$b1 + coeff$b2)/delta + (coeff$d1 - coeff$d2)/(3*delta^2))/2
    pi_2 = tau_sqr_1*pi/delta +  (1 + (coeff$b1 + coeff$b2)/delta - (coeff$d1 - coeff$d2)/(3*delta^2))/2
    a = coeff$a
    b =  coeff$b1 + pi*tau_sqr_1
    c = coeff$c1
    d = coeff$d1 - 3*pi*tau_sqr_1^2
    k = 1 - pi
    e1 = (k*(d*pi_1 - 3*delta^2*pi_1*pi_2 - a*c) - (3*delta*pi_1 + b + 2*delta*pi_2)*(b*pi_1 - a^2 - delta*pi_1*pi_2))*(k*(b*pi_1 + 2*delta*pi_1*pi_2 - a^2 - b*pi_2 + delta*pi_2^2) + 2*a^2*pi_2)
    e2 = (k*(c*pi_1 - a*b + a*delta*pi_2) - a*(b*pi_1 - a^2 - delta*pi_1*pi_2))*(k*(c*pi_1 - a*b  -  2*a*delta*pi_2 - c*pi_2) + 2*a*pi_2*(3*delta*pi_1 +b + 2*delta*pi_2))
    return( e1 - e2)
  }
  pi_0 = rootSolve::uniroot.all(f, lower = 0, upper =1)
  a = coeff$a
  b =  coeff$b1 + pi_0*tau_sqr_1
  c = coeff$c1
  d = coeff$d1 - 3*pi_0*tau_sqr_1^2
  k = 1 - pi_0
  pi_2 = tau_sqr_1*pi_0/delta +  (1 + (coeff$b1 + coeff$b2)/delta - (coeff$d1 - coeff$d2)/(3*delta^2))/2
  pi_1 = 1 - pi_2 - pi_0
  mu_2 = (k*(d*pi_1 - 3*delta^2*pi_1*pi_2 - a*c) - (3*delta*pi_1 + b + 2*delta*pi_2)*(b*pi_1 - a^2 - delta*pi_1*pi_2))/(k*(c*pi_1 - a*b  -  2*a*delta*pi_2 - c*pi_2) + 2*a*pi_2*(3*delta*pi_1 +b + 2*delta*pi_2))
  mu_1 = (d - mu_2*c - 3*delta*pi_2*mu_2^2 - 3*delta^2*pi_2)/(c - mu_2*b - 2*delta*pi_2*mu_2)
  #mu_1 = (c - mu_2 *b - 2*delta*pi_2*mu_2)/(b - mu_2*a - delta*pi_2)
  index = which(pi_1 < 1 & pi_2<1 & pi_1 >=0 & pi_2 >=0 & pi_0 <u &pi_0 >l)
  return(c(pi_0[index],pi_1[index],pi_2[index],mu_1[index] + nu_0,mu_2[index] + nu_0))
}



Parameter_estimation_fund_nu <- function(Hi_hat,a_i_n2,p,tau_sqr,nu_0,u,l){
  Hi_hat = Hi_hat - nu_0
  a = sum(Hi_hat)/p
  b = sum(Hi_hat^2 - tau_sqr - a_i_n2)/p
  c = sum(Hi_hat^3 - 3*(tau_sqr + a_i_n2)*Hi_hat)/p
  d = sum(Hi_hat^4 - 3*(tau_sqr + a_i_n2)^2 - 6*(tau_sqr + a_i_n2)*(Hi_hat^2 - tau_sqr - a_i_n2))/p
  c_3 = 2*tau_sqr^3
  c_2 = - 3*tau_sqr^3 - d*tau_sqr
  c_1 = 3*(a^2 - b)*tau_sqr^2 + (-3*b^2 + 2*a*c + d)*tau_sqr + c^2 - b*d
  c_0 = - b^3 + 2*a*b*c - a^2*d - c^2 + b*d
  roots = Re(polyroot(c(c_0,c_1,c_2,c_3)))
  pi_0 = roots[which(roots > l & roots < u)] #propotion of zeros should be high
  A = (b + tau_sqr*pi_0)^2 - a*c
  B = a*(d - 3*tau_sqr^2*pi_0) - c*(b + tau_sqr*pi_0)
  C = c^2 - (d - 3*tau_sqr^2*pi_0) *(b + tau_sqr*pi_0)
  in_0 = which(B^2 - 4*A*C >= 0)
  if(length(in_0) == 0){return(NULL)}
  pi_0 = pi_0[in_0]
  delta_sqrt = sqrt(B^2 - 4*A*C)
  u_1 = (- B - delta_sqrt)/(2*A)
  u_2 = (- B + delta_sqrt)/(2*A)
  pi_1 = (b + pi_0*tau_sqr - u_2*a)/(u_1*(u_1 - u_2))
  pi_2 = (b + pi_0*tau_sqr - u_1*a)/(u_2*(u_2 - u_1))
  index = which(pi_1 < 1 & pi_2<1 & pi_1 >=0 & pi_2 >=0 )
  return(c(pi_0[index],pi_1[index],pi_2[index],u_1[index] + nu_0,u_2[index] + nu_0))}



TVD<-function(x,f1,f2){
  n<-length(x)
  delta = (x[2:n]-x[1:(n-1)])
  value = (pmin(f1,f2)[2:n]+pmin(f1,f2)[1:(n-1)])/2
  inter<- min(delta%*%value,1)
  return(1-inter)
}


Terms_for_pdf_0.1.1 <- function(Sigma, eig,tol){
  ev = eig
  p =  length( which(ev$values >= tol))
  lambda_p = ev$values[p]
  dim_W = length( which(ev$values - lambda_p > 1e-6))
  B = ev$vectors[,1:dim_W] %*% diag(sqrt(ev$values[1:dim_W] - lambda_p),nrow = dim_W,ncol = dim_W)
  table = list(B = B ,dim_W = dim_W,lambda_p = lambda_p)
  return(table)
}



Terms_for_pdf_1 <- function(Z,pi_0,mu_0, tau_sqr,B,lambda_p,r_W){
  Z = matrix(Z)
  p = nrow(Z)
  if(p != nrow(B)){stop('dim does not match')}
  set_1 = c(rep(NULL,p))
  N_1 = c(rep(NULL,p))
  N_2 = c(rep(NULL,p))
  for (j in 1:p){
    N_1[j] = dnorm(Z[j], mean = B[j,] %*% r_W, sd = sqrt(lambda_p))
    # Zi|W ~ N(biW,lambda_p) when ui = 0 for i = 1,...,p
    N_2[j] = dnorm(Z[j], mean = mu_0 + B[j,] %*% r_W, sd = sqrt(lambda_p + tau_sqr))
    # Zi|W ~ N(mu_0 + biW,lambda_p + tao_sqr) when ui ~ N(mu_0,tao_sqr) for i = 1,...,p
    set_1[j] = pi_0*N_1[j] + (1 - pi_0)*N_2[j]
    # Zi|W ~ pi_0*N(biW,lambda_p) + (1-pi_0)*N(mu_0 + biW,lambda_p + tao_sqr)
    # since ui ~ pi_0*delta_0 + (1-pi_0)*N(mu_0,tao_sqr)
  }
  return(list(N_1 = N_1, set_1 = set_1))
}





Terms_for_pdf_less_0_nu_3 <- function(Z,pi_0,pi_1,pi_2,nu_0,mu_1,mu_2, tau_sqr_1,tau_sqr_2,B,lambda_p,r_W){
  Z = matrix(Z)
  p = nrow(Z)
  if(p != nrow(B)){stop('dim does not match')}
  set_1 = c(rep(NULL,p))
  set_2 = c(rep(NULL,p))
  N_0 = c(rep(NULL,p))
  N_1 = c(rep(NULL,p))
  N_2 = c(rep(NULL,p))
  N_11 = c(rep(NULL,p))
  N_21 = c(rep(NULL,p))
  sig_sqr_1 = (tau_sqr_1*lambda_p)/(lambda_p + tau_sqr_1)
  sig_sqr_2 = (tau_sqr_2*lambda_p)/(lambda_p + tau_sqr_2)
  c_1 = 1/(sqrt(2*pi*(lambda_p +tau_sqr_1)))
  c_2 = 1/(sqrt(2*pi*(lambda_p +tau_sqr_2)))
  if(nu_0 > 0){
    for (j in 1:p){
      BjW = B[j,] %*% r_W
      N_0[j] = dnorm(Z[j], mean = nu_0 + BjW, sd = sqrt(lambda_p))
      N_1[j] = dnorm(Z[j], mean = mu_1 + BjW, sd = sqrt(lambda_p +tau_sqr_1))
      N_2[j] = dnorm(Z[j], mean = mu_2 + BjW, sd = sqrt(lambda_p + tau_sqr_2))
      set_1[j] = pi_0*N_0[j] + pi_1*N_1[j] + pi_2*N_2[j]
      mu_0 = mu_1
      beta= (tau_sqr_1*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_1)
      N_11[j] =  c_1*pnorm(0, mean = beta, sd = sqrt(sig_sqr_1), lower.tail =  TRUE)*exp(beta^2/(2*sig_sqr_1)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_1))
      mu_0 = mu_2
      beta= (tau_sqr_2*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_2)
      N_21[j] =  c_2*pnorm(0, mean = beta, sd = sqrt(sig_sqr_2), lower.tail =  TRUE)*exp(beta^2/(2*sig_sqr_2)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_2))
      set_2[j] = pi_1*N_11[j] + pi_2* N_21[j]
    }}
  else{for (j in 1:p){
    BjW = B[j,] %*% r_W
    N_0[j] = dnorm(Z[j], mean = nu_0 + BjW, sd = sqrt(lambda_p))
    N_1[j] = dnorm(Z[j], mean = mu_1 + BjW, sd = sqrt(lambda_p +tau_sqr_1))
    N_2[j] = dnorm(Z[j], mean = mu_2 + BjW, sd = sqrt(lambda_p + tau_sqr_2))
    set_1[j] = pi_0*N_0[j] + pi_1*N_1[j] + pi_2*N_2[j]
    mu_0 = mu_1
    beta= (tau_sqr_1*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_1)
    N_11[j] =  c_1*pnorm(0, mean = beta, sd = sqrt(sig_sqr_1), lower.tail =  TRUE)*exp(beta^2/(2*sig_sqr_1)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_1))
    mu_0 = mu_2
    beta= (tau_sqr_2*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_2)
    N_21[j] =  c_2*pnorm(0, mean = beta, sd = sqrt(sig_sqr_2), lower.tail =  TRUE)*exp(beta^2/(2*sig_sqr_2)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_2))
    set_2[j] = pi_0*N_0[j] + pi_1*N_11[j] + pi_2*N_21[j]
  }}
  return(list( set_1 = set_1, set_2 = set_2))
}



Terms_for_pdf_greater_0_nu_3 <- function(Z,pi_0,pi_1,pi_2,nu_0, mu_1,mu_2, tau_sqr_1, tau_sqr_2 ,B,lambda_p,r_W){
  Z = matrix(Z)
  p = nrow(Z)
  if(p != nrow(B)){stop('dim does not match')}
  set_1 = c(rep(NULL,p))
  set_2 = c(rep(NULL,p))
  N_0 = c(rep(NULL,p))
  N_1 = c(rep(NULL,p))
  N_2 = c(rep(NULL,p))
  N_11 = c(rep(NULL,p))
  N_21 = c(rep(NULL,p))
  sig_sqr_1 = (tau_sqr_1*lambda_p)/(lambda_p + tau_sqr_1)
  sig_sqr_2 = (tau_sqr_2*lambda_p)/(lambda_p + tau_sqr_2)
  c_1 = 1/(sqrt(2*pi*(lambda_p +tau_sqr_1)))
  c_2 = 1/(sqrt(2*pi*(lambda_p +tau_sqr_2)))
  if(nu_0 < 0){
    for (j in 1:p){
      BjW = B[j,] %*% r_W
      N_0[j] = dnorm(Z[j], mean = nu_0 + BjW, sd = sqrt(lambda_p))
      N_1[j] = dnorm(Z[j], mean = mu_1 + BjW, sd = sqrt(lambda_p +tau_sqr_1))
      N_2[j] = dnorm(Z[j], mean = mu_2 + BjW, sd = sqrt(lambda_p + tau_sqr_2))
      set_1[j] = pi_0*N_0[j] + pi_1*N_1[j] + pi_2*N_2[j]
      mu_0 = mu_1
      beta= (tau_sqr_1*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_1)
      N_11[j] =  c_1*pnorm(0, mean = beta, sd = sqrt(sig_sqr_1), lower.tail =  FALSE)*exp(beta^2/(2*sig_sqr_1)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_1))
      mu_0 = mu_2
      beta= (tau_sqr_2*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_2)
      N_21[j] =  c_2*pnorm(0, mean = beta, sd = sqrt(sig_sqr_2), lower.tail =  FALSE)*exp(beta^2/(2*sig_sqr_2)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_2))
      set_2[j] = pi_1*N_11[j] + pi_2* N_21[j]
    }}
  else{for (j in 1:p){
    BjW = B[j,] %*% r_W
    N_0[j] = dnorm(Z[j], mean = nu_0 + BjW, sd = sqrt(lambda_p))
    N_1[j] = dnorm(Z[j], mean = mu_1 + BjW, sd = sqrt(lambda_p +tau_sqr_1))
    N_2[j] = dnorm(Z[j], mean = mu_2 + BjW, sd = sqrt(lambda_p + tau_sqr_2))
    set_1[j] = pi_0*N_0[j] + pi_1*N_1[j] + pi_2*N_2[j]
    mu_0 = mu_1
    beta= (tau_sqr_1*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_1)
    N_11[j] =  c_1*pnorm(0, mean = beta, sd = sqrt(sig_sqr_1), lower.tail =  FALSE)*exp(beta^2/(2*sig_sqr_1)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_1))
    mu_0 = mu_2
    beta= (tau_sqr_2*(Z[j] - BjW) + mu_0*lambda_p)/(lambda_p +tau_sqr_2)
    N_21[j] =  c_2*pnorm(0, mean = beta, sd = sqrt(sig_sqr_2), lower.tail =  FALSE)*exp(beta^2/(2*sig_sqr_2)-(1/2)*((Z[j] - BjW)^2/lambda_p + mu_0^2/tau_sqr_2))
    set_2[j] = pi_0*N_0[j] + pi_1*N_11[j] + pi_2* N_21[j]
  }}
  return(list( set_1 = set_1, set_2 = set_2))
}




posterior_prob_fund_nu_2 <- function(Z,Sigma,n,pi_0,pi_1,pi_2,nu_0,mu_1,mu_2,terms, tau_sqr_1,tau_sqr_2,p){
  B = terms$B
  dim_W = terms$dim_W
  lambda_p = terms$lambda_p
  log_fz_set = c(rep(NULL,n))
  log_fz_ui = matrix(c(rep(0,n*p)),nrow = p, ncol = n)
  for (j in 1:n){
    r_W = c(rnorm(dim_W,mean = 0, sd = 1))
    terms_1 = Terms_for_pdf_less_0_nu_3(Z,pi_0,pi_1,pi_2,nu_0,mu_1,mu_2, tau_sqr_1,tau_sqr_2,B,lambda_p,r_W)
    set_1 = terms_1$set_1
    set_2 = terms_1$set_2
    log_set_1  = log(set_1)
    log_set_2 = log(set_2)
    log_fz_set[j] = sum(log_set_1)
    log_fz_ui[,j] =  log_set_2 -log_set_1 + log_fz_set[j]
  }
  adj = max(log_fz_set)
  fz_sum = sum(exp(log_fz_set - adj))
  fz_ui_sum = rowSums(exp(log_fz_ui - adj))
  return(list(prob = fz_ui_sum/fz_sum))}



posterior_prob_fund_nu_3 <- function(Z,Sigma,n,pi_0,pi_1,pi_2,nu_0,mu_1,mu_2,terms, tau_sqr_1,tau_sqr_2,p){
  B = terms$B
  dim_W = terms$dim_W
  lambda_p = terms$lambda_p
  log_fz_set = c(rep(NULL,n))
  log_fz_ui = matrix(c(rep(0,n*p)),nrow = p, ncol = n)
  for (j in 1:n){
    r_W = c(rnorm(dim_W,mean = 0, sd = 1))
    terms_1 = Terms_for_pdf_greater_0_nu_3(Z,pi_0,pi_1,pi_2,nu_0,mu_1,mu_2, tau_sqr_1,tau_sqr_2,B,lambda_p,r_W)
    set_1 = terms_1$set_1
    set_2 = terms_1$set_2
    log_set_1  = log(set_1)
    log_set_2 = log(set_2)
    log_fz_set[j] = sum(log_set_1)
    log_fz_ui[,j] =  log_set_2 -log_set_1 + log_fz_set[j]
  }
  adj = max(log_fz_set)
  fz_sum = sum(exp(log_fz_set - adj))
  fz_ui_sum = rowSums(exp(log_fz_ui - adj))
  return(list(prob = fz_ui_sum/fz_sum))}


######################################################


#' AEB
#'
#' Estimate the parameters in the three-part mixture
#' @param Z a vector of test statistics
#' @param Sigma covariance matrix
#' @param eig eig value information
#' @param eig_tol the smallest eigen value used in the calulation
#' @param set_nu  a search region for nu_0
#' @param set1 a search region for tau_sqr_1
#' @param set2 a search region for tau_sqr_2
#' @param setp a search region for proportion
#' @examples
#' p = 500
#'  n_col = 10
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#' Z = rnorm(p,0,1) #this is just an example for testing the algorithm.
#' #Not true test statistics with respect to Sigma
#' best_set =  AEB(Z,Sigma)
#' print(c(best_set$pi_0, best_set$pi_1, best_set$pi_2, best_set$mu_1, best_set$mu_2))
#'\donttest{
#' library(MASS)
#' ######################################
#' #construct a test statistic vector Z
#' p = 1000
#' n_col = 4
#' pi_0 = 0.6
#' pi_1 = 0.2
#' pi_2 = 0.2
#' nu_0 = 0
#' mu_1 = -1.5
#' mu_2 = 1.5
#' tau_sqr_1 = 0.1
#' tau_sqr_2 = 0.1
#'
#'
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#'
#'
#'
#'
#' b = rmultinom(p, size = 1, prob = c(pi_0,pi_1,pi_2))
#' ui = b[1,]*nu_0 + b[2,]*rnorm(p, mean = mu_1,
#'      sd = sqrt(tau_sqr_1)) + b[3,]*rnorm(p, mean = mu_2,
#'       sd = sqrt(tau_sqr_2)) # actual situation
#' Z = mvrnorm(n = 1,ui, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#'
#' best_set =AEB(Z,Sigma)
#' print(c(best_set$pi_0, best_set$pi_1, best_set$pi_2, best_set$mu_1, best_set$mu_2))
#'}
#' @return The return of the function is a list
#' in the form of list(nu_0, tau_sqr_1, tau_sqr_2, pi_0, pi_1, pi_2, mu_1, mu_2, Z_hat).
#'
#' nu_0, tau_sqr_1, tau_sqr_2: The best combination of \eqn{(\nu_0, \tau_1^2, \tau_2^2)} that minimize the total variance from the regions \eqn{(D_{\nu_0}, D_{\tau_1^2}, D_{\tau_2^2})}.
#'
#' pi_0, pi_1, pi_2, mu_1, mu_2: The estimates of parameters \eqn{\pi_0, \pi_1, \pi_2, \mu_1, \mu_2}.
#'
#' Z_hat: A vector of simulated data  base on the parameter estimates.
#' @details Estimate the parameters in the three-part mixture \eqn{Z|\mu ~ N_p(\mu,\Sigma )}
#' where \eqn{\mu_i ~ \pi_0 \delta_ {\nu_0} + \pi_1  N(\mu_1, \tau_1^2)  + \pi_2  N(\mu_2, \tau_2^2), i = 1, ..., p}
#' @export
AEB <- function(Z,Sigma, eig = eigs_sym(Sigma,min(400,length(Z)), which = "LM"),eig_tol = 1,set_nu = c(0),set1 = c(0:10)*0.01 + 0.01,set2 = c(0:10)*0.01 + 0.01,setp =  c(1:7)*0.1){
  terms_1000 = Terms_for_pdf_2.1.1(Sigma,eig,eig_tol) #terms for estimating parameters
  B = terms_1000$B
  Br = terms_1000$Br
  a_i_n2 = terms_1000$a_i_n2
  dim_rW = terms_1000$nc - dim(terms_1000$B)[2]
  rW = c(rnorm(dim_rW,mean = 0, sd = 1))
  BW_res =  terms_1000$Br%*%rW
  vari_l = 1
  p = length(c(Z))
  if(min(eig$values) <0){stop("Sigma is not non-negative definite")}
  if(p != dim(Sigma)[1] | p != dim(Sigma)[1] ){stop("The dimensions of Z and Sigma are not match")}
  for(proportion in setp){
    parameter = Parameter_estimation_0(Z,terms_1000$B,p,proportion)
    Hi_hat = parameter$Hi_hat
    for(nu_0 in set_nu){
      for(tau_sqr_1 in set1){
        for(tau_sqr_2 in set2){
          if(tau_sqr_1 == tau_sqr_2 ){ Parameter_hat_set = Parameter_estimation_fund_nu(Hi_hat,a_i_n2,p,tau_sqr_1,nu_0,1,0)}else{
            coeff = Parameter_estimation_fund_1_c(Hi_hat,a_i_n2,p,tau_sqr_1, tau_sqr_2,nu_0)
            Parameter_hat_set = Parameter_estimation_fund_c(coeff, tau_sqr_1,tau_sqr_2,0,1,nu_0)
          }
          n_root = length(Parameter_hat_set)/5
          if(n_root != 0){
            for( i in 1:n_root){
              Parameter_hat = c(Parameter_hat_set[i], Parameter_hat_set[i + n_root],Parameter_hat_set[i + 2*n_root],Parameter_hat_set[i + 3*n_root],Parameter_hat_set[i + 4*n_root])
              if(abs(Parameter_hat[4]) <10){
                num = c( rmultinom(1, size = p, prob= c(Parameter_hat[1],Parameter_hat[2],Parameter_hat[3])))
                ui_set = c(rep(nu_0,num[1]),rnorm(num[2], mean = Parameter_hat[4], sd = sqrt(tau_sqr_1)),rnorm(num[3], mean = Parameter_hat[5], sd = sqrt(tau_sqr_2)))
                Z_gener = ui_set + BW_res
                Z_hat = Z_gener + parameter$BW_hat
                A <- density(Z, from = -5, to = 5, width = 1.5, kernel = "gaussian")
                B <- density(Z_hat, from = -5, to = 5, width = 1.5, kernel = "gaussian")
                vari = TVD(A$x,A$y,B$y)
                if(vari < vari_l){
                  vari_l = vari
                  best_set = c(proportion,tau_sqr_1,tau_sqr_2,nu_0,Parameter_hat)
                  Z_hat_best = Z_hat

                }}
            }

          }}}}}

  result = list(nu_0 = best_set[4],
                tau_sqr_1 = best_set[2] ,
                tau_sqr_2 = best_set[3],
                pi_0 = best_set[5],
                pi_1= best_set[6],
                pi_2 = best_set[7],
                mu_1 = best_set[8],
                mu_2 = best_set[9],
                Z_hat_best = Z_hat_best )

  return(result)}





#d-value calculation
#necessary part for calculate d-values

#' d_value
#'
#' Calculating the estimates for \eqn{P(\mu_i \le 0 | Z)}
#' @param Z a vector of test statistics
#' @param Sigma covariance matrix
#' @param eig eig value information
#' @param eig_value the smallest eigen value used in the calulation
#' @param sim_size simulation size
#' @param best_set a list of parameters (list(nu_0 = ..., tau_sqr_1 = ..., tau_sqr_2 = ..., pi_0 = ..., pi_1= ..., pi_2 = ..., mu_1 = ..., mu_2 = ...)) or returns from AEB
#' @return a vector of estimates for \eqn{P(\mu_i \le 0 | Z), i = 1, ..., p}
#' @examples
#' p = 500
#' n_col = 10
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#' Z = rnorm(p,0,1) #this is just an example for testing the algorithm.
#' #Not true test statistics with respect to Sigma
#'  d_value(Z,Sigma,sim_size = 5)
#' #To save time, we set the simulation size to be 10, but the default value is much better.
#'\donttest{
#' library(MASS)
#' ######################################
#' #construct a test statistic vector Z
#' p = 1000
#' n_col = 4
#' pi_0 = 0.6
#' pi_1 = 0.2
#' pi_2 = 0.2
#' nu_0 = 0
#' mu_1 = -1.5
#' mu_2 = 1.5
#' tau_sqr_1 = 0.1
#' tau_sqr_2 = 0.1
#'
#'
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#'
#'
#'
#'
#' b = rmultinom(p, size = 1, prob = c(pi_0,pi_1,pi_2))
#' ui = b[1,]*nu_0 + b[2,]*rnorm(p, mean = mu_1,
#'      sd = sqrt(tau_sqr_1)) + b[3,]*rnorm(p, mean = mu_2,
#'       sd = sqrt(tau_sqr_2)) # actual situation
#' Z = mvrnorm(n = 1,ui, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#'
#' d_value(Z,Sigma)
#'
#'}
#' @export
d_value <-function(Z,Sigma, best_set =  AEB(Z,Sigma), eig = eigs_sym(Sigma,min(400,length(Z)), which = "LM"), sim_size = 3000,eig_value = 0.35){
  p = length(c(Z))
  if(p != dim(Sigma)[1] | p != dim(Sigma)[1] ){stop("The dimensions of Z and Sigma are not match")}
  if(min(eig$values) <0){stop("Sigma is not non-negative definite")}
  terms = Terms_for_pdf_0.1.1(Sigma,eig,eig_value)
  prob_1 = posterior_prob_fund_nu_2(Z,Sigma,sim_size,best_set$pi_0,best_set$pi_1,best_set$pi_2,best_set$nu_0,best_set$mu_1,best_set$mu_2,terms,best_set$tau_sqr_1, best_set$tau_sqr_2,p)
  prob_1 = prob_1$prob
  return(prob_1)}




#' l_value
#'
#' Calculating the estimates for \eqn{P(\mu_i \ge 0 | Z)}
#' @param Z a vector of test statistics
#' @param Sigma covariance matrix
#' @param eig eig value information
#' @param sim_size simulation size
#' @param eig_value the smallest eigen value used in the calulation
#' @param best_set a list of parameters (list(nu_0 = ..., tau_sqr_1 = ..., tau_sqr_2 = ..., pi_0 = ..., pi_1= ..., pi_2 = ..., mu_1 = ..., mu_2 = ...)) or returns from Fund_parameter_estimation
#' @examples
#' p = 500
#'  n_col = 10
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#' Z = rnorm(p,0,1) #this is just an example for testing the algorithm.
#' #Not true test statistics with respect to Sigma
#'  l_value(Z,Sigma,sim_size = 5)
#'  #To save time, we set the simulation size to be 10, but the default value is much better.
#'\donttest{
#' library(MASS)
#' ######################################
#' #construct a test statistic vector Z
#' p = 1000
#' n_col = 4
#' pi_0 = 0.6
#' pi_1 = 0.2
#' pi_2 = 0.2
#' nu_0 = 0
#' mu_1 = -1.5
#' mu_2 = 1.5
#' tau_sqr_1 = 0.1
#' tau_sqr_2 = 0.1
#'
#'
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#'
#'
#'
#'
#' b = rmultinom(p, size = 1, prob = c(pi_0,pi_1,pi_2))
#' ui = b[1,]*nu_0 + b[2,]*rnorm(p, mean = mu_1,
#'      sd = sqrt(tau_sqr_1)) + b[3,]*rnorm(p, mean = mu_2,
#'       sd = sqrt(tau_sqr_2)) # actual situation
#' Z = mvrnorm(n = 1,ui, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#'
#' l_value(Z,Sigma)
#'
#'}
#' @return a vector of estimates for \eqn{P(\mu_i \ge 0 | Z)}
#' @export
l_value <-function(Z,Sigma, best_set =AEB(Z,Sigma), eig = eigs_sym(Sigma,min(400,length(Z)), which = "LM"), sim_size = 3000,eig_value = 0.35){
  p = length(c(Z))
  if(p != dim(Sigma)[1] | p != dim(Sigma)[1] ){stop("The dimensions of Z and Sigma are not match")}
  if(min(eig$values) <0){stop("Sigma is not non-negative definite")}
  terms = Terms_for_pdf_0.1.1(Sigma,eig,eig_value)
  prob_1 = posterior_prob_fund_nu_3(Z,Sigma,sim_size,best_set$pi_0,best_set$pi_1,best_set$pi_2,best_set$nu_0,best_set$mu_1,best_set$mu_2,terms,best_set$tau_sqr_1, best_set$tau_sqr_2,p)
  prob_1 = prob_1$prob
  return(prob_1)}


#' Optimal_procedure_3
#'
#' decision process
#' @param prob_0 d-values or l-values
#' @param alpha significance level
#' @examples
#' prob = runif(100,0,1) #assume this is the posterior probability vector
#' level = 0.3 #the significance level
#' Optimal_procedure_3(prob,level)
#' \donttest{
#' library(MASS)
#' ######################################
#' #construct a test statistic vector Z
#' p = 1000
#' n_col = 4
#' pi_0 = 0.6
#' pi_1 = 0.2
#' pi_2 = 0.2
#' nu_0 = 0
#' mu_1 = -1.5
#' mu_2 = 1.5
#' tau_sqr_1 = 0.1
#' tau_sqr_2 = 0.1
#'
#'
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#'
#'
#'
#'
#' b = rmultinom(p, size = 1, prob = c(pi_0,pi_1,pi_2))
#' ui = b[1,]*nu_0 + b[2,]*rnorm(p, mean = mu_1,
#'      sd = sqrt(tau_sqr_1)) + b[3,]*rnorm(p, mean = mu_2,
#'       sd = sqrt(tau_sqr_2)) # actual situation
#' Z = mvrnorm(n = 1,ui, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#'
#' prob_p = d_value(Z,Sigma)
#' #decision
#' level = 0.1 #significance level
#' decision_p = Optimal_procedure_3(prob_p,level)
#' decision_p$cj
#' head(decision_p$ai)
#'}
#' @return
#' ai: a vector of decisions. (1 indicates rejection)
#'
#' cj: The number of rejections
#'
#' FDR_hat: The estimated false discovery rate (FDR).
#'
#' FNR_hat: The estimated false non-discovery rate (FNR).
#' @export
Optimal_procedure_3 <- function(prob_0,alpha){
  p = length(prob_0)
  optimal_set = data.frame(i = c(1:p),prob_0 = prob_0)
  optimal_set = optimal_set[order(optimal_set$prob_0),]
  sum_prob_1 = sum(1 - prob_0)
  sum_prob_0 = 0
  for (cj in 1:p ){
    sum_prob_0 = sum_prob_0 + optimal_set$prob_0[cj]
    FDR_AVG = sum_prob_0/cj
    if(FDR_AVG > alpha){
      cj = cj - 1
      break}
    FDR_hat = FDR_AVG
  }
  a_set =  c(rep(0,p))
  if (cj != 0){
    index = optimal_set$i[1:cj]
    a_set[index] = 1}
  FNR_hat = 1 - sum( optimal_set$prob_0[(cj+1): p] )/(p - cj)
  if(cj == 0){FDR_hat  = 0}
  return(list(ai = a_set, cj = cj, FDR_hat = FDR_hat,FNR_hat = FNR_hat))
}


#' FDP_compute
#'
#' False discovery proportion and False non-discovery proportion computation
#' @param decision returns from the function Optimal_procedure_3
#' @param ui true mean vector
#' @param positive TRUE/FALSE valued. TRUE: H0: ui no greater than 0.  FALSE: H0: ui no less than 0.
#' @examples
#' ui = rnorm(10,0,1) #assume this is true parameter
#' decision = rbinom(10,1, 0.4) #assume this is decision vector
#' FDP_compute(decision,ui,TRUE)
#' \donttest{
#' library(MASS)
#' ######################################
#' #construct a test statistic vector Z
#' p = 1000
#' n_col = 4
#' pi_0 = 0.6
#' pi_1 = 0.2
#' pi_2 = 0.2
#' nu_0 = 0
#' mu_1 = -1.5
#' mu_2 = 1.5
#' tau_sqr_1 = 0.1
#' tau_sqr_2 = 0.1
#'
#'
#' A = matrix(rnorm(p*n_col,0,1), nrow = p, ncol = n_col, byrow = TRUE)
#' Sigma = A %*% t(A) +diag(p)
#' Sigma = cov2cor(Sigma) #covariance matrix
#'
#'
#'
#'
#' b = rmultinom(p, size = 1, prob = c(pi_0,pi_1,pi_2))
#' ui = b[1,]*nu_0 + b[2,]*rnorm(p, mean = mu_1,
#'      sd = sqrt(tau_sqr_1)) + b[3,]*rnorm(p, mean = mu_2,
#'       sd = sqrt(tau_sqr_2)) # actual situation
#' Z = mvrnorm(n = 1,ui, Sigma, tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
#'
#' prob_p = d_value(Z,Sigma)
#' #decision
#' level = 0.1 #significance level
#' decision_p = Optimal_procedure_3(prob_p,level)
#' FDP_compute(decision_p$ai,ui,TRUE)
#'}
#' @return False discovery proportion (FDP) and False non-discovery proportion (FNP)
#' @export
FDP_compute <- function(decision,ui,positive){
  if(positive){
    FDP = 1 - sum(which(decision == 1) %in% which(ui > 0))/length(which(decision == 1))
    FNP =  sum(which(decision == 0) %in% which(ui > 0))/length(which(decision == 0))}else{
      FDP = 1 - sum(which(decision == 1) %in% which(ui < 0))/length(which(decision == 1))
      FNP =  sum(which(decision == 0) %in% which(ui < 0))/length(which(decision == 0))}
  if(sum(decision == 1) ==0 ){FDP = 0}
  if(sum(decision == 0) ==0 ){FNP = 0}
  return(list(FDP = FDP,FNP = FNP))}


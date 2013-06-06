#
#  Code developed for the paper "Regularized Robust Portfolio Estimation"
#
#  Copyright 2012-2013 Theodoros Evgeniou, Massimiliano Pontil,
#                      Diomidis Spinellis, Rafal Swiderski, and
#                      Nick Nassuphis.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#


######################
## the SRL method ####
######################


SRL <- function(portfolio.dataset=NULL,signal.dataset,all.lambdas,mean.mom,use.mean){

  if(abs(mean.mom) != 1) stop("The third argument must be 1 for trend or -1 for mean-reversion.")

  n <- ncol(signal.dataset)
  tol <- 1.e-8

  s_t_1=make.lagged.x(signal.dataset,1)
  s_t=make.lagged.y(signal.dataset,1)

  if(!use.mean){
    s_t <- scale(s_t)
    s_t_1 <- scale(s_t_1)
    Gamma  = cov(s_t) + 0.000000001*diag(rep(1,ncol(s_t)))
    Theta  = cov(s_t,s_t_1)
  }
  else{
    s_t <- apply(s_t,2,function(column){column/sqrt(sum(column*column)/length(column))})
    s_t_1 <- apply(s_t_1,2,function(column){column/sqrt(sum(column*column)/length(column))})
    Gamma <- t(s_t) %*% s_t/nrow(s_t)
    Theta = t(s_t) %*% s_t_1/nrow(s_t)
  }
  Theta = mean.mom * Theta # take the minus theta when doing mean reversion
  S <- 0.5*(Theta + t(Theta))
  eigenS <- eigen(S,symmetric=T)
  P <- eigenS$vectors %*% apply(diag(eigenS$values),1:2,function(x){max(x,0)}) %*% t(eigenS$vectors)
  N <- P-S
  norm.N <- norm(N,"2")

  SRL.for.given.x0 <- function(eps,x0){

    eta <- function(x){as.double((t(x) %*% N %*% x + eps*(sum(abs(x))^2) - t(x) %*% P %*% x)/(t(x) %*% Gamma %*% x))}

    arg_of_prox <- function(x){
      if(eta(x) > 0) return((diag(n)+(1/norm.N)*(S+eta(x)*Gamma)) %*% x)
      else return((diag(n)+(1/norm(N-eta(x)*Gamma,"2"))*(S+eta(x)*Gamma)) %*% x)}

    prox <- function(z,last_eta){
      if(last_eta > 0) kappa = eps/(2*norm.N)
      else kappa = eps/(2*norm(N-last_eta*Gamma,"2"))

      arg_of_argmin <- function(x){0.5*sum((x-z)^2)+kappa*sum(abs(x))^2}
      sum_of_lambdas_minus_1 <- function(rho){sum(sapply(rho*abs(z)-2*kappa,function(w){max(w,0)}))-1}

      # Find the value of rho for which the lambdas sum up to 1
      z_nonzero <- which(z!=0)
      L <- 0
      R <- (1/length(z_nonzero)+2*kappa)/min(abs(z[z_nonzero]))
      M <- (L+R)/2
      while(abs(sum_of_lambdas_minus_1(M)) > tol ){
        if(sum_of_lambdas_minus_1(M) > 0){R=M; M=(L+R)/2}
        else{L=M;M=(L+R)/2}
      }
      the_correct_rho <- M

      lambda <- sapply(the_correct_rho*abs(z)-2*kappa,function(w){max(w,0)})
      lambda=lambda/sum(abs(lambda))
      return(lambda*z/(2*kappa+lambda))
    }

    # SPARSE RATIO LEARNER
    last_x <- x0
    for (k in 1:10000){
      x <- prox(arg_of_prox(last_x),eta(last_x))
      stopifnot(eta(last_x) - eta(x) >= -tol)
      if(abs((eta(last_x)-eta(x))) <= tol) break
      last_x <- x
    }
    return(x/sum(abs(x)))
  }

  eta <- function(x,eps){as.double((t(x) %*% N %*% x + eps*(sum(abs(x))^2) - t(x) %*% P %*% x)/(t(x) %*% Gamma %*% x))}

  ##################################################
  ####### continuation method left to right ########
  ##################################################

  # compute the non-sparse portfolio by solving the generalised eigenvalue problem
  triple.product <- solve(msqr(Gamma)) %*% S %*% solve(msqr(Gamma))
  max.EV.index <- which.max(eigen(triple.product)$values)
  nonsparse.max.autocor.portf <- solve(msqr(Gamma)) %*% eigen(triple.product)$vectors[,max.EV.index]
  nonsparse.max.autocor.portf <- nonsparse.max.autocor.portf/sum(abs(nonsparse.max.autocor.portf))

  x.matrix.L2R = matrix(nrow=n,ncol=length(all.lambdas))
  last.x = nonsparse.max.autocor.portf
  for(k in 1:length(all.lambdas)){
    x <- SRL.for.given.x0(all.lambdas[k],last.x)
    x.matrix.L2R[,k] <- x/sum(abs(x))
    last.x = x
  }

  ##################################################
  ####### continuation method right to left ########
  ##################################################


  highest.autocor.stock = which.max(diag(Theta)[1:n]/diag(Gamma)[1:n])
  highest.autocor.stock.portf <- rep(0,n)
  highest.autocor.stock.portf[highest.autocor.stock] <- 1

  x.matrix.R2L = matrix(nrow=n,ncol=length(all.lambdas))
  last.x = highest.autocor.stock.portf
  for(k in length(all.lambdas):1){
    x <- SRL.for.given.x0(all.lambdas[k],last.x)
    x.matrix.R2L[,k] <- x/sum(abs(x))
    last.x = x
  }

  #############################################
  #### choose the better of L2R and R2L #######
  #############################################

  x.matrix.L2R.R2L = matrix(nrow=n,ncol=length(all.lambdas))
  for(k in 1:length(all.lambdas)){
    if(eta(x.matrix.L2R[,k],all.lambdas[k]) < eta(x.matrix.R2L[,k],all.lambdas[k])){
      x.matrix.L2R.R2L[,k] = x.matrix.L2R[,k]
    }
    else{
      x.matrix.L2R.R2L[,k] = x.matrix.R2L[,k]
    }
  }
  #####################################################################
  ##### scale the portfolio by the sd or square root of 2nd moment ####
  #####################################################################

  s_t=make.lagged.y(signal.dataset,1)

  if(use.mean){
    Gamma  = cov(s_t) + 0.000000001*diag(rep(1,ncol(signal.dataset)))
  }
  else{
    # Massi's method
    Gamma <- t(s_t) %*% s_t/nrow(s_t)

    # Rafal's method
    #Gamma  = cov(s_t) + apply(s_t,2,mean) %*% t(apply(s_t,2,mean)) + 0.000000001*diag(rep(1,ncol(s_t)))
  }

  x.matrix.L2R.R2L <- solve(diag(sqrt(diag(Gamma)))) %*% x.matrix.L2R.R2L


  portfolio <- x.matrix.L2R.R2L
  portfolio=portfolio*(abs(portfolio)>0.00000001)
  portfolio = apply(portfolio,2,function(column){column/sum(abs(column))})
  rownames(portfolio) <- colnames(signal.dataset)

  list(portfolios=portfolio,signals=portfolio)
  # the convention is that we consider the sign of the signal in PostProcess
}

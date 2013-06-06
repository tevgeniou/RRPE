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


##################################
#### robust CCA method ###########
##################################

robust.CCA <- function(portfolio.dataset,signal.dataset,all.lambdas,mean.mom=NULL,use.mean,q){

  if (use.mean)
    use.mean=TRUE
  if (!use.mean)
    use.mean=FALSE

  if(nrow(portfolio.dataset) != nrow(signal.dataset)) stop("portfolio.dataset and signal.dataset must have the same number of rows.")

  signal.dataset=make.lagged.x(signal.dataset,1)
  portfolio.dataset=make.lagged.y(portfolio.dataset,1)

  if(q!=1 & q!=2) stop("q must be either 1 or 2.")
  if(class(use.mean)!="logical") stop("use.mean must be TRUE or FALSE.")
  m <- ncol(portfolio.dataset)
  n <- ncol(signal.dataset)
  tol <- 1.e-8
  s_t = portfolio.dataset
  r_t = signal.dataset

  #######################
  ## scale the data #####
  #######################

  if(use.mean == TRUE){
    r_t = scale(r_t)
    s_t = scale(s_t)
  }else{
    r_t <- apply(r_t,2,function(column){column/sqrt(sum(column*column)/(length(column)-1))})
    s_t <- apply(s_t,2,function(column){column/sqrt(sum(column*column)/(length(column)-1))})
  }

  ###########################################
  ##### define Theta, Gamma and Sigma #######
  ###########################################
  Theta = t(r_t)%*%s_t/(nrow(r_t)-1)
  if(use.mean == TRUE){
    Gamma  = cov(r_t) + 0.000000001*diag(rep(1,ncol(r_t)))
    Sigma  = cov(s_t) + 0.000000001*diag(rep(1,ncol(s_t)))
  } else{
    Gamma = t(r_t) %*% r_t/(nrow(r_t)-1)
    Sigma = t(s_t) %*% s_t/(nrow(s_t)-1)
  }

  ##############################
  ##### define A and B #########
  ##############################

  A = matrix(0,nrow=nrow(Theta)+ncol(Theta),ncol=nrow(Theta)+ncol(Theta))
  A[(nrow(Theta)+1):(nrow(Theta)+ncol(Theta)),1:nrow(Theta)] = t(Theta)
  A[1:nrow(Theta),(nrow(Theta)+1):(nrow(Theta)+ncol(Theta))] = Theta

  B = matrix(0,nrow=nrow(Gamma)+nrow(Sigma),ncol=ncol(Gamma)+ncol(Sigma))
  B[1:nrow(Gamma),1:ncol(Gamma)] = Gamma
  B[(nrow(Gamma)+1):(nrow(Gamma)+nrow(Sigma)),(ncol(Gamma)+1):(ncol(Gamma)+ncol(Sigma))] = Sigma

  ##############################
  ##### define P and N #########
  ##############################

  eigenA <- eigen(-A,symmetric=T)
  P <- eigenA$vectors %*% apply(diag(eigenA$values),1:2,function(x){max(x,0)}) %*% t(eigenA$vectors)
  N <- P+A
  #######################
  ### define eta ########
  #######################
  eta1 = function(z,eps){
    z.x = z[1:nrow(Gamma)]
    z.y = z[(nrow(Gamma)+1):length(z)]
    as.double((-t(z.x)%*%Theta%*%z.y+eps*L1.norm(z.x)*L1.norm(z.y))/(sqrt(t(z.x)%*%Gamma%*%z.x)*sqrt(t(z.y)%*%Sigma%*%z.y)))}
  eta2 = function(z,eps){
    z.x = z[1:nrow(Gamma)]
    z.y = z[(nrow(Gamma)+1):length(z)]
    as.double((-t(z.x)%*%Theta%*%z.y+eps*L2.norm(z.x)*L2.norm(z.y))/(sqrt(t(z.x)%*%Gamma%*%z.x)*sqrt(t(z.y)%*%Sigma%*%z.y)))}
  if(q==1){eta=eta1}else{eta=eta2}
  ##############################################
  ### robust CCA for given starting point ######
  ##############################################
  robust.CCA.for.given.z0 <- function(eps,z0){
    arg_of_prox <- function(z){
      #browser()
      z.x = z[1:n]
      z.y = z[(n+1):length(z)]
      if(q==1){u = c(2*L1.norm(z.x)*sign(z.x),2*L1.norm(z.y)*sign(z.y))}else{u=2*z}
      chain.rule.prefactor = 2*(sqrt(t(z.x)%*%Gamma%*%z.x)+sqrt(t(z.y)%*%Sigma%*%z.y))
      v = chain.rule.prefactor*c(Gamma%*%z.x/sqrt(as.double(t(z.x)%*%Gamma%*%z.x)),Sigma%*%z.y/sqrt(as.double(t(z.y)%*%Sigma%*%z.y)))
      return(z-((-A+abs(eta(z,eps))*B)%*%z-eps*u/2-(eta(z,eps)>0)*eta(z,eps)*v/2)/norm(P+abs(eta(z,eps))*B,"2"))}
    prox1 <- function(w,last_eta){
      #browser()
      r=eps*0.5/norm(P+abs(last_eta)*B,"2")
      arg_of_argmin <- function(u){0.5*sum((u-w)^2)+r*sum(abs(u))^2}
      sum_of_lambdas_minus_1 <- function(rho){sum(sapply(rho*abs(w)-2*r,function(w){max(w,0)}))-1}
      w_nonzero <- which(w!=0)
      L <- 0
      R <- (1/length(w_nonzero)+2*r)/min(abs(w[w_nonzero]))
      M <- (L+R)/2
      while(abs(sum_of_lambdas_minus_1(M)) > tol ){
        if(sum_of_lambdas_minus_1(M) > 0){R=M; M=(L+R)/2}
        else{L=M;M=(L+R)/2}
      }
      the_correct_rho <- M
      lambda <- sapply(the_correct_rho*abs(w)-2*r,function(w){max(w,0)})
      lambda=lambda/sum(abs(lambda))
      return(lambda*w/(2*r+lambda))
    }
    prox2 <- function(w,last_eta){
      #browser()
      r=eps*0.5/norm(P+abs(last_eta)*B,"2")
      arg_of_argmin <- function(u){0.5*sum((u-w)^2)+r*(L2.norm(u[1:n])+L2.norm(u[(n+1):(n+m)]))^2}
      sum_of_lambdas_minus_1 <- function(rho){max(rho*L2.norm(w[1:n])-2*r,0)+max(rho*L2.norm(w[(n+1):(n+m)])-2*r,0)-1}
      L <- 0
      R <- (0.5+2*r)/min(L2.norm(w[1:n]),L2.norm(w[(n+1):(n+m)]))
      M <- (L+R)/2
      while(abs(sum_of_lambdas_minus_1(M)) > tol ){
        if(sum_of_lambdas_minus_1(M) > 0){R=M; M=(L+R)/2}
        else{L=M;M=(L+R)/2}
      }
      the_correct_rho <- M
      lambda1 <- max(the_correct_rho*L2.norm(w[1:n])-2*r,0)
      lambda2 <- max(the_correct_rho*L2.norm(w[(n+1):(n+m)])-2*r,0)
      return(c(lambda1*w[1:n]/(2*r+lambda1),lambda2*w[(n+1):(n+m)]/(2*r+lambda2)))
    }

    if(q==1){prox=prox1}else{prox=prox2}
    eta.plot = vector(length=5000)
    last_z <- as.vector(z0)
    last_z.x = last_z[1:n]
    last_z.x = last_z.x/sqrt(t(last_z.x)%*%Gamma%*%last_z.x)
    last_z.y = last_z[(n+1):length(last_z)]
    last_z.y = last_z.y/sqrt(t(last_z.y)%*%Sigma%*%last_z.y)
    last_z = c(last_z.x,last_z.y)
    for (k in 1:5000){
      #browser()
      z <- prox(arg_of_prox(last_z),eta(last_z,eps))
      #stopifnot(eta(last_z,eps) - eta(z,eps) >= -tol)
      z.x = z[1:n]
      z.x = z.x/sqrt(t(z.x)%*%Gamma%*%z.x)
      z.y = z[(n+1):length(z)]
      z.y = z.y/sqrt(t(z.y)%*%Sigma%*%z.y)
      z = c(z.x,z.y)
      eta.plot[k] = eta(z,eps)
      if(eta(last_z,eps) - eta(z,eps) < -tol) cat("eta just increased by ",abs(eta(last_z,eps) - eta(z,eps)),"    eta:",eta(z,eps),"    epsilon:",eps,"    iteration:",k,"\n")
      if(abs(eta(last_z,eps)-eta(z,eps)) <= tol){
        eta.plot=eta.plot[1:k];
        #plot(eta.plot,main=paste("eta"),xlab="",ylab="",type="l");
        break}
      last_z <- z
      #if(k==5000){plot(eta.plot,main=paste("eta"),type="l",xlab="",ylab="")}
    }
    z[1:n] = z[1:n]/L1.norm(z[1:n])
    z[(n+1):length(z)] = z[(n+1):length(z)]/L1.norm(z[(n+1):length(z)])
    return(z)
  }

  ##################################################
  ####### continuation method left to right ########
  ##################################################

  sig_t = msqr(solve(B))
  A.tilde = sig_t %*% A %*% sig_t
  max.EV.index = which.max(eigen(A.tilde)$values)
  z.nonsparse = sig_t %*% eigen(A.tilde)$vectors[,max.EV.index]
  z.nonsparse[1:n] = z.nonsparse[1:n]/L1.norm(z.nonsparse[1:n])
  z.nonsparse[(n+1):(n+m)] = z.nonsparse[(n+1):(n+m)]/L1.norm(z.nonsparse[(n+1):(n+m)])

  z.matrix.L2R = matrix(NA,nrow=n+m,ncol=length(all.lambdas))
  last.z = z.nonsparse
  for(j in 1:length(all.lambdas)){
    #cat(j," ")
    z <- robust.CCA.for.given.z0(all.lambdas[j],last.z)
    z.matrix.L2R[,j] = z
    last.z = z
  }

  ##################################################
  ####### continuation method right to left ########
  ##################################################

  if(q==1){
    most.autocor.pair = arrayInd(which.max(Theta),dim(Theta))
    x.one.sparse = rep(0,n)
    x.one.sparse[most.autocor.pair[1]]=1
    y.one.sparse = rep(0,m)
    y.one.sparse[most.autocor.pair[2]]=1
    z.sparse = c(x.one.sparse,y.one.sparse)

    z.matrix.R2L = matrix(NA,nrow=n+m,ncol=length(all.lambdas))
    last.z = z.sparse
    for(j in length(all.lambdas):1){
      z <- robust.CCA.for.given.z0(all.lambdas[j],last.z)
      z.matrix.R2L[,j] = z
      last.z = z
    }
  }else{
    sig_t = msqr(solve(Gamma))
    Gamma.inverse = solve(Gamma)
    min.EV.index = which.min(eigen(Gamma.inverse)$values)
    x.Gamma = sig_t %*% eigen(Gamma.inverse)$vectors[,min.EV.index]
    x.Gamma = x.Gamma/L1.norm(x.Gamma)

    sig_t = msqr(solve(Sigma))
    Sigma.inverse = solve(Sigma)
    min.EV.index = which.min(eigen(Sigma.inverse)$values)
    y.Sigma = sig_t %*% eigen(Sigma.inverse)$vectors[,min.EV.index]
    y.Sigma = y.Sigma/L1.norm(y.Sigma)

    if(eta(c(x.Gamma,y.Sigma),0) > 0) y.Sigma = -y.Sigma

    z.matrix.R2L = matrix(NA,nrow=n+m,ncol=length(all.lambdas))
    last.z = c(x.Gamma,y.Sigma)
    for(j in length(all.lambdas):1){
      #cat(j," ")
      z <- robust.CCA.for.given.z0(all.lambdas[j],last.z)
      z.matrix.R2L[,j] = z
      last.z = z
    }
  }

  ##############################################
  ###### select the better of L2R and R2L ######
  ##############################################

  z.matrix = matrix(NA,nrow=n+m,ncol=length(all.lambdas))
  for(k in 1:length(all.lambdas)){
    if(eta(z.matrix.L2R[,k],all.lambdas[k]) < eta(z.matrix.R2L[,k],all.lambdas[k])){
      z.matrix[,k] = z.matrix.L2R[,k]
    }
    else{
      z.matrix[,k] = z.matrix.R2L[,k]
    }
  }
  #####################################################################
  ##### scale the portfolio by the sd or square root of 2nd moment ####
  #####################################################################

  s_t=portfolio.dataset
  r_t=signal.dataset

  if(use.mean == TRUE){
    Gamma  = cov(r_t) + 0.000000001*diag(rep(1,ncol(r_t)))
    Sigma  = cov(s_t) + 0.000000001*diag(rep(1,ncol(s_t)))
  }
  else{
    Gamma = t(r_t) %*% r_t/(nrow(r_t)-1)
    Sigma = t(s_t) %*% s_t/(nrow(s_t)-1)
  }

  x.matrix <- z.matrix[1:n,]
  y.matrix <- z.matrix[(n+1):(n+m),]
  x.matrix <- solve(diag(sqrt(diag(Gamma)))) %*% x.matrix
  y.matrix <- solve(diag(sqrt(diag(Sigma)))) %*% y.matrix

  portfolio.x = x.matrix
  portfolio.y = y.matrix
  portfolio.x=portfolio.x*(abs(portfolio.x)>0.00000001)
  portfolio.y=portfolio.y*(abs(portfolio.y)>0.00000001)
  portfolio.x = apply(portfolio.x,2,function(column){column/L1.norm(column)})
  portfolio.y = apply(portfolio.y,2,function(column){column/L1.norm(column)})
  rownames(portfolio.x) <- colnames(signal.dataset)
  rownames(portfolio.y) <- colnames(portfolio.dataset)

  list(portfolios=portfolio.y,signals=portfolio.x)
} # end of robust CCA method


CCA1 <- function(portfolio.dataset,signal.dataset,all.lambdas,mean.mom,use.mean)
  robust.CCA(portfolio.dataset,signal.dataset,all.lambdas,mean.mom=NULL,use.mean,1)

CCA2 <- function(portfolio.dataset,signal.dataset,all.lambdas,mean.mom,use.mean)
  robust.CCA(portfolio.dataset,signal.dataset,all.lambdas,mean.mom=NULL,use.mean,2)

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



###################################################################################################################
### The L2 learner
###################################################################################################################

# This returns a list of two matrices, one for the portfolios and one for the signals, each matrix having one column per lambda

# this is for 1 lambda
max.lag1.portfolio <- function(portfolio.dataset=NULL,signal.dataset,lambda,mean.mom,use.mean){

  ################## NOTE ##############
  mean.mom=-mean.mom # to make this like in SRL
  ######################################

  s_t_1=make.lagged.x(signal.dataset,1)
  s_t=make.lagged.y(signal.dataset,1)

  if (!use.mean){
    Theta= (nrow(signal.dataset)-1)*t(s_t_1)%*%s_t # just to scale lambda to a reasonable range
    Gamma  = t(s_t)%*%s_t/(nrow(signal.dataset)-1)+small_number*diag(rep(1,ncol(s_t)))
  } else {
    Theta= (nrow(signal.dataset)-1)*cov(s_t_1,s_t) # just to scale lambda to a reasonable range
    Gamma  = cov(s_t)+small_number*diag(rep(1,ncol(s_t)))
  }

  Theta=mean.mom*0.5*(Theta+t(Theta)) # take the minus theta when doing mean reversion
  Theta=Theta+lambda*diag(rep(1,ncol(Theta)))
  sig_t=msqr(solve(Gamma))
  sandwich=sig_t%*%Theta%*%sig_t

  ptfraw=eigen(sandwich, symmetric=TRUE)$vectors
  ptfraw=sig_t%*%ptfraw

  if (lambda == large_number){
    ptfraw=eigen(Gamma)$vectors
    ptfraw=ptfraw[,ncol(ptfraw):1]
  }
  ptf = apply(ptfraw,2,function(x){ tot<-sum(abs(x)); if(tot>0) {x/tot} else x })

  rownames(ptf)<-colnames(s_t)
  colnames(ptf)<-paste("e",1:ncol(ptf),sep=".")
  return(ptf)
}

# now for many lambdas (which is what we call to be consistent with SRL)
# assumes the last lambda is "infinity" and the lambdas as sorted from lowest to highest value
max.lag1.portfolio.manylambdas <- function(portfolio.dataset=NULL,signal.dataset,all.lambdas,mean.mom,use.mean){

  s_t_1=make.lagged.x(signal.dataset,1)
  s_t=make.lagged.y(signal.dataset,1)

  ####### we need this to calculate nu fast
  if (!use.mean){
    Theta= (nrow(signal.dataset)-1)*t(s_t_1)%*%s_t # just to scale lambda to a reasonable range
    Gamma  = t(s_t)%*%s_t/(nrow(signal.dataset)-1)+small_number*diag(rep(1,ncol(s_t)))
  } else {
    Theta= (nrow(signal.dataset)-1)*cov(s_t_1,s_t) # just to scale lambda to a reasonable range
    Gamma  = cov(s_t)+small_number*diag(rep(1,ncol(s_t)))
  }
  Theta=0.5*(Theta+t(Theta))
  ##########

  all.portf = matrix(nrow=ncol(signal.dataset),ncol=length(all.lambdas))
  lambda.index=length(all.lambdas)
  all.portf[,lambda.index]=max.lag1.portfolio(portfolio.dataset,signal.dataset,all.lambdas[lambda.index],mean.mom,use.mean)[,ncol(signal.dataset)]
  for(lambda.index in (length(all.lambdas)-1):1){
    train_pnl_old=(signal.dataset%*%all.portf[,lambda.index+1,drop=F])
    portf.sd = sd(train_pnl_old)
    P = train_pnl_old/portf.sd
    dollars=(1-exp(-P))/(1+exp(-P))
    train_pnl_old=-mean.mom*head(dollars,-1)*tail(train_pnl_old,-1)

    x=all.portf[,lambda.index+1,drop=F]
    nu_old=as.double((t(x) %*% Theta %*% x)/(t(x) %*% Gamma %*% x))

    all_eigens = max.lag1.portfolio(portfolio.dataset,signal.dataset,all.lambdas[lambda.index],mean.mom,use.mean)
    all_eigens=apply(all_eigens, 2, function(r){
      res=r
      if (sum(r*all.portf[,length(all.lambdas),drop=F])<0)
        res=-r
      res})
    similar.eigen=which.min(apply(all_eigens, 2, function(r){
      train_pnl_new=(signal.dataset%*%matrix(r,ncol=1))
      portf.sd = sd(train_pnl_new)
      P = train_pnl_new/portf.sd
      dollars=(1-exp(-P))/(1+exp(-P))
      train_pnl_new=-mean.mom*head(dollars,-1)*tail(train_pnl_new,-1)

      x=matrix(r,ncol=1)
      nu_new=as.double((t(x) %*% Theta %*% x)/(t(x) %*% Gamma %*% x))

      #Linfinity.norm(nu_new-nu_old)
      Linfinity.norm(train_pnl_new-train_pnl_old)
      }))

    all.portf[,lambda.index]=all_eigens[,similar.eigen]
  }
  list(portfolios=all.portf,signals=all.portf)
}

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



large_number=10^10; small_number=0.00001

################################################################################################
# CREATE THE SECTORS NOW
################################################################################################
# Sectors from http://www.nasdaq.com/screening/industries.aspx (sorted by market Capitalization, and using only those with market
# cap more than 10billion in April 2013)
financials=intersect(colnames(dataset),c("HBC","JPM","WFC","BAC","C","RY","LFC","WBK","MTU","ITUB","TD","SAN","GS","BNS","AXP","BBD","UBS","USB","BCS","LYG","AIG","SMFG","MFG","BBVA","MS","DB","BLK","BMO","MET","PUK","CS","ING","PNC","COF","CM","BK","RBS","TRV","HDB","BEN","ACE","BSBR","MFC","PRU","STT","IBN","AFL","CB","ALL","SCHW","BBT","NMR","MMC","DFS","CME","TROW","AON","SLF","L","AV","STI","FITB","BCH","PGR","AMP","MTB","KB","NTRS","IX","IVZ","AEG","BAP","RF","ICE","AMTD","CTRX","HIG"))
energy=intersect(colnames(dataset),c("XOM","PTR","GE","CVX","BBL","BP","TOT","PBR","SLB","SNP","CEO","SI","E","STO","OXY","COP","SU","APC","EMR","PSX","HAL","IMO","EOG","ENB","APA","CNQ","SSL","NOV","PHG","MPC","VLO","CVE","DVN","MRO","HES","CMI","NBL","RIG","BHI","PC","PAA","CAM","PXD","CLR","ESV","ECA","NXY","CHK","TLM","SWN","RRC","MUR","FTI","COG","HFC","MMP","DO","CIE","WFT","CXO","NE"))
#"XOM" may need to be removed as it is the largest company in the world and behaves special
healthcare=intersect(colnames(dataset),c("MMM","AMGN","CVS","GILD","UNH","MDT","ESRX","WAG","BIIB","BAX","COV","SYK","MCK","ISRG","LUX","WLP","CI","BDX","AET","HCA","CAH","ZMH","STJ","HUM","DVA","ABC","LIFE","SNN","BSX","FMS","EW"))
tech=intersect(colnames(dataset),c("AAPL","GOOG","MSFT","IBM","ORCL","QCOM","CSCO","INTC","TSM","SAP","DCM","EMC","FB","TXN","VMW","HPQ","INFY","ADP","ITW","ETN","BIDU","CRM","YHOO","DELL","CTSH","WIT","ADBE","INTU","BRCM","MSI","AMAT","LNKD","KYO","SYMC","CERN","NOK","ATVI","OMC","ADI","CTXS","NTAP","DOV","STX","SNDK","WDC","CA","ALTR","JNPR","FISV","CHKP","PNR","RHT","TDC","RAX"))
################################################################################################

###################################################################################################################
## Some Helper Functions
###################################################################################################################

sharpe<-function(x)ifelse(sum(abs(x))>0,16*mean(x)/sd(x),0)

msqr <-cmpfun(function(x){
  x=x+0.000000001*diag(rep(1,ncol(x)))
  e=eigen(x,symmetric=TRUE)
  v=e$vectors
  d=if(length(e$values)==1) {as.matrix(sqrt(e$values))} else {sqrt(diag(abs(e$values)))}
  mr=v %*% d %*% t(v)
  mr
})

make.lagged.x<-cmpfun(function(s,lags){
  p<-nrow(s)-lags
  Reduce(cbind,lapply(1:lags,function(i)tail(head(s,-i),p)))
})

make.lagged.y<-cmpfun(function(s,lags) tail(s,-lags))

L2.norm <- function(x){sqrt(sum(x*x))}
L1.norm <- function(x){sum(abs(x))}
Linfinity.norm <- function(x){max(abs(x))}

nu.function <- function(portfolio,signal,portfolio.dataset,signal.dataset,use.cca,use.mean){
  s_t_1=make.lagged.x(signal.dataset,1)
  s_t=make.lagged.y(portfolio.dataset,1)
  if (!use.mean){
    Theta= t(s_t_1)%*%s_t/(nrow(portfolio.dataset)-1)
    Gamma  = t(s_t)%*%s_t/(nrow(portfolio.dataset)-1)+small_number*diag(rep(1,ncol(s_t)))
    Sigma = t(s_t_1) %*% s_t_1/(nrow(s_t_1)-1)
  }else {
    Theta= cov(s_t_1,s_t)
    Gamma  = cov(s_t)+small_number*diag(rep(1,ncol(s_t)))
    Sigma  = cov(s_t_1)+small_number*diag(rep(1,ncol(s_t_1)))
  }
  if (!use.cca){
    res=as.double((t(signal) %*% Theta %*% portfolio)/(t(signal) %*% Gamma %*% portfolio))
  } else {
    res =  as.double((t(signal)%*%Theta%*%portfolio)/(sqrt(t(signal)%*%Sigma%*%signal)*sqrt(t(portfolio)%*%Gamma%*%portfolio)))
  }
  res
}

###########################################################################################
# Post Processing function
###########################################################################################

get_stats<-function(portfolios,signals,portfolio.dataset,signal.dataset,get.nu,use.cca,use.mean){
  performance_iter=Reduce(cbind,lapply(1:ncol(signals), function(i) {
    if (get.nu){
      res=nu.function(matrix(portfolios[,i],ncol=1),matrix(signals[,i],ncol=1),portfolio.dataset,signal.dataset,use.cca,use.mean)
    } else {
      test_ptf=portfolio.dataset%*%matrix(portfolios[,i],ncol=1)
      test_sign=signal.dataset%*%matrix(signals[,i],ncol=1)
      res=sign(head(test_sign,-1))*tail(test_ptf,-1)
    }
  }))
  performance_iter
}


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


rm(list=ls())
library(compiler); library(parallel)

load("NIPS2013.data")
source("Helpers.R")

#################################################################################
# THE CHOICES TO MAKE
Run.L2=1; Run.SRL=0; Run.CCA1=0; Run.CCA2=0 # choose which method to run
train=1000; validation=250; test=nrow(dataset)-train-validation # this should be the same as in RunExperiments.R
#################################################################################


savefile.ini="Results"
if (Run.L2){
  savefile.ini=paste(savefile.ini,"L2",sep="_")
  epsilon=sort(0.1^seq(-6,6,by=0.005))
  epsilon=c(0,epsilon,large_number) # just make sure lambda="infinity" is at the very end
}
if (Run.SRL){
  savefile.ini=paste(savefile.ini,"SRL",sep="_")
  epsilon=sort(0.1^seq(0.05,5,by=0.05))
}
if (Run.CCA1){
  savefile.ini=paste(savefile.ini,"CCA1",sep="_")
  epsilon = sort(0.1^seq(-2,5,by=0.02))
}
if (Run.CCA2){
  savefile.ini=paste(savefile.ini,"CCA2",sep="_")
  epsilon = sort(0.1^seq(-3,3,by=0.02))
}

log10epsilon=log10(epsilon)
use.cca=Run.CCA1+Run.CCA2
fix.sector=0
if (use.cca*fix.sector){
  savefile.ini=paste(savefile.ini,"financials",sep="ccasignal_")
}

for (lag1.sign in c(-1,1)){   # -1 is for mean reversion, 1 is for momentum
  for (sectornow in names(all_sectors)){

    cat("\n\n***** Sector now: ",sectornow)
    if (lag1.sign==1)
      cat(" Momentum **********")
    if (lag1.sign==-1)
      cat(" Mean Reversion **********")

    sector=all_sectors[[which(names(all_sectors)==sectornow)]]
    dir.create(savefile.ini, showWarnings = FALSE)
    resfile=paste(paste(savefile.ini,sectornow,sep="/"),"data",sep="_");  load(resfile)
    cat("\nData File Used:",resfile)

    if (lag1.sign==1)
      plotfile.ini=paste(paste(savefile.ini,sectornow,sep="/"),"momentum",sep="_")
    if (lag1.sign==-1)
      plotfile.ini=paste(paste(savefile.ini,sectornow,sep="/"),"meanrevert",sep="_")

    train.portfolio.dataset=dataset[1:train,sector]
    validation.portfolio.dataset=dataset[(train+1):(train+validation),sector]
    test.portfolio.dataset = dataset[(train+validation+1):(nrow(dataset)),sector]

    train.signal.dataset=dataset[1:train,sector]
    validation.signal.dataset=dataset[(train+1):(train+validation),sector]
    test.signal.dataset = dataset[(train+validation+1):(nrow(dataset)),sector]

    if (use.cca*fix.sector){
      train.signal.dataset=dataset[1:train,all_sectors$financials]
      validation.signal.dataset=dataset[(train+1):(train+validation),all_sectors$financials]
      test.signal.dataset=dataset[(train+validation+1):(nrow(dataset)),all_sectors$financials]
    }

    ######################################################################################
    ######################################################################################
    for (use.mean in 1:2){ # havethe choice to subtract the mean or not
      cat("\n\nUse Mean now:",use.mean-1)
      plotfile=paste(plotfile.ini,use.mean-1,sep="_")
      all.solutions=SECTOR_RES[[use.mean]]
      if (lag1.sign==-1){
        solution=all.solutions$meanreversion
      } else {
        solution=all.solutions$momentum
      }
      portfolios=solution$portfolios
      signals=solution$signals

      validation_series=lag1.sign*100*get_stats(portfolios,signals,validation.portfolio.dataset,validation.signal.dataset,get.nu=0,use.cca,use.mean)
      test_series=lag1.sign*100*get_stats(portfolios,signals,test.portfolio.dataset,test.signal.dataset,get.nu=0,use.cca,use.mean)

      validation_pnl=apply(validation_series,2,sum)
      validation_sharpe=apply(validation_series,2,sharpe)
      test_pnl=apply(test_series,2,sum)
      test_sharpe=apply(test_series,2,sharpe)
      train_nu=lag1.sign*drop(get_stats(portfolios,signals,train.portfolio.dataset,train.signal.dataset,get.nu=1,use.cca,use.mean))
      validation_nu=lag1.sign*drop(get_stats(portfolios,signals,validation.portfolio.dataset,validation.signal.dataset,get.nu=1,use.cca,use.mean))
      test_nu=lag1.sign*drop(get_stats(portfolios,signals,test.portfolio.dataset,test.signal.dataset,get.nu=1,use.cca,use.mean))

      lambda.chosen=which.max(validation_pnl)
      lambda.hindsight=which.max(test_pnl)
      market=apply(test.portfolio.dataset,1,mean)
      market.lag1=100*lag1.sign*sign(head(market,-1))*tail(market,-1)
      market=100*tail(market,-1)

      all.pnls=cbind(market,market.lag1,test_series[,1],apply(test_series[,lambda.chosen,drop=F],1,mean),
                     test_series[,ncol(test_series)],test_series[,lambda.hindsight])
      cat("\nDays Market, MarketLag1, Lambda0, Chosen,Gamma,hindsight")
      cat("\nSUM ALL:",nrow(all.pnls),apply(all.pnls,2,sum))
      cat("\nSHARPE ALL:",nrow(all.pnls),apply(all.pnls,2,sharpe))

      if (sectornow=="financials"){
        plot(cumsum(all.pnls[,1]),type="l", main="market all",cex.main=1.5,lwd=2.8)
        dev.copy(png,filename=paste(plotfile,"market.png",sep="_")); dev.off ();
        plot(cumsum(all.pnls[,2]),type="l", main="market.lag1 all",cex.main=1.5,lwd=2.8)
        dev.copy(png,filename=paste(plotfile,"marketlag1.png",sep="_")); dev.off ();
        plot(cumsum(all.pnls[,3]),type="l", main="lambda0 all",cex.main=1.5,lwd=2.8)
        dev.copy(png,filename=paste(plotfile,"lambda0.png",sep="_")); dev.off ();
        plot(cumsum(all.pnls[,4]),type="l", main="algo all",cex.main=1.5,lwd=2.8)
        dev.copy(png,filename=paste(plotfile,"algo.png",sep="_")); dev.off ();
        plot(cumsum(all.pnls[,5]),type="l", main="lambda infinity all",cex.main=1.5,lwd=2.8)
        dev.copy(png,filename=paste(plotfile,"lambdainfinity.png",sep="_")); dev.off ();
        plot(cumsum(all.pnls[,6]),type="l", main="lambda hindsight all",cex.main=1.5,lwd=2.8)
        dev.copy(png,filename=paste(plotfile,"lambdahindsight.png",sep="_")); dev.off ();
      }

      plot(log10epsilon,test_nu,type="l",main="Test Lag-1 Autocorrelation",cex.main=1.5,lwd=2.8,ylab="Test Lag-1 Autocorrelation",xlab=expression(paste(log(epsilon))),cex.lab=1.2)
      dev.copy(png,filename=paste(plotfile,"testnu.png",sep="_")); dev.off (); dev.off ();
      plot(log10epsilon,test_pnl,type="l",main="Test Cumulative Returns",cex.main=1.5,lwd=2.8,ylab="Test Cumulative Returns",xlab=expression(paste(log(epsilon))),cex.lab=1.2)
      dev.copy(png,filename=paste(plotfile,"testpnl.png",sep="_")); dev.off (); dev.off ();
      plot(log10epsilon,test_sharpe,type="l",main="Test Sharpe",cex.main=1.5,lwd=2.8,ylab="Test Sharpe",xlab=expression(paste(log(epsilon))),cex.lab=1.2)
      dev.copy(png,filename=paste(plotfile,"testsharpe.png",sep="_")); dev.off (); dev.off ();

      plot(log10epsilon,train_nu,type="l",main="Train Lag-1 Autocorrelation",cex.main=1.5,lwd=2.8,ylab="Train Lag-1 Autocorrelation",xlab=expression(paste(log(epsilon))),cex.lab=1.2)
      dev.copy(png,filename=paste(plotfile,"trainnu.png",sep="_")); dev.off ();   dev.off ();
      plot(log10epsilon,validation_nu,type="l",main="Validation Lag-1 Autocorrelation",cex.main=1.5,lwd=2.8,ylab="Validation Lag-1 Autocorrelation",xlab=expression(paste(log(epsilon))),cex.lab=1.2)
      dev.copy(png,filename=paste(plotfile,"valnu.png",sep="_")); dev.off ();   dev.off ();
      plot(log10epsilon,validation_pnl,type="l",main="Validation Cumulative Returns",cex.main=1.5,lwd=2.8,ylab="Validation Cumulative Returns",xlab=expression(paste(log(epsilon))),cex.lab=1.2)
      dev.copy(png,filename=paste(plotfile,"valpnl.png",sep="_")); dev.off ();   dev.off ();
      plot(log10epsilon,validation_sharpe,type="l",main="Validation Sharpe",cex.main=1.5,lwd=2.8,ylab="Validation Sharpe",xlab=expression(paste(log(epsilon))),cex.lab=1.2)
      dev.copy(png,filename=paste(plotfile,"valsharpe.png",sep="_")); dev.off (); dev.off ();
    }
  }
}

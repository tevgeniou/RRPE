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


rm(list=ls()); library(compiler); library(parallel)


load("NIPS2013.data")
savefile.ini="Results"

source("Helpers.R")
source("L2learner.R")
source("SRlearner.R")
source("CCAlearner.R")

#################################################################################
#################################################################################
# THE CHOICES TO MAKE

Run.L2=0; Run.SRL=0; Run.CCA1=1; Run.CCA2=0 # choose which method to run
train=1000; validation=250; test=nrow(dataset)-train-validation
use.cca=Run.CCA1+Run.CCA2
# Need to define the all.lambdas below as well as fix nu for CCA in Helpers.R

#################################################################################
#################################################################################


if (Run.L2){
  cat("\n*********** RUNNING THE L2 LEARNER ************\n\n")
  Learner.Used=max.lag1.portfolio.manylambdas
  savefile.ini=paste(savefile.ini,"L2",sep="_");
  all.lambdas=sort(0.1^seq(-6,6,by=0.005))
  all.lambdas=c(0,all.lambdas,large_number) # just make sure lambda="infinity" is at the very end
}
if (Run.SRL){
  cat("\n*********** RUNNING THE SRL LEARNER ************\n\n")
  Learner.Used=SRL
  savefile.ini=paste(savefile.ini,"SRL",sep="_");
  all.lambdas=sort(0.1^seq(0.05,5,by=0.05))
}
if (Run.CCA1){
  cat("\n*********** RUNNING THE CCA1 LEARNER ************\n\n")
  Learner.Used=CCA1
  savefile.ini=paste(savefile.ini,"CCA1",sep="_");
  all.lambdas = sort(0.1^seq(-2,5,by=0.02))
}
if (Run.CCA2){
  cat("\n*********** RUNNING THE CCA2 LEARNER ************\n\n")
  Learner.Used=CCA2
  savefile.ini=paste(savefile.ini,"CCA2",sep="_");
  all.lambdas = sort(0.1^seq(-3,3,by=0.02))
}

############################################################################################
fix.sector=0
if (use.cca*fix.sector){
  savefile.ini=paste(savefile.ini,"financials",sep="ccasignal_")
}

ALL=mclapply(1:length(all_sectors),function(s){
  sector=all_sectors[[s]]
  dir.create(savefile.ini, showWarnings = FALSE)
  savefile=paste(savefile.ini,names(all_sectors)[s],sep="/");
  savefile=paste(savefile,"data",sep="_")

  portfolio.dataset=NULL
  if (use.cca)
    portfolio.dataset=dataset[1:train,sector]
  signal.dataset=dataset[1:train,sector]
  if (use.cca*fix.sector)
    signal.dataset=dataset[1:train,all_sectors$financials]

  SECTOR_RES=mclapply(0:1,function(use.mean){
    cat("\n######################\n","use.mean=",use.mean,"sector=",names(all_sectors)[s],"\n#####################\n")
    # momentum first, mean reverrsion next
    momentum=Learner.Used(portfolio.dataset,signal.dataset,all.lambdas,1,use.mean)
    meanreversion=momentum
    if (!use.cca)
      meanreversion=Learner.Used(portfolio.dataset,signal.dataset,all.lambdas,-1,use.mean)
    list(momentum=momentum,meanreversion=meanreversion)
  })
  save(SECTOR_RES,file=savefile)
})


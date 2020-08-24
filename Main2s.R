#Load packages
library(tidyverse)
library(stringr)
library(caret)
library(tm)
library(dplyr)
library(readxl)
library(e1071)
library(ranger)
library(RcppEigen)
library(ggplot2)
library(devtools)
library(varhandle)
library(MASS)
library(numDeriv)
library(matrixcalc)
library(pracma)
library(qlcMatrix)

###########################################
###### PART 1 INITIALIZE MODEL ############
###########################################

#Set exogenous paramters E class
sigmaE=1.5
betaE=.9932
probE=matrix(1,2,2)
probE[1,1]=.7
probE[1,2]=1-probE[1,1]
probE[2,2]=.95
probE[2,1]=1-probE[2,2]

nshocks=length(probE[,1])
techE=matrix(1,1,2)
techE[1,1]=.1
techE[1,2]=1
ncols=length(techE)

dist_E<-invdist(probE)
E_bar<-dist_E %*% t(techE)

#Set exogenous paramters N class
sigmaN=1.5
betaN=.9932
probN=matrix(1,2,2)
probN[1,1]=.5
probN[1,2]=1-probN[1,1]
probN[2,2]=.5
probN[2,1]=1-probN[2,2]

techN=matrix(1,1,2)
techN[1,1]=.1
techN[1,2]=1

dist_N<-invdist(probN)
N_bar<- dist_N %*% t(techN)

#Construct Grid for N Class
agridmin_N=0
agridmax_N=2
acutpoint=0

nassets1=100
nassets2=100
  
stepsize1=(acutpoint-agridmin_N)/(nassets1-1)
stepsize2=(agridmax_N-acutpoint)/nassets2

agrid1=seq(agridmin_N,acutpoint,by=stepsize1)
agrid2=seq(acutpoint+stepsize2,agridmax_N,by=stepsize2)
agrid_N=c(agrid1,agrid2)
nassetsN=length(agrid_N)

#Construct Grid for E Class
agridmin_E=-2
agridmax_E=2
acutpoint=0

nassets1=100
nassests2=100

stepsize1=(acutpoint-agridmin_E)/(nassets1-1)
stepsize2=(agridmax_E-acutpoint)/nassets2

agrid1=seq(agridmin_E,acutpoint,by=stepsize1)
agrid2=seq(acutpoint+stepsize2,agridmax_E,by=stepsize2)
agrid_E=c(agrid1,agrid2)
nassetsE=length(agrid_E)


#vectorize parameters for (1) discrete shock (2) agrid
ntotalgrid_E=nshocks*nassetsE
ntotalgrid_N=nshocks*nassetsN

prob_vec_E=matrix(rep( t(probE),nassetsE),ncol=2,byrow=TRUE)
tech_vec_E<-t(matrix( rep( techE, nassetsE ), 1 , byrow=TRUE))
agrid_vec_E<-kronecker(agrid_E,matrix(1,nshocks,1))

prob_vec_N=matrix(rep( t(probN),nassetsN),ncol=2,byrow=TRUE)
tech_vec_N<-t(matrix( rep( techN, nassetsN ), 1 , byrow=TRUE))
agrid_vec_N<-kronecker(agrid_N,matrix(1,nshocks,1))

v_N=matrix(0,nshocks,nassetsN)
dec_N=matrix(0,nshocks,nassetsN)

v_E=matrix(0,nshocks,nassetsE)
dec_E=matrix(0,nshocks,nassetsE)


maxiter_v=1000
maxiter_out=250
tol_v=.01
tol_p=.01
tol_R=.01
tol_a=.01

RPath=matrix(0,maxiter_out,1)
aPath=matrix(0,maxiter_out,1)

relax=.50

R_low=-.10
R_high=1/betaE - 1
R=relax*(R_low)+(1-relax)*R_high
niter_out=1

#Outer Loop
for(niter_out in 1:maxiter_out){

print(niter_out)
wage=1
a_bar=0

total_vec_N=(1+R)*agrid_vec_N+wage*tech_vec_N
aprime_mat_N=repmat(agrid_N,ntotalgrid_N,1)
total_mat_N=repmat(total_vec_N,1,nassetsN)
consum_mat_N=total_mat_N-aprime_mat_N

util_mat_N=replace(consum_mat_N,consum_mat_N<0,-Inf)
util_mat_N=(util_mat_N^(1-sigmaN))/(1-sigmaN)
util_mat_N=replace(util_mat_N,is.nan(util_mat_N),-Inf)

total_vec_E=(1+R)*agrid_vec_E+wage*tech_vec_E
aprime_mat_E=repmat(agrid_E,ntotalgrid_E,1)
total_mat_E=repmat(total_vec_E,1,nassetsE)
consum_mat_E=total_mat_E-aprime_mat_E

util_mat_E=replace(consum_mat_E,consum_mat_E<0,-Inf)
util_mat_E=(util_mat_E^(1-sigmaE))/(1-sigmaE)
util_mat_E=replace(util_mat_E,is.nan(util_mat_E),-Inf)
  

#Inner Loop ATTEMPT ONE: TWO INNERS (ONE FOR EACH CLASS)

for(i in 1:maxiter_v){
  
  expect_value_mat_N<-repmat(betaN*probN %*% v_N,nassetsN,1)
  totalval_mat_N=util_mat_N+expect_value_mat_N
  tv_N<-apply(X=totalval_mat_N,MARGIN=1,FUN=max)
  indx_tdec_N<-max.col(totalval_mat_N)
  tdec_N=agrid_N[indx_tdec_N]
  
  metric_v=max(abs(tv_N-as.vector(v_N)))/max(abs(as.vector(v_N)))
  metric_p=max(abs(tdec_N-as.vector(dec_N)))
  
  if((metric_v<tol_v)&(metric_p<tol_p))
  {break}else
  {v_N=matrix(tv_N,nshocks,nassetsN)
    dec_N=matrix(tdec_N,nshocks,nassetsN)}
}

for(i in 1:maxiter_v){
 
  expect_value_mat_E<-repmat(betaE*probE %*% v_E,nassetsE,1)
  totalval_mat_E=util_mat_E+expect_value_mat_E
  tv_E<-apply(X=totalval_mat_E,MARGIN=1,FUN=max)
  indx_tdec_E<-max.col(totalval_mat_E)
  tdec_E=agrid_E[indx_tdec_E]
  
  metric_v=max(abs(tv_E-as.vector(v_E)))/max(abs(as.vector(v_E)))
  metric_p=max(abs(tdec_E-as.vector(dec_E)))
  
  if((metric_v<tol_v)&(metric_p<tol_p))
  {break}else
  {v_E=matrix(tv_E,nshocks,nassetsE)
  dec_E=matrix(tdec_E,nshocks,nassetsE)}
}

#Compute Invarient distribution of assets and shocks -- class N
trans=matrix(0,ntotalgrid_N,ntotalgrid_N)

for(k in 1:ntotalgrid_N){
  idx_start=(indx_tdec_N[k]-1)*nshocks+1
  idx_end=idx_start+nshocks-1
  trans[k,idx_start:idx_end]=prob_vec_N[k,]}

dist_vec_N=invdist(trans)
dist_2d_N=dist_vec_N[1,]
dist_2d_N=matrix(dist_2d_N,nrow=2,byrow = FALSE)
dist_agrid_N=colSums(dist_2d_N)


#Compute Invarient distribution of assets and shocks -- class E
trans=matrix(0,ntotalgrid_E,ntotalgrid_E)

for(k in 1:ntotalgrid_E){
  idx_start=(indx_tdec_E[k]-1)*nshocks+1
  idx_end=idx_start+nshocks-1
  trans[k,idx_start:idx_end]=prob_vec_E[k,]}

dist_vec_E=invdist(trans)
dist_2d_E=dist_vec_E[1,]
dist_2d_E=matrix(dist_2d_E,nrow=2,byrow = FALSE)
dist_agrid_E=colSums(dist_2d_E)

#Equilibrium Updating
na_bar=(dist_agrid_E %*% agrid_E) + (dist_agrid_N %*% agrid_N)
a_netsupply=na_bar-a_bar
metric_a=a_netsupply

if(metric_a>0)
  {R_high=R}else{if(metric_a<0)
  {R_low=R}else{R_high=R
   R_low=R}}

R_new = relax*R_low+(1-relax)*R_high
metric_R=abs(R_new-R)

RPath[niter_out]=R
aPath[niter_out]=a_netsupply
if((abs(metric_a)<tol_a)&(abs(metric_R)<tol_R)){break}else{R=R_new}
}







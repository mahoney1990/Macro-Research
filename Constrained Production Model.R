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

#Set exogenous paramters H class
betaH=.95
betaL=0.95
alphaH=0.7
alphaL=0.5
A_H=2
A_L=1
sigmaH=1.5
sigmaL=1.2

#Construct Grid for H Class
agridmin_N=-1
agridmax_N=1
acutpoint=0

nassets1=100
nassets2=100

stepsize1=(acutpoint-agridmin_N)/(nassets1-1)
stepsize2=(agridmax_N-acutpoint)/nassets2

agrid1=seq(agridmin_N,acutpoint,by=stepsize1)
agrid2=seq(acutpoint+stepsize2,agridmax_N,by=stepsize2)
agrid_H=c(agrid1,agrid2)
agrid_L=c(agrid1,agrid2)
nassetsL=length(agrid_L)
nassetsH=length(agrid_H)


#Simple case: no disutility to labor supply
K=.5 #Constraint on L sector production
w_H=alphaH*A_H
w_L=alphaL*A_L*((K/A_L)^((alphaL-1)/alphaL))

#Initialize Value Function Iteration
r=.01

Y_H=A_H
Y_L=A_L*((K/A_L)^((alphaL-1)/alphaL))
Y=Y_L+Y_H

a_bar=0
c_bar=Y

v_H=matrix(0,1,nassetsH)
v_L=matrix(0,1,nassetsL)
dec_H=matrix(0,1,nassetsH)
dec_L=matrix(0,1,nassetsL)

tol_v=.001
tol_p=.001

maxiter_v=1000

#Outer loop (iterate on R to find fixed point in assets)

#Construct Utility Matrix for H type
total_vec_H=(1+r)*agrid_H+w_H
aprime_mat_H=repmat(agrid_H,nassetsH,1)
total_mat_H=repmat(total_vec_H,nassetsH,1)
consum_mat_H=total_mat_H-t(aprime_mat_H)

util_mat_H<-replace(consum_mat_H, consum_mat_H<0,-Inf)
util_mat_H=(util_mat_H^(1-sigmaH))/(1-sigmaH)
util_mat_H=replace(util_mat_H, is.nan(util_mat_H),-Inf)

#Construct Utility Matrix for L type
total_vec_L=(1+r)*agrid_L+w_L
aprime_mat_L=repmat(agrid_L,nassetsL,1)
total_mat_L=repmat(total_vec_L,nassetsL,1)
consum_mat_L=(total_mat_L)-t(aprime_mat_L)

util_mat_L<-replace(consum_mat_L, consum_mat_L<0,-Inf)
util_mat_L=(util_mat_L^(1-sigmaL))/(1-sigmaL)
util_mat_L=replace(util_mat_L, is.nan(util_mat_L),-Inf)


#Inner loop H 
for(i in 1:maxiter_v){
  print(i)
value_mat_H=repmat(betaH*v_H,nassetsH,1)
totalval_mat_H=util_mat_H+value_mat_H
tv_H<-apply(X=totalval_mat_H,MARGIN = 2,FUN=max)
indx_tdec_H=max.col(t(total_mat_H))
tdec_H=agrid_H[indx_tdec_H]

metric_v=max(abs(tv_H-as.vector(v_H)))/max(abs(as.vector(v_H)))
metric_p=max(abs(tdec_H-as.vector(dec_H)))

if((metric_v<tol_v)&(metric_p<tol_p))
{break}else
{v_H=matrix(tv_H,1,nassetsH)
    dec_H=matrix(tdec_H,1,nassetsH)}
}


value_mat_L=repmat(betaL*v_L,nassetsL,1)
totalval_mat_L=util_mat_L+value_mat_L
tv_L<-apply(totalval_mat_L,MARGIN = 1,FUN=max)
  


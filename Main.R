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

#Set exogenous paramters
sigma=1.5
beta=.9932
prob=matrix(1,2,2)
prob[1,1]=.7
prob[1,2]=1-prob[1,1]
prob[2,2]=.95
prob[2,1]=1-prob[2,2]

nshocks=length(prob[,1])
tech=matrix(1,1,2)
tech[1,1]=.1
tech[1,2]=1
ncol=length(tech)
nShock=size(prob,1)
#Calculate the Invariant Distribtuion of the Markov Process

dist_s<-invdist(prob)
e_bar  = dist_s %*% t(tech)

#Construct Grids
agridMin = -1.0 # Borrowing constraint                       
agridMax =  2.0 # Maximum wealth                      
aCutPoint = 0.0 # Partition Point

nasset1=10
nasset2=10

stepsize1=(aCutPoint-agridMin)/(nasset1-1)
stepsize2=(agridMax-aCutPoint)/nasset2

agrid1=seq(agridMin, aCutPoint, by=stepsize1)
agrid2=seq(aCutPoint+stepsize2,agridMax, by=stepsize2)

agrid=c(agrid1,agrid2) # Grid of potential asset values
nassets=length(agrid) # size of grid

#vectorize parameters with order (1)discrete shock;(2)agrid;
nTotalGrid = nshocks*nassets

prob_vec<-matrix( rep( t(prob) , nassets ) , ncol=2 , byrow = TRUE )
tech_vec<-t(matrix( rep( tech, nassets ), 1 , byrow=TRUE))
agrid_vec<-kronecker(agrid,matrix(1,nshocks,1))
agridMin_vec <-  matrix( rep( agridMin, nTotalGrid) , 1, byrow=TRUE)
agridMax_vec <-  matrix( rep( agridMax, nTotalGrid) , 1, byrow=TRUE)

nAssetInvDist = 1001;
StepSizeInvDist = (agridMax-agridMin)/(nAssetInvDist-1) 
agridInvDist = seq(agridMin, agridMax, by=StepSizeInvDist)

nTotalGridInvDist = nShock*nAssetInvDist;          

probInvDist_vec <- matrix( rep( prob , nAssetInvDist ) , ncol , byrow = TRUE )
techInvDist_vec <- t(matrix( rep( tech, nAssetInvDist ), 1 , byrow=TRUE))
agridInvDist_vec<- kronecker(agridInvDist,matrix(1,nshocks,1))
agridMin_vec <-  matrix( rep( agridMin, nTotalGridInvDist) , 1, byrow=TRUE)
agridMax_vec <-  matrix( rep( agridMax, nTotalGridInvDist) , 1, byrow=TRUE)

######################################################
####### PART TWO -- VALUE FUNCTION ITERATION #########
######################################################

maxiter_out= 100
maxiter_v = 2000
maxiter_h = 100
maxiter_optim = 100

tol_a=.0001
tol_R=.0001
tol_v=.000001
tol_p=.000001
tol_optim=.000001

RPath = matrix(0,maxiter_out,1)
aPath = matrix(0,maxiter_out,1)

relax=.50

R_low=-.1
R_high= 1/beta - 1
R = relax*R_low + (1-relax)*R_high

wage=1
a_bar=0


v=matrix(0,nshocks,nassets)

dec=matrix(0,nshocks,nassets)

#Outer loop
for(niter_out in 1:maxiter_out ){
  print(niter_out)
wage=1
a_bar=0

total_vec=(1+R)*agrid_vec+wage*tech_vec
aprime_mat=repmat(agrid,nTotalGrid,1)
total_mat=repmat(total_vec,1,nassets)
consum_mat=total_mat-aprime_mat

util_mat <- replace(consum_mat, consum_mat < 0, -Inf)
util_mat=(util_mat^(1-sigma))/(1-sigma)
util_mat=replace(util_mat, is.nan(util_mat), -Inf)

#Inner loop
for(i in 1:maxiter_v){
  
  expect_val_mat=repmat(beta*prob %*% v,nassets,1) #First Error
  totalval_mat=util_mat+expect_val_mat
  tv<-apply(X=totalval_mat, MARGIN=1, FUN=max)
  indx_tdec=max.col(totalval_mat)
  tdec=agrid[indx_tdec]
  
  metric_v=max((abs(tv-as.vector(v)))/max(abs(as.vector(v)))) #Second error
  metric_p=max(abs(tdec-as.vector(dec))) #Note differences
  
  if((metric_v<tol_v)&(metric_p<tol_p))
  {break}else
  {v=matrix(tv,nshocks,nassets)
   dec=matrix(tdec,nshocks,nassets)}
}

#Compuete invariant density of assets and shocks
trans=matrix(0, nTotalGrid, nTotalGrid)

for(i in 1:nTotalGrid)
  {idx_start = (indx_tdec[i] - 1)*nShock + 1
   idx_end = idx_start + nShock - 1 
   trans[i,idx_start:idx_end] = prob_vec[i,]}

dist_vec=invdist(trans)
dist_2d=dist_vec[1,]
dist_2d=matrix(dist_2d,nrow=2,byrow=FALSE)
dist_agrid=colSums(dist_2d)

#Equlibrium updating
na_bar = dist_agrid %*% agrid       #aggregate asset in n-th iteration
a_netSupply = na_bar - a_bar  # net supply of assets
metric_a = a_netSupply     


if(metric_a > 0)        
  {R_high = R}else{if(metric_a < 0)        
  {R_low  = R}else{R_high = R
  R_low  = R}}
         
R_new  = relax * R_low + (1 - relax) * R_high
metric_R = abs(R_new - R)

RPath[niter_out] = R
aPath[niter_out] = a_netSupply
if((abs(metric_a) < tol_a)&(abs(metric_R) < tol_R)){break}else{R = R_new}
}




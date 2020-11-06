#Part One data smoothing
rm(list = ls()) #Clear environment
library(readxl)
#Make sure plm package is not loaded, plm lag function replaces base R lag fn
if(any(search()=="package:plm")){
detach("package:plm", unload = TRUE)}



setwd("C:/Users/Makkel McDougal/Desktop/Macro Research")

data<-read_xlsx("Covid Data Backup.xlsx")
data<-subset(data, Date>"2020-01-29")
data<-data[!(data$State=="Alaska" |data$State=="Maine"|data$State=="Montana"|data$State=="Puerto Rico"|data$State=="Vermont"|data$State=="Hawaii"| data$State=="West Virginia"|data$State=="Wyoming"|data$State=="Tennessee"|data$State=="South Dakota"|data$State=="North Dakota"|data$State=="Nebraska"|data$State=="Missouri"),]
#data<-subset(data, LD_Start!=as.Date("2100-01-01")) # Drop States with no lockdown 
#This is for data manipulation pruposes, will fix later

library(mFilter)
library(lubridate)
library(dplyr)
library(zoo)

gamma=.20 #Biological parameters -- transition rate
theta=.10 #Biological parameters -- recovery rate
window=99 #Choose window


st_list<-unique(data$State)

beta_df <- data.frame(matrix(ncol = length(st_list), nrow = window-1))
R_df <- data.frame(matrix(ncol = length(st_list), nrow = window-1))

colnames(beta_df)<-st_list
#beta_df$Date<-unique(data$Date)

colnames(R_df)<-st_list
#R_df$Date<-unique(data$Date)

start=matrix(0,3,length(st_list))
start[1,]=st_list

end=matrix(0,3,length(st_list))
end[1,]=st_list
i=21

for(i in 1:length(st_list)){
  print(i)
  ts<-subset(data,State==st_list[i])$DeathMA #Note change here
  
  start[2,i]=as.numeric(min(which(ts>0))+5)
  start[3,i]=as.Date(data$Date[as.numeric(start[2,i])])
  
  end[2,i]=as.numeric(start[2,i])+window
  end[3,i]=as.Date(data$Date[as.numeric(end[2,i])])
  
  ts<-ts[start[2,i]:end[2,i]]
  hp<-hpfilter(as.ts(ts),(250))
  
  smoothed_trend_ts<-as.vector(hp$trend)
  #smoothed_trend_ts[smoothed_trend_ts<0]=gamma
  d_t3=smoothed_trend_ts[3:(window+1)]
  d_t2=lag(smoothed_trend_ts,1)[3:(window+1)]
  d_t1=lag(smoothed_trend_ts,2)[3:(window+1)]
  d_t=lag(smoothed_trend_ts,3)[3:(window+1)]
  
  
  delta_3=hpfilter(as.ts(d_t3-d_t2),50)
  delta_3=delta_3$trend
  delta_2=hpfilter(as.ts(d_t2-d_t1),50)
  delta_2=delta_2$trend
  #Key identifying equation from SIRD(2020)
  beta=as.ts(gamma+((1/theta)*(delta_3-delta_2)+delta_2)/((1/theta)*delta_2+d_t1))
  
  beta[beta<0] <- gamma
  beta_df[[i]]<-beta
  R_df[[i]]<-beta*(1/gamma)}

#Form some plots
for(i in 1:length(st_list)){
  print(i)
  plot(R_df[[i]],ylab=st_list[i])
}


#Part Two: Match observations to dates/lockdown/betas

#Construct matrix of cutoff values
cutoff<-matrix(0,length(st_list),7)
cutoff[,1]<-st_list
i=3
for(i in 1:length(st_list)){

st_data<-subset(data,State==st_list[i])
lower=as.Date(as.numeric(start[3,i]))
upper=as.Date(as.numeric(end[3,i]))

first_death<-as.numeric(min(which(st_data$Daily_Deaths>0)))
data$Date[56]


LD_start<-st_data$LD_Start[1]
LD_end<-st_data$LD_End[1]


if(LD_start>lower){
  span=seq(lower, upper, by="days")}else{
  span=seq(LD_start,LD_start+(60*60*24*window),by="days")}

if(LD_start==as.Date("2100-01-01")){idx=7
  idx_l=1}else{
idx<-which(span==LD_start)
idx_l =max(idx-7,1)}

cutoff[i,2]<-mean(R_df[[i]][idx_l:idx],na.rm=TRUE)
cutoff[i,3]<-as.Date(as.numeric(lower))
cutoff[i,4]<-as.Date(as.numeric(upper))
cutoff[i,5]<-as.Date(LD_start) #Change Here
cutoff[i,6]<-as.Date(LD_end)
cutoff[i,7]<-idx
}

state_info<-data.frame(cutoff)
colnames(state_info)<-c("State","Cutoff","lower","upper","LD_Start","LD_End","Index")


library(ggplot2)
library(reshape2)

t<-seq(1,window-1)
beta_df$time=t
R_df$time=t



#Build Regression Data Matrix
roll=7

weekly_beta<-matrix(NA,length(beta_df$Alabama),length(st_list))

for(i in 1:length(st_list)){
  weekly_beta[roll:(window-1),i]<-rollmean(beta_df[,i],7)
}

weekly_beta<-data.frame(weekly_beta)

reg_mat=matrix(0,length(t)*length(st_list),7)
i=1
for(i in 1:length(st_list)){
  l=(i-1)*length(t)+1
  u=i*length(t)
  reg_mat[l:u,1]=as.numeric(t)
  reg_mat[l:u,2]=st_list[i]
  reg_mat[l:u,3]=weekly_beta[,i]
  reg_mat[l:u,4]=beta_df[[i]]
  
  lower=as.Date(as.numeric(start[3,i]))+2
  upper=as.Date(as.numeric(end[3,i]))
  span=seq(lower, upper, by="days")
  reg_mat[l:u,5]=span
  reg_mat[l:u,5]=as.Date(as.numeric(reg_mat[l:u,5]))
  
  
  reg_mat[l:u,6]=as.Date(as.numeric(state_info$LD_Start[i]))
  reg_mat[l:u,7]=as.Date(as.numeric(state_info$LD_End[i]))}

  daily_data_df<-data.frame(reg_mat)
  colnames(daily_data_df)<-c("time","State","BetaMA","Beta","Date","LDstart","LDend")
  
  daily_data_df$Date<-as.Date(as.numeric(daily_data_df$Date))
  daily_data_df$LDstart<-as.Date(as.numeric(daily_data_df$LDstart))
  daily_data_df$LDend<-as.Date(as.numeric(daily_data_df$LDend))
  
  #reg_mat[,7]=ifelse(reg_mat[,4]>reg_mat[,5],1,0)
  #reg_mat[,8]=ifelse(reg_mat[,4]>reg_mat[,6],2,0)

 
  daily_data_df$Date<-as.POSIXct.Date(as.numeric(daily_data_df$Date))
  
  #Form Daily Panel
  daily_panel <- merge(daily_data_df, data, by = c("Date","State"), sort = F, all.x = T)
  daily_panel$Phase<-as.factor(daily_panel$Phase)
  
  #Import Weekly Data Unemployment Data
  weekly_data<-read_xlsx("Unemployment Claims.xlsx")
  weekly_data<-subset(weekly_data, Date> "2020-01-31")
  
  #Merge weekly observations into daily panel
  #We need to interpolate missing daily values
  daily_panel<-merge(daily_panel,weekly_data,by=c("State","Date"),all.x =TRUE )
  
  #Interpolate Labor Market Data to Daily
  for(j in 1:length(st_list)){
    index<-which(daily_panel$State==st_list[j])
    start<-min(which(daily_panel$Lconstraint[index]>0))
    start<-index[start]
    stop<-max(which(daily_panel$Lconstraint[index]>0))
    stop<-index[stop]
    daily_panel$Lconstraint[start:stop]<- na.approx(daily_panel$Lconstraint[start:stop])}
  
  #Remove Extraneous Data
  daily_panel<-subset(daily_panel, select = c("Date","time","State","Beta","Lconstraint","Phase.x"))
  
  #Run full panel regression
  #X<-panel_model(daily_panel)
  #beta_reg<-X$beta_reg
  #L_reg<-X$L_reg
  
  #Run SIRD simulation
  #Generate counterfactual dynamic paths for beta
  #CF_data<-read_xlsx("CF Test.xlsx")
  #CF_data$time<-as.factor(CF_data$time)
  #CF_data$Phase<-as.factor(CF_data$Phase)
  #CF<-predict(beta_reg,CF_data)
  
  #SIRD function goes here
  
  
  
  #daily_panel=na.omit(daily_panel) not yet!

  
  
  
  colnames(weekly_beta)=st_list
  weekly_beta<-na.omit(weekly_beta)
  library(gridExtra)
  



  
  

#Part Three Compute Weekly averages & match with unemployment data
#Form Weekly Panel
weekly_data<-read_xlsx("Unemployment Claims.xlsx")
weekly_data<-subset(weekly_data, Date> "2020-01-31")

panel <- merge(daily_data_df, weekly_data, by = c("Date","State"), sort = F, all.x = T)
panel<-na.omit(panel)
panel<-as.data.frame(panel)
panel$Stay<-ifelse(panel$Date>panel$LDstart & panel$Date<panel$LDend, 1, 0 )
panel$Post<-ifelse(panel$Date>panel$LDstart & panel$Date>panel$LDend, 1, 0 )
panel$Beta<-as.numeric(panel$Beta)
panel$State<-as.factor(panel$State)

#Part Four Estimation
#Run Both Smoothing scripts first 
library(plm)


panel<-pdata.frame(panel,index=c("State","Date")) 
lead_beta=lead(panel$Beta,1)
panel$lead_beta<-lead_beta 
panel<-na.omit(panel) 
panel$BetaMA<-as.numeric(panel$BetaMA)


for(i in 1:length(st_list)){
idx<-which(panel$State==st_list[i])
span<-seq(1:length(idx))
panel$time[idx]<-span}

idx<-which(panel$Phase==0)
panel$Phase==as.numeric(panel$Phase)
panel$Phase[idx]=4
panel$Phase=as.factor(panel$Phase)

beta_reg<-glm(BetaMA ~ Phase+State+time,data=panel)
L_reg<-glm(Lconstraint~Phase+time+State,data=panel)











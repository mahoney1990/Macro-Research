library(readxl)
setwd("C:/Users/Makkel McDougal/Desktop/Macro Research")

data<-read_xlsx("Covid Data Backup.xlsx")
data<-subset(data, Date>"2020-01-29")
data<-subset(data, LD_Start!=as.Date("1900-01-01")) # Drop States with no lockdown 
#This is for data manipulation pruposes, will fix later

library(mFilter)
library(lubridate)
library(dplyr)
library(zoo)

gamma=.20 #Biological parameters -- transition rate
theta=.10 #Biological parameters -- recovery rate
window=89 #Choose window


st_list<-unique(data$State)

beta_df <- data.frame(matrix(ncol = length(st_list), nrow = window+1))
R_df <- data.frame(matrix(ncol = length(st_list), nrow = window+1))

colnames(beta_df)<-st_list
#beta_df$Date<-unique(data$Date)

colnames(R_df)<-st_list
#R_df$Date<-unique(data$Date)

start=matrix(0,3,length(st_list))
start[1,]=st_list

end=matrix(0,3,length(st_list))
end[1,]=st_list

for(i in 1:length(st_list)){
  print(i)
  ts<-subset(data,State==st_list[i])$DeathMA
 
  start[2,i]=as.numeric(min(which(ts>0))-7)
  start[3,i]=as.Date(data$Date[as.numeric(start[2,i])])
  
  end[2,i]=as.numeric(start[2,i])+window
  end[3,i]=as.Date(data$Date[as.numeric(end[2,i])])
  
  ts<-ts[start[2,i]:end[2,i]]
  hp<-hpfilter(as.ts(ts),120)
  
  smoothed_trend_ts<-as.vector(hp$trend)
  d_t3=smoothed_trend_ts
  d_t2=lag(smoothed_trend_ts,1)
  d_t1=lag(smoothed_trend_ts,2)
  d_t=lag(smoothed_trend_ts,3)
  
  delta_3=d_t3-d_t2
  delta_2=d_t2-d_t1
  
  #Key identifying equation from SIRD(2020)
  beta=as.ts(gamma+((1/theta)*(delta_3-delta_2)+delta_2)/((1/theta)*delta_2+d_t1))
  beta[beta<0] <- 0
  beta_df[[i]]<-beta
  R_df[[i]]<-beta*(1/gamma)
}


#Part Two: Match observations to dates/lockdown/betas

#Construct matrix of cutoff values
cutoff<-matrix(0,length(st_list),7)
cutoff[,1]<-st_list

for(i in 1:length(st_list)){

st_data<-subset(data,State==st_list[i])
lower=as.Date(as.numeric(start[3,i]))
upper=as.Date(as.numeric(end[3,i]))

first_death<-as.numeric(min(which(st_data$Daily_Deaths>0)))
data$Date[57]


LD_start<-st_data$LD_Start[1]
LD_end<-st_data$LD_End[1]

if(LD_start>lower){
  span=seq(lower, upper, by="days")}else{
  span=seq(LD_start,LD_start+(60*60*24*window),by="days")}


idx<-which(span==LD_start)
idx_l=max(idx-7,1)

cutoff[i,2]<-mean(R_df[[i]][idx_l:idx],na.rm=TRUE)
cutoff[i,3]<-as.Date(as.numeric(lower))
cutoff[i,4]<-as.Date(as.numeric(upper))
cutoff[i,5]<-as.Date(as.numeric(LD_start))
cutoff[i,6]<-as.Date(as.numeric(LD_end))
cutoff[i,7]<-idx
}

state_info<-data.frame(cutoff)
colnames(state_info)<-c("State","Cutoff","lower","upper","LD_Start","LD_End","Index")


library(ggplot2)
library(reshape2)

t<-seq(1,window+1)
beta_df$time=t
R_df$time=t



#Build Regression Data Matrix -- Fuck it, just brute force it with a loop
roll=7

weekly_beta<-matrix(NA,length(beta_df$Alabama),length(st_list))

for(i in 1:length(st_list)){
  weekly_beta[roll:(window+1),i]<-rollmean(beta_df[,i],7)
}

weekly_beta<-data.frame(weekly_beta)

reg_mat=matrix(0,length(t)*length(st_list),8)

for(i in 1:length(st_list)){
  l=(i-1)*length(t)+1
  u=i*length(t)
  reg_mat[l:u,1]=as.numeric(t)
  reg_mat[l:u,2]=st_list[i]
  reg_mat[l:u,3]=weekly_beta[,i]
  
  
  lower=as.Date(as.numeric(start[3,i]))
  upper=as.Date(as.numeric(end[3,i]))
  span=seq(lower, upper, by="days")
  reg_mat[l:u,4]=span
  reg_mat[l:u,4]=as.Date(as.numeric(reg_mat[l:u,4]))
  
  reg_mat[l:u,5]=as.Date(as.numeric(state_info$LD_Start[i]))
  reg_mat[l:u,6]=as.Date(as.numeric(state_info$LD_End[i]))}

  daily_data_df<-data.frame(reg_mat)
  colnames(daily_data_df)<-c("time","State","Beta","Date","LDstart","LDend","Stay","Post")
  
  daily_data_df$Date<-as.Date(as.numeric(daily_data_df$Date))
  daily_data_df$LDstart<-as.Date(as.numeric(daily_data_df$LDstart))
  daily_data_df$LDend<-as.Date(as.numeric(daily_data_df$LDend))
  
  reg_mat[,7]=ifelse(reg_mat[,4]>reg_mat[,5],1,0)
  reg_mat[,8]=ifelse(reg_mat[,4]>reg_mat[,6],2,0)

 
  daily_data_df$Date<-as.POSIXct.Date(as.numeric(daily_data_df$Date))
  
  #Form some plots
  plot(R_df$California)
  plot(R_df$Illinois)
  plot(R_df$Louisiana)
  plot(R_df$Washington)
  plot(R_df$Oregon)
  
  
  
  ggplot() + 
    geom_line(data = R_df, aes(x = time, y = Washington), color = "red") +
    geom_line(data = R_df, aes(x = time, y = Oregon), color = "blue") +
    xlab('Time') +
    ylab('Reproduction Rate')
  
  

#Part Three Compute Weekly averages & match with unemployment data
weekly_data<-read_xlsx("Unemployment Claims.xlsx")
weekly_data<-subset(weekly_data, Date> "2020-01-31")

panel <- merge(daily_data_df, weekly_data, by = c("Date","State"), sort = F, all.x = T)
panel<-na.omit(panel)
panel<-as.data.fram(panel)
panel$LD=lag()



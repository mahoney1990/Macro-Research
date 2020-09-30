library(readxl)
setwd("C:/Users/Makkel McDougal/Desktop/Macro Research")

US_data<-read_xlsx("UStotals.xlsx")


library(mFilter)
library(lubridate)

ts<-subset(US_data,Date>"2020-02-29")$DeathMA
hp<-hpfilter(as.ts(ts),50)

gamma=.20
theta=.10

smoothed_trend_ts<-as.vector(hp$trend)
d_t3=smoothed_trend_ts
d_t2=lag(smoothed_trend_ts,1)
d_t1=lag(smoothed_trend_ts,2)
d_t=lag(smoothed_trend_ts,3)

delta_3=d_t3-d_t2
delta_2=d_t2-d_t1

beta=as.ts(gamma+((1/theta)*(delta_3-delta_2)+delta_2)/((1/theta)*delta_2+d_t1))
beta[beta<0] <- 0

plot(beta)

mean(US_data$Date)





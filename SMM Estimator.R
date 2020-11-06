SMM<-function(w){
 
  V_L=colSums(L_st^2)
  V_d=w*colSums(d_st^2)
  
  loss=-V_L+V_d
  
  min_idx=which.min(loss)
  
  diff=sum((as.numeric(Observed_Phase)-as.numeric(Valid_Seq[,min_idx]))^2)
  
  return(diff)}

output=matrix(0,length(st_list),4)
output[,1]=st_list
colnames(output)=c("State","W_hat","Sequence","Index")

w0=.5

for(j in 1:5){
print(st_list[j])
idx=which(panel$State==st_list[j])
Observed_Phase=panel$Phase[idx]
#Observed_Phase=as.numeric(Observed_Phase)-1

L_st<-L_predictions[,,j]
d_st<-d_predictions[,,j]

#if(length(d_st[,1])==14){d_st<-d_st[-1,]}

eta=.9999

indx=min(which(data$State==st_list[j]))
pop=data$Pop[indx]/1000

for(i in 1:M){
  L_st[i,]=eta^(i-1)*L_st[i,]
  d_st[i,]=eta^(i-1)*d_st[i,]*pop}

#which.min(loss)
est<-optim(w0,SMM,method="SANN",control=list(maxit=25000))
output[j,2]=as.numeric(est$par)
output[j,3]=as.numeric(est$value)

V_L=colSums(L_st^2)
V_d=as.numeric(output[j,2])*colSums(d_st^2)
loss=-V_L+V_d
output[j,4]=which.min(loss)

}

as.numeric(panel$Phase[which(panel$State=="Florida")])

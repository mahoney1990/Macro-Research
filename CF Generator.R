#CF data generation
N<-length(which(panel$State=="Alabama" & panel$Lconstraint>0))


span=which(panel$State==st_list[1])
CF$Phase=panel$Phase[span]

start<-which(daily_panel$State==CF$State[1] & daily_panel$Lconstraint>0)
CF$Date<-daily_panel$Date[start]

col_size<-3^N
#Construct Matrix of all possible phase-sequences
CF_matrix<-matrix(0,nrow=N,ncol=col_size)


vec <- c(1, 2,3)
lst <- lapply(numeric(N), function(x) vec)
X<-t(as.matrix(expand.grid(lst)))
K<-length(X[1,])
Z<-matrix(0,1,K)

for(j in 1:K){
  print(j)
  idx<-seq(from=min(which(X[,j]<3)), to=N)
  if(idx==13){next}
  Y<-X[,j]
  Y<-Y[idx]

  for(i in 1:(length(idx)-1)){
    tf<-(Y[i]>Y[i+1])
    tf<-as.logical(tf)
    if(tf==TRUE){Z[j]<-1}
  }  
}

Z<-which(Z>0)
Valid_Seq<-X[,-Z]

#Import Valid 3dim sequences
Valid_Seq<-read.csv("Valid Sequences 3 dim.csv")

#Import Valid 4dim sequences
Valid_Seq<-read.csv("Valid_Seq_4dim.csv")
Valid_Seq<-Valid_Seq[,-1]

K<-length(Valid_Seq[1,])

for(i in 1:K){
idx<-which(Valid_Seq[,i]==4)
Valid_Seq[idx,i]=0
Valid_Seq[,i]<-as.factor(Valid_Seq[,i])
}



#Generate predictions across time,phase sequences, and states. 
N<-length(which(panel$State=="Alabama" & panel$Lconstraint>0))

M<-length(Valid_Seq[,1])
beta_state<-matrix(0,M,K)
L_state<-matrix(0,M,K)

beta_predictions<-array(0,dim=c(M,K,length(st_list)))
L_predictions<-array(0,dim=c(M,K,length(st_list)))

for(j in 1:length(st_list)){
  print(st_list[j])
  for(i in 1:K){
    CF<-data.frame(matrix(0,nrow=N,ncol=3))
    colnames(CF)<-c("time","State","Phase")

    CF$time<-seq(1:N)
    CF$time<-as.factor(CF$time)
    CF$State=st_list[j]

    CF$Phase=as.factor(Valid_Seq[,i])

    beta_hat<-predict(beta_reg,CF)
    beta_hat<-pmax(beta_hat,0)
    beta_state[,i]<-beta_hat

    L_state[,i]<-predict(L_reg,CF)}

  beta_predictions[,,j]<-beta_state
  L_predictions[,,j]<-L_state
}

#Generate Counterfactual Pandemic Dynamics
#Determine starting infected values
  #Do this later --sounds boring

D_predictions<-array(0,dim=c(M+1,K,length(st_list)))
I_predictions<-array(0,dim=c(M+1,K,length(st_list)))
S_predictions<-array(0,dim=c(M+1,K,length(st_list)))
R_predictions<-array(0,dim=c(M+1,K,length(st_list)))
C_predictions<-array(0,dim=c(M+1,K,length(st_list)))
d_predictions<-array(0,dim=c(M+1,K,length(st_list)))


for(l in 1:length(st_list)){
  for(j in 1:K){

    I=1
    S=1000
    b1=beta_predictions[,j,l]

    gamma=.20
    theta=.10
    delta=.01
    T=length(b1)+1

    t=matrix(0,T,1)
    S_vec=matrix(0,T,1)
    S_vec[1]=S
  
    N=S+I
  
    I_vec=matrix(0,T,1)
    I_vec[1]=I
    R_vec=matrix(0,T,1)
    D_vec=matrix(0,T,1)
    C_vec=matrix(0,T,1)
    d_vec=matrix(0,T,1)
    
    for(i in 2:length(S_vec)){
      S_vec[i]=max(S_vec[i-1]-(b1[i-1]/N)*S_vec[i-1]*I_vec[i-1],0)
      I_vec[i]=min(max(I_vec[i-1]+(b1[i-1]/N)*S_vec[i-1]*I_vec[i-1]-gamma*I_vec[i-1],0),S)
      R_vec[i]=min(max(R_vec[i-1]+gamma*I_vec[i-1]-theta*R_vec[i-1],0),S)
      D_vec[i]=min(D_vec[i-1]+delta*theta*R_vec[i-1],S)
      C_vec[i]=min(C_vec[i-1]+(1-delta)*theta*R_vec[i-1],S)
      d_vec[i]=D_vec[i]-D_vec[i-1]}

    S_predictions[,j,l]=S_vec
    I_predictions[,j,l]=I_vec
    R_predictions[,j,l]=R_vec
    C_predictions[,j,l]=C_vec
    D_predictions[,j,l]=D_vec*7 #Correct adjustment for weekly data?
    d_predictions[,j,l]=d_vec*7  
  }
}

#Extract each state's observed lockdown sequence
idx=which(panel$State==st_list[27])
Observed_Phase=panel$Phase[idx]
Observed_Phase=as.numeric(Observed_Phase)-1

which(Observed_Phase)

L_alabama<-L_predictions[,,27]
d_alabama<-d_predictions[,,27]
d_alabama<-d_alabama[-1,]
eta=.9999

for(i in 1:13){
  L_alabama[i,]=eta^(i-1)*L_alabama[i,]
  d_alabama[i,]=eta^(i-1)*d_alabama[i,]*4900}

w=.25

V_L=colSums(L_alabama^2)
V_d=w*colSums(d_alabama^2)

loss=-V_L+V_d
min_idx=which.min(loss)
print(Valid_Seq[,min_idx])
diff=sum(abs(Observed_Phase-Valid_Seq[,min_idx]))/13
#diff=sum(((Observed_Phase-Valid_Seq[,min_idx])^2))/13

#521
for (i in 1:length(Valid_Seq[1,])){
  dec=i
  log=(Observed_Phase==Valid_Seq[,i])
  if(all(log)==TRUE){break}
}

dec



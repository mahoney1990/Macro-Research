library(ggplot2)

#Infection Dynamics -- No Feedback
S=1000
I=1
b1=.001
b2=.5

t=matrix(0,100,1)
S_vec=matrix(0,100,1)
S_vec[1]=S

I_vec=matrix(0,100,1)
I_vec[1]=I

R_vec=matrix(0,100,1)

for(i in 2:length(S_vec)){
  S_vec[i]=max(S_vec[i-1]-b1*S_vec[i-1]*I_vec[i-1],0)
  I_vec[i]=min(max(I_vec[i-1]+b1*S_vec[i-1]*I_vec[i-1]-b2*I_vec[i-1],0),S)
  R_vec[i]=min(max(R_vec[i-1]+b2*I_vec[i-1],0),S)
}

plot_data=data.frame("I" = I_vec,"S" = S_vec, "R"=R_vec,"t"=seq(1:100))


ggplot() + 
  geom_line(data = plot_data, aes(x = t, y = I), color = "red") +
  geom_line(data = plot_data, aes(x = t, y = R), color = "blue") +
  geom_line(data = plot_data, aes(x = t, y = S), color = "green") +
  xlab('Time') +
  ylab('Number')

#Infection Dynamics -- With Feedback
S=1000
I=1
b1=.001
b2=.1
alpha=.6
t=matrix(0,100,1)
S_vec=matrix(0,100,1)
S_vec[1]=S

I_vec=matrix(0,100,1)
I_vec[1]=I

R_vec=matrix(0,100,1)
Y_vec=matrix(0,100,1)
A=1
b_mid=.5*b1
b_high=0
z1=.1*S
z2=.05*S
k1=1
k2=.5
for(i in 2:length(S_vec)){
  if(I_vec[i-1]<z1){
  S_vec[i]=max(S_vec[i-1]-b1*S_vec[i-1]*I_vec[i-1],0)
  I_vec[i]=min(max(I_vec[i-1]+b1*S_vec[i-1]*I_vec[i-1]-b2*I_vec[i-1],0),S)
  R_vec[i]=min(max(R_vec[i-1]+b2*I_vec[i-1],0),S)
  K=k1
  Y_vec[i]=A*(K^alpha)
  w_vec[i]=alpha*A*K^(alpha-1)}
  else{print(i)
    S_vec[i]=max(S_vec[i-1]-b_high*S_vec[i-1]*I_vec[i-1],0)
    I_vec[i]=min(max(I_vec[i-1]+b_high*S_vec[i-1]*I_vec[i-1]-b2*I_vec[i-1],0),S)
    R_vec[i]=min(max(R_vec[i-1]+b2*I_vec[i-1],0),S)
    K=k2
    Y_vec[i]=A*(K^alpha)
    w_vec[i]=alpha*A*K^(alpha-1)}}

plot_data=data.frame("I" = I_vec,"S" = S_vec, "R"=R_vec,"t"=seq(1:100))


P2<-ggplot() + 
  geom_line(data = plot_data, aes(x = t, y = I), color = "red") +
  geom_line(data = plot_data, aes(x = t, y = R), color = "blue") +
  geom_line(data = plot_data, aes(x = t, y = S), color = "green") +
  xlab('Time') +
  ylab('Number')

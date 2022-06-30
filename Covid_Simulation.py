

# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 10:38:03 2021

@author: mahon
"""

import os
import numpy as np
from quantecon.optimize import brent_max, brentq
from interpolation import interp
import matplotlib.pyplot as plt

from quantecon import MarkovChain
from scipy.optimize import minimize_scalar, bisect
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import random
from numpy.random import binomial
from numpy.random import multinomial
import pandas as pd
from numpy import asarray
from numpy import savetxt
from numpy import loadtxt

class IFP:
    
    def __init__(self,
                 r=0.0001,
                 beta=.9999,
                 gamma=.5,
                 grid_max=10,
                 grid_size=501,
                 uy_I=-.50,
                 u_H=0,
                 uo_I=-1,
                 uo_D=-5,
                 uy_D=-5):

        self.R=(1+r)
        self.beta, self.gamma = beta, gamma
        self.uo_I,self.uy_I,self.u_H,self.uo_D,self.uy_D,=uo_I,uy_I,u_H,uo_D,uy_D
        self.asset_grid = np.linspace(0,grid_max,grid_size)
        
        assert self.R * self.beta < 1, "Model Unstable"
        
    def u_y(self,c,sick,dead):
        return ((c)**(1-self.gamma))/(1-self.gamma) + self.uy_I*sick + self.u_H*(1-sick)+self.uy_D*dead

    
    
    def u_o(self,c,sick,dead):
        return ((c)**(1-self.gamma))/(1-self.gamma) + self.uo_I*sick + self.u_H*(1-sick)+self.uo_D*dead


    def u_prime(self,c):
        return (c+.00000000001)**(-self.gamma)
    
def euler_diff(c,a,z,sigma_vals,ifp):
    """Parameters:
    ----------
    c : consumption choice
    a : asset holdings -- given at period t
    z : stochastic state -- high or low income
    sigma_vals : policy function represented by a matrix
    ifp : instance of IFP class
    Returns : diff """
    
    R, beta, gamma = ifp.R, ifp.beta, ifp.gamma
    agrid, u_prime = ifp.asset_grid, ifp.u_prime
    
    #Interpolate consumptuion function sigma    
    def sigma(a,z):
        return interp(agrid, sigma_vals[:, z], a)
    
    expect=0.0
    
    
    for z_hat in range(n):
        expect += u_prime(sigma(R*(a)-c+y[z_hat],z_hat)) * P[z,z_hat]

    return (u_prime(c)-max(beta*R*expect, u_prime(a)))**2



def K(sigma,ifp):
    """
    Iterate on objective function for all a,z values
    Parameters:
    ----------
    sigma : policy function represented by a matrix
    ifp : instance of IFP class
    Returns : sigma fn approximation """
    
    sigma_new=np.empty_like(sigma)
    for i,a in enumerate(ifp.asset_grid):
        for z in range(n):
            result = minimize_scalar(euler_diff, 
                             bounds=(0,a+y[z]),
                             args=(a,z,sigma,ifp),
                             method='bounded')
            sigma_new[i,z]=result.x
    return sigma_new 

print_skip=10

def solve_model_time_iter(model,sigma,tol=.001,max_iter=1000,verbose=True, print_skip=10):
     error = tol +1
     i=0
     while i<max_iter and error>tol:
         sigma_new=K(sigma,model)
         error=np.max(np.abs(sigma-sigma_new))
         i += 1
         if verbose and i % print_skip ==0:
             print(f"Error at Iteration #"+str(i)+":"+str(error))
         sigma=sigma_new
         
     if i == max_iter:
        print("Convergence Failed, check functions.")
        
     if verbose and i < max_iter:
        print("Convergence Achieved in "+str(i)+" iterations.")
        
     return sigma
util_dict={}
util_dict_old={}

#%%
#Current run -- young and old, baseline, lower penalties, fixed asset dist
ifp=IFP()
os.chdir(r'C:\Users\mahon\Documents\Python Scripts')
mdata=pd.read_csv("State_Markov_Sequences_v2.csv")
ddata=pd.read_csv("Old_pop.csv")
st_list=pd.read_csv("st_list.csv")
st_list=list(st_list["x"])


T=32
burn_in=3

c_mat_ag=np.zeros((T,39))
s_mat_ag=np.zeros((T,39))
a_mat_ag=np.zeros((T,39))
U_ag_mat=np.zeros((T,39))
consutil_mat_ag=np.zeros((T,39))

it=0
run_type="Young"
exercise="State"




for it in range(39):
    policy={}
    
    time=list(range(T))
    
    state=str(st_list[it])
    print(state)
    
    #Define policy sequences for young households
    
    L_prob=np.array(mdata[state+"_L_prob"])
    L_prob1=L_prob+(1-L_prob)/2
    L_prob2=L_prob-(1-L_prob)/2
    U_prob=np.array(mdata[state+"_U_prob"])
    I_prob=np.array(mdata[state+"_I_prob"])
    
    I_term=I_prob[15]
    I_seq=(np.linspace(I_term,0,T-len(I_prob)-burn_in))
    H_seq=1-I_seq
    
    P_death_young=.001
    D_seq=np.linspace(P_death_young,0,T-len(I_prob)-burn_in)
    
    P_recover_young=.75
    R_seq=np.linspace(P_recover_young,1,T-len(I_prob)-burn_in)
    
    L1_term=L_prob[15]
    L2_term=L_prob2[15]
    
    L1_seq=(np.linspace(L1_term,1,T-len(I_prob)-burn_in))
    L2_seq=(np.linspace(L2_term,1,T-len(I_prob)-burn_in))
    
    
    ad=burn_in+16
    #1-L_prob[i-burn_in]-I_prob[i-burn_in]
    for i in range(T):
        if i < burn_in:
            policy[i]=((.99,.01,0,0),(.98,.02,0,0),(1,0,0,0),(1,0,0,0))
        elif i < ad:
            policy[i]=((L_prob[i-burn_in],max(1-L_prob[i-burn_in]-I_prob[i-burn_in],0),I_prob[i-burn_in],0),
                       (L_prob2[i-burn_in],max(1-L_prob2[i-burn_in]-I_prob[i-burn_in],0),I_prob[i-burn_in],0),
                       (0,P_recover_young,1-P_death_young-P_recover_young,P_death_young),
                       (0,0,0,1))
            
        else:
            policy[i]=((L1_seq[i-ad],max(1-I_seq[i-ad]-L1_seq[i-ad],0),I_seq[i-ad],0),
                       (L2_seq[i-ad],max(1-L2_seq[i-ad]-I_seq[i-ad],0),I_seq[i-ad],0),
                       (0,P_recover_young,1-P_death_young-P_recover_young,P_death_young),
                       (0,0,0,1))    
    
    U=list(range(T))
    policy_seq=list(range(T))
    agrid=ifp.asset_grid
    a_size=len(agrid)


    z_size=len(policy[0])
    n=len(policy[0])

    #Initializing guess for consumption fn sigma -- consume all savings
    sigma_init=np.repeat(agrid.reshape(a_size,1), z_size, axis=1 )
    sigma_vals=sigma_init

    y=[1,.25,.25,0]
    n_periods=len(policy_seq)
    z_int=0
    z_vec=np.empty_like(policy_seq)
    
    #Simulate Household sequence -- young
    
    
    sigma_dict={}
    
 
    for k in range(T):
        print(k)
        P=np.array(policy[k])
        if k==0:
            sigma_dict[k]=solve_model_time_iter(ifp, sigma_init,tol=.001)
        else:
            sigma_dict[k]=solve_model_time_iter(ifp, sigma_dict[k-1],tol=.001)
    
    random.seed(1337)
    n_sims=100000
    n_dist=5
    step=int(n_sims/n_dist)
        
    import numpy
    from scipy import interpolate
    
    
    a_list=numpy.zeros((T,n_sims))
    c_list=numpy.zeros((T,n_sims))
    u_list=numpy.zeros((T,n_sims))
    consutil_list=numpy.zeros((T,n_sims))
    z_dict={}
    
    for j in range(n_dist):
        a_int=agrid[50+50*j]
        a_list[0,(j*(step)):((step)*(j+1)-1)]=a_int
    
    a_list[0,0:50000]=1
    a_list[0,50001:60000]=1.5
    a_list[0,60001:70000]=2
    a_list[0,70001:80000]=3
    a_list[0,80001:90000]=4     
    a_list[0,90001:95000]=5
    a_list[0,95001:99999]=10    
    
    print_skip=1000
    
    #Load pre calculated draws
    #z_mat=loadtxt(run_type+"_"+state+"_draws.csv" ,delimiter=",")
    
    for j in range(n_sims):
        
        if j % print_skip ==0:
            print(f"Simulating Household Responses. Progress:"+str(j/n_sims*100)+"%")
             
        z_vec=numpy.zeros((T,1))
        z_vec=z_vec.astype(int)
    
        
        for k in range(n_periods):
            if k==0:
                 z_vec[k]=z_int
            else:
                z_vec[k]=z_new
                
            pol=policy_seq[k]
            P=np.array(policy[pol])    
            
            roll=multinomial(1,[P[z_vec[k][0],0],P[z_vec[k][0],1],P[z_vec[k][0],2],P[z_vec[k][0],3]])
            z_new=np.where(roll==1)[0][0]
            
       
        z_dict[j]=z_vec
        
        Util=np.empty(len(policy_seq))
        cons=np.empty(len(policy_seq))
        income=np.empty(len(policy_seq))
        a=np.empty(len(policy_seq))
        
        for k in range(T-1):
            sigma=sigma_dict[policy_seq[k]]
            z=z_vec[k][0]
            f=interp1d(agrid, sigma[:, z],kind='quadratic',fill_value="extrapolate")
            c_list[k,j]=max(f(a_list[k,j]),0.0000000000000001)
            consutil_list[k,j]=(ifp.beta**k)*(1/ifp.gamma)*c_list[k,j]**(1-ifp.gamma)
            income[k]=y[z]
            if z==2:
                sick=1
                dead=0
            elif z == 3:
                dead=1
                sick=0
            else:
                sick=0
                dead=0
            u_list[k,j]=ifp.u_y(c_list[k,j],sick,dead)
            if k<T:
                a_list[k+1,j]=a_list[k,j]+income[k]-c_list[k,j]
    
    c_mat_ag[:,it]=np.sum(c_list,axis=1)/n_sims
    a_mat_ag[:,it]=np.sum(a_list,axis=1)/n_sims
    consutil_mat_ag[:,it]=np.sum(consutil_list,axis=1)/n_sims
    
    for j in range(T-1):
        if j<T:
            s_mat_ag[j,it]=a_mat_ag[j+1,it]-a_mat_ag[j,it]
    
    
    utility=0
    
    
    E_util=np.zeros(T)
    E_cons=np.zeros(T)
    E_assets=np.zeros(T-1)
    
    #AGGREGATION
    for j in range(T-1):
        utility=0
        consum=0
        Assets=0
        for i in range(n_sims):
            if(z_dict[i][j][0]==2):
                sick=1
                dead=0
            elif(z_dict[i][j][0]==3):
                sick=0
                dead=1
            else:
                sick=0
                dead=0
            if (c_list[j,i]==0):
                break
            utility+=(ifp.beta**j)*ifp.u_y(c_list[j,i], sick,dead)
            consum+=c_list[j,i]/n_sims
            Assets+=a_list[j,i]/n_sims
        
        E_util[j]=utility/n_sims
        E_cons[j]=consum
        E_assets[j]=Assets
        U_ag_mat[j,it]=E_util[j]
    
    
    savings=E_assets[2:(T-1)]-E_assets[1:(T-2)]
    t = np.linspace(1, len(savings), len(savings))
    
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(t,savings,"-b", label="Savings") 
    plt.plot(t,E_cons[1:30],"-r", label="Consumption")
    plt.legend(loc="upper left")
    plt.title(str(state)+" Aggregated "+"Savings/Consumption")
    plt.savefig(str(st_list[it])+'_young_v2.png')
    
    
        
    
        
    
    print(str(state)+"Expected Utility:")
    print(sum(E_util))    
    
    util_dict[it]=E_util    
    
    pd.DataFrame(consutil_mat_ag).to_csv("consutil_mat_state.csv")
    pd.DataFrame(c_mat_ag).to_csv("c_mat_state.csv")
    pd.DataFrame(s_mat_ag).to_csv("s_mat_ntl_state.csv")
    pd.DataFrame(a_mat_ag).to_csv("a_mat_ntl_state.csv")
    
Exp_Utilit=np.zeros(39)
for i in range(39):
    Exp_Utilit[i]=sum(util_dict[i])
pd.DataFrame(Exp_Utilit).to_csv("Young_Util_state.csv")
with open(run_type+"_"+exercise+"_util_state.csv", 'w') as f:
    for key in util_dict.keys():
        f.write("%s,%s\n"%(key,util_dict[key]))




####OLD RUN
ddata=pd.read_csv("Old_pop.csv")
st_list=pd.read_csv("st_list.csv")
st_list=list(st_list["x"])
sigma_save={}


run_type="Old"
exercise="National"
it=0
T=32
burn_in=3

time=list(range(T))

c_old_mat_ag=np.zeros((T,39))
s_old_mat_ag=np.zeros((T,39))
a_old_mat_ag=np.zeros((T,39))
U_old_ag_mat=np.zeros((T,39))
ConsUtil_old_mat_ag=np.zeros((T,39))


#Current run -- old, national, low penalties
for it in range(39):
    policy={}

    
    
    
    state=str(st_list[it])
    print(state)
    
    #Define policy sequences for young households
    
    L_prob=np.array(mdata[state+"_L_prob"])
    L_prob1=L_prob+(1-L_prob)/2
    L_prob2=L_prob-(1-L_prob)/2
    U_prob=np.array(mdata[state+"_U_prob"])
    I_prob=np.array(mdata[state+"_I_prob"])
    
    I_term=I_prob[15]
    I_seq=np.linspace(I_term,0,T-len(I_prob)-burn_in)/2
    H_seq=1-I_seq
    
    P_death_old=.025
    D_seq=np.linspace(P_death_old,0,T-len(I_prob)-burn_in)
    
    P_recover_old=.20
    R_seq=np.linspace(P_recover_old,1,T-len(I_prob)-burn_in)
    
    L1_term=L_prob[15]
    L2_term=L_prob2[15]
    
    L1_seq=np.sqrt(np.linspace(L1_term,1,T-len(I_prob)-burn_in))
    L2_seq=np.sqrt(np.linspace(L2_term,1,T-len(I_prob)-burn_in))
    
    
    ad=burn_in+16
    for i in range(T):
        if i < burn_in:
            policy[i]=((1,0,0),(1,0,0),(1,0,0))
        elif i < ad:
            policy[i]=((1-I_prob[i-burn_in],I_prob[i-burn_in],0),
                       (P_recover_old,1-P_death_old-P_recover_old,P_death_old),
                       (0,0,1))
            
        else:
            policy[i]=((1-I_seq[i-ad],I_seq[i-ad],0),
                       (P_recover_old,1-P_death_old-P_recover_old,P_death_old),
                       (0,0,1))    
    
    U=list(range(T))
    policy_seq=list(range(T))
    agrid=ifp.asset_grid
    a_size=len(agrid)


    z_size=len(policy[0])
    n=len(policy[0])

    #Initializing guess for consumption fn sigma -- consume all savings
    sigma_init=np.repeat(agrid.reshape(a_size,1), z_size, axis=1 )
    sigma_vals=sigma_init

    y=[.25,.25,0]
    n_periods=len(policy_seq)
    z_int=0
    z_vec=np.empty_like(policy_seq)
    
    #Simulate Household sequence -- young
    
    
    sigma_dict={}
    

    for k in range(T):
        P=np.array(policy[k])
        if k==0:
            sigma_dict[k]=solve_model_time_iter(ifp, sigma_init,tol=.001)
        else:
            sigma_dict[k]=solve_model_time_iter(ifp, sigma_dict[k-1],tol=.001)
    
    sigma_save[it]=sigma_dict
    
    
    random.seed(1337)
    n_sims=100000
    n_dist=5
    step=int(n_sims/n_dist)
        
    import numpy
    from scipy import interpolate
    
    
    a_list=numpy.zeros((T,n_sims))
    c_list=numpy.zeros((T,n_sims))
    u_list=numpy.zeros((T,n_sims))
    consutil_list=numpy.zeros((T,n_sims))
    z_dict={}
    

    a_list[0,0:50000]=3
    a_list[0,50001:70000]=4    
    a_list[0,70001:99999]=5  
    print_skip=2000
    

    
    for j in range(n_sims):
        
        if j % print_skip ==0:
            print(f"Simulating Household Responses. Progress:"+str(j/n_sims*100)+"%")
             
        z_vec=numpy.zeros((T,1))
        z_vec=z_vec.astype(int)
    
        
        for k in range(n_periods):
            if k==0:
                z_vec[k]=z_int
            else:
                z_vec[k]=z_new
                
            pol=policy_seq[k]
            P=np.array(policy[pol])    
            
            roll=multinomial(1,[P[z_vec[k][0],0],P[z_vec[k][0],1],P[z_vec[k][0],2]])
            z_new=np.where(roll==1)[0][0]
            
       
        z_dict[j]=z_vec
        
        Util=np.empty(len(policy_seq))
        cons=np.empty(len(policy_seq))
        income=np.empty(len(policy_seq))
        a=np.empty(len(policy_seq))
        
        for k in range(T-1):
            sigma=sigma_dict[policy_seq[k]]
            z=z_vec[k][0]
            f=interp1d(agrid, sigma[:, z],kind='quadratic',fill_value="extrapolate")
            c_list[k,j]=f(a_list[k,j])
            consutil_list[k,j]=(ifp.beta**k)*(1/ifp.gamma)*c_list[k,j]**(1-ifp.gamma)
            
            income[k]=y[z]
            if z==2:
                sick=1
                dead=0
            elif z==3:
                sick=0
                dead=1
            else:
                dead=0
                sick=0
            u_list[k,j]=ifp.u_o(c_list[k,j],sick,dead)
            if k<T:
                a_list[k+1,j]=a_list[k,j]+income[k]-c_list[k,j]
    
    ConsUtil_old_mat_ag[:,it]=np.sum(u_list,axis=1)/n_sims
    c_old_mat_ag[:,it]=np.sum(c_list,axis=1)/n_sims
    a_old_mat_ag[:,it]=np.sum(a_list,axis=1)/n_sims
    
    for j in range(T-1):
        if j<T:
            s_old_mat_ag[j,it]=a_old_mat_ag[j+1,it]-a_old_mat_ag[j,it]
    
        
    utility=0
    E_util=np.zeros(T)
    E_cons=np.zeros(T)
    E_assets=np.zeros(T)
    
    print("Aggregating...")
    E_util=np.zeros(T)
    for j in range(0,T-1):
        utility=0
        consum=0
        Assets=0
        for i in range(n_sims):
            if(z_dict[i][j]==1):
                sick=1
                dead=0
            elif(z_dict[i][j]==2):
                sick=0
                dead=1
            else:
                sick=0
                dead=0
            if (c_list[j,i]==0):
                break
            utility+=(ifp.beta**j)*ifp.u_o(c_list[j,i], sick,dead)
            consum+=c_list[j,i]/n_sims
            Assets+=a_list[j,i]/n_sims
        
        E_util[j]=utility/n_sims
        E_cons[j]=consum
        print(str(E_cons[j]))
        E_assets[j]=Assets
        U_old_ag_mat[j,it]=E_util[j]
        
        
    savings=E_assets[2:(T-1)]-E_assets[1:(T-2)]
    t = np.linspace(1, len(savings), len(savings))
    
    fig = plt.figure()
    ax = plt.axes()
    plt.plot(t,savings,"-b", label="Savings") 
    plt.plot(t,E_cons[2:31],"-r", label="Consumption")
    plt.legend(loc="upper left")
    plt.title(str(state)+" Aggregated "+"Savings/Consumption")
    plt.savefig(str(st_list[it])+'_old_state_v2.png')
    
    
    print(str(state)+"Expected Utility:")
    print(sum(E_util))    
    
    util_dict_old[it]=E_util 
    pd.DataFrame(c_old_mat_ag).to_csv("c_old_mat_state.csv")
    pd.DataFrame(s_old_mat_ag).to_csv("s_old_mat_state.csv")
    pd.DataFrame(ConsUtil_old_mat_ag).to_csv("consutil_old_state.csv")

Exp_Utilit_Old=np.zeros(39)
for i in range(39):
    Exp_Utilit_Old[i]=sum(util_dict_old[i])

pd.DataFrame(Exp_Utilit_Old).to_csv("Old_Baseline_Util_state.csv")
with open(run_type+"_"+exercise+'_util_state.csv', 'w') as f:
    for key in util_dict_old.keys():
        f.write("%s,%s\n"%(key,util_dict_old[key]))


age_data=pd.read_csv('age_data.csv')
old_prop=list(age_data['Weight'])

baseline_rf=56.211

Cons_Util=[]
Final_Util=[]
difference=[]
CEU=[]
for i in range(39):
    Cons_Util.append((1-.165)*sum(consutil_mat_ag[:,i])+.165*sum(ConsUtil_old_mat_ag[:,i]))
    Final_Util.append((1-.165)*Exp_Utilit[i]+.165*Exp_Utilit_Old[i])
    difference.append(Final_Util[i]-Cons_Util[i])
    CEU.append(((baseline_rf-difference[i])/Cons_Util[i])**(1/(1-ifp.gamma))-1)







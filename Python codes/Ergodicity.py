# -*- coding: utf-8 -*-
"""
Created on Wed Nov 01 16:17:48 2017

@author: Swaroop
"""

import numpy as np
N=1000
X=np.zeros((N,N),dtype=float)
E_Ensemble=np.zeros(N,dtype=float)
E_time=np.zeros(N,dtype=float)
Error=np.zeros(N,dtype=float)
for i in range(N):
    
        X[i]=np.random.randn(N)
        E_time[i]=np.sum(X[i])/float(N)
        
#print(X)
print(E_time)        

for j in range(N):
    sum=0.0
    for k in range(N):
        
        sum = sum + X[k][j]
        E_Ensemble[j]=sum/float(N) 
    
print(E_Ensemble)

#Error=(E_time - E_Ensemble).sum()

#print(Error)    
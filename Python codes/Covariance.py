# -*- coding: utf-8 -*-
"""
Created on Wed Nov 01 17:02:03 2017

@author: Swaroop
"""

import numpy as np
N=1000
X1=np.random.randn(N)
X2=np.random.randn(N)
X=np.matrix(np.reshape(X1,(1,-1)))
Y=np.matrix(np.reshape(X2,(-1,1)))
Z=np.zeros(N)
print(np.shape(X),np.shape(Y))

for i in range(N):
    Z[i]=Y[i,0]*X[0,i]
print(len(Z))    
COV=np.sum(Z)/float(N) 
print(COV)   
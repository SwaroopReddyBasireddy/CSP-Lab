# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 11:28:58 2017

@author: Swaroop
"""

import numpy as np

h=([1+1j,2+2j,3+3j,4+4j,5+5j])
N=100
#x=np.zeros((1,5),dtype=complex)
x=([2+2j,3+3j,4+4j,5+5j,6+6j])
#x=(np.random.normal(0,1,N)+1j*np.random.normal(0,1,N))
m=len(x)
n=len(h)
a=np.zeros((n,m),dtype=complex)

for i in range(n):
       for j in range(m):
           a[i][j] += h[i]*x[j]
           
#print(a) 
y=np.zeros(n+m-1,dtype=complex) 

for k in range(n+m):
    
    for i in range(n):
        for j in range(m):
            if (i+j)==k:
                y[k] += a[i][j]
                
print(y,np.shape(h),np.shape(y))                
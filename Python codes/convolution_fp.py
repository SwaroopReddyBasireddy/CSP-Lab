# -*- coding: utf-8 -*-
"""
Created on Fri Nov 10 09:50:19 2017

@author: Swaroop
"""

import numpy as np


h=([1,0.8,0.2])
x=([-1,0.01,0.99])
h_fp=np.zeros(len(h),dtype=np.int16)
x_fp=np.zeros(len(x),dtype=np.int16)
m=len(x)
n=len(h)
a=np.zeros((n,m))
a_fp=np.zeros((n,m))

for i in range(n):
    h_fp[i]= np.int16(h[i]*(2**14))
for j in range(m):
    
    x_fp[j]=np.int16(x[j]*(2**13))

#print(h_fp,np.size(h_fp))

for i in range(n):
       for j in range(m):
           a_fp[i][j] += h_fp[i]*x_fp[j]
           a[i][j] += h[i]*x[j]           
#print(a) 
y=np.zeros(n+m-1)
y_fp=np.zeros(n+m-1,dtype=np.int16)
y1=np.zeros(n+m-1,dtype=float) 

for k in range(n+m):
    
    for i in range(n):
        for j in range(m):
            if (i+j)==k:
                y[k] += a[i][j]
                y_fp[k] += a_fp[i][j]
                y1[k] =y_fp[k]/(2**13)
                
print(y)
print(y1)                
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 12:00:31 2017

@author: Swaroop
"""
import numpy as np

h=([1,0.8,0.2])
x=([-1,0.01,0.99])

m=len(x)
n=len(h)
a=np.zeros((n,m))

for i in range(n):
       for j in range(m):
           a[i][j] += h[i]*x[j]
           
print(a) 
y=np.zeros(n+m-1) 

for k in range(n+m):
    
    for i in range(n):
        for j in range(m):
            if (i+j)==k:
                y[k] += a[i][j]
                
print(y)            
N=3
z=np.zeros(N)
if N<n+m-1:
   for i in range(n+m-1):
       if i<N:
           z[i]=y[i]
       else:
           z[i-N]+=y[i]
else:
  for i in range(N):
      if i<(n+m-1):
          z[i]=y[i]
      else:
          z[i]=0
          
print(z)          
      

          
           
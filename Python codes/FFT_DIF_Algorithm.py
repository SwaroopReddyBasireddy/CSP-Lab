# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 10:20:51 2017

@author: Swaroop
"""

import numpy as np
import math

def fft_stage(N,stage,x):
   w=twiddle_factor(N/2**(stage-1)) 
   #print stage
   y=np.zeros(N,dtype=complex)
   for j in range(2**(stage-1)):
       
     for i in range(0,N/(2**stage),1):
        m=2*j*N/(2**(stage)) 
        #m=N/(2**(stage-1))
        #print(j)
        y[i+m]=x[i+m]+x[i+m+N/(2**stage)]
        if (2**stage) != N:
           y[i+m+N/(2**stage)]=(x[i+m]-x[i+m+N/(2**stage)])*w[i] 
        else:
           y[i+m+N/(2**stage)]=(x[i+m]-x[i+m+N/(2**stage)])
        
   return(y)    
def twiddle_factor(N):
    w=np.zeros(N,dtype=complex)
    for i in range(N):
        c=2*np.pi/N
        w[i]=np.round((np.cos(c*i)+1j*np.sin(c*i)),3)
    return w 

def swap(x,N):
    for i in range(1,N/2,2):
        c=x[i]
        x[i]=x[i-1+N/2]
        x[i-1+N/2]=c
    return x

    
N=8

no_stages=int(math.log(N,2))

x=([1,2,3,4,5,6,7,8])
z=np.zeros((no_stages,N),dtype=complex)

#stage=np.linspace(1,no_stages,1)
for i in range(no_stages):
    if i==0:
        z[i]=fft_stage(N,i+1,x)
    else:    
        
      z[i]=fft_stage(N,i+1,z[i-1])
      
    
k=swap(z[i],N)            

#X=z[no_stages]

for i in range(no_stages):
    
    print(z[i])
    
print(np.fft.fft(x))

print(k)   
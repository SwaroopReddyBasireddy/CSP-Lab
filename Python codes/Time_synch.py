# -*- coding: utf-8 -*-
"""
Created on Wed Nov 08 15:25:57 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt

kk=np.sqrt(0.5)*np.array(([0, 0, 1+1j, 0, 0, 0, -1+1j, 0, 0, 0, 1+1j, 0, 0, 0, -1-1j, 0, 0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, 0,0, 0, 0, -1-1j, 0, 0, 0, -1-1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0, 0, 0, 1+1j, 0,0]))
x=len(kk)
a=kk[len(kk)-27:len(kk)-1]
b=np.zeros(11)
c=kk[0:26]

k1=np.concatenate(([0],a),axis=0)
k2=np.concatenate((k1,b),axis=0)
kk1=np.concatenate((k2,c),axis=0)
#print(len(kk1))
kk2=np.sqrt(64)*np.fft.ifft(kk1)
B1=kk2[len(kk2)-10:len(kk2)]
Tx=np.concatenate((B1,kk2))
B2=Tx[len(B1):len(B1)+16]
#print(Tx)
#print(len(Tx))

snrlen=15
EbNodB=np.arange(0,snrlen,1)
P=np.zeros(snrlen,dtype=float)
Y=np.zeros(len(B2),dtype=complex)
for i in range(snrlen):
    N=10000
    count=0
    for j in range(N):
        h=(np.random.randn(1) + 1j*np.random.randn(1))/np.sqrt(2.0)
        sigma=10**(-EbNodB[i]/20.0)
        
        noise=(np.random.normal(0,sigma,len(Tx)) + 1j*np.random.normal(0,sigma,len(Tx)))/np.sqrt(2.0)
        y=np.convolve(Tx,h)+noise
        
        delta=3
        R=np.zeros(len(B2),dtype=float)
        r=np.zeros(len(B2),dtype=complex)
        
        for ii in range(len(B2)):
           y1=y[delta+ii:delta+ii+len(B2)]
           X=np.matrix(np.reshape(y1,(1,-1)))
           for k in range(len(B2)):
               Y[k]=np.conj(B2[k])
           Y1=np.matrix(np.reshape(Y,(-1,1)))
           for jj in range(len(B1)):
             r[ii] += Y1[jj,0]*X[0,jj]
           R[ii]=(np.abs(r[ii]))**2.0
        max_position=np.argmax(R)
        #print(max_position)
           
        if (delta+max_position == len(B1)):
             count +=1
    P[i]=count/float(N) 
        
print(P)        
        
plt.semilogy(EbNodB,P,'b') 
plt.title('channel synchronization')
plt.show()        
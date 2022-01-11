# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 02:37:18 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special


N=int(2e4)
snrlen=11
Eb_N0_dB=np.arange(0,snrlen,1)
#Eb_N0_dB=np.matrix((Eb_N0_dB))
s=np.zeros(int(N/2),np.complex)
noise=np.zeros(int(N/2),np.complex)
y=np.zeros(int(N/2),np.complex)
nErr_soft=np.zeros((1,snrlen))
simBer_soft=np.zeros((1,snrlen))


for i in range(snrlen):
    
    ip=np.random.randint(2,size=N)
    ip=np.array(ip)
    ipM=np.matrix(np.reshape(ip,(-1,2)))
    #print((ipM))
    r=np.zeros((np.shape(ipM)))
    s=-2*ipM+1
    m=([1,1j])
    Sm=s*np.matrix.transpose(np.matrix(m))/np.sqrt(2.0)
    #print(Sm)
    sigma=np.sqrt(10**(-Eb_N0_dB[i]/10.0))
    noise_real=np.random.normal(0,sigma,int(N/2))
    noise_img=np.random.normal(0,sigma,int(N/2))
    noise=np.vectorize(complex)(noise_real,noise_img)/np.sqrt(2.0)
    #n=np.matrix(np.reshape(np.matrix(noise),(-1,2)))
    
    y=np.matrix.transpose(Sm)+np.matrix(noise)
    #print(np.shape(Sm),np.shape(noise))        
    #print(np.shape(y)) 
    #print(y,np.shape(y),cmath.phase(y[0,2000])) 
           
            
    
    #y=np.array((y))
    for k in range(int(N/2)):
        if (np.real(y[0,k])>0.0 and np.imag(y[0,k])>0.0):
            r[k]=([0,0])
        elif (np.real(y[0,k])<0.0 and np.imag(y[0,k])>0.0): 
            r[k]=([1,0])
        elif (np.real(y[0,k])<0.0 and np.imag(y[0,k])<0.0):
             r[k]=([1,1])
        else:
           r[k]=([0,1])
           
           #print(r)  
    rM=np.matrix(np.reshape(r,(1,N)))
    print(np.size(rM))
    print(N)
    nErr_soft[0,i]=np.count_nonzero(ip-rM)
print(nErr_soft)

simBer_soft=nErr_soft[0]/float(N)

theoryBer=0.5*special.erfc(np.sqrt(0.5*(10**(0.1*Eb_N0_dB))))
print(Eb_N0_dB)
print(simBer_soft)
print(theoryBer)
print(np.shape(simBer_soft),np.shape(Eb_N0_dB),np.shape(theoryBer))
plt.plot(Eb_N0_dB,theoryBer,'r',Eb_N0_dB, simBer_soft,'b')
plt.legend(['theory','Practical'],loc=1)
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for QPSK in AWGN')
#plt.xtricks([0,1,2,3,4,5,6,7,8,9,10,11])
plt.grid()
plt.show()                
    
    

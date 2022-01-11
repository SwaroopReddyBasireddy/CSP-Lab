# -*- coding: utf-8 -*-
"""
Created on Sun Oct 08 03:03:27 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
N=int(1e6)
snrlen=16
Eb_N0_dB=np.arange(0,snrlen,1)
M=4.0

nErr_soft=np.zeros((1,snrlen))
simBer_soft=np.zeros((1,snrlen))


for i in range(snrlen):
    ip=np.random.randint(0,2,N)
    ipM=np.matrix(np.reshape(ip,(-1,2)))
    #print(ipM)
    s=-2*ipM+1
    #print(s)
    sM=s*np.matrix.transpose(np.matrix(([2,1])))/np.sqrt(5.0)
    #print(sM)
    
    sigma=np.sqrt(10**((-Eb_N0_dB[i])/10.0))
    noise_real=np.random.normal(0,sigma,N/2)
    noise_img=np.random.normal(0,sigma,N/2)
    noise=np.vectorize(complex)(noise_real,noise_img)/np.sqrt(2.0)
    y=np.matrix.transpose(sM) + np.matrix(noise)
    #r=np.matrix(r)
    r=np.zeros(np.shape(ipM))
    #print(np.shape(r),np.shape([0,1]))
    for k in range(N/2):
        if (np.real(y[0,k])>0.0 and np.real(y[0,k])<=2.0/np.sqrt(5.0)):
            r[k]=[0,1]
        elif (np.real(y[0,k])>2.0/np.sqrt(5.0)): 
            r[k]=[0,0]
        elif (np.real(y[0,k])> -2.0/np.sqrt(5.0) and np.real(y[0,k])<0.0):
             r[k]=[1,0]
        else:
           r[k]=[1,1]
           
    #print(r) 
    rM=np.matrix(np.reshape(r,(1,N)))
    nErr_soft[0,i]=np.count_nonzero(ipM-r)
    
simBer_soft = nErr_soft[0]/float(N/2)    
   

theoryBer=2*(M-1)*0.5*scipy.special.erfc(np.sqrt(0.2*(10**(0.1*Eb_N0_dB))))/M
#print(Eb_N0_dB)
#print(simBer_soft)
#print(theoryBer)
#print(np.shape(simBer_soft),np.shape(Eb_N0_dB),np.shape(theoryBer)) 
plt.plot(Eb_N0_dB,theoryBer,'r',Eb_N0_dB, simBer_soft,'b')
plt.legend(['theory','Practical'],loc=1)
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for 4-PAM in AWGN')
#plt.xtricks([0,1,2,3,4,5,6,7,8,9,10,11])
plt.grid()
plt.show()
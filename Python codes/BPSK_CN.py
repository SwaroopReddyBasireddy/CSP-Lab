# -*- coding: utf-8 -*-
"""
Created on Wed Oct 04 02:43:00 2017

@author: Swaroop
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Sep 11 04:53:28 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import special


N=int(1e6)
snrlen=11
Eb_N0_dB=np.arange(0,snrlen)
s=np.zeros(N)
r=np.zeros(N)
nErr_soft=np.zeros((1,snrlen))

for i in range(snrlen):
    ip=np.random.randint(2,size=N)
    s=2*ip-1
    sigma=np.sqrt((10**(-Eb_N0_dB[i]/10.0)))
    noise_real=np.random.normal(0,sigma,size=N)
    noise_img=np.random.normal(0,sigma,size=N)
    noise=(noise_real+1j*noise_img)/np.sqrt(2.0)
    y=s+noise
    print y
    cipSoftM=np.reshape(np.real(y),(1,N))
    for k in range(N):
        if y[k]>0:
            r[k]=1
        else:
            r[k]=0
            
    nErr_soft[0,i]=np.count_nonzero(ip-r)
theoryBer=0.5*special.erfc(np.sqrt((10**(Eb_N0_dB/10.0))))
#t h e o r e t i c a l be r uncoded AWGN
simBersoft=nErr_soft/float(N)
plt.plot(Eb_N0_dB,theoryBer,'r',Eb_N0_dB,simBersoft[0],'b')
plt.legend(['theory','Practical'],loc=1)
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for BPSK in AWGN')
#plt.xtricks([0,1,2,3,4,5,6,7,8,9,10,11])
plt.grid()
plt.show()

        
    
    

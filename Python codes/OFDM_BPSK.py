# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 02:14:56 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy

n=2**20
N=256
snrlen=11
Eb_N0_dB=np.arange(0,snrlen)
s_IFFT=np.zeros((n/N,N),dtype=complex)
y=np.zeros((n/N,N),dtype=complex)
r_FFT=np.zeros((n/N,N),dtype=complex)
nErr_soft=np.zeros((1,snrlen))
noise=np.zeros((n/N,N),dtype=complex)
rM=np.zeros((1,n))

for i in range(snrlen):
    
    ip=np.random.randint(2,size=n)
    ip=np.array(ip)
    s=2*ip-1
    ipM=np.matrix(np.reshape(s,(-1,N)))
    #print(np.shape((ipM)),np.shape((s_IFFT)))
    for j in range(n/N):
        s_IFFT[j:]=np.sqrt(float(N))*np.fft.ifft((ipM[j]))
    #print((ipM))
    #print((s_IFFT))
    sigma=np.sqrt((10**(-Eb_N0_dB[i]/10.0)))
    noise_real=np.random.normal(0,sigma,n)
    noise_img=np.random.normal(0,sigma,n)
    #noise=np.vectorize(complex)(noise_real,noise_img)/np.sqrt(2.0)
    noise=(noise_real+1j*noise_img)/np.sqrt(2.0)
    noise=np.array(noise)
    noise_matrix=np.matrix(np.reshape(noise,(-1,N)))
    #print((np.shape(n)))
    y=s_IFFT+noise_matrix
    #print(np.shape((y)),(y),np.shape((r_FFT)))
    for k in range(n/N):
        r_FFT[k:]=np.fft.fft(y[k])
    #print((r_FFT),np.shape((r_FFT)))
    r=np.matrix(np.reshape(np.real(r_FFT),(1,n)))
    #r=np.array(r)
    #print(r)
    for k in range(n):
        if (r[0,k]>0):
            rM[0,k]=1
        else:
            rM[0,k]=0
    #print(rM)
    nErr_soft[0,i]=np.count_nonzero(ip-rM[0])
print((nErr_soft[0]))    
theoryBer=0.5*scipy.special.erfc(np.sqrt((10**(Eb_N0_dB/10.0))))
#t h e o r e t i c a l be r uncoded AWGN
simBersoft=nErr_soft/float(n)
plt.plot(Eb_N0_dB,theoryBer,'r',Eb_N0_dB,simBersoft[0],'b')
plt.legend(['theory','Practical'],loc=1)
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for OFDM with BPSK Modulation  in AWGN')
#plt.xtricks([0,1,2,3,4,5,6,7,8,9,10,11])
plt.grid()
plt.show()

               
    
   
     
    
#print((ipM[5]),s_IFFT[5],np.shape((ipM[5])),np.shape((s_IFFT[5])))    
    
    
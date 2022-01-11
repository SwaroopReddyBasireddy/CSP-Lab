# -*- coding: utf-8 -*-
"""
Created on Fri Oct 20 09:29:38 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy

n=2**20
N=256
L=16
snrlen=11
Eb_N0_dB=np.arange(0,snrlen)
s_IFFT=np.zeros((n/N,N),dtype=complex)
s_OFDM=np.zeros((n/N,N+L-1),dtype=complex)
#y=np.zeros((n/N,N),dtype=complex)
y_Rx=np.zeros((n/N,N),dtype=complex)
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
        CP_add=s_IFFT[j,N-L+1:]
        s_CP=np.concatenate((CP_add,s_IFFT[j]))
        s_OFDM[j]=s_CP
    #print(np.shape(s_IFFT),np.shape(s_OFDM))
    #print(s_IFFT)
    #print(s_OFDM)
    s_Tx=np.matrix(np.reshape(s_OFDM,(1,-1)))
    sigma=np.sqrt((10**(-Eb_N0_dB[i]/10.0)))
    noise_real=np.random.normal(0,sigma,n+(n/N)*(L-1))
    noise_img=np.random.normal(0,sigma,n+(n/N)*(L-1))
    #noise=np.vectorize(complex)(noise_real,noise_img)/np.sqrt(2.0)
    noise=(noise_real+1j*noise_img)/np.sqrt(2.0)
    noise=np.array(noise)
    #noise_matrix=np.matrix(np.reshape(noise,(-1,N)))
    #print((np.shape(n)))
    y=s_Tx+noise
    #print(y)
    y_OFDM=np.matrix(np.reshape(y,(-1,N+L-1)))
    #print(np.shape(y_OFDM))
    for k in range(n/N):
        CP_rem=y_OFDM[k,L-1:]
        y_Rx[k]=CP_rem
        r_FFT[k:]=np.fft.fft(y_Rx[k])
        
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
plt.title('BER for OFDM_CP with BPSK Modulation  in AWGN')
#plt.xtricks([0,1,2,3,4,5,6,7,8,9,10,11])
plt.grid()
plt.show()              
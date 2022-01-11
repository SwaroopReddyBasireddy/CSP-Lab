# -*- coding: utf-8 -*-
"""
Created on Fri Sep 08 21:42:36 2017

@author: Swaroop
"""

import numpy as np
import scipy
import matplotlib.pyplot as plt

#Block Length

n=7

# Systematic bits
k=4

# Code Rate
R=float(k)/n

#number of bits
N=int(2e6)

#SNR length
snrlen=11

#SNR for uncoded system
Eb_No_dB=np.arange(0,snrlen)

#SNR for coded system
Ec_No_dB=Eb_No_dB-10*np.log10(1/R)


#parity Matrix
h=np.matrix([[1,0,1],[1,1,1],[1,1,0],[0,1,1]])

#Gene rator Mat r i x f o r Encoding
g=np.column_stack((np.eye(4),h))

# P a r i t y Check Mat r i x f o r Decoding
ht=(np.row_stack((h,np.eye(3)))).T

#Codebook
c_vec=np.zeros((2**k,n))

# nErr hard=np . z e r o s ( ( 1 , s n r l e n ) 
nErr_soft=np.zeros((1,snrlen))

for kk in range(2**k):
    m_vec=np.matrix(map(int,np.binary_repr(kk,width=k)))
    c_vec[kk,:]=(m_vec*g)%2
    
for yy in range(len(Eb_No_dB)):
    # t r a s m i t t e r
    ip = np.random.randint(2 ,size=N) 
    # g e n e r a t i n g 0 ,1 wi t h equap r o b a b i l i t y
    ip=np.array(ip)
    #hamming coding ( 7 , 4 )
    ipM=np.matrix(np.reshape(ip,(-1,4)))
    ipC=(ipM*g)%2
    cip=np.reshape(ipC,(1,(N/4)*7))
    # mo d u l a t i o n
    s=2*cip-1  #BPSK mo d u l a t i o n 0 −> −1; 1 −> 0
    # channel −AWGN
    sigma1=np.sqrt(0.5*(10**(-Ec_No_dB[yy]/10.0)))
    noise = np.random.normal(0,sigma1,np.shape(cip))
    y=s+noise
    y=np.array(y)
    # r e c e i v e r
    # S o f t d e c i s i o n Hamming de code r
    cipSoftM=np.reshape(np.real(y),(-1,n))
    c_vec=np.matrix(c_vec)
    corr=cipSoftM*(2*c_vec.T-1)
    idx=corr.argmax(axis=1)
    ipHatsoft=[]
    for i1 in range(np.shape(idx)[0]):
        aa=list(np.binary_repr(idx[i1,0],width=k))
        for j1 in range(k):
            ipHatsoft.append(int(aa[j1]))
    ipHatsoft=np.array(ipHatsoft)
    # c o u n t i n g t h e e r r o r s
    nErr_soft[0,yy]=np.count_nonzero(ip-ipHatsoft)
    
theoryBer=0.5*scipy.special.erfc(np.sqrt(10**(Eb_No_dB/10.0)))
#t h e o r e t i c a l be r uncoded AWGN
simBersoft=nErr_soft/float(N)
plt.plot(Eb_No_dB,theoryBer,'b',Eb_No_dB,simBersoft[0],'g')
plt.legend(['theory-Uncoded','coded-soft'],loc=1)
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for BPSK in AWGN with Hamming (7,4) code')
plt.xtricks([0,1,2,3,4,5,6,7,8,9,10,11])
plt.grid()
plt.show()


    
    






# -*- coding: utf-8 -*-
"""
Created on Mon Dec 04 19:10:58 2017

@author: Swaroop
"""

from numpy import sqrt
from numpy.random import rand, randn
import numpy as np
import matplotlib.pyplot as plt

n_cyc=10
n_fft=64
n_sym=10000
n_bits=n_fft*n_sym
snrlen=11
n_tap=5
Eb=1.0
delta=3

Eb_No_dB=np.arange(0,snrlen,1)
EbNo=np.zeros(snrlen,dtype=float)
ber=np.zeros(snrlen,dtype=float)
Y=np.zeros(n_fft+n_cyc,dtype=complex)
R=np.zeros(n_fft+n_cyc,dtype=float)
rec=np.zeros(n_fft+n_cyc,dtype=complex)
z=np.zeros(n_fft+n_cyc,dtype=complex)
h=np.zeros((n_sym,n_tap),dtype=complex)

s_in=2*(rand(n_bits)>=0.5)-1

r=np.zeros((n_sym,n_fft+n_cyc),dtype=complex)
for i in range(snrlen):
  for n in range (0, n_sym):
      s_ofdm=sqrt(n_fft)*np.fft.ifft(s_in[n*n_fft:(n+1)*n_fft])
      s_ofdm_cyc=np.concatenate((s_ofdm[n_fft-n_cyc:n_fft], s_ofdm), axis=0)
      ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(n_tap)+1j*randn(n_tap));
      h[n]=ht
      r_ofdm=np.convolve(s_ofdm_cyc,ht)
      r_ofdm_cyc=(r_ofdm[0:n_fft+n_cyc])
      wn=randn(n_fft+n_cyc)+1j*randn(n_fft+n_cyc)
      EbNo[i]=10**(Eb_No_dB[i]/10.0)
      r_ofdm_cyc=r_ofdm_cyc+sqrt(Eb/(2*EbNo[i]))*wn
      r[n]=r_ofdm_cyc
  y=np.reshape(r,(1,-1)) 
  y=np.array(y)[0]  
  x=np.reshape(h,(1,-1))
  #y=np.array(y)[0]
  x=np.array(x)[0]  
  #print(np.shape(y))
  #print(np.shape(r))
  
  pre=y[n_cyc:n_cyc+16]
  for k in range(len(pre)):
     Y[k]=np.conj(pre[k])
  Y1=np.matrix(np.reshape(Y,(-1,1)))
  #for k in range(n_fft):
   #   Y[k]=np.conj(pre[k])
  for ii in range(n_cyc):
      y1=y[ii+delta:delta+16+ii]
      X=np.matrix(np.reshape(y1,(1,-1)))
      for jj in range(16):
           z[ii] += Y1[jj,0]*X[0,jj]
      R[ii]=(np.abs(z[ii]))**2.0 
           
           
              
  delta1=delta+np.argmax(R)
  #print(delta1)
  n0=n_cyc-np.argmax(R)
  s_est=np.zeros(n_bits)
  
  #h=np.concatenate((np.zeros(delta1),ht))
  #H=np.fft.fft(h,n_fft)

  for m in range(0,n_sym):
      h1=x[m*n_tap:(m+1)*n_tap]
      h2=np.concatenate((np.zeros(n0),h1))
      #for k in range(n_tap):
       #   h1[k]=h1[n_tap-1-k]
      H1=np.fft.fft(h2,n_fft)
      #if m==0:
         
       #  p=delta1
        # s_est[m*n_fft:(m+1)*n_fft]=(1/sqrt(n_fft))*np.real(np.fft.fft(y[p:p+n_fft],n_fft)/H)
         
      #else:
      p=n_cyc
      s_est[m*n_fft:(m+1)*n_fft]=(1/sqrt(n_fft))*np.real(np.fft.fft(y[m*n_fft+p:(m+1)*n_fft+p],n_fft)/H1)         
  
  s_out=2*(s_est>=0)-1  
  errors=(s_in!=s_out).sum() 
  ber[i]=1.0*errors/n_bits

print('Bit Error Rate', ber)  

plt.plot(Eb_No_dB,ber,'b')
plt.legend(['theory','Practical'],loc=1)
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for OFDM with BPSK Modulation  in AWGN')
plt.grid()
plt.show()     
  
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 28 00:13:41 2017

@author: Swaroop
"""

from numpy import sqrt
from numpy.random import rand, randn
import numpy as np
import matplotlib.pyplot as plt


n_cyc=32
n_fft=128
n_sym=10000
n_bits=n_fft*n_sym
snrlen=11
n_tap=5
Eb=1.0
EbNodB=np.arange(0,snrlen,1)
EbNo=np.zeros(snrlen,dtype=float)
ber=np.zeros(snrlen,dtype=float)
s_in=2*(rand(n_bits)>=0.5)-1
s_est=np.zeros(n_bits)
for i in range(snrlen):
  for n in range (0, n_sym):
    s_ofdm=sqrt(n_fft)*np.fft.ifft(s_in[n*n_fft:(n+1)*n_fft])
    s_ofdm_cyc=np.concatenate((s_ofdm[n_fft-n_cyc:n_fft], s_ofdm), axis=0)
    ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(n_tap)+1j*randn(n_tap));
    Hf=np.fft.fft(ht,n_fft)
    r_ofdm=np.convolve(s_ofdm_cyc,ht)
    r_ofdm_cyc=(r_ofdm[0:n_fft+n_cyc])
    wn=randn(n_fft+n_cyc)+1j*randn(n_fft+n_cyc)
    EbNo[i]=10**(EbNodB[i]/10.0)
    r_ofdm_cyc=r_ofdm_cyc+sqrt(Eb/(2*EbNo[i]))*wn
    r_ofdm_cropped=r_ofdm_cyc[n_cyc:n_fft+n_cyc]
    s_est[n*n_fft:(n+1)*n_fft]=(1/sqrt(n_fft))*np.real((np.fft.fft(r_ofdm_cropped))/Hf)
    
  s_out=2*(s_est>=0)-1  
  errors=(s_in!=s_out).sum()
  ber[i]=1.0*errors/n_bits

#print('Total number of OFDM symbols', n_sym)
#print('Total number of bits', n_bits)
#print('Total number of bits in error', errors)
#print('Energy per bit to noise PSD(dB)', EbNodB)
print('Bit Error Rate', ber)
#theoryBer=0.5*scipy.special.erfc(np.sqrt((10**(EbNodB/10.0))))

plt.plot(EbNodB,ber,'b')
plt.legend(['theory','Practical'],loc=1)
plt.yscale('log')
plt.ylabel('BitErrorRate')
plt.xlabel('Eb/N0 in dB')
plt.title('BER for OFDM with BPSK Modulation  in AWGN')
#plt.xtricks([0,1,2,3,4,5,6,7,8,9,10,11])
plt.grid()
plt.show()
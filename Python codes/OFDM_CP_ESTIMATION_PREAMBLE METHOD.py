# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 22:14:28 2017

@author: Theresh
"""
from numpy import sqrt
import numpy as np
from numpy.random import rand,randn
import matplotlib.pyplot as plt


n_cyc=10
n_fft=64
n_sym=10000
n_bits=n_fft*n_sym
snrlen=11
n_tap=5
tau=3
Eb_No_dB=np.arange(0,snrlen,1)
EbNo=np.zeros(snrlen,dtype=float)
ber=np.zeros(snrlen,dtype=float)
z=np.zeros(16,dtype=complex)
c=np.zeros(16,dtype=complex)
C=np.zeros(16,dtype=float)
ber=np.zeros(snrlen,dtype=float)
for i in range(snrlen):
    errors=0
    for j in range(n_sym):
        errors_sym=0
        s_in=2*(rand(n_fft)>=0.5)-1
        s_ofdm=np.sqrt(n_fft)*np.fft.ifft(s_in)
        s_ofdm_cyc=np.concatenate((s_ofdm[n_fft-n_cyc:n_fft], s_ofdm), axis=0)
        ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(n_tap)+1j*randn(n_tap));
        r_ofdm=np.convolve(s_ofdm_cyc,ht)
        r_ofdm_cyc=(r_ofdm[0:n_fft+n_cyc])
        wn=randn(n_fft+n_cyc)+1j*randn(n_fft+n_cyc)
        EbNo[i]=10**(-Eb_No_dB[i]/10.0)
        r_ofdm_cyc=r_ofdm_cyc+sqrt(0.5*EbNo[i])*wn
        pre_conj=np.conj(r_ofdm_cyc[n_cyc:n_cyc+16])
        pre_conj=np.matrix(np.reshape(pre_conj,(-1,1)))
        for k in range(n_cyc):
            z=r_ofdm_cyc[tau+k:tau+k+16]
            X=np.matrix(np.reshape(z,(1,-1)))
            for m in range(16):
                c[k]+=pre_conj[m,0]*X[0,m]
            C[k]=(np.abs(c[k]))**2.0
        tau1=np.argmax(C)
        peak=tau+tau1-1
        n0=n_cyc-tau1
        if peak>n_cyc:
            tau1=n_cyc
            n0=n_cyc-tau1
        h=np.concatenate((np.zeros(n0),ht)) 
        H=np.fft.fft(h,n_fft)
        r=r_ofdm_cyc[tau1:64+tau1]
        R=np.fft.fft(r)
        s_est=(1/sqrt(n_fft))*np.real(R/H)
        s_out=2*(s_est>=0)-1 
        for m in range(len(s_in)):
            if (s_out[m]!=s_in[m]):
                errors_sym+=1
        errors +=errors_sym      
    ber[i]=1.0*errors/n_bits
plt.semilogy(Eb_No_dB,ber)
plt.title('Time synchronized based OFDM Detection')         
plt.xlabel("Eb_No_dB----->")
plt.ylabel("BER------->")
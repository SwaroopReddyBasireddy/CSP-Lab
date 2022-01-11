# -*- coding: utf-8 -*-
"""
Created on Wed Dec 06 09:25:47 2017

@author: 
    
"""

from numpy import sqrt
from numpy.random import rand, randn
import numpy as np
import matplotlib.pyplot as plt

n_cyc=10
n_fft=64
n_sym=10000
n_bits=n_fft*n_sym
snrlen=25
n_tap=5

Eb_No_dB=np.arange(0,snrlen,1)
EbNo=np.zeros(snrlen,dtype=float)
ber=np.zeros(snrlen,dtype=float)
Y=np.zeros(n_fft+n_cyc,dtype=complex)
R=np.zeros(n_fft+n_cyc,dtype=float)
rec=np.zeros(n_fft+n_cyc,dtype=complex)
z=np.zeros(n_fft+n_cyc,dtype=complex)
h=np.zeros((n_sym,n_tap),dtype=complex)



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

C=np.zeros(len(B2),dtype=float)
c=np.zeros(len(B2),dtype=complex)
tau=3
for i in range(snrlen):
    errors=0
    for j in range(n_sym):
        errors_sym=0
        ht=(1/sqrt(2))*(1/sqrt(n_tap))*(randn(n_tap)+1j*randn(n_tap));
        if j==0:
            s_in=2*(np.real(kk)>0)-1
            s_ofdm_cyc=Tx
        else:
            s_in=2*(rand(n_fft)>=0.5)-1
            s_ofdm=np.sqrt(n_fft)*np.fft.ifft(s_in)
            s_ofdm_cyc=np.concatenate((s_ofdm[n_fft-n_cyc:n_fft], s_ofdm), axis=0)
            #s_ofdm_cyc=np.concatenate((s_ofdm[n_fft-n_cyc:n_fft], s_ofdm), axis=0)
        r_ofdm=np.convolve(s_ofdm_cyc,ht)
        r_ofdm_cyc=(r_ofdm[0:n_fft+n_cyc])
        wn=randn(n_fft+n_cyc)+1j*randn(n_fft+n_cyc)
        EbNo[i]=10**(-Eb_No_dB[i]/10.0)
        r_ofdm_cyc=r_ofdm_cyc+sqrt(0.5*EbNo[i])*wn
        if j==0:
            x=B2
            for k in range(len(B2)):
               Y[k]=np.conj(x[k])
            Y1=np.matrix(np.reshape(Y,(-1,1)))
            
                       
            for ii in range(n_cyc):
                z=r_ofdm_cyc[tau+ii:tau+ii+len(x)]
                X=np.matrix(np.reshape(z,(1,-1)))
                for jj in range(16):
                  c[ii] += Y1[jj,0]*X[0,jj]
                C[ii]=(np.abs(c[ii]))**2.0
            tau1=np.argmax(C)    
            peak=tau+tau1-1
            n0=n_cyc-tau1
            if peak>10:
                tau1=n_cyc
                
                n0=n_cyc-tau1
            #print(peak)
            #print(n0)   
            
        h=np.concatenate((np.zeros(n0),ht)) 
        H=np.fft.fft(h,n_fft)
        r=r_ofdm_cyc[tau1:64+tau1]
        R=np.fft.fft(r)
        s_est=(1/sqrt(n_fft))*np.real(R/H)
        s_out=2*(s_est>=0)-1 
        for m in range(len(s_in)):
            if (s_out[m]!=s_in[m]):
               errors_sym+=1
        #print(errors_sym)       
        errors +=errors_sym
    ber[i]=1.0*errors/n_bits    
         
        #errors_sym=(s_in!=s_out).sum()
    
        
print(ber) 
plt.semilogy(Eb_No_dB,ber)       
plt.title('Time synchronized based OFDM Detection')         
plt.xlabel("Eb_No_dB")
plt.ylabel("BER")            
            
            
        


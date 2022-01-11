# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 15:26:30 2017

@author: Swaroop
"""

import numpy as np
from numpy import random as random
from matplotlib import pyplot as plt
from scipy.special import erfc
N=1000000
x=random.randint(0,2,N)
mx_fp=np.int16((2*x-1)*(2**13))
mx=2*x-1
ebnodb=np.arange(0,11,1)
e=np.size(ebnodb)
errors=0
ber=np.zeros(e)
#pe=np.zeros(e)
for t in range(e):
    sigma=(10.0**(-ebnodb[t]/20.0))
    noise1=np.random.normal(0,sigma,N)
    noise=np.int16((np.divide(noise1,np.sqrt(2)))*(2**13))
    tx=mx_fp+noise
    rx1=(tx>=0)
    rx=2*rx1-1
    errors = (mx!=rx).sum()
    ber[t]=1.0*errors/N  
pe=0.5*erfc(np.sqrt(10.0**(ebnodb/10.0)))
print (ber,pe)
plt.semilogy(ebnodb,pe,'r',linewidth=2)
plt.semilogy(ebnodb,ber,'-s')
plt.grid(True)
plt.legend(('analytical','simulation'))
plt.xlabel('---Eb/No (dB)--->')
plt.ylabel('------BER---->')
plt.title('BPSK BER')
plt.show()
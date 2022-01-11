# -*- coding: utf-8 -*-
"""
Created on Wed Nov 01 16:07:07 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt

N=100000

Y=np.random.randn(N)
Y=np.array(Y)
print(Y,np.shape(Y))
#X=np.hist(Y)

sum=0.0
for i in range(N):
  sum = sum + Y[i]
Mean=sum/float(N)

print(Mean)
Var_sum=0.0
for j in range(N):
    Var_sum += (Y[j]-Mean)**2.0
    
Variance=Var_sum/float(N)

plt.hist(Y)   

print('Mean=', Mean,'Vaiance=',Variance) 

Mean_error=0.0-Mean
Var_error=1.0-Variance

print('Mean_error=',Mean_error,'Var_error=',Var_error)
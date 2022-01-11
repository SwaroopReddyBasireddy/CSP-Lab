
import numpy as np
from matplotlib import pyplot as plt
t=np.linspace(-1,1,201)
a=np.int16(np.linspace(-1,1,201)*(2**13))   ## a=([-1:0.01:1])
b=([1,0.8,0.3])
m=len(a)
n=len(b)
c=np.zeros((m,n))
for i in range(m):
    for j in range(n):
        c[i][j]=a[i]*b[j]
print c
y=np.zeros(n+m-1)
for k in range(n+m-1):
    for i in range(m):
        for j in range(n):
            if i+j==k:
                y[k]+=c[i][j]
print"The linear convolution is", np.round((y/(2**13)),2)
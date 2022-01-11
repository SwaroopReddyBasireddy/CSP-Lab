import numpy as np
import matplotlib.pyplot as plt
from scipy import special


N=int(1e6)
snrlen= 15
Eb_N0_dB=np.arange(0,snrlen)
ber1 = list()
ber2 = list()
for i in range(snrlen):
	s_in = 2*np.random.randint(2,size=N) - 1
	x= np.random.normal(0,1,N)
	y=np.random.normal(0,1,N)
	h = np.sqrt((x**2 + y**2)/2.0)
	sigma = np.sqrt(10**(-0.1*Eb_N0_dB[i]))
	noise = np.random.normal(0,sigma,N)
	y1 = h * s_in + noise
	y2 = s_in + noise
				
	x_hat1=2*(y1>=0)-1
	x_hat2=2*(y2>=0)-1
	errors1=(x_hat1!=s_in).sum()
	errors2=(x_hat2!=s_in).sum()
	ber1.append(errors1/float(N))
	ber2.append(errors2/float(N))
SNR = 10**(0.1*Eb_N0_dB)	
ber1_theory = 0.5 *(1 - np.sqrt(SNR/(2+SNR)))
ber2_theory = 0.5*special.erfc(np.sqrt(SNR/2.0))
plt.semilogy(Eb_N0_dB,ber1,color='g')
plt.semilogy(Eb_N0_dB,ber2,color='b')
plt.semilogy(Eb_N0_dB,ber1_theory,'*',color='r')
plt.semilogy(Eb_N0_dB,ber2_theory,'+',color='r')

plt.grid()
plt.xlabel('$ SNR $')
plt.ylabel('$ BER $')
plt.legend(['Rayleigh fading channel','AWGN channel','Rayleigh fading_theory','AWGN_theory'],loc='best')
plt.savefig('AWGN Vs Rayleigh_fading.eps')
plt.show()

	

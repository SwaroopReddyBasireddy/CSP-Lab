# -*- coding: utf-8 -*-
"""
Created on Sat Nov 04 17:40:10 2017

@author: Swaroop
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy

beta = 0.20 # roll off factor

Tsample = 1.0 # sampling period, should at least twice the rate of the symbol

oversampling_rate = 8 # oversampling of the bit stream, this gives samples per symbol
# must be at least 2X the bit rate

Tsymbol = oversampling_rate * Tsample # pulse duration should be at least 2 * Ts
span = 50 # number of symbols to span, must be even
n = span*oversampling_rate # length of the filter = samples per symbol * symbol span

# t_step must be from -span/2 to +span/2 symbols.
# each symbol has 'sps' number of samples per second.
t_step = Tsample * np.linspace(-n/2,n/2,n+1) # n+1 to include 0 time

BW = (1 + beta) / Tsymbol
a = np.zeros_like(t_step)

for item in list(enumerate(t_step)):
    i,t = item 
    # t is n*Ts
    if (1-(2.0*beta*t/Tsymbol)**2) == 0:
        a[i] = np.pi/4 * np.sinc(t/Tsymbol)
        print 'i = %d' % i
    elif t == 0:
        a[i] = np.cos(beta * np.pi * t / Tsymbol)/ (1-(2.0*beta*t/Tsymbol)**2)
        print 't = 0 captured'
        print 'i = %d' % i 

    else:
        numerator = np.sinc( np.pi * t/Tsymbol )*np.cos( np.pi*beta*t/Tsymbol )
        denominator = (1.0 - (2.0*beta*t/Tsymbol)**2)
        a[i] =  numerator / denominator

#a = a/sum(a) # normalize total power

plot_filter = 0
if plot_filter == 1:

    w,h = scipy.signal.freqz(a)
    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.title('Digital filter (raised cosine) frequency response')
    ax1 = fig.add_subplot(211)
    plt.plot(w/np.pi, 20*np.log10(abs(h)),'b')
    #plt.plot(w/np.pi, abs(h),'b')
    plt.ylabel('Amplitude (dB)', color = 'b')
    plt.xlabel(r'Normalized Frequency ($\pi$ rad/sample)')

    ax2 = ax1.twinx()
    angles = np.unwrap(np.angle(h))
    plt.plot(w/np.pi, angles, 'g')
    plt.ylabel('Angle (radians)', color = 'g')
    plt.grid()
    plt.axis('tight')
    plt.show()


    plt.subplot(2,1,2)
    plt.stem(a)
    plt.show()
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 09:55:32 2018

@author: dori
"""

import refractive
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

c = 299792458.0

temp = 273.15
frequencies = 1.0e9*np.linspace(1,250,1000)
wavelenghts = c/frequencies
wavenumbers = 2.0*np.pi/wavelenghts

ni = refractive.ice.n(temp,frequencies)
nw = refractive.water.n(temp,frequencies)

ds = 1.0e-6*np.array([10,20,30,40])
f, (ax1) = plt.subplots(1,1)
ax1t = ax1.twinx()
for d in ds:
    print(d)
    mkdi = wavenumbers*d*np.abs(ni)
    mkdw = wavenumbers*d*np.abs(nw)
    ax1.plot(frequencies*1.0e-9,mkdi,label=str(int(round(1.0e6*d))))
    ax1t.plot(frequencies*1.0e-9,mkdw,'--',label=str(int(round(1.0e6*d))))
ax1.set_ylim(0,0.6)
ax1t.set_ylim(0,0.6)
ax1.grid()
ax1.legend(loc=2,title='d  [$\mu$m]')
#ax1t.legend(loc=6,title='water dipole')
ax1t.tick_params( axis='y',
                 which='both',
                 left='off',
                 right='off',
                 labelright='off',
                 labelleft='off')
ax1.set_xlabel('Frequency    [GHz]')
ax1.set_ylabel('|m|kd')
f.savefig('mkd.pdf')
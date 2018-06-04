# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 09:55:32 2018

@author: dori
"""

import refractiveIndex
import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

from matplotlib.lines import Line2D

c = 299792458.0

temp = 273.15
frequencies = 1.0e9*np.linspace(1,250,1000)
wavelenghts = c/frequencies
wavenumbers = 2.0*np.pi/wavelenghts

ni = refractiveIndex.ice.n(temp,frequencies)
nw = refractiveIndex.water.n(temp,frequencies)
delei = refractiveIndex.skin_depth(frequency=frequencies,refractive_index=ni)
delew = refractiveIndex.skin_depth(frequency=frequencies,refractive_index=nw)

fig, ((ax00,ax01),(ax10,ax11)) = plt.subplots(2,2,figsize=(8,6))

ds = 1.0e-6*np.array([10,20,30,40,50])
ax00t = ax00.twinx()
#c='#1f77b4''#ff7f0e'
ax00.vlines([2.8,5.6,9.6,13.6,35.6,94.0,220.0],ymin=1,ymax=10,colors='g')
ax00.text(2.2,8,'S',color='g')
ax00.text(4.6,8,'C',color='g')
ax00.text(7.8,8,'X',color='g')
ax00.text(14,8,'Ku',color='g')
ax00.text(26,8,'Ka',color='g')
ax00.text(70,8,'W',color='g')
ax00.text(170,8,'G',color='g')
ax00.plot(frequencies*1e-9,ni.real,label='ice')
ax00.plot(frequencies*1e-9,nw.real,label='water')
ax00t.plot(frequencies*1e-9,ni.imag,'--')
ax00t.plot(frequencies*1e-9,nw.imag,'--')
ax00.set_xscale('log')
ax00t.set_xscale('log')
ax00t.set_yscale('log')

ax01.plot(frequencies*1e-9,delei,label='ice')
ax01.plot(frequencies*1e-9,delew,label='water')
ax01.set_xscale('log')
ax01.set_yscale('log')
ax00.set_ylabel('Real(m)')
ax00t.set_ylabel('Imag(m)')
ax01.set_ylabel('skin depth   $\delta_e$   [m]')
ax00.grid()
ax01.legend()
custom_lines = [Line2D([0], [0], color='k', lw=3),
                Line2D([0], [0], linestyle=':', color='k', lw=3)]
ax00.legend(loc=6)
ax00t.legend(custom_lines,['Real(m)','Imag(m)'],loc=7)

ax01.grid()
ax00.set_xlabel('Frequency    [GHz]')
ax01.set_xlabel('Frequency    [GHz]')

ax10t = ax10.twinx()
ax11t = ax11.twinx()
for d in ds:
    print(d)
    mkdi = wavenumbers*d*np.abs(ni)
    mkdw = wavenumbers*d*np.abs(nw)
    magw = refractiveIndex.magnetic2electric_ratio(size=d,frequency=frequencies,refractive_index=nw)
    magi = refractiveIndex.magnetic2electric_ratio(size=d,frequency=frequencies,refractive_index=ni)
    ax10.plot(frequencies*1.0e-9,mkdi,label=str(int(round(1.0e6*d))))
    ax10t.plot(frequencies*1.0e-9,mkdw,'--',label=str(int(round(1.0e6*d))))
    ax11.plot(frequencies*1.0e-9,magi,label=str(int(round(1.0e6*d))))
    ax11t.plot(frequencies*1.0e-9,magw,'--',label=str(int(round(1.0e6*d))))
ax10.set_ylim(0,0.6)
ax10t.set_ylim(0,0.6)
ax11.set_ylim(0,0.03)
ax11t.set_ylim(0,0.03)

ax10.grid()
ax11.grid()
ax10.legend(loc=2,title='d  [$\mu$m]')
ax11.legend(loc=2,title='d  [$\mu$m]')
#ax1t.legend(loc=6,title='water dipole')
ax10t.tick_params( axis='y',
                 which='both',
                 left='off',
                 right='off',
                 labelright='off',
                 labelleft='off')
ax11t.tick_params( axis='y',
                 which='both',
                 left='off',
                 right='off',
                 labelright='off',
                 labelleft='off')
ax10.set_xlabel('Frequency    [GHz]')
ax11.set_xlabel('Frequency    [GHz]')
ax10.set_ylabel('|m|kd')
ax11.set_ylabel('$C^m_{abs}/C^e_{abs}$')
ax10t.legend(custom_lines,['ice','water'],loc=9)
ax11t.legend(custom_lines,['ice','water'],loc=9)
fig.tight_layout()
fig.savefig('mkd.pdf')
fig.savefig('mkd.png',dpi=300)
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  7 10:46:58 2018

@author: dori
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from scipy import integrate
import gc
import resource

sys.path.append('/work/DBs/scattnlay')
from scattnlay import scattnlay
from scattnlay import fieldnlay

import refractiveIndex as ref

import matplotlib as mpl
from cycler import cycler
mpl.rcParams['axes.prop_cycle'] = cycler(u'color', [u'#d62728', u'#ff7f0e', u'#2ca02c', u'#1f77b4',  u'#9467bd'])
from collections import OrderedDict

c = 299792458. # m/s
size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c*1.e-6/f # expected GHz returns mm

freqs=OrderedDict([('X',9.6),('Ku',13.6),('Ka',35.6),('W',94),('G',220)])
thickness=np.array([0.0,0.001,0.01,0.05,0.1])
sizes = np.linspace(2,20,1000)

temp = 273.15
Ncomput = len(sizes)#1
Nlayers = 2

fig,axs = plt.subplots(3,5,figsize=(29,21))

i = 0
for lay in thickness:
    #f,ax = plt.subplots(1,1)
    ax = axs[0][i]
    for fk in freqs.keys():
        f = freqs[fk]
        wl = fr2lam(f)
        mi = ref.ice.n(temp,f*1e9)
        mw = ref.water.n(temp,f*1e9)
        x = np.ndarray((Ncomput,Nlayers),dtype=np.float64)
        m = np.ndarray((Ncomput,Nlayers),dtype=np.complex128)
        outer_size = sizes
        outer_x = size2x(outer_size,wl)
        inner_size = outer_size - lay
        inner_x = size2x(inner_size,wl)
        x[:,0] = inner_x
        x[:,1] = outer_x
        m[:,0] = mi
        m[:,1] = mw
        #thetas = deg2rad(mueller.index.values)
        thetas = None
        terms, MQe, MQs, MQa, MQb, MQp, Mg, Mssa, S1, S2 = scattnlay(x,m)
        ax.plot(sizes,MQb,label=fk)
        ax.set_title('water layer thickness '+str(lay))
    ax.legend(title='frequency')
    ax.set_xlabel('sphere outer diameter')
    ax.set_ylabel('backscattering efficiency')
    ax.grid()
    i = i+1

i = 0
for fk in freqs.keys():
    #f,ax = plt.subplots(1,1)
    ax = axs[1][i]
    for lay in thickness:
        f = freqs[fk]
        wl = fr2lam(f)
        mi = ref.ice.n(temp,f*1e9)
        mw = ref.water.n(temp,f*1e9)
        x = np.ndarray((Ncomput,Nlayers),dtype=np.float64)
        m = np.ndarray((Ncomput,Nlayers),dtype=np.complex128)
        outer_size = sizes
        outer_x = size2x(outer_size,wl)
        inner_size = outer_size - lay
        inner_x = size2x(inner_size,wl)
        x[:,0] = inner_x
        x[:,1] = outer_x
        m[:,0] = mi
        m[:,1] = mw
        #thetas = deg2rad(mueller.index.values)
        thetas = None
        terms, MQe, MQs, MQa, MQb, MQp, Mg, Mssa, S1, S2 = scattnlay(x,m)
        ax.plot(sizes,MQb,label=str(lay))
        ax.set_title(fk)
    ax.legend(title='water thickness')
    ax.set_xlabel('sphere outer diameter')
    ax.set_ylabel('backscattering efficiency')
    ax.grid()
    i = i+1

i = 0
ratios = np.array([0.99,0.95,0.90,0.85,0.80])
for fk in freqs.keys():
    #f,ax = plt.subplots(1,1)
    ax = axs[2][i]
    for r in ratios:
        f = freqs[fk]
        wl = fr2lam(f)
        mi = ref.ice.n(temp,f*1e9)
        mw = ref.water.n(temp,f*1e9)
        x = np.ndarray((Ncomput,Nlayers),dtype=np.float64)
        m = np.ndarray((Ncomput,Nlayers),dtype=np.complex128)
        outer_size = sizes
        outer_x = size2x(outer_size,wl)
        inner_size = r*outer_size
        inner_x = size2x(inner_size,wl)
        x[:,0] = inner_x
        x[:,1] = outer_x
        m[:,0] = mi
        m[:,1] = mw
        #thetas = deg2rad(mueller.index.values)
        thetas = None
        terms, MQe, MQs, MQa, MQb, MQp, Mg, Mssa, S1, S2 = scattnlay(x,m)
        ax.plot(sizes,MQb,label=str(r))
        ax.set_title(fk)
    ax.legend(title='diameters ratio')
    ax.set_xlabel('sphere outer diameter')
    ax.set_ylabel('backscattering efficiency')
    ax.grid()
    i = i+1
fig.savefig('multiplot_sphere.pdf')
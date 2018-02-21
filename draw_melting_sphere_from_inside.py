# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 09:47:56 2018

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

size = 1.0
melt_fracs = np.arange(0.2,1.0,0.2)
cx = cy = size*0.5

f,axes = plt.subplots(1,len(melt_fracs))

for mf,ax in zip(melt_fracs,axes):
    fr = 1-mf
    insize = fr*size
    outcircle = plt.Circle((cx,cy),size*0.5,color='xkcd:azure', fill=True,alpha=0.7)
    incircle = plt.Circle((cx,cy),insize*0.5,color='xkcd:blue grey', fill=True)
    ax.add_artist(outcircle)
    ax.add_artist(incircle)
    outring = plt.Circle((cx,cy),size*0.5,color='k', fill=False)
    inring = plt.Circle((cx,cy),insize*0.5,color='k', fill=False)
    ax.add_artist(outring)
    ax.add_artist(inring)
    ax.set_xlabel('f='+str(mf))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect(1.0)
f.savefig('melting_sphere.png',dpi=600,bbox_inches='tight')
f.savefig('melting_sphere.pdf',bbox_inches='tight')


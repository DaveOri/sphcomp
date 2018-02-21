# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:46:45 2018

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

size = 1.0 # size of the domain
cx = cy = size*0.5
Dratio = 0.6
insize = Dratio*size

fig, ax = plt.subplots(figsize=(6,6))
outcircle = plt.Circle((cx,cy),size*0.5,color='k', fill=False)
incircle = plt.Circle((cx,cy),insize*0.5,color='k', fill=False)

ax.add_artist(outcircle)
ax.add_artist(incircle)

jagged=1
resolved=2

resolution = 8 # how many dipoles per side
d = size/resolution
halfd = 0.5*d 
ccoordinates = np.arange(halfd,size,d)
ptsx, ptsy = np.meshgrid(ccoordinates,ccoordinates)
ptsx = ptsx.flatten()
ptsy = ptsy.flatten()
rad2 = ((ptsx-cx)**2+(ptsy-cy)**2)

if not (resolved-1 or jagged-1):
    rboolin = (rad2 < (0.5*insize)**2)
    ptinx, ptiny = ptsx[rboolin],ptsy[rboolin]
    ax.scatter(ptinx,ptiny,c='xkcd:blue')
    
    rboolout = (rad2 < (0.5*size)**2) * (rad2 > (0.5*insize)**2)
    ptoutx, ptouty = ptsx[rboolout],ptsy[rboolout]
    ax.scatter(ptoutx,ptouty,c='xkcd:bright orange')
    for x,y in zip(ptinx,ptiny):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:turquoise')
        ax.add_patch(p)
    
    for x,y in zip(ptoutx,ptouty):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:orange')
        ax.add_patch(p)
    
    ax.set_xlim([0,size])
    ax.set_ylim([0,size])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect(1.0)
    fig.savefig('grid_layered.pdf',bbox_inches='tight')
    fig.savefig('grid_layered.png',dpi=300,bbox_inches='tight')

if jagged-1:
    dj = d/jagged
    halfdj = 0.5*dj
    for cxj,cyj in zip(ptoutx,ptouty):
        ccoojx = np.arange(cxj-halfd+halfdj,cxj+halfd,dj)
        ccoojy = np.arange(cyj-halfd+halfdj,cyj+halfd,dj)
        ptsxj, ptsyj = np.meshgrid(ccoojx,ccoojy)
        ptsxj = ptsxj.flatten()
        ptsyj = ptsyj.flatten()
        ax.scatter(ptsxj,ptsyj,c='xkcd:bright orange',s=10)

        for x,y in zip(ptsxj,ptsyj):
            p = patches.Rectangle((x-halfdj,y-halfdj),dj,dj,fill=True,alpha=0.2,color='xkcd:orange')
            ax.add_patch(p)
        
    for cxj,cyj in zip(ptinx,ptiny):
        ccoojx = np.arange(cxj-halfd+halfdj,cxj+halfd,dj)
        ccoojy = np.arange(cyj-halfd+halfdj,cyj+halfd,dj)
        ptsxj, ptsyj = np.meshgrid(ccoojx,ccoojy)
        ptsxj = ptsxj.flatten()
        ptsyj = ptsyj.flatten()
        ax.scatter(ptsxj,ptsyj,c='xkcd:blue',s=10)
        
        for x,y in zip(ptsxj,ptsyj):
            p = patches.Rectangle((x-halfdj,y-halfdj),dj,dj,fill=True,alpha=0.2,color='xkcd:turquoise')
            ax.add_patch(p)

    ax.set_xlim([0,size])
    ax.set_ylim([0,size])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect(1.0)
    fig.savefig('grid_layered_jagged.pdf',bbox_inches='tight')
    fig.savefig('grid_layered_jagged.png',dpi=300,bbox_inches='tight')

if resolved-1:
    resolution = resolution*resolved # how many dipoles per side
    d = size/resolution
    halfd = 0.5*d 
    ccoordinates = np.arange(halfd,size,d)
    ptsx, ptsy = np.meshgrid(ccoordinates,ccoordinates)
    ptsx = ptsx.flatten()
    ptsy = ptsy.flatten()
    rad2 = ((ptsx-cx)**2+(ptsy-cy)**2)
    
    rboolin = (rad2 < (0.5*insize)**2)
    ptinx, ptiny = ptsx[rboolin],ptsy[rboolin]
    ax.scatter(ptinx,ptiny,c='xkcd:blue',s=10)

    rboolout = (rad2 < (0.5*size)**2) * (rad2 > (0.5*insize)**2)
    ptoutx, ptouty = ptsx[rboolout],ptsy[rboolout]
    ax.scatter(ptoutx,ptouty,c='xkcd:bright orange',s=10)
    
    for x,y in zip(ptinx,ptiny):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:turquoise')
        ax.add_patch(p)

    for x,y in zip(ptoutx,ptouty):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:orange')
        ax.add_patch(p)

    ax.set_xlim([0,size])
    ax.set_ylim([0,size])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_aspect(1.0)
    fig.savefig('grid_layered_hires.pdf',bbox_inches='tight')
    fig.savefig('grid_layered_hires.png',dpi=300,bbox_inches='tight')

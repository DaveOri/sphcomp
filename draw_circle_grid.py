# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:46:45 2018

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec

plt.close('all')

size = 1.0 # size of the domain
cx = cy = size*0.5
Dratio = 0.6
insize = Dratio*size

#fig, ax = plt.subplots(figsize=(6,6))
#fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,6))
fig = plt.figure(figsize=(20,6))
gs = gridspec.GridSpec(2,12)
axa = plt.subplot(gs[0,0:3]) #plt.subplot(241)
axb = plt.subplot(gs[0,3:6])#plt.subplot(242)
axc = plt.subplot(gs[0,6:9])#plt.subplot(243)
axd = plt.subplot(gs[0,9:12])#plt.subplot(244)
ax1 = plt.subplot(gs[1,0:4])#plt.subplot(234)
ax2 = plt.subplot(gs[1,4:8])#plt.subplot(235)
ax3 = plt.subplot(gs[1,8:12])#plt.subplot(236)
outcircle = plt.Circle((cx,cy),size*0.5,color='k', fill=False)
incircle = plt.Circle((cx,cy),insize*0.5,color='k', fill=False)

ax1.add_artist(outcircle)
ax1.add_artist(incircle)

jagged=2
resolved=2

resolution = 8 # how many dipoles per side
d = size/resolution
halfd = 0.5*d 
ccoordinates = np.arange(halfd,size,d)
ptsx, ptsy = np.meshgrid(ccoordinates,ccoordinates)
ptsx = ptsx.flatten()
ptsy = ptsy.flatten()
rad2 = ((ptsx-cx)**2+(ptsy-cy)**2)

if 1:#not (resolved-1 or jagged-1):
    rboolin = (rad2 < (0.5*insize)**2)
    ptinx, ptiny = ptsx[rboolin],ptsy[rboolin]
    ax1.scatter(ptinx,ptiny,c='xkcd:blue grey')
    
    rboolout = (rad2 < (0.5*size)**2) * (rad2 > (0.5*insize)**2)
    ptoutx, ptouty = ptsx[rboolout],ptsy[rboolout]
    ax1.scatter(ptoutx,ptouty,c='xkcd:electric blue')
    for x,y in zip(ptinx,ptiny):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:blue grey')
        ax1.add_patch(p)
    
    for x,y in zip(ptoutx,ptouty):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:azure')
        ax1.add_patch(p)
    
    ax1.set_xlim([0,size])
    ax1.set_ylim([0,size])
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_aspect(1.0)
    #fig.savefig('grid_layered.pdf',bbox_inches='tight')
    #fig.savefig('grid_layered.png',dpi=300,bbox_inches='tight')

if jagged-1:
    outcircle = plt.Circle((cx,cy),size*0.5,color='k', fill=False)
    incircle = plt.Circle((cx,cy),insize*0.5,color='k', fill=False)
    ax2.add_artist(outcircle)
    ax2.add_artist(incircle)
    dj = d/jagged
    halfdj = 0.5*dj
    for cxj,cyj in zip(ptoutx,ptouty):
        ccoojx = np.arange(cxj-halfd+halfdj,cxj+halfd,dj)
        ccoojy = np.arange(cyj-halfd+halfdj,cyj+halfd,dj)
        ptsxj, ptsyj = np.meshgrid(ccoojx,ccoojy)
        ptsxj = ptsxj.flatten()
        ptsyj = ptsyj.flatten()
        ax2.scatter(ptsxj,ptsyj,c='xkcd:electric blue',s=10)

        for x,y in zip(ptsxj,ptsyj):
            p = patches.Rectangle((x-halfdj,y-halfdj),dj,dj,fill=True,alpha=0.2,color='xkcd:azure')
            ax2.add_patch(p)
        
    for cxj,cyj in zip(ptinx,ptiny):
        ccoojx = np.arange(cxj-halfd+halfdj,cxj+halfd,dj)
        ccoojy = np.arange(cyj-halfd+halfdj,cyj+halfd,dj)
        ptsxj, ptsyj = np.meshgrid(ccoojx,ccoojy)
        ptsxj = ptsxj.flatten()
        ptsyj = ptsyj.flatten()
        ax2.scatter(ptsxj,ptsyj,c='xkcd:blue grey',s=10)
        
        for x,y in zip(ptsxj,ptsyj):
            p = patches.Rectangle((x-halfdj,y-halfdj),dj,dj,fill=True,alpha=0.2,color='xkcd:blue grey')
            ax2.add_patch(p)

    ax2.set_xlim([0,size])
    ax2.set_ylim([0,size])
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_aspect(1.0)
    #fig.savefig('grid_layered_jagged.pdf',bbox_inches='tight')
    #fig.savefig('grid_layered_jagged.png',dpi=300,bbox_inches='tight')

if resolved-1:
    outcircle = plt.Circle((cx,cy),size*0.5,color='k', fill=False)
    incircle = plt.Circle((cx,cy),insize*0.5,color='k', fill=False)
    ax3.add_artist(outcircle)
    ax3.add_artist(incircle)
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
    ax3.scatter(ptinx,ptiny,c='xkcd:blue grey',s=10)

    rboolout = (rad2 < (0.5*size)**2) * (rad2 > (0.5*insize)**2)
    ptoutx, ptouty = ptsx[rboolout],ptsy[rboolout]
    ax3.scatter(ptoutx,ptouty,c='xkcd:electric blue',s=10)
    
    for x,y in zip(ptinx,ptiny):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:blue grey')
        ax3.add_patch(p)

    for x,y in zip(ptoutx,ptouty):
        p = patches.Rectangle((x-halfd,y-halfd),d,d,fill=True,alpha=0.2,color='xkcd:azure')
        ax3.add_patch(p)

    ax3.set_xlim([0,size])
    ax3.set_ylim([0,size])
    ax3.set_xticks([])
    ax3.set_yticks([])
    ax3.set_aspect(1.0)
    #fig.savefig('grid_layered_hires.pdf',bbox_inches='tight')
    #fig.savefig('grid_layered_hires.png',dpi=300,bbox_inches='tight')
ax1.set_xlabel('original')
ax2.set_xlabel('half dipole length')
ax3.set_xlabel('double resolution')

size = 1.0
melt_fracs = np.arange(0.2,1.0,0.2)
cx = cy = size*0.5

axes = [axa,axb,axc,axd]

for mf,ax in zip(melt_fracs,axes):
    fr = 1-mf
    insize = size*fr**(1.0/3.0)
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
ax1.text(0.05,0.9,'a)',fontweight='bold')
ax2.text(0.05,0.9,'b)',fontweight='bold')
ax3.text(0.05,0.9,'c)',fontweight='bold')
gs.tight_layout(fig,rect=[0.5,0.0,1.0,1.0], h_pad=0.1)
fig.savefig('combo_grid_refinement.pdf',bbox_inches='tight')
fig.savefig('combo_grid_refinement.png',dpi=600,bbox_inches='tight')
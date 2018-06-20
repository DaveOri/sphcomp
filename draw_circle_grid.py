# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 11:46:45 2018

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.colors as clrs

from mpl_toolkits.mplot3d import Axes3D

plt.close('all')

size = 1.0 # size of the domain
cx = cy = size*0.5
Dratio = 0.6
insize = Dratio*size

#fig, ax = plt.subplots(figsize=(6,6))
#fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(10,6))

#fig = plt.figure(figsize=(20,6))
#gs = gridspec.GridSpec(2,12)
fig = plt.figure(figsize=(15,6.75)) #20,9
gs = gridspec.GridSpec(3,12)

axa = plt.subplot(gs[0,0:3]) #plt.subplot(241)
axb = plt.subplot(gs[0,3:6])#plt.subplot(242)
axc = plt.subplot(gs[0,6:9])#plt.subplot(243)
axd = plt.subplot(gs[0,9:12])#plt.subplot(244)
ax1 = plt.subplot(gs[1,0:4])#plt.subplot(234)
ax2 = plt.subplot(gs[1,4:8])#plt.subplot(235)
ax3 = plt.subplot(gs[1,8:12])#plt.subplot(236)
axm1= plt.subplot(gs[2,1:5],projection='3d')
#axm1.set_aspect('equal')
axm2= plt.subplot(gs[2,7:11],projection='3d')
#axm2.set_aspect('equal')

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

s=8
if 1:#not (resolved-1 or jagged-1):
    rboolin = (rad2 < (0.5*insize)**2)
    ptinx, ptiny = ptsx[rboolin],ptsy[rboolin]
    ax1.scatter(ptinx,ptiny,c='xkcd:blue grey',s=s)
    
    rboolout = (rad2 < (0.5*size)**2) * (rad2 > (0.5*insize)**2)
    ptoutx, ptouty = ptsx[rboolout],ptsy[rboolout]
    ax1.scatter(ptoutx,ptouty,c='xkcd:electric blue',s=s)
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
s=4
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
        ax2.scatter(ptsxj,ptsyj,c='xkcd:electric blue',s=s)

        for x,y in zip(ptsxj,ptsyj):
            p = patches.Rectangle((x-halfdj,y-halfdj),dj,dj,fill=True,alpha=0.2,color='xkcd:azure')
            ax2.add_patch(p)
        
    for cxj,cyj in zip(ptinx,ptiny):
        ccoojx = np.arange(cxj-halfd+halfdj,cxj+halfd,dj)
        ccoojy = np.arange(cyj-halfd+halfdj,cyj+halfd,dj)
        ptsxj, ptsyj = np.meshgrid(ccoojx,ccoojy)
        ptsxj = ptsxj.flatten()
        ptsyj = ptsyj.flatten()
        ax2.scatter(ptsxj,ptsyj,c='xkcd:blue grey',s=s)
        
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
    ax3.scatter(ptinx,ptiny,c='xkcd:blue grey',s=s)

    rboolout = (rad2 < (0.5*size)**2) * (rad2 > (0.5*insize)**2)
    ptoutx, ptouty = ptsx[rboolout],ptsy[rboolout]
    ax3.scatter(ptoutx,ptouty,c='xkcd:electric blue',s=s)
    
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

#axm1.set_aspect(1.0)
#axm2.set_aspect(1.0)

import sys
sys.path.append('/work/DBs/scattDB/scattDB/')
import shape
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

def cuboid_data(o, size=(1,1,1)):
    X = [[[0, 1, 0], [0, 0, 0], [1, 0, 0], [1, 1, 0]],
         [[0, 0, 0], [0, 0, 1], [1, 0, 1], [1, 0, 0]],
         [[1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1]],
         [[0, 0, 1], [0, 0, 0], [0, 1, 0], [0, 1, 1]],
         [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
         [[0, 1, 1], [0, 0, 1], [1, 0, 1], [1, 1, 1]]]
    X = np.array(X).astype(float)
    for i in range(3):
        X[:,:,i] *= size[i]
    X += np.array(o)
    return X

def plotCubeAt(positions,sizes=None,colors=None, **kwargs):
    if not isinstance(colors,(list,np.ndarray)): colors=["C0"]*len(positions)
    if not isinstance(sizes,(list,np.ndarray)): sizes=[(1,1,1)]*len(positions)
    g = []
    for p,s,c in zip(positions,sizes,colors):
        g.append( cuboid_data(p, size=s) )
    return Poly3DCollection(np.concatenate(g),  
                            facecolors=np.repeat(colors,6, axis=0), **kwargs)

iceshape = shape.shp1200.shape[shape.shp1200.shape['CX']==1]
meltshape = shape.shp1200.shape[shape.shp1200.shape['CX']==2]

xs = iceshape.X.values
ys = iceshape.Y.values
zs = iceshape.Z.values

grey = clrs.to_rgb('xkcd:bluegrey')
azure = clrs.to_rgb('xkcd:electric blue')

positions = np.c_[xs,ys,zs]
pc11 = plotCubeAt(positions, colors=np.tile(grey,(len(xs),1)))
axm1.add_collection3d(pc11)
pc12 = plotCubeAt(positions, colors=np.tile(grey,(len(xs),1)),edgecolor='k')
axm2.add_collection3d(pc12)

xs = meltshape.X.values
ys = meltshape.Y.values
zs = meltshape.Z.values
positions = np.c_[xs,ys,zs]
pc21 = plotCubeAt(positions, colors=np.tile(azure,(len(xs),1)))
axm1.add_collection3d(pc21)
pc22 = plotCubeAt(positions, colors=np.tile(azure,(len(xs),1)),edgecolor='k')
axm2.add_collection3d(pc22)

axm1.view_init(30,55)
axm2.view_init(30,55)
axm1.set_xlim([xs.min(),xs.max()])
axm1.set_ylim([ys.min(),ys.max()])
axm1.set_zlim([zs.min(),zs.max()])
axm2.set_xlim([xs.min(),xs.max()])
axm2.set_ylim([ys.min(),ys.max()])
axm2.set_zlim([zs.min(),zs.max()])

axm2.set_xlim([-10,0])
axm2.set_ylim([10,20])
axm2.set_zlim([-20,-10])

Rect = plt.Rectangle((1,1),1,1)
axm1.add_artist(Rect)

axm2.axis('off')
axm1.axis('off')
#axm1.set_xlabel('X')
#axm1.set_xticks([])
#axm1.set_yticks([])
#axm1.set_zticks([])
#axm2.set_xticks([])
#axm2.set_yticks([])
#axm2.set_zticks([])

gs.tight_layout(fig,rect=[0.5,0.0,1.0,1.0], h_pad=0.1)
last = fig.add_axes([0.63,0.1,0.04,0.05])
last.axis('off')
last.add_patch(patches.Rectangle((0,0),1,1,fill=False,edgecolor='red',facecolor='none',linewidth=5))
#last.plot((0,1),(0,1),'r')
#last.patch.set_alpha(0)

fig.savefig('combo_grid_refinement.pdf',bbox_inches='tight')
fig.savefig('combo_grid_refinement.png',dpi=600,bbox_inches='tight')
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

import matplotlib as mpl
from cycler import cycler
#mpl.rcParams['axes.prop_cycle'] = cycler(u'color', [u'#d62728', u'#ff7f0e', u'#2ca02c', u'#1f77b4',  u'#9467bd'])
mpl.rcParams['axes.prop_cycle'] = cycler(u'color', [u'#d62728', u'#ff7f0e', u'#2ca02c',u'#75bbfd', u'#1f77b4',  u'#9467bd',u'#ed0dd9'])

from collections import OrderedDict

c = 299792458. # m/s
size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c*1.e-9/f # expected GHz

#freqs={'X':9.6,'Ku':13.6,'Ka':35.6,'W':94,'G':220}
#freqs=OrderedDict([('X',9.6),('Ku',13.6),('Ka',35.6),('W',94),('G',220)])
freqs=OrderedDict([('S',2.8),('C',5.6),('X',9.6),('Ku',13.6),('Ka',35.6),('W',94),('G',220)])

part_size = '5'
mfrac='0_98'

vlin=np.array([50,25,16,12.5])

def ycoord(bounds,alpha):
    return bounds[0]+alpha*(bounds[1]-bounds[0])
def xcooadd(bounds,alpha):
    return alpha*(bounds[1]-bounds[0])
    
def get_line(lines,string):
    return [x for x in lines if string in x][0]

def get_log_numbers(logfilename):
    logfile = open(logfilename,'r')
    lines = logfile.readlines()
    lam_dda = float(get_line(lines,'lambda:').split()[-1])
    dpl = float(get_line(lines,'Dipoles/lambda:').split()[-1])
    Ndipoles = int(get_line(lines,'Total number of occupied dipoles:').split()[-1])
    N1 = int(get_line(lines,'  per domain: 1.').split()[-1])
    N2 = Ndipoles - N1 # this parsing will become a problem when 3 layers will be added
    dda_d = lam_dda/dpl
    n1 = complex(get_line(lines,'refractive index: 1.').split()[-1].replace('i','j'))
    n2 = complex(get_line(lines,'                  2.').split()[-1].replace('i','j'))
    logfile.close()
    return lam_dda, dpl, Ndipoles,N1,N2,dda_d,n1,n2

def get_cross_sections(CrossSecFileName):
    CrossSecfile = open(CrossSecFileName,'r')
    lines = CrossSecfile.readlines()
    Cext = float(get_line(lines,'Cext').split()[-1])
    Qext = float(get_line(lines,'Qext').split()[-1])
    area = Cext/Qext
    Cabs = float(get_line(lines,'Cabs').split()[-1])
    Qabs = float(get_line(lines,'Qabs').split()[-1])
    Csca = Cext - Cabs
    Qsca = Csca/area
    CrossSecfile.close()
    return Cext, Qext, Csca, Qsca, Cabs, Qabs, area

def Ampl2Mueller(S1,S2):
    P11 = 0.5*(S1*S1.conj()+S2*S2.conj()).real
    P12 = 0.5*(S2*S2.conj()-S1*S1.conj()).real
    P33 = 0.5*(S1*S2.conj()+S2*S1.conj()).real
    P34 = (complex(0.5,0.5)*(S1*S2.conj()-S2*S1.conj())).real
    return P11, P12, P33, P34

def plotMueller(angles,data,tags,title,figname=None,normalized=False,logscale=False,ax=None):
    save = False
    if ax is None:
        plt.figure()
        ax = plt.gca()
        save = True
    for x,y,tag in zip(angles,data,tags):
        if normalized:
            cosx = np.cos(np.pi*x/180.0)
            intP = -integrate.trapz(y,cosx)
            y = y/intP
        ax.plot(x,y,label=tag)
    if logscale:
        ax.set_yscale('log')
    if save:
        ax.legend()
        ax.grid()
        ax.set_title(title)
        ax.set_xlabel('Scattering angle')
        ax.set_ylabel(title)
        plt.tight_layout()
        plt.savefig(figname,dpi=300)
        plt.close()
        plt.clf()

deg2rad = lambda angles: np.pi*angles/180.0
rad2deg = lambda angles: 180.0*angles/np.pi
moment  = lambda  x,y,k: integrate.trapz(y*x**k,x)

def plot_field(field,savepath,what='|E|^2',name='intensity',radius=0):
    if what=='|E|^2':
        vmin = 0.0
        vmax = 0.1
    else:
        vmin = -0.3
        vmax = 0.3
    
    intFieldX = field[field.x == min(abs(field.x))]
    intFieldY = field[field.y == min(abs(field.y))]
    intFieldZ = field[field.z == min(abs(field.z))]
    X = sorted(intFieldY.x.drop_duplicates())
    Y = sorted(intFieldZ.y.drop_duplicates())
    Z = sorted(intFieldX.z.drop_duplicates())
    
    circle1 = plt.Circle((0, 0), radius, color='k', fill=False)
    circle2 = plt.Circle((0, 0), radius, color='k', fill=False)
    circle3 = plt.Circle((0, 0), radius, color='k', fill=False)
    
    yi = intFieldX.y.apply(Y.index)
    zi = intFieldX.z.apply(Z.index)
    xv, yv = np.meshgrid(Y, Z)
    zv = np.nan*xv
    zv[zi,yi] = intFieldX[what]
    plt.figure(figsize=(8,8),dpi=300)
    #plt.contourf(1000*xv,1000*yv,1000*zv,cmap='jet')
    plt.pcolormesh(1000*xv,1000*yv,zv,cmap='jet',vmin=vmin,vmax=vmax)
    plt.xlabel('Y')
    plt.ylabel('Z')
    plt.colorbar()
    plt.gca().add_artist(circle1)
    plt.gca().set_aspect(1.0)
    plt.tight_layout()
    plt.savefig(savepath+'X'+name+'.png')
    
    xi = intFieldY.x.apply(X.index)
    zi = intFieldY.z.apply(Z.index)
    xv, yv = np.meshgrid(X, Z)
    zv = np.nan*xv
    zv[zi,xi] = intFieldY[what]
    plt.figure(figsize=(8,8),dpi=300)
    #plt.contourf(1000*xv,1000*yv,1000*zv,cmap='jet')
    plt.pcolormesh(1000*xv,1000*yv,zv,cmap='jet',vmin=vmin,vmax=vmax)
    plt.xlabel('X')
    plt.ylabel('Z')
    plt.colorbar()
    plt.gca().add_artist(circle2)
    plt.gca().set_aspect(1.0)
    plt.tight_layout()
    plt.savefig(savepath+'Y'+name+'.png')
    
    xi = intFieldZ.x.apply(X.index)
    yi = intFieldZ.y.apply(Y.index)
    xv, yv = np.meshgrid(X, Y)
    zv = np.nan*xv
    zv[yi,xi] = intFieldZ[what]
    plt.figure(figsize=(8,8),dpi=300)
    #plt.contourf(1000*xv,1000*yv,1000*zv,cmap='jet')
    plt.pcolormesh(1000*xv,1000*yv,zv,cmap='jet',vmin=vmin,vmax=vmax)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar()
    plt.gca().add_artist(circle3)
    plt.gca().set_aspect(1.0)
    plt.tight_layout()
    plt.savefig(savepath+'Z'+name+'.png')
    plt.close('all')
    plt.clf()
    del xv,yv,zv,xi,zi,intFieldX,intFieldY,intFieldZ,X,Y,Z

def plot_field_Mie(Xf,Yf,Zf,vlim,xlabel,ylabel,savename,radius=0.0):
    vmin,vmax = vlim
    
    circle1 = plt.Circle((0, 0), radius[0], color='k', fill=False)
    circle2 = plt.Circle((0, 0), radius[1], color='k', fill=False)
    
    plt.figure(figsize=(8,8),dpi=300)
    plt.pcolormesh(Xf,Yf,Zf,vmin=vmin,vmax=vmax,cmap='jet')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.colorbar()
    plt.gca().add_artist(circle1)
    plt.gca().add_artist(circle2)
    plt.gca().set_aspect(1.0)
    plt.tight_layout()
    plt.savefig(savename)
    plt.close()
    plt.clf()

def compute_plot_field_Mie(x,m,folder,plane='X'):
    factor = 1.0
    radii = [factor*x[0,0],factor*x[0,1]]
    npts=1000
    scan = np.linspace(-factor*x[0, -1], factor*x[0, -1], npts)
    coord1, coord2 = np.meshgrid(scan, scan)
    coord1.resize(npts*npts)
    coord2.resize(npts*npts)
    coord0 = np.zeros(npts*npts, dtype = np.float64)
            
    if plane == 'X':
        coord = np.vstack((coord0, coord1, coord2)).transpose()
        xlabel, ylabel = 'Y', 'Z'
    elif plane == 'Y':
        coord = np.vstack((coord1, coord0, coord2)).transpose()
        xlabel, ylabel = 'X', 'Z'
    elif plane == 'Z':
        coord = np.vstack((coord1, coord2, coord0)).transpose()
        xlabel, ylabel = 'X', 'Y'

    vlim = [-0.3,0.3]
    terms, ME, MH = fieldnlay(x, m, coord)
    
    Mplot = ME[0,:,0].reshape(npts,npts).real
    savename = folder+'/Mie'+plane+'E0.r'+'.png'
    plot_field_Mie(coord1.reshape((npts,npts)),coord2.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    Mplot = ME[0,:,0].reshape(npts,npts).imag
    savename = folder+'/Mie'+plane+'E0.i'+'.png'
    plot_field_Mie(coord1.reshape((npts,npts)),coord2.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    #Mplot = MH[0,:,0].reshape(npts,npts).real
    #savename = folder+'/Mie'+plane+'H0.r'+'.png'
    #plot_field_Mie(coordX.reshape((npts,npts)),coordZ.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    #Mplot = MH[0,:,0].reshape(npts,npts).imag
    #savename = folder+'/Mie'+plane+'H0.i'+'.png'
    #plot_field_Mie(coordX.reshape((npts,npts)),coordZ.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    
    Mplot = ME[0,:,1].reshape(npts,npts).real
    savename = particle_folder+'/Mie'+plane+'E1.r'+'.png'
    plot_field_Mie(coord1.reshape((npts,npts)),coord2.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    Mplot = ME[0,:,1].reshape(npts,npts).imag
    savename = particle_folder+'/Mie'+plane+'E1.i'+'.png'
    plot_field_Mie(coord1.reshape((npts,npts)),coord2.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    #Mplot = MH[0,:,1].reshape(npts,npts).real
    #savename = particle_folder+'/Mie'+plane+'H1.r'+'.png'
    #plot_field_Mie(coordX.reshape((npts,npts)),coordZ.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    #Mplot = MH[0,:,1].reshape(npts,npts).imag
    #savename = particle_folder+'/Mie'+plane+'H1.i'+'.png'
    #plot_field_Mie(coordX.reshape((npts,npts)),coordZ.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)

    Mplot = ME[0,:,2].reshape(npts,npts).real
    savename = particle_folder+'/Mie'+plane+'E2.r'+'.png'
    plot_field_Mie(coord1.reshape((npts,npts)),coord2.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    Mplot = ME[0,:,2].reshape(npts,npts).imag
    savename = particle_folder+'/Mie'+plane+'E2.i'+'.png'
    plot_field_Mie(coord1.reshape((npts,npts)),coord2.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    #Mplot = MH[0,:,2].reshape(npts,npts).real
    #savename = particle_folder+'/Mie'+plane+'H2.r'+'.png'
    #plot_field_Mie(coordX.reshape((npts,npts)),coordZ.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    #Mplot = MH[0,:,2].reshape(npts,npts).imag
    #savename = particle_folder+'/Mie'+plane+'H2.i'+'.png'
    #plot_field_Mie(coordX.reshape((npts,npts)),coordZ.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)
    
    Mplot = (ME[0,:,0]*ME[0,:,0].conj()+ME[0,:,1]*ME[0,:,1].conj()+ME[0,:,2]*ME[0,:,2].conj()).reshape(npts,npts).real
    savename = particle_folder+'/Mie'+plane+'intensity'+'.png'
    vlim = [0.0,0.1]
    plot_field_Mie(coord1.reshape((npts,npts)),coord2.reshape((npts,npts)),Mplot,vlim,xlabel,ylabel,savename,radius=radii)

MIEdict = OrderedDict()
DDAdict = OrderedDict()

Pdda = OrderedDict()
Pmie = OrderedDict()

data_folder = '/data/optimice/scattering_databases/melting_sphere/test_resolution/'+str(part_size)+'mm/'+mfrac
for freq_str in freqs.keys():#[0:1]:
    f = freqs[freq_str]
    lam = fr2lam(f)
    k2 = 4.*(np.pi/lam)**2
    print(f,lam)
    particles_folders = glob(data_folder+'/'+freq_str+'/*')
    DDA = pd.DataFrame(index=range(len(particles_folders)),columns=['Dratio','Qext','Qabs','Qsca','Qbk','g','ssa','dipole_spacing'] )
    MIE = pd.DataFrame(index=range(len(particles_folders)),columns=['Dratio','Qext','Qabs','Qsca','Qbk','g','ssa','dipole_spacing'] )
    plt.plot()
    ax = plt.gca()

    P11dda = pd.DataFrame(index=np.linspace(0.0,180.0,721),columns=np.arange(1,len(particles_folders)+1,1))
    P11mie = pd.DataFrame(index=np.linspace(0.0,180.0,721),columns=np.arange(1,len(particles_folders)+1,1))

    particles_folders = sorted(particles_folders)#[0:10] ##
    for particle_folder,i in zip(particles_folders,range(len(particles_folders))):
        print(particle_folder)
        #dipoles = pd.read_csv(particle_folder+'/coated.geom',sep=' ',header=4,names=['X','Y','Z','M'])
        mueller = pd.read_csv(particle_folder+'/mueller',sep=' ',index_col='theta')

        logfilename = particle_folder+'/log'
        lam_dda, dpl, Ndipoles,N1,N2,dda_d,n1,n2 = get_log_numbers(logfilename)
        #print(lam_dda, dpl, Ndipoles,N1,N2,dda_d,n1,n2)

        CrossSecFileName = particle_folder+'/CrossSec-Y'
        Cext, Qext, Csca, Qsca, Cabs, Qabs, area = get_cross_sections(CrossSecFileName)

        back = mueller.loc[180]
        Cbck = 2.0*np.pi*(back.s11+back.s22+2.0*back.s12)/k2
        Qbck = Cbck/area
        volume_ratio = float(N1)/float(Ndipoles)

        Ncomput = 1
        Nlayers = 2
        x = np.ndarray((Ncomput,Nlayers),dtype=np.float64)
        m = np.ndarray((Ncomput,Nlayers),dtype=np.complex128)
        
        #outer_size = 0.5*0.001*float(part_size)
        outer_size=dda_d*np.cbrt((6.0*Ndipoles/np.pi))*0.5 
        outer_x = size2x(outer_size,lam_dda)
        #inner_size = outer_size*float(mfrac[2:])*0.01
        inner_size=dda_d*np.cbrt((6.0*N2/np.pi))*0.5
        inner_x = size2x(inner_size,lam_dda)
        print(outer_size,inner_size)
        x[:,0] = inner_x
        x[:,1] = outer_x
        m[:,0] = n2
        m[:,1] = n1
        #print(x)
        #print(m)
        thetas = deg2rad(mueller.index.values)
        g = moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,1)/moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,0)
        DDA.loc[i] = inner_size/outer_size, Qext, Qabs, Qsca, Qbck, g, Qsca/Qext, dda_d
        
        terms, MQe, MQs, MQa, MQb, MQp, Mg, Mssa, S1, S2 = scattnlay(x,m,theta=thetas)
        MIE.loc[i] = inner_size/outer_size, MQe[0], MQa[0], MQs[0], MQb[0], Mg[0], Mssa[0], dda_d
        P11, P12, P33, P34 = Ampl2Mueller(S1[0],S2[0])
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s11.values,P11],tags=['DDA','MIE'],title='P11',figname=particle_folder+'/P11.png',logscale=True)
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s12.values,P12],tags=['DDA','MIE'],title='P12',figname=particle_folder+'/P12.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s33.values,P33],tags=['DDA','MIE'],title='P33',figname=particle_folder+'/P33.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s34.values,P34],tags=['DDA','MIE'],title='P34',figname=particle_folder+'/P34.png')
        #print(inner_size/outer_size,(outer_size-inner_size)*1.0e6,MQe[0]/Qext,MQs[0]/Qsca,MQa[0]/Qabs,MQb[0]/Qbck,Mg[0]/g)
        P11dda.loc[:,i+1] = mueller.s11
        P11mie.loc[:,i+1] = P11
#        try:
#            print('intfield')
#            intField = pd.read_csv(particle_folder+'/IntField-y.dat',sep=' ')
#            plot_field(intField,savepath=particle_folder+'/',what='|E|^2',name='intensity',radius=inner_size*1000)
#            plot_field(intField,savepath=particle_folder+'/',what='Ex.r',name='Exr',radius=inner_size*1000)
#            plot_field(intField,savepath=particle_folder+'/',what='Ex.i',name='Exr',radius=inner_size*1000)
#            plot_field(intField,savepath=particle_folder+'/',what='Ey.r',name='Eyr',radius=inner_size*1000)
#            plot_field(intField,savepath=particle_folder+'/',what='Ey.i',name='Eyi',radius=inner_size*1000)
#            plot_field(intField,savepath=particle_folder+'/',what='Ez.r',name='Ezr',radius=inner_size*1000)
#            plot_field(intField,savepath=particle_folder+'/',what='Ez.i',name='Ezi',radius=inner_size*1000)
#    
#            print('fieldnalay')
#            compute_plot_field_Mie(x,m,particle_folder,plane='X')
#            compute_plot_field_Mie(x,m,particle_folder,plane='Y')
#            compute_plot_field_Mie(x,m,particle_folder,plane='Z')
#        except:
#            print('cannot do internal fields')
#            pass
        #del intField, P11, P12, P33, P34, S1, S2, mueller, thetas
        print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        gc.collect()
        print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    DDA.sort_values('dipole_spacing',inplace=True)
    MIE.sort_values('dipole_spacing',inplace=True)
    DDAdict[freq_str] = DDA
    MIEdict[freq_str] = MIE
    Pdda[freq_str] = P11dda
    Pmie[freq_str] = P11mie

#%%

def plot_comparison(dda_data,mie_data,quantity,folder,ax=None):
    if ax is None:
        plt.figure()
        ax = plt.gca()
        figname = folder + '/' + quantity + 'relative_diff.png'
        for f in dda_data.keys():
            DFdda = dda_data[f]
            DFmie = mie_data[f]
            #lay_thick_mie = Dout*(1.0 - DFmie['Dratio'])
            #lay_thick_dda = Dout*(1.0 - DFdda['Dratio'])
            ax.plot(1.0e6*DFmie['dipole_spacing'],DFmie[quantity],label=f)
            ax.scatter(1.0e6*DFdda['dipole_spacing'],DFdda[quantity].values,marker='.')
        ax.set_xlabel('dipole spacing [um]')
        ax.set_ylabel(quantity)
        ax.grid()
        ax.legend(loc=2)
        print(figname)
        plt.savefig(figname,dpi=300)
        plt.close()
        plt.clf()
    else:
        for f in dda_data.keys():
            DFmie = mie_data[f]
            #lay_thick_mie = Dout*(1.0 - DFmie['Dratio'])
            ax.plot(1.0e6*DFmie['dipole_spacing'],DFmie[quantity],label=f)
        ax.set_color_cycle(None)
        for f in dda_data.keys():
            DFdda = dda_data[f]
            #lay_thick_dda = Dout*(1.0 - DFdda['Dratio'])
            ax.plot(1.0e6*DFdda['dipole_spacing'],DFdda[quantity].values,marker='.',linewidth=0)
        i=0
        for vl in vlin:
            vline2d = ax.axvline(vl,ls='--',c='k')
            dataline = vline2d.get_data()
            i = i+1
            print(ax.get_ybound())
            ax.text(dataline[0][0]-xcooadd(ax.get_xbound(),0.05),ycoord(ax.get_ybound(),0.2),str(i))
        ax.grid()

def plot_difference(dda_data,mie_data,quantity,folder,ax=None):
    if ax is None:
        plt.figure()
        ax = plt.gca()
        figname = folder + '/' + quantity + 'diff.png'
        for f in dda_data.keys():
            DFdda = dda_data[f]
            DFmie = mie_data[f]
            ax.plot(1.0e6*DFmie['dipole_spacing'],DFmie[quantity]-DFdda[quantity].values,label=f)
        ax.set_xlabel('dipole spacing [um]')
        ax.set_ylabel(quantity + ' absolute difference')
        ax.set_xscale('log')
        ax.grid()
        ax.legend(loc=2)
        print(figname)
        plt.savefig(figname,dpi=300)
        plt.close()
        plt.clf()
    else:
        for f in dda_data.keys():
            DFdda = dda_data[f]
            DFmie = mie_data[f]
            ax.plot(1.0e6*DFmie['dipole_spacing'],(DFmie[quantity]-DFdda[quantity].values),label=f)
        i = 0
        for vl in vlin:
                    vline2d = ax.axvline(vl,ls='--',c='k')
                    dataline = vline2d.get_data()
                    i = i+1
                    print(ax.get_ybound())
                    ax.text(dataline[0][0]-xcooadd(ax.get_xbound(),0.05),ycoord(ax.get_ybound(),0.2),str(i))
        ax.grid()

    
def plot_relative_difference(dda_data,mie_data,quantity,folder,ax=None):
    if ax is None:
        plt.figure()
        ax = plt.gca()
        figname = folder + '/' + quantity + 'relative_diff.png'
        for f in dda_data.keys():
            DFdda = dda_data[f]
            DFmie = mie_data[f]
            ax.plot(1.0e6*DFmie['dipole_spacing'],(DFmie[quantity]-DFdda[quantity].values)/DFmie[quantity].values,label=f)
        ax.set_xlabel('dipole spacing [um]')
        ax.set_ylabel(quantity + ' relative difference')
        ax.set_xscale('log')
        ax.grid()
        ax.legend(loc=2)
        print(figname)
        plt.savefig(figname,dpi=300)
        plt.close()
        plt.clf()
    else:
        for f in dda_data.keys():
            DFdda = dda_data[f]
            DFmie = mie_data[f]
            ax.plot(1.0e6*DFmie['dipole_spacing'],100*(DFmie[quantity]-DFdda[quantity].values)/DFmie[quantity].values,label=f)
        i=0
        for vl in vlin:
            vline2d = ax.axvline(vl,ls='--',c='k')
            dataline = vline2d.get_data()
            i = i+1
            print(ax.get_ybound())
            ax.text(dataline[0][0]-xcooadd(ax.get_xbound(),0.05),ycoord(ax.get_ybound(),0.8),str(i))
        ax.grid()

plot_comparison(DDAdict,MIEdict,quantity='Qext',folder=data_folder)
plot_comparison(DDAdict,MIEdict,quantity='Qsca',folder=data_folder)
plot_comparison(DDAdict,MIEdict,quantity='Qabs',folder=data_folder)
plot_comparison(DDAdict,MIEdict,quantity='Qbk',folder=data_folder)
plot_comparison(DDAdict,MIEdict,quantity='g',folder=data_folder)
plot_comparison(DDAdict,MIEdict,quantity='ssa',folder=data_folder)

plot_difference(DDAdict,MIEdict,quantity='Qext',folder=data_folder)
plot_difference(DDAdict,MIEdict,quantity='Qsca',folder=data_folder)
plot_difference(DDAdict,MIEdict,quantity='Qabs',folder=data_folder)
plot_difference(DDAdict,MIEdict,quantity='Qbk',folder=data_folder)
plot_difference(DDAdict,MIEdict,quantity='g',folder=data_folder)
plot_difference(DDAdict,MIEdict,quantity='ssa',folder=data_folder)

plot_relative_difference(DDAdict,MIEdict,quantity='Qext',folder=data_folder)
plot_relative_difference(DDAdict,MIEdict,quantity='Qsca',folder=data_folder)
plot_relative_difference(DDAdict,MIEdict,quantity='Qabs',folder=data_folder)
plot_relative_difference(DDAdict,MIEdict,quantity='Qbk',folder=data_folder)
plot_relative_difference(DDAdict,MIEdict,quantity='g',folder=data_folder)
plot_relative_difference(DDAdict,MIEdict,quantity='ssa',folder=data_folder)


f,((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,sharex=True,figsize=(7,7))
plot_comparison(DDAdict,MIEdict,quantity='Qsca',folder=data_folder,ax=ax1)
plot_comparison(DDAdict,MIEdict,quantity='Qabs',folder=data_folder,ax=ax3)
plot_comparison(DDAdict,MIEdict,quantity='Qbk',folder=data_folder,ax=ax5)
plot_comparison(DDAdict,MIEdict,quantity='g',folder=data_folder,ax=ax7)

plot_relative_difference(DDAdict,MIEdict,quantity='Qsca',folder=data_folder,ax=ax2)
plot_relative_difference(DDAdict,MIEdict,quantity='Qabs',folder=data_folder,ax=ax4)
plot_relative_difference(DDAdict,MIEdict,quantity='Qbk',folder=data_folder,ax=ax6)
plot_difference(DDAdict,MIEdict,quantity='g',folder=data_folder,ax=ax8)

ax2.legend(loc=4,ncol=1)
ax1.set_ylabel('Q$_{sca}$')
ax2.set_ylabel('$\Delta$Q$_{sca}$/$Q_{sca}^{mie}$     [%]')
ax3.set_ylabel('Q$_{abs}$')
ax4.set_ylabel('$\Delta$Q$_{abs}$/$Q_{abs}^{mie}$     [%]')
ax5.set_ylabel('Q$_{bk}$')
ax6.set_ylabel('$\Delta$Q$_{bk}$/$Q_{bk}^{mie}$      [%]',labelpad=0)
ax7.set_ylabel('g')
ax8.set_ylabel('$\Delta$ g')
ax8.set_xlabel('dipole spacing [um]')
ax7.set_xlabel('dipole spacing [um]')
#f.savefig(data_folder + '/'+'8_panel.png',dpi=300)
#f.savefig(data_folder + '/'+'8_panel.pdf',dpi=300)
f.suptitle(part_size+'mm sphere 50um water refine discretization and shape',y=0.999999)
f.tight_layout()
f.savefig(data_folder + '/'+'resolution_8panel.png',dpi=300,bbox_inches='tight')
f.savefig(data_folder + '/'+'resolution_8panel.pdf',dpi=300,bbox_inches='tight')

def plot_phase3(ax,MIE,DDA,band,op='='):
    angles = MIE[band].index.values
    mie = MIE[band].loc[:,1].values
    coarse = DDA[band].loc[:,3].values
    refined = DDA[band].loc[:,1].values
    if op == '=':
        plotMueller(angles=[angles,angles,angles],
                    data=[mie,coarse,refined],logscale=True,
                    tags=['MIE','DDA coarse','DDA refined'], title='P11',ax=ax)
    elif op == '-':
        plotMueller(angles=[angles,angles,angles],#logscale=True,
                    #data=[abs(mie-coarse),abs(mie-refined),abs(coarse-refined)],
                    data=[(mie-coarse)/mie,(mie-refined)/mie,(coarse-refined)/mie],
                    tags=['total','shape','discretization'], title='error',ax=ax)


#f, ((ax11,ax12),(ax21,ax22),(ax31,ax32),(ax41,ax42),(ax51,ax52)) = plt.subplots(5,2,sharex=True,figsize=(9.5,9.5))
f, ((ax11,ax12),(ax21,ax22),(ax31,ax32)) = plt.subplots(3,2,sharex=True,figsize=(9.5,9.5))
plot_phase3(ax11,Pmie,Pdda,'S','=')
plot_phase3(ax21,Pmie,Pdda,'X','=')
plot_phase3(ax31,Pmie,Pdda,'G','=')
#plot_phase3(ax41,Pmie,Pdda,'W','=')
#plot_phase3(ax51,Pmie,Pdda,'G','=')
plot_phase3(ax12,Pmie,Pdda,'S','-')
plot_phase3(ax22,Pmie,Pdda,'X','-')
plot_phase3(ax32,Pmie,Pdda,'G','-')
#plot_phase3(ax42,Pmie,Pdda,'W','-')
#plot_phase3(ax52,Pmie,Pdda,'G','-')
ax11.legend()
ax12.legend()
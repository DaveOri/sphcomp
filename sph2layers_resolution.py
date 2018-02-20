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
#sys.path.append('/home/dori/pymiecoated')
#from pymiecoated import Mie

c = 299792458. # m/s
size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c*1.e-9/f # expected GHz

freqs={'X':9.6,'Ku':13.6,'Ka':35.6,'W':94,'G':220}

part_size = '10'
#from sys import argv
#scriptname, part_size = argv

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

def plotMueller(angles,data,tags,title,figname):
    plt.figure()
    ax = plt.gca()
    for x,y,tag in zip(angles,data,tags):
        ax.plot(x,y,label=tag)
    ax.legend()
    ax.grid()
    ax.set_title(title)
    ax.set_xlabel('Scattering angle')
    ax.set_ylabel(title)
    ax.legend()
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

MIEdict = {}
DDAdict = {}

data_folder = '/data/optimice/scattering_databases/melting_sphere/test_resolution/'+str(part_size)+'mm/0_99'
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

    particles_folders = sorted(particles_folders)#[0:10] ##
    for particle_folder,i in zip(particles_folders,range(len(particles_folders))):
        print(particle_folder)
        #dipoles = pd.read_csv(particle_folder+'/coated.geom',sep=' ',header=4,names=['X','Y','Z','M'])
        mueller = pd.read_csv(particle_folder+'/mueller',sep=' ',index_col='theta')

        logfilename = particle_folder+'/log'
        lam_dda, dpl, Ndipoles,N1,N2,dda_d,n1,n2 = get_log_numbers(logfilename)
        print(lam_dda, dpl, Ndipoles,N1,N2,dda_d,n1,n2)

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
        
        outer_size = dda_d*np.cbrt((6.0*Ndipoles/np.pi))*0.5 
        outer_x = size2x(outer_size,lam_dda)
        inner_size = dda_d*np.cbrt((6.0*N2/np.pi))*0.5
        inner_x = size2x(inner_size,lam_dda)
        x[:,0] = inner_x
        x[:,1] = outer_x
        m[:,0] = n2
        m[:,1] = n1
        print(x)
        print(m)
        thetas = deg2rad(mueller.index.values)
        g = moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,1)/moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,0)
        DDA.loc[i] = inner_size/outer_size, Qext, Qabs, Qsca, Qbck, g, Qsca/Qext, dda_d
        
        terms, MQe, MQs, MQa, MQb, MQp, Mg, Mssa, S1, S2 = scattnlay(x,m,theta=thetas)
        MIE.loc[i] = inner_size/outer_size, MQe[0], MQa[0], MQs[0], MQb[0], Mg[0], Mssa[0], dda_d
        P11, P12, P33, P34 = Ampl2Mueller(S1[0],S2[0])
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s11.values,P11],tags=['DDA','MIE'],title='P11',figname=particle_folder+'/P11.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s12.values,P12],tags=['DDA','MIE'],title='P12',figname=particle_folder+'/P12.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s33.values,P33],tags=['DDA','MIE'],title='P33',figname=particle_folder+'/P33.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s34.values,P34],tags=['DDA','MIE'],title='P34',figname=particle_folder+'/P34.png')
        #print(inner_size/outer_size,(outer_size-inner_size)*1.0e6,MQe[0]/Qext,MQs[0]/Qsca,MQa[0]/Qabs,MQb[0]/Qbck,Mg[0]/g)
        
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
#            pass
        #del intField, P11, P12, P33, P34, S1, S2, mueller, thetas
        print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        gc.collect()
        print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
    DDA.sort_values('dipole_spacing',inplace=True)
    MIE.sort_values('dipole_spacing',inplace=True)
    DDAdict[freq_str] = DDA
    MIEdict[freq_str] = MIE

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
        ax.grid()
        ax.set_ylim([-20,20])

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


f,((ax1,ax2),(ax3,ax4),(ax5,ax6),(ax7,ax8)) = plt.subplots(4,2,sharex=True,figsize=(14,12))
plot_comparison(DDAdict,MIEdict,quantity='Qsca',folder=data_folder,ax=ax1)
plot_comparison(DDAdict,MIEdict,quantity='Qabs',folder=data_folder,ax=ax3)
plot_comparison(DDAdict,MIEdict,quantity='Qbk',folder=data_folder,ax=ax5)
plot_comparison(DDAdict,MIEdict,quantity='g',folder=data_folder,ax=ax7)

plot_relative_difference(DDAdict,MIEdict,quantity='Qsca',folder=data_folder,ax=ax2)
plot_relative_difference(DDAdict,MIEdict,quantity='Qabs',folder=data_folder,ax=ax4)
plot_relative_difference(DDAdict,MIEdict,quantity='Qbk',folder=data_folder,ax=ax6)
plot_relative_difference(DDAdict,MIEdict,quantity='g',folder=data_folder,ax=ax8)

ax2.legend(ncol=5)
ax1.set_ylabel('Q$_{sca}$')
ax2.set_ylabel('$\Delta$Q$_{sca}$     [%]')
ax3.set_ylabel('Q$_{abs}$')
ax4.set_ylabel('$\Delta$Q$_{abs}$     [%]')
ax5.set_ylabel('Q$_{bk}$')
ax6.set_ylabel('$\Delta$Q$_{bk}$     [%]')
ax7.set_ylabel('g')
ax8.set_ylabel('$\Delta$ g     [%]')
ax8.set_xlabel('dipole spacing [um]')
ax7.set_xlabel('dipole spacing [um]')
f.savefig(data_folder + '/'+'8_panel.png',dpi=300)
f.savefig(data_folder + '/'+'8_panel.pdf',dpi=300)
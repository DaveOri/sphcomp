import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from scipy import integrate
#from matplotlib.mlab import griddata

#sys.path.append('/work/DBs/scattnlay')
#from scattnlay import scattnlay
#from scattnlay import fieldnlay
#sys.path.append('/home/dori/pymiecoated')
#from pymiecoated import Mie

c = 299792458. # m/s
size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c*1.e-9/f # expected GHz

freqs={'X':9.6,'Ku':13.6,'Ka':35.6,'W':94}
part_size = '10'

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
    P22 = 0.5*(S2*S2.conj()-S1*S1.conj()).real
    P33 = 0.5*(S1*S2.conj()+S2*S1.conj()).real
    P34 = (complex(0.5,0.5)*(S1*S2.conj()-S2*S1.conj())).real
    return P11, P22, P33, P34

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
    plt.savefig(figname)
    plt.close()

deg2rad = lambda angles: np.pi*angles/180.0
rad2deg = lambda angles: 180.0*angles/np.pi
moment  = lambda  x,y,k: integrate.trapz(y*x**k,x)

def plot_field(field,savepath,what='|E|^2',name='intensity',radius=0):
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
    zv[yi,zi] = intFieldX[what]
    plt.figure(figsize=(8,8),dpi=300)
    plt.contourf(1000*xv,1000*yv,1000*zv,cmap='jet')
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
    zv[xi,zi] = intFieldY[what]
    plt.figure(figsize=(8,8),dpi=300)
    plt.contourf(1000*xv,1000*yv,1000*zv,cmap='jet')
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
    zv[xi,yi] = intFieldZ[what]
    plt.figure(figsize=(8,8),dpi=300)
    plt.contourf(1000*xv,1000*yv,1000*zv,cmap='jet')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar()
    plt.gca().add_artist(circle3)
    plt.gca().set_aspect(1.0)
    plt.tight_layout()
    plt.savefig(savepath+'Z'+name+'.png')
    plt.close('all')
    

data_folder = './'+str(part_size)+'mm'
for freq_str in freqs.keys():
    f = freqs[freq_str]
    lam = fr2lam(f)
    k2 = 4.*(np.pi/lam)**2
    print(f,lam)
    particles_folders = glob(data_folder+'/'+freq_str+'/*')
    DDA = pd.DataFrame(index=range(len(particles_folders)),columns=['Dratio','Qext','Qabs','Qsca','Qbk','g','ssa'] )
    MIE = pd.DataFrame(index=range(len(particles_folders)),columns=['Dratio','Qext','Qabs','Qsca','Qbk','g','ssa'] )
    plt.plot()
    ax = plt.gca()

#    particles_folders = particles_folders[0:5]
    for particle_folder,i in zip(particles_folders,range(len(particles_folders))):
        print(particle_folder)
        dipoles = pd.read_csv(particle_folder+'/coated.geom',sep=' ',header=4,names=['X','Y','Z','M'])
        mueller = pd.read_csv(particle_folder+'/mueller',sep=' ',index_col='theta')

        logfilename = particle_folder+'/log'
        lam_dda, dpl, Ndipoles,N1,N2,dda_d,n1,n2 = get_log_numbers(logfilename)

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
        outer_x = size2x(outer_size,lam)
        inner_size = dda_d*np.cbrt((6.0*N2/np.pi))*0.5
        inner_x = size2x(inner_size,lam)
        x[:,0] = inner_x
        x[:,1] = outer_x
        m[:,0] = n2
        m[:,1] = n1
        thetas = deg2rad(mueller.index.values)
        g = moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,1)/moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,0)
        DDA.loc[i] = inner_size/outer_size, Qext, Qsca, Qabs, Qbck, g, Qsca/Qext
        
        #terms, MQe, MQs, MQa, MQb, MQp, Mg, Mssa, S1, S2 = = scattnlay(x,m,theta=thetas)
        #MIE.loc[i] = inner_size/outer_size, MQe, MQs, MQa, MQb, Mg, Mssa
        #P11, P12, P33, P34 = Ampl2Mueller(S1,S2)
        #plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s11.values,P11],tags=['DDA','MIE'],title='P11',figname=particle_folder+'/P11.png')
        #plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s12.values,P12],tags=['DDA','MIE'],title='P12',figname=particle_folder+'/P12.png')
        #plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s33.values,P33],tags=['DDA','MIE'],title='P33',figname=particle_folder+'/P33.png')
        #plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s34.values,P34],tags=['DDA','MIE'],title='P34',figname=particle_folder+'/P34.png')
        #print(inner_size/outer_size,(outer_size-inner_size)*1.0e6,MQe/Qext,MQs/Qsca,MQa/Qabs,MQb/Qbck,g/Mg)
        
        try:
            intField = pd.read_csv(particle_folder+'/IntField-Y',sep=' ')
            plot_field(intField,savepath=particle_folder+'/',what='|E|^2',name='intensity',radius=inner_size*1000)
            plot_field(intField,savepath=particle_folder+'/',what='Ex.r',name='Exr',radius=inner_size*1000)
            plot_field(intField,savepath=particle_folder+'/',what='Ex.i',name='Exr',radius=inner_size*1000)
            plot_field(intField,savepath=particle_folder+'/',what='Ey.r',name='Eyr',radius=inner_size*1000)
            plot_field(intField,savepath=particle_folder+'/',what='Ey.i',name='Eyi',radius=inner_size*1000)
            plot_field(intField,savepath=particle_folder+'/',what='Ez.r',name='Ezr',radius=inner_size*1000)
            plot_field(intField,savepath=particle_folder+'/',what='Ez.i',name='Ezi',radius=inner_size*1000)
            #factor=2.5
            #scan = np.linspace(-factor*x[0, 2], factor*x[0, 2], npts)
            #coordX, coordZ = np.meshgrid(scan, scan)
            #coordX.resize(npts*npts)
            #coordZ.resize(npts*npts)
            #coordY = np.zeros(npts*npts, dtype = np.float64)
            #coord = np.vstack((coordX, coordY, coordZ)).transpose()
            #terms, ME, MH = fieldnlay(x, m, coord)
        except:
            pass

    DDA.sort('Dratio',inplace=True)
#    MIE.sort('Dratio',inplace=True)
#    ax.plot(MIE.Dratio,MIE.Qabs)
    ax.scatter(DDA.Dratio,DDA.Qabs)
    ax.grid()
#    plt.show()

                


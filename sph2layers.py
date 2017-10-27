import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from glob import glob
from scipy import integrate

sys.path.append('/work/DBs/scattnlay')
sys.path.append('/home/dori/pymiecoated')
from scattnlay import scattnlay
#from pymiecoated import Mie

c = 299792458. # m/s
size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c*1.e-9/f # expected GHz

freqs={'X':9.6,'Ku':13.6,'Ka':35.6,'W':94}
part_size = '4'

def get_line(lines,string):
    return [x for x in lines if string in x][0]

def get_log_numbers(logfilename):
    logfile = open(logfilename,'r')
    lines = logfile.readlines()
    #lam_lines = get_line(lines,'lambda:')
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

    particles_folders = particles_folders[0:5]
    for particle_folder,i in zip(particles_folders,range(len(particles_folders))):
        print(particle_folder)
        dipoles = pd.read_csv(particle_folder+'/coated.geom',sep=' ',header=4,names=['X','Y','Z','M'])
        mueller = pd.read_csv(particle_folder+'/mueller',sep=' ',index_col='theta')
#        intField = pd.read_csv(particle_folder+'/IntField-Y',sep=' ')
#        print(intField.iloc[0])

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
        #betas = np.linspace(0,2*np.pi,360) # compute S at this scattering angles, NOW [radians]
        outer_size = dda_d*np.cbrt((6.0*Ndipoles/np.pi))*0.5 
        outer_x = size2x(outer_size,lam)
        inner_size = dda_d*np.cbrt((6.0*N2/np.pi))*0.5
    	inner_x = size2x(inner_size,lam)
        x[:,0] = inner_x
        x[:,1] = outer_x
        m[:,0] = n2
        m[:,1] = n1
        thetas = deg2rad(mueller.index.values)
    	results = scattnlay(x,m,theta=thetas)
        g = moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,1)/moment(np.cos(deg2rad(mueller.index.values)),mueller.s11.values,0)
        DDA.loc[i] = inner_size/outer_size, Qext, Qsca, Qabs, Qbck, g, Qsca/Qext
        MIE.loc[i] = inner_size/outer_size, results[1][0], results[2][0], results[3][0], results[4][0], results[6][0], results[7][0]
        S1 = results[8][0]
        S2 = results[9][0]
        P11, P12, P33, P34 = Ampl2Mueller(S1,S2)
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s11.values,P11],tags=['DDA','MIE'],title='P11',figname=particle_folder+'/P11.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s12.values,P12],tags=['DDA','MIE'],title='P12',figname=particle_folder+'/P12.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s33.values,P33],tags=['DDA','MIE'],title='P33',figname=particle_folder+'/P33.png')
        plotMueller(angles=[mueller.index.values,rad2deg(thetas)],data=[mueller.s34.values,P34],tags=['DDA','MIE'],title='P34',figname=particle_folder+'/P34.png')
        print(inner_size/outer_size,(outer_size-inner_size)*1.0e6,results[1][0]/Qext,results[2][0]/Qsca,results[3][0]/Qabs,results[4][0]/Qbck,g/results[6][0])
        #print(results)
        print(len(results))

    DDA.sort('Dratio',inplace=True)
    MIE.sort('Dratio',inplace=True)
    ax.plot(MIE.Dratio,MIE.Qabs)
    ax.scatter(DDA.Dratio,DDA.Qabs)
    ax.grid()
#    plt.show()

                


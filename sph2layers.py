import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
import subprocess
from copy import deepcopy
import pandas as pd
from glob import glob

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
    return Cext, Qext, Csca, Qsca, Cabs, Qabs, area

data_folder = './'+str(part_size)+'mm'
for freq_str in freqs.keys():
    f = freqs[freq_str]
    lam = fr2lam(f)
    k2 = 4.*(np.pi/lam)**2
    print(f,lam)
    particles_folders = glob(data_folder+'/'+freq_str+'/*')
    DDA = pd.DataFrame(index=range(len(particles_folders)),columns=['Dratio','Qext','Qabs','Qsca','Qbk'] )
    MIE = pd.DataFrame(index=range(len(particles_folders)),columns=['Dratio','Qext','Qabs','Qsca','Qbk'] )
    plt.plot()
    ax = plt.gca()

    particles_folders = particles_folders[0:1]
    for particle_folder,i in zip(particles_folders,range(len(particles_folders))):
        dipoles = pd.read_csv(particle_folder+'/coated.geom',sep=' ',header=4,names=['X','Y','Z','M'])
        mueller = pd.read_csv(particle_folder+'/mueller',sep=' ',index_col='theta')
        intField = pd.read_csv(particle_folder+'/IntField-Y',sep=' ')
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
    	results = scattnlay(x,m)
        DDA.loc[i,'Dratio'] = inner_size/outer_size
        DDA.loc[i,'Qext']   = Qext
        DDA.loc[i,'Qsca']   = Qsca
        DDA.loc[i,'Qabs']   = Qabs
        DDA.loc[i,'Qbck']   = Qbck
        MIE.loc[i,'Dratio'] = inner_size/outer_size
        MIE.loc[i,'Qext']   = results[1][0]
        MIE.loc[i,'Qsca']   = results[2][0]
        MIE.loc[i,'Qabs']   = results[3][0]
        MIE.loc[i,'Qbck']   = results[4][0]
        print(inner_size/outer_size,(outer_size-inner_size)*1.0e6,results[1][0]/Qext,results[2][0]/Qsca,results[3][0]/Qabs,results[4][0]/Qbck)
    DDA.sort('Dratio',inplace=True)
    MIE.sort('Dratio',inplace=True)
    ax.plot(MIE.Dratio,MIE.Qabs)
    ax.scatter(DDA.Dratio,DDA.Qabs)
    ax.grid()
#    plt.show()

                


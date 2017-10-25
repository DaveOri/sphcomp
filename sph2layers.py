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
#from scattnlay import scattnlay
#from pymiecoated import Mie

c = 299792458. # m/s
size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c*1.e-9/f # expected GHz

freqs={'X':9.6,'Ku':13.6,'Ka':35.6,'W':94}
part_size = '4'

def get_line(lines,string):
    return [x for x in lines if string in x][0]

data_folder = './'+str(part_size)+'mm'
for freq_str in freqs.keys()[0]:
    f = freqs[freq_str]
    lam = fr2lam(f)
    print(f,lam)
    particles_folders = glob(data_folder+'/'+freq_str+'/*')
    for particle_folder in particles_folders[0:1]:
        print(particle_folder)
        dipoles = pd.read_csv(particle_folder+'/coated.geom',sep=' ',header=4,names=['X','Y','Z','M'])
        print(dipoles.iloc[0])
        mueller = pd.read_csv(particle_folder+'/mueller',sep=' ',index_col='theta')
        print(mueller.iloc[0])
#        intField = pd.read_csv(particle_folder+'/IntField-Y',sep=' ')
#        print(intField.iloc[0])
        logfile = open(particle_folder+'/log','r')
        lines = logfile.readlines()
        lam_lines = get_line(lines,'lambda:')
        lam_dda = float(get_line(lines,'lambda:').split()[-1])
        dpl = float(get_line(lines,'Dipoles/lambda:').split()[-1])
        Ndipoles = int(get_line(lines,'Total number of occupied dipoles:').split()[-1])
        N1 = int(get_line(lines,'  per domain: 1.').split()[-1])
        N2 = Ndipoles - N1 # this parsing will become a problem when 3 layers will be added
        print(N1,N2)
        dda_d = lam_dda/dpl
        n1 = complex(get_line(lines,'refractive index: 1.').split()[-1].replace('i','j'))
        n2 = complex(get_line(lines,'                  2.').split()[-1].replace('i','j'))
        print(n1,n2)
        print(dda_d,lam_dda,lam)
        CrossSecfile = open(particle_folder+'/CrossSec-Y','r')
        lines = CrossSecfile.readlines()
        Cext = float(get_line(lines,'Cext').split()[-1])
        Cabs = float(get_line(lines,'Cabs').split()[-1])
        Csca = Cext - Cabs
        print(Csca, Cabs, Cext)

        


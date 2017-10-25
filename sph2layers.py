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

data_folder = './'+str(part_size)+'mm'
for freq_str in freqs.keys():
    data_folders = glob(data_folder+'/'+freq_str+'/*')
    print(data_folders)

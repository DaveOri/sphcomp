import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd
from pytmatrix import tmatrix, tmatrix_aux, scatter, refractive, radar
#sys.path.append('/work/DBs/scattnlay')
#from scattnlay import scattnlay

#from pytmatrix.test import test_tmatrix
#test_tmatrix.run_tests()

c = 299792458. # m/s
size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c/f

sizes = np.linspace(1,10,101) # millimeters
cols = ['x','lam','mr','mi','Cs','Cb']
data = pd.DataFrame(index=sizes,columns=cols)

lam = tmatrix_aux.wl_X
for size in sizes:
	m = refractive.m_w_0C[lam]
	x = size2x(size,lam)
	scatterer = tmatrix.Scatterer(radius=size,wavelength=lam, m=m, axis_ratio=1.0)
	Cs = scatter.sca_xsect(scatterer)
	Cb = radar.radar_xsect(scatterer)
	data.loc[size] = [x,lam,m.real,m.imag,Cs,Cb]
data.index.name = 'size'
data.to_csv('Xband.csv')

lam = tmatrix_aux.wl_Ka
for size in sizes:
	m = refractive.m_w_0C[lam]
	x = size2x(size,lam)
	scatterer = tmatrix.Scatterer(radius=size,wavelength=lam, m=m, axis_ratio=1.0)
	Cs = scatter.sca_xsect(scatterer)
	Cb = radar.radar_xsect(scatterer)
	data.loc[size] = [x,lam,m.real,m.imag,Cs,Cb]
data.index.name = 'size'
data.to_csv('Kaband.csv')

lam = tmatrix_aux.wl_W
for size in sizes:
	m = refractive.m_w_0C[lam]
	x = size2x(size,lam)
	scatterer = tmatrix.Scatterer(radius=size,wavelength=lam, m=m, axis_ratio=1.0)
	Cs = scatter.sca_xsect(scatterer)
	Cb = radar.radar_xsect(scatterer)
	data.loc[size] = [x,lam,m.real,m.imag,Cs,Cb]
data.index.name = 'size'
data.to_csv('Wband.csv')



import numpy as np
import sys
import matplotlib.pyplot as plt
import pandas as pd

c = 299792458000000. # micron/s
refld = '/work/DBs/melted_aggregate_scaled_reff_Ku_Ka_W_89_165_183/melt3a_aggregate2_010_20091210_222748_halfscale_f000001_AEFF_1000_ROT535_13.4/'
icefile = refld + 'ice_266_WB08.tab'
watfile = refld + 'water_273.tab'

ice = pd.read_csv(icefile,header=2,delim_whitespace=True)
iceidx = 1e-9*c/ice.LAMBDA 
ice.index = iceidx
ice.index.name = 'freq'
ice.loc[9.6] = np.nan
ice.loc[13.6] = np.nan
ice.loc[35.6] = np.nan
ice.loc[94.0] = np.nan
ice = ice.sort_index().interpolate('cubic')

wat = pd.read_csv(watfile,header=2,delim_whitespace=True)
watidx = 1e-9*c/wat.LAMBDA 
wat.index = watidx
wat.index.name = 'freq'
wat.loc[9.6] = np.nan
wat.loc[13.6] = np.nan
wat.loc[35.6] = np.nan
wat.loc[94.0] = np.nan
wat = wat.sort_index().interpolate('cubic')

#print(ice, wat)

sys.path.append('/work/DBs/scattnlay')
sys.path.append('/home/dori/pymiecoated')
from scattnlay import scattnlay
from pymiecoated import Mie

mie = Mie(x=1.5,m=complex(3.0,0.5))
print(mie.qsca())

size2x = lambda s,l: 2.*np.pi*s/l
fr2lam = lambda f: c*1.e-6/f

plt.figure()
ax = plt.gca()
for freq in [9.6,13.6,35.6,94]:
	Nlayers = 2 # Number of layers (to be consistent across the computations)
	Ncomput = 1000 # Number of computations / (x,m) pairs
	part_size = 0.005 # meters

	mi = complex(ice.loc[freq,'Re(N)'],ice.loc[freq,'Im(N)'])
	mw = complex(wat.loc[freq,'Re(N)'],wat.loc[freq,'Im(N)'])

	x = np.ndarray((Ncomput,Nlayers),dtype=np.float64)
	m = np.ndarray((Ncomput,Nlayers),dtype=np.complex128)
	betas = np.linspace(0,2*np.pi,360) # compute S at this scattering angles, NOW [radians]

	part_x = size2x(part_size,fr2lam(freq*1e9))
	rad_ratio = np.linspace(0.1,0.99,Ncomput)
	xl = rad_ratio*part_x
	x[:,0] = xl         #[5.0]#,2.0]
	x[:,1] = part_x  #[10.0]#,2.0]
	m[:,0] = mi
	m[:,1] = mw

	print(freq,m[0,0],m[0,1])
	results = scattnlay(x,m)

	ax.plot(rad_ratio,results[2],label=freq)
#	resj = results
#	for i in range(len(xl)):
#		mie = Mie(x=xl[i],m=mi,y=part_x,m2=mw)
#		resj[1][i] = mie.qext()
#		resj[2][i] = mie.qsca()
#		resj[3][i] = mie.qabs()
#		resj[4][i] = mie.qb()
#	ax.plot(rad_ratio,resj[2],'-+',label=str(freq)+' j')

ax.legend()
ax.grid()
#ax.set_ylim([0,1])
plt.show()





#S1 = results[8][0]
#S2 = results[9][0]
#F11 = (S1*S1.conjugate()+S2*S2.conjugate()).real
#ax = plt.subplot(111, projection='polar')
#ax.set_rscale('log')
#ax.plot(betas, F11)
##ax.set_rmax(2)
##ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
#ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
#ax.grid(True)
#ax.set_title("Phase function", va='bottom')
#plt.show()
#plt.close()

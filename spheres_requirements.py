# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:30:34 2017

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')

def memory(size,jobs,d=0.02):
    #d = 0.02 # millimeters
    Nr=(np.pi*size**3.0/6.0)/d**3.0
    nx=size/d
    N=nx**3.0
    return [((288+384*jobs/nx+192/jobs)*N+463*Nr)*1e-9,nx/jobs]

jobss = [8**n for n in range(5) ]
sizes = np.arange(1,31,1)

f, ((ax1,ax2),(ax3,ax4),(ax5,ax6)) = plt.subplots(3,2,sharex=True)
for jobs in jobss:
    mems, nxperjob = memory(sizes,jobs,d=0.01)
    ax1.plot(sizes,mems,label=str(jobs))
    ax2.plot(sizes,mems/jobs,label=str(jobs))
    mems, nxperjob = memory(sizes,jobs,d=0.02)
    ax3.plot(sizes,mems,label=str(jobs))
    ax4.plot(sizes,mems/jobs,label=str(jobs))
    mems, nxperjob = memory(sizes,jobs,d=0.04)
    ax5.plot(sizes,mems,label=str(jobs))
    ax6.plot(sizes,mems/jobs,label=str(jobs))
    print(mems[-1])

ax1.legend(title='N procs')
ax1.grid()
ax1.set_title('Memory     [GB]')
ax1.set_ylabel('d = 10 $\mu$m')

ax2.set_yscale('log')
ax2.set_title('Memory per job     [GB]')
ax2.grid()

ax3.set_ylabel('d = 20 $\mu$m')
ax3.grid()

ax4.set_yscale('log')
ax4.grid()

ax5.set_ylabel('d = 40 $\mu$m')
ax5.set_xlabel('Particle dimension [mm]')
ax5.grid()

ax6.set_yscale('log')
ax6.set_xlabel('Particle dimension [mm]')
ax6.grid()

plt.tight_layout()
plt.savefig('DDA_requirements.pdf')

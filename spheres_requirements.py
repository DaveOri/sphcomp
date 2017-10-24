# -*- coding: utf-8 -*-
"""
Created on Thu Sep 28 16:30:34 2017

@author: dori
"""

import numpy as np
import matplotlib.pyplot as plt

def memory(size,jobs):
    d = 0.02 # millimeters
    Nr=(np.pi*size**3.0/6.0)/d**3.0
    nx=size/d
    N=nx**3.0
    return [((288+384*jobs/nx+192/jobs)*N+463*Nr)*1e-9,nx/jobs]

jobss = [2**n for n in range(10) ]
sizes = np.arange(1,21,1)

plt.figure()
ax1 = plt.gca()
plt.figure()
ax2 = plt.gca()
for jobs in jobss:
    mems, nxperjob = memory(sizes,jobs)
    ax1.plot(sizes,mems,label=str(jobs))
    ax2.plot(sizes,mems/jobs,label=str(jobs))
    print(mems[-1])

ax1.legend()
ax1.grid()
ax2.set_yscale('log')
ax2.legend()
ax2.grid()
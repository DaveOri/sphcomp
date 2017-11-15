# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 18:43:26 2017

@author: dori
"""

from sys import argv
from glob import glob
tol = 25.e-6
    
#filename = '/work/DBs/sphcomp/10mm/W/91/IntField-Y'
#filenamew = 'IntField-y.dat'

#script, filename, filenamew = argv
    
#f = open(filename)
#firstline = f.readline()
#tol = 25.e-6
#
#fw = open(filenamew,'w+')
#fw.write(firstline)
#
#for line in f:
#    x,y,z = line.split()[0:3]
#    x = abs(float(x))
#    y = abs(float(y))
#    z = abs(float(z))
#    if (x<tol) + (y<tol) + (z<tol):
#        fw.write(line)


files = glob('4mm/*/*/IntField-Y')
for filename in files:
    filenamew = filename[:-10] + 'IntField-y.dat'
    f = open(filename)
    firstline = f.readline()
    fw = open(filenamew,'w+')
    fw.write(firstline)
    for line in f:
        x,y,z = line.split()[0:3]
        x = abs(float(x))
        y = abs(float(y))
        z = abs(float(z))
        if (x<tol) + (y<tol) + (z<tol):
            fw.write(line)
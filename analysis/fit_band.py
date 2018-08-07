#!/usr/bin/python3
# -*- coding: UTF-8 -*-

from models import kagome
import numpy as np


## Get the distances between vacancies in order to get the correct lattice
#  parameters and reciprocal vectors



k = (0,0,0)
_,pos,latt,_ = kagome.cell(1)
a1,a2 = latt
e0 = 0.0
t1 = -1
t2 = -0.3
t3 = 0.0
H = kagome.hamil(k,a1,a2,e0,t1,t2,t3,t4=0,t5=0)





fname = '../../../Documents/data/artificial_lattices/kagome/OUTS/1orb/simple/n45_l2/nv3_na0/d0/alpha0.0/e-0.104/dfct.bands'
X,Y,_ = np.loadtxt(fname,unpack=True)

Xx = np.unique(X)
Ximp = []
Yimp = []
for x in Xx:
   for iy in Y[X==x][np.argsort(np.abs(Y[X==x]))][0:2]:
      Ximp.append( x )
      Yimp.append( iy )
Ximp=np.array(Ximp)
Yimp=np.array(Yimp)


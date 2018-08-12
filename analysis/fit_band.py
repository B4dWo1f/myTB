#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import geometry as geo
import matplotlib.pyplot as plt
from models import kagome,graphene,triangular
import numpy as np

R = '../../../Documents/data/artificial_lattices/'
fol = R+'triangular/OUTS/1orb/simple/n45_l2/nv1_na0/dNone/alpha0.0/e-0.2/'
fol = R+'graphene/OUTS/1orb/simple/n60_l2/nv2_na0/d0/alpha0.0/e-0.2/'
fol = R+'kagome/OUTS/1orb/simple/n60_l2/nv3_na0/dNone/alpha0.0/e-0.2/'


models = {'kagome':kagome, 'graphene':graphene,'triangular':triangular}

model = models[fol.split('/')[-10]]
Nbands = int(fol.split('/')[-5].split('_')[0].replace('nv',''))
dfct_pos = fol + 'base_dfct.xyz'
dfct_bands = fol + 'dfct.bands'

## Get the distances between vacancies in order to get the correct lattice
#  parameters and reciprocal vectors
try: lines = open(dfct_pos).readlines()
except FileNotFoundError: lines = open(dfct_pos.replace('base_','')).readlines()
lines[1] = lines[1].split('#')[0] # clean comments
vecs = lines[1].replace(' ','').lstrip().rstrip().split('][')
vecs = [v.replace('[','').replace(']','') for v in vecs]
latt = [np.array(list(map(float,v.split(',')))) for v in vecs]
a1,a2 = latt
_,pos,latt,_ = model.cell(1,a=np.linalg.norm(a1)/np.sqrt(3))
# position and lattice of the meta-crystal

try: X,Y,_ = np.loadtxt(dfct_bands,unpack=True)
except FileNotFoundError:
   M = np.load(dfct_bands+'.npy')
   X = M[:,0]
   Y = M[:,1]

## High-symmetry points to define the K-path
# Reciprocal vecctors
recip = geo.reciprocal(latt)
G = 0*recip[0] + 0*recip[1]
K = (2*recip[0]+recip[1])/3.
Kp = (recip[0]+2*recip[1])/3.
M = (recip[0]+recip[1])/2.

## K-path
points = [G,K,Kp,G]
rec = geo.recorrido(points,100)


def myerror(rec,E0,t1,t2,t3):
   y = []
   for k in rec:
      for e in reversed(np.linalg.eigvalsh(model.hamil(k,E0,t1,t2,t3,latt))):
         y.append(e)
   y = np.array(y)
   #error = np.sum((y-Yimp) * (y-Yimp))
   #return error
   return y


Xx = np.unique(X)
Ximp = []
Yimp = []
for x in Xx:
   for iy in Y[X==x][np.argsort(np.abs(Y[X==x]))][0:Nbands]:
      Ximp.append( x )
      Yimp.append( iy )
Ximp=np.array(Ximp)
Yimp=np.array(Yimp)
cent = np.mean(Yimp)


from scipy.optimize import curve_fit
E0,t1,t2,t3 = -0.03,-0.01,0.,0.0
E0,t1,t2,t3 = cent, -0.00330413354607, 0.000188387934477, 0.0002
a,b = curve_fit(myerror,rec,Yimp,p0=(E0,t1,t2,t3))
E0a,t1a,t2a,t3a = a
#E0a,t1a,t2a,t3a = E0,t1,t2,t3   # ignore fit


print('%s, %s, %s, %s'%(E0a,t1a,t2a,t3a))

## Calculate bands with the fitted parameters (and initial guess)
x,y,ya = [],[],[]
for i in range(len(rec)):
   k = rec[i]
   for e in np.linalg.eigvalsh(model.hamil(k,E0,t1,t2,t3,latt)):
      y.append(e)
      x.append(i)
   for e in np.linalg.eigvalsh(model.hamil(k,E0a,t1a,t2a,t3a,latt)):
      ya.append(e)




print(len(Yimp))
print(len(ya))

fig, ax = plt.subplots()
ax.scatter(Ximp,Yimp,label='data')
ax.scatter(x,ya,label='fit')
m = min([np.min(y),np.min(Yimp)])
M = max([np.max(y),np.max(Yimp)])
ax.set_ylim([m,M])
ax.set_xlim([0,300])
plt.show()


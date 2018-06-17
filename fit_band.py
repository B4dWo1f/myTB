#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import geometry as geo

import sys
try: fname = sys.argv[1]
except IndexError:
   print('File not specified')
   exit()


X,Y,_ = np.loadtxt(fname,unpack=True)

Xx = np.unique(X)
Yimp = []
for x in Xx:
   Yimp.append( Y[X==x][np.argmin(np.abs(Y[X==x]))] )
Yimp=np.array(Yimp)

r3 = np.sqrt(3)
pos = [np.array([0,0,0])]
latt = [np.array([65.1,40.0103736548,0.0]),
        np.array([67.2,-36.3730669589,0.0])]
a1,a2 = latt
#latt = [np.array([a,0,0]),
#        np.array([a/2,r3*a/2,0])]

recip = geo.reciprocal(latt)



EG = -0.0522227718136
EK = -0.0352771046745
EM = -0.0306864668879

## Check!!!
b1,b2 = recip
G = 0*b1 + 0*b2
K = (2*b1+b2)/3.
Kp = (b1+2*b2)/3.
M = (b1+b2)/2.
points = [G,K,Kp,G]

y0 = Y[np.argmin(np.abs(Y[X==0]))]
rec = geo.recorrido(points,100)
E0=-0.0375
t = (y0 - E0)/6


def Hk(k,E0,t1,t2,t3,recip):
   a1,a2 = recip
   return E0 + 2*t1*( np.cos(np.dot(k,a1)) + \
                      np.cos(np.dot(k,a2)) +\
                      np.cos(np.dot(k,a1-a2)) )  +\
               2*t2*( np.cos(np.dot(k,a1+a2)) +\
                      np.cos(np.dot(k,2*a1-a2)) +\
                      np.cos(np.dot(k,-a1+2*a2)) )+\
               2*t3*( np.cos(np.dot(k,2*a1)) + \
                      np.cos(np.dot(k,2*a2)) +\
                      np.cos(np.dot(k,2*(a1-a2))) )

def myerror(rec,E0,t1,t2,t3):
   a1,a2 = recip
   y = np.array([Hk(k,E0,t1,t2,t3,latt).real for k in rec])
   error = np.sum((y-Yimp) * (y-Yimp))
   return error

from scipy.optimize import curve_fit
E0,t1,t2,t3 = -0.036846, -0.00212014, -0.0006598, 0.000144978
a,b = curve_fit(myerror,rec,Yimp,p0=(E0,t1,t2,t3))
E0,t1,t2,t3 = a

print('%s, %s, %s, %s'%(E0,t1,t2,t3))

x = list(range(len(rec)))
y = np.array([Hk(k,E0,t1,t2,t3,latt).real for k in rec])


fig, ax = plt.subplots()
ax.scatter(X,Y)
ax.scatter(x,y)
ax.set_ylim([-0.06,-0.025])
ax.set_xlim([0,300])
ax.grid()
plt.show()

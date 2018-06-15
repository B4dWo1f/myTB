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

def Hk(k,E0,t1,t2,recip):
   a1,a2 = recip
   return E0 + t1*( np.exp(1j*np.dot(k,a1)) + np.exp(-1j*np.dot(k,a1)) + \
                  np.exp(1j*np.dot(k,a2)) + np.exp(-1j*np.dot(k,a2)) +\
                  np.exp(1j*np.dot(k,a1-a2)) + np.exp(-1j*np.dot(k,a1-a2))) +\
               t2*( np.exp(1j*np.dot(k,a1+a2)) + np.exp(-1j*np.dot(k,a1+a2))+\
                  np.exp(1j*np.dot(k,2*a1-a2)) + np.exp(-1j*np.dot(k,2*a1-a2))+\
                  np.exp(1j*np.dot(k,2*a2-a1)) + np.exp(-1j*np.dot(k,2*a2-a1)) )


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



def myerror(rec,E0,t1,t2):
   a1,a2 = recip
   def hk(k):
      aux = E0 + t1*( np.exp(1j*np.dot(k,a1)) + np.exp(-1j*np.dot(k,a1)) + \
                  np.exp(1j*np.dot(k,a2)) + np.exp(-1j*np.dot(k,a2)) +\
                  np.exp(1j*np.dot(k,a1-a2)) + np.exp(-1j*np.dot(k,a1-a2))) +\
               t2*( np.exp(1j*np.dot(k,a1+a2)) + np.exp(-1j*np.dot(k,a1+a2))+\
                  np.exp(1j*np.dot(k,2*a1-a2)) + np.exp(-1j*np.dot(k,2*a1-a2))+\
                  np.exp(1j*np.dot(k,2*a2-a1)) + np.exp(-1j*np.dot(k,2*a2-a1)) )
      return aux.real
   y = np.array([Hk(k,E0,t1,t2,latt).real for k in rec])
   error = np.sum((y-Yimp) * (y-Yimp))
   return error

from scipy.optimize import curve_fit
E0,t1,t2 = -0.0375, -0.002, 0.0
E0,t1,t2 = -0.0368451757658, -0.00206159311621, -0.000659068631154
E0,t1,t2 = -0.0368448939822, -0.00206165908683, -0.000659074828917
a,b = curve_fit(myerror,rec,Yimp,p0=(E0,t1,t2))
E0,t1,t2 = a


#params,Err = [],[]
#E0s,t1s,t2s = -0.0375, -0.00204081632653, -0.000505050505051
#dE0 = E0s*0.1
#dt1 = t1s*0.1
#dt2 = t2s*0.1
#E  = np.linspace(E0s-dE0, E0s+dE0,3)          #-0.06,-0.025,10)
#T1 = np.linspace(t1s-dt1, t1s+dt1,50)          #-0.1,0.1,50)
#T2 = np.linspace(t2s-dt2, t2s+dt2,100)          #-0.01,0.01,100)
#for E0 in [-0.0375]: #E:
#   for t1 in T1:
#      for t2 in T2:
#         y = np.array([Hk(k,E0,t1,t2,latt).real for k in rec])
#         error = np.sum((y-Yimp) * (y-Yimp))
#         Err.append( error )
#         params.append( (E0,t1,t2) )
#
#ind = np.argmin(Err)
#E0,t1,t2 = params[ind]
print('%s, %s, %s'%(E0,t1,t2))

#E0 = -0.0375
#t1 = -0.00245379530227
y = np.array([Hk(k,E0,t1,t2,latt).real for k in rec])



x = list(range(len(rec)))


fig, ax = plt.subplots()
ax.scatter(X,Y)
ax.scatter(x,y)
ax.set_ylim([-0.06,-0.025])
ax.set_xlim([0,300])
ax.grid()
plt.show()

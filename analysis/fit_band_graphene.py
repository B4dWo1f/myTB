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
Ximp = []
Yimp = []
for x in Xx:
   for iy in Y[X==x][np.argsort(np.abs(Y[X==x]))][0:2]:
      Ximp.append( x )
      Yimp.append( iy )
Ximp=np.array(Ximp)
Yimp=np.array(Yimp)
#fig, ax = plt.subplots()
#ax.scatter(X,Y)
#ax.scatter(Ximp,Yimp)


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
def Hk(k,E0,t1,t2,recip):
   a1,a2 = recip
   M = E0 + t1*(1+ np.exp(1j*np.dot(k,a1)) + np.exp(-1j*np.dot(k,a1)) +\
                   np.exp(1j*np.dot(k,a2)) + np.exp(-1j*np.dot(k,a2))   )
   H0 = np.matrix([[E0,t1],[t1,E0]])
   V1 = np.matrix([[t2,0],[t1,t2]])
   V2 = np.matrix([[t2,0],[t1,t2]])
   V1m2= np.matrix([[t2,0],[0,t2]])
   H = H0 +\
       V1*np.exp(1j*np.dot(k,a1)) + V1.H*np.exp(-1j*np.dot(k,a1)) +\
       V2*np.exp(1j*np.dot(k,a2)) + V2.H*np.exp(-1j*np.dot(k,a2)) +\
       V1m2*np.exp(1j*np.dot(k,a1-a2)) + V1m2.H*np.exp(-1j*np.dot(k,a1-a2))
   return np.linalg.eigvalsh(H)


def myerror(rec,E0,t1,t2):
   a1,a2 = recip
   y = []
   for k in rec:
      for e in Hk(k,E0,t1,t2,latt):
         y.append(e)
   y = np.array(y)
   #error = np.sum((y-Yimp) * (y-Yimp))
   #return error
   return y

## Greaphene check
#E0=0
#t1=2.7
#t2=0
#x,y = [],[]
#for i in range(len(rec)):
#   k=rec[i]
#   for e in Hk(k,E0,t1,t2,latt):
#      y.append(e)
#      x.append(i)
#import matplotlib.pyplot as plt
#fig, ax = plt.subplots()
#ax.scatter(x,y)
#plt.show()

#N=20
#print('writting')
#f = open('NNfit/data.train','w')
#for E0 in np.linspace(-0.04,-0.02,N):
#   for t1 in np.linspace(-0.02,0.02,N):
#      for t2 in np.linspace(-0.02,0.02,N):
#         y = []
#         for k in rec:
#            for e in Hk(k,E0,t1,t2,latt):
#               y.append(e)
#         for iy in y:
#            f.write(str(iy)+'  ')
#         f.write('%s  %s  %s\n'%(E0,t1,t2))
#f.close()
#print('written')



from scipy.optimize import curve_fit
E0 = -0.036
t1 = -0.0108
t2 = t1*1.1
E0,t1,t2 = -0.03946126898954034, -0.01, 0.006
E0,t1,t2 = -0.03945970556669155, 0.01, -0.02
#E0,t1,t2 = -0.040633, 0.01, -0.0017853856263088663
E0,t1,t2 = -0.037, -0.004, 0.0001
try: a,b = curve_fit(myerror,rec,Yimp,p0=(E0,t1,t2))
except RuntimeError:
   print('No Converged')
   a = E0,t1,t2
E0a,t1a,t2a = a

print('%s, %s, %s'%(E0a,t1a,t2a))


#x = list(range(len(rec)))
#y = np.array([Hk(k,E0,t1,t2,t3,latt).real for k in rec])
x = []
y = []
ya = []
for i in range(len(rec)):
   k = rec[i]
   for e in Hk(k,E0,t1,t2,latt):
      y.append(e)
      x.append(i)
   for e in Hk(k,E0a,t1a,t2a,latt):
      ya.append(e)

my_y = np.array(y)
error = np.sum( (my_y - Yimp)*(my_y - Yimp) )
print(error)



fig, ax = plt.subplots()
ax.scatter(X,Y)
ax.scatter(x,y)
#ax.scatter(x,ya)
m = min([np.min(y),np.min(Yimp)])
M = max([np.max(y),np.max(Yimp)])
ax.set_ylim([m,M])
ax.set_xlim([0,300])
plt.show()

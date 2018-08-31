#!/usr/bin/python3
# -*- coding: UTF-8 -*-


import numpy as np
import mwe_exchange as ex

fol = '../../../Documents/data_test/OUTS/1orb/ac/n45_l2/nv6_na0/d80.0/alpha0.0/e-0.2/'


#def model(x,a1,a2,a3):
def model(x,E,t,t2,t3):
   #E = -0.0365225358108
   #t1 = 0.00311256676964
   #t2 = 2.6742787152e-05
   #t3 = -0.000362200345163
   a1 = -0.0143634638305
   a2 = -0.0143634624306
   a3 = -0.00830741564187
   H = np.matrix([[E ,t ,t2,t3,t2,t ],
                  [t ,E ,t ,t2,t3,t2],
                  [t2,t ,E ,t ,t2,t3],
                  [t3,t2,t ,E ,t ,t2],
                  [t2,t3,t2,t ,E ,t ],
                  [t ,t2,t3,t2,t ,E]])
   A = np.matrix([[0 ,a1,0 ,0 ,0 ,a3],
                  [a1,0 ,a2,0 ,0 ,0 ],
                  [0 ,a2,0 ,a3,0 ,0 ],
                  [0 ,0 ,a3,0 ,a1,0 ],
                  [0 ,0 ,0 ,a1,0 ,a2],
                  [a3,0 ,0 ,0 ,a2,0 ]])
   return np.linalg.eigvalsh(H)

S = ex.Spectrum(fol)
Es = S.E_ingap

from scipy.optimize import curve_fit
E = np.mean(Es)
t = abs(np.max(Es) - np.min(Es))
t2 = t/10
t3 = t/15
print('Initial:',E,t,t2,t3)
a,b = curve_fit(model,[],Es,p0=(E,t,t2,t3))
print('Fit:',*a)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(Es,'o-')
ax.plot(model([],*a),'o-')
plt.show()

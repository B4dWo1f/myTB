#!/usr/bin/python3
# -*- coding: UTF-8 -*-


import numpy as np

def blue(J,D,tRL,tLR,UR,UL,e1,e2):
   """
   order: {u,u},{d,d},{u,d},{d,u},{ud,},{,ud}
   """
   M = np.matrix([[-J/4,  0 ,  0 ,  0 ,   0  ,   0  ],
                  [  0 ,-J/4,  0 ,  0 ,   0  ,   0  ],
                  [  0 ,  0 , J/4,-J/2,  tRL ,  tLR ],
                  [  0 ,  0 ,-J/2, J/4, -tRL , -tLR ],
                  [  0 ,  0 , tRL,-tRL,UR-J/4,   D  ],
                  [  0 ,  0 , tLR,-tLR,  D   ,UL-J/4]])

   K = np.matrix([[e1+e2,  0  ,  0  ,  0  ,  0  ,  0  ],
                  [  0  ,e1+e2,  0  ,  0  ,  0  ,  0  ],
                  [  0  ,  0  ,e1+e2,  0  ,  0  ,  0  ],
                  [  0  ,  0  ,  0  ,e1+e2,  0  ,  0  ],
                  [  0  ,  0  ,  0  ,  0  ,e1+e1,  0  ],
                  [  0  ,  0  ,  0  ,  0  ,  0  ,e2+e2]])
   return M+K


import sys
try:
   fname = sys.argv[1:]
   d,e = map(float,fname)
except IndexError:
   print('Distance or electric file not specified')
   exit()



Ds,E,JF,D,TRL,TLR,UR,UL,E1,E2 = np.loadtxt('datos.dat',unpack=True)


dist_d = np.argmin(np.abs(Ds-d))
d = Ds[dist_d]
E = E[Ds==d]
JF = JF[Ds==d]
D = D[Ds==d]
TRL = TRL[Ds==d]
TLR = TLR[Ds==d]
UR = UR[Ds==d]
UL = UL[Ds==d]
E1 = E1[Ds==d]
E2 = E2[Ds==d]
Ds = Ds[Ds==d]

dist_e = np.argmin(np.abs(E-e))
e = E[dist_e]
Ds = Ds[E==e]
JF = JF[E==e]
D = D[E==e]
TRL = TRL[E==e]
TLR = TLR[E==e]
UR = UR[E==e]
UL = UL[E==e]
E1 = E1[E==e]
E2 = E2[E==e]
E = E[E==e]

print(d,e)

import matplotlib.pyplot as plt
fig, ax = plt.subplots()
for i in range(len(E)):
   Hb = blue(JF[i],D[i],TRL[i],TLR[i],UR[i],UL[i],E1[i],E2[i])
   es = np.linalg.eigvalsh(Hb)
   ax.plot(es,'o',label='E=%s'%(E[i]))
ax.legend(loc=4)
plt.show()

#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 15




## Basis
class basis(object):
   """
     Minimum dummy class to store indices of atom and orbital
   """
   def __init__(self,fname='dfct.basis'):
      self.ind,self.n = np.loadtxt(fname,usecols=(0,1),unpack=True,dtype=int)

def plot_spectrum(E,V,pos,Ba=[]): #inds=[]):   #,ax=None):
   def state_space(v,pos,base,S=100):
      """
      Plot real space distribution of a given wave-function
        v: vector to plot, with shape: (N,) 
        X,Y,Z: x,y,z coordinates of each of the atoms
      """
      X = pos[:,0]
      Y = pos[:,1]
      Z = np.where(pos[:,2]>0,1,-1)
      V = v * np.conj(v)
      V = np.bincount(Ba.n, weights=V)
      s = 1000/np.sqrt(len(X))
      fig, ax = plt.subplots()
      ax.scatter(X,Y,c='k',s=s,edgecolors='none',alpha=0.3,zorder=0)
      ax.scatter(X,Y,c=Z,s=900*s*V,edgecolors='none',cmap='rainbow')
      ax.set_xlim([min(X)-1,max(X)+1])
      ax.set_ylim([min(Y)-1,max(Y)+1])
      ax.set_aspect('equal')
      fig.tight_layout()
      plt.show()

   def picker_wrapper(maxd=0.075):
      def picker(line, Mevent):
         """
         find the points within a certain distance from the mouseclick in
         data coords and attach some extra attributes, pickx and picky
         which are the data points that were picked
         """
         if Mevent.xdata is None:
            return False, dict()
         x0 = line.get_xdata()
         y0 = line.get_ydata()
         d = np.sqrt((x0 - Mevent.xdata)**2. + (y0 - Mevent.ydata)**2.)
         ind = np.nonzero(np.less_equal(d, maxd))
         if len(ind):
            pickx = np.take(x0, ind)
            picky = np.take(y0, ind)
            props = dict(ind=ind, pickx=pickx, picky=picky)
            return True, props
         else:
            return False, dict()
      return picker

   def onpick_wrapper(energies,vecs,pos,base,plot):
      def onpick(event):
         try:
            ind = event.pickx[0,0]
            v = vecs[ind,:]
            print('Eigenvalue:',ind,'  Energy: %seV'%(energies[ind]))
            plot(v,pos,base)
            with open('/tmp/ingap.txt','a') as f:
               f.write(str(ind)+'\n')
         except IndexError: pass   #print('no data selected')
      return onpick

   my_onpick = onpick_wrapper(E,V,pos,Ba,state_space)
   my_picker = picker_wrapper()

   fig, ax = plt.subplots()
   X = range(len(E))
   line, = ax.plot(X,E, 'o', picker=my_picker)
   ax.set_xlim([-1,len(X)+1])
   ax.grid()

   fig.canvas.mpl_connect('pick_event', my_onpick)
   plt.show()


def get_exchanges(vL,vR,U=2.7):
   """
     Compute all the effective exchange couplings
   """
   JF, UL, UR, D, tLR, tRL = 0., 0., 0., 0., 0., 0.
   for iv in range(len(vL)):
      JF  += U * ( np.conj(vL[iv])*vL[iv] ) * ( np.conj(vR[iv])*vR[iv] )
      # Renormalized Hubbard ---> U*IPR
      UL  += U * ( np.conj(vL[iv])*vL[iv] ) * ( np.conj(vL[iv])*vL[iv] )
      UR  += U * ( np.conj(vR[iv])*vR[iv] ) * ( np.conj(vR[iv])*vR[iv] )
      D   += U * ( np.conj(vL[iv]) * vR[iv] )**2
      tLR += U *  np.conj(vL[iv])*vL[iv] * np.conj(vL[iv])*vR[iv]
      tRL += U *  np.conj(vR[iv])*vR[iv] * np.conj(vR[iv])*vL[iv]
   return UL, UR, 2*JF, D, tLR, tRL


class Spectrum(object):
   def __init__(self,fname):
      """
        fname is the folder containing the output of the calculations.
        It must containthe following files:
         - pris_spectrum.npy           - dfct_spectrum.npy
         - base_pris.xyz               - base_dfct.xyz
         - pris.basis (optional)       - dfct.basis
      """
      # Useful files
      self.f = fname
      pris = fname + 'pris_spectrum.npy'
      dfct = fname + 'dfct_spectrum.npy'
      pris_pos = fname + 'base_pris.xyz'
      dfct_pos = fname + 'base_dfct.xyz'
      dfct_basis = fname + 'dfct.basis'
      log = fname + 'main.log'

      try:  # Basis for multiorbital systems
         self.Ba = basis(dfct_basis)
      except:
         print('No basis')
         self.Ba = []
      #Spectrum
      M = np.load(pris)
      Ep = M[:,0].real
      M = np.load(dfct)
      self.E = M[:,0].real
      self.V = M[:,1:]
      self.pos = np.loadtxt(dfct_pos,skiprows=2,usecols=(1,2,3))
   def plot(self):
      plot_spectrum(self.E, self.V, self.pos, self.Ba)
   def select_ingap(self):
      # Workaround to select states
      tmp = '/tmp/ingap.txt'
      f = open(tmp,'w')
      f.close()
      plot_spectrum(self.E, self.V, self.pos, self.Ba)
      self.inds = np.loadtxt(tmp,dtype=int)
      f = open(tmp,'w')
      f.close()
   def analyze_ingap(self):
      try: self.inds
      except AttributeError: self.select_ingap()
      inds = self.inds
      print('Analyzing states',self.inds)
      self.E_ingap = self.E[inds]
      self.V_ingap = self.V[inds]
   def get_blue_parameters(self):
      UL, UR, JF, D, tLR, tRL = get_exchanges(*self.V_ingap)
      e1,e2 = self.E_ingap
      print(JF,D,tRL,tLR,UR,UL,e1,e2)
      H = blue(JF,D,tRL,tLR,UR,UL,e1,e2)
      es,v = np.linalg.eigh(H)
      v = v.transpose()
      return es,v


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

if __name__ == '__main__':
   import sys
   fol = sys.argv[1]
   print(fol)
   A = Spectrum(fol)
   #A.select_ingap()
   A.analyze_ingap()
   es,v = A.get_blue_parameters()
   for i in range(len(es)):
      print(i,es[i],np.round(v[i,:],4))
   import matplotlib.pyplot as plt
   plt.close('all')
   fig, ax = plt.subplots()
   ax.plot(es,'o')
   ax.grid()
   plt.show()

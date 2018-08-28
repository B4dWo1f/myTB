#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl



def plot_spectrum(E,V,pos,Ba=[],Ep=[],inds=[],title=''):   #,ax=None):
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

   def sync_spectrums(Ep,Ed):
      """
        Unfinished attempt to align in the X axis the pristine and defected
        spectrum
      """
      Pcond = np.max(Ep[Ep<0])
      Pvale = np.min(Ep[Ep>0])

      Dcond = Ed[np.argmin(np.abs(Ed-Pcond))]
      Dvale = Ed[np.argmin(np.abs(Ed-Pvale))]
      InGap = Ed[Ed>Dcond]
      InGap = InGap[InGap<Dvale]
      Ng = len(InGap)

      Xd = np.array(range(len(Ed)))
      indp, = np.where(Ep == Pcond)
      indd, = np.where(Ed == Dcond)
      if indp == indd:
         print('extra state in E>0')
         Xp = []
         for i in range(len(Ed)):
            if Ep[i] <= Pcond: Xp.append(i)
            else: Xp.append(i+Ng)
      else:
         print('extra state in E<0')
         Xp = []
         for i in range(len(Ed)):
            if Ep[i] <= Pcond: Xp.append(i-Ng+1)
            else: Xp.append(i+Ng-1)
      return np.array(Xp), Xd

   Xp,X = sync_spectrums(Ep,E)
   inds = np.array(inds)
   fig, ax = plt.subplots()
   line, = ax.plot(X,E, 'o', picker=my_picker,label='dfct')
   if len(Ep) != 0: ax.plot(Xp+0.1,Ep, 'o',label='pris')
   if len(inds) != 0: ax.scatter(X[inds],E[inds],c='r',s=100,zorder=0)
   ax.set_xlim([-1,len(X)+1])
   ax.legend(loc=2)
   if len(title) != 0: ax.set_title(title)

   fig.canvas.mpl_connect('pick_event', my_onpick)
   plt.show()


def get_geo_stuff(fol):
   """
     Scrap the log file for geometry information
     TODO: Move the geometri info to an independent output file
   """
   log = fol + 'main.log'
   ## Get Vacancies indices
   vacs = os.popen('grep " Changing onsite of atom:" %s'%(log)).readlines()
   vacs = [v.split(':')[-1] for v in vacs]
   vacs = list(set([int(v.lstrip().rstrip()) for v in vacs]))
   Nv = len(vacs)
   ## Get distance
   try:
      com = 'grep " Requested-dist/Real-dist: " %s'%(log)
      dis = os.popen(com).read().split()[-1].split('/')[-1]
   except IndexError:
      #print('Distance not found')
      dis = 0.0
   dis = float(dis)
   ## Get angle
   try:
      com = 'grep " Requested-angle/Real-angle: " %s'%(log)
      ang = os.popen(com).read().split()[-1].split('/')[-1]
   except IndexError:
      #print('Angle not found')
      ang = 0.0
   ang = float(ang)
   return Nv, dis, ang, vacs

class InGap(Exception):
   def __init__(self, value):
      self.value = value
   def __str__(self):
      return repr(self.value)

def get_ingap(E,V,Ba,pos,vacs,cond,vale):
   # Check in-gap
   v = E[E>vale]
   v= v[v<cond]
   if v.shape[0] == len(vacs):
      inds = np.array([np.where(E==iv)[0][0] for iv in v])
      return inds,False
   else:
      inds = []
      for r0 in pos[vacs]:
         ls = []
         for i in range(V.shape[0]):
            v = V[i,:]
            l = get_lc(Ba,pos,v,r0)
            ls.append(l)
         inds.append( np.argmin(ls) )
      return inds, True

def get_cond_vale(E,Ep):
   """
     In the pristine case, vale is the last occupied (E<0) state and
     cond is the first empty state (E>0).
     In the defective case cond and vale are the states closest to the
     pristine case
   """
   condP = np.min(Ep[Ep>0])
   valeP = np.max(Ep[Ep<0])
   cond = E[np.argmin(np.abs(E-condP))]
   vale = E[np.argmin(np.abs(E-valeP))]
   cond = E[np.argmin(np.abs(E-condP))]
   return condP,valeP, cond, vale

def get_lc(Ba,pos,v,r0,lim=0.9):
   """ Calculate the confinement length """
   v = np.conj(v)*v
   v = v.real   # by  construction v is real
   vp = np.bincount(Ba.n, weights=v)
   dist_aux = np.array([np.linalg.norm(r-r0) for r in pos])
   inds_dist = np.argsort(dist_aux)
   v_sorted = vp[inds_dist]  #np.conj(v[inds_dist]) * v[inds_dist]
   dist_aux = dist_aux[inds_dist]
   aux = 0.0
   for i in range(len(v_sorted)):
      aux += v_sorted[i]
      if aux >= lim: return dist_aux[i]
   return dist_aux[i]


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


def mean_val(v,M):
   """ Returns the expected value: <v|M|v> """
   return np.sum(np.conj(v)*M*v)


## Basis
class basis(object):
   """
     Minimum dummy class to store indices of atom and orbital
   """
   def __init__(self,fname='dfct.basis'):
      self.ind,self.n = np.loadtxt(fname,usecols=(0,1),unpack=True,dtype=int)


class Spectrum(object):
   def __init__(self,fname,nv=None,slct=False):
      """
        fname is the folder containing the output of the calculations.
        It must containthe following files:
         - pris_spectrum.npy           - dfct_spectrum.npy
         - base_pris.xyz               - base_dfct.xyz
         - pris.basis (optional)       - dfct.basis
        slct: In case of non-clear in-gap states, display the spectrum to
              manually select the in-gap states
      """
      if fname[-1] != '/': fname += '/'
      # Useful files
      self.f = fname
      pris = fname + 'pris.spectrum.npy'
      dfct = fname + 'dfct.spectrum.npy'
      pris_pos = fname + 'pris.xyz'
      dfct_pos = fname + 'dfct.xyz'
      dfct_basis = fname + 'dfct.basis'
      log = fname + 'main.log'

      try:  # Basis for multiorbital systems
         self.Ba = basis(dfct_basis)
      except:
         print('No basis')
         self.Ba = []

      self.elec = float(fname.split('/')[-2][1:]) # TODO Fix this!!!

      # Spectrum
      M = np.load(pris)
      self.Ep = M[:,0].real
      M = np.load(dfct)
      self.E = M[:,0].real
      self.V = M[:,1:]
      self.condP,self.valeP, self.cond,self.vale = get_cond_vale(self.E,self.Ep)
      self.gapP = self.condP - self.valeP
      self.gap = self.cond - self.vale

      self.pos = np.loadtxt(dfct_pos,skiprows=2,usecols=(1,2,3,4))
      self.sub = self.pos[:,3]
      self.pos = self.pos[:,:-1]
      # Requires hexagon with 1 side parallel to X or Y
      X = self.pos[:,0]
      Y = self.pos[:,1]
      Z = self.pos[:,2]
      lx = np.max(X) - np.min(X)   # diameter
      ly = np.max(Y) - np.min(Y)
      self.a = min((lx/2.,ly/2.))
      self.diameter = max((lx,ly))
      self.ap = self.a * np.sqrt(3)/2.

      # Geo stuff
      if nv == None:
         self.Nv, self.dist, self.alpha, self.vacs = get_geo_stuff(fname)
      else:
         _, self.dist, self.alpha, self.vacs = get_geo_stuff(fname)
         self.Nv = nv
      if self.Nv == 2:
         r1 = self.pos[self.vacs][0]
         r2 = self.pos[self.vacs][1]
         r = r1-r2
         self.alpha = np.degrees(np.arctan(r[1]/r[0]))
      self.analyze_ingap(slct)
   def plot(self):
      plot_spectrum(self.E, self.V, self.pos, self.Ba, self.Ep,self.inds)
   def select_ingap(self):
      # Workaround to select states
      tmp = '/tmp/ingap.txt'
      f = open(tmp,'w')
      f.close()
      plot_spectrum(self.E, self.V, self.pos, self.Ba, self.Ep,title='Select in-gap')
      self.inds = np.unique(np.sort(np.loadtxt(tmp,dtype=int)))
      f = open(tmp,'w')
      f.close()
   def analyze_ingap(self,select=False):
      try: self.inds
      except AttributeError:
         inds,warn = get_ingap(self.E,self.V,self.Ba,self.pos,self.vacs,
                                                          self.cond,self.vale)
         #try: inds = get_ingap(self.E,self.V,self.Ba,self.pos,self.vacs,
         #                                                 self.cond,self.vale)
         #except InGap:
         #   #print('raised and addressed')
         #   #inds =  self.select_ingap()
         #   inds = get_ingap(self.E,self.V,self.Ba,self.pos,self.vacs,
         #                                                 self.cond,self.vale)
         if warn:
            self.warning = True
            if select: self.select_ingap()
      self.inds = inds
      #inds = self.inds
      #print('Analyzing states',self.inds)
      self.E_ingap = self.E[inds]
      self.V_ingap = self.V[inds]
      # Properties of the in-gap states
      SUB = self.sub[self.Ba.n]
      self.SP = [mean_val(v,SUB) for v in self.V_ingap]
      LAY = np.where(self.pos[:,2]>0,1,-1)[self.Ba.n]
      self.SL = [mean_val(v,LAY) for v in self.V_ingap]
      norm = np.linalg.norm
      self.lc = []
      for v in self.V_ingap:
         r0 = self.pos[self.vacs[0]]
         dists = [norm(self.pos[i]-r0)**2 for i in range(len(self.pos))]
         dists = np.array(dists)
         self.lc.append( mean_val(v,dists)/self.diameter )
      self.lc90 = [ get_lc(self.Ba,self.pos,v,r0) for v in self.V_ingap]
   def get_blue_parameters(self):
      UL, UR, JF, D, tLR, tRL = get_exchanges(*self.V_ingap)
      e1,e2 = self.E_ingap
      return JF,D,tRL,tLR,UR,UL,e1,e2
   def __str__(self):
      msg = 'Analysis of the system: %s\n'%(self.f)
      msg += 'Nc: %s'%(self.pos.shape[0])
      msg += '   Nv: %s'%(self.Nv)
      msg += '   elec: %s\n'%(self.elec)
      msg += 'Gap:\n'
      msg += '  -Pris: %s\n'%(self.gapP)
      msg += '  -Dfct: %s\n'%(self.gap)
      if self.Nv == 2: msg += 'dist: %s   alpha: %s\n'%(self.dist,self.alpha)
      msg += 'Geometry:\n'
      msg += 'Hexagon  ->  a: %s  ;  a\': %s\n'%(self.a, self.ap)
      msg += 'Confinement Length:'
      if self.Nv == 1:
         msg +='\n   lc: %.4f\n'%(self.lc[0])
      elif self.Nv == 2:
         msg +='\n   Left: %.4f   Right: %4f\n'%(self.lc[0],self.lc[1])
         aux = np.abs(self.E_ingap[1]-self.E_ingap[0])
         msg += '\nIn-gap splitting: %s (eV)\n'%(aux)
         fac = 1000   #XXX meV!!!
         msg += 'Effective Exchange couplings (meV)\n'
         msg += '  UL: %s    ;    UR: %s\n'%(self.UL*fac, self.UR*fac)
         msg += '  J_f: %s    ;    D: %s\n'%(self.J_f*fac,self.D*fac)
         msg += '  tLR: %s    ;    tRL: %s\n'%(self.tLR*fac, self.tRL*fac)
         msg += '  J_af: %s\n'%(self.J_af*fac)
      msg += '\n ------ Pristine ------           ------ Defected ------\n'
      ic = 0 # In-gap counter
      for i in range(len(self.Ep)):
         ep = self.Ep[i]
         ed = self.E[i]
         space = '                '
         if ep in [self.valeP, self.condP]: msg += '    * %+.7f'%(ep)
         else: msg += '      %+.7f'%(ep)
         msg += space
         if ed in [self.vale, self.cond]: msg += '  * %+.7f\n'%(ed)
         elif ed in self.E_ingap:
               msg += ' -> %+.7f'%(ed)
               msg +='   SP:%.4f  SL: %.4f\n'%(self.SP[ic],self.SL[ic])
               #msg +='   SP:%.4f  SL: %.4f\n'%(self.SP_eig[ic],self.SL_eig[ic])
               ic += 1
         else: msg += '    %+.7f\n'%(ed)
      return msg


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
   A = Spectrum(fol,slct=True)
   A.plot()
   print(A)
   exit()
   A.select_ingap()
   #A.analyze_ingap()
   JF,D,tRL,tLR,UR,UL,e1,e2 = A.get_blue_parameters()
   H = blue(JF,D,tRL,tLR,UR,UL,e1,e2)
   es,v = np.linalg.eigh(H)
   v = v.transpose()

   for i in range(len(es)):
      print(i,es[i],np.round(v[i,:],4))
   import matplotlib.pyplot as plt
   plt.close('all')
   fig, ax = plt.subplots()
   ax.plot(es,'o')
   ax.grid()
   plt.show()

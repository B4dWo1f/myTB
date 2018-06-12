#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['font.size'] = 15


def find_nearest(array,value):
   """
     Returns the index of the element of an array closest to a certain value
   """
   return (np.abs(array-value)).argmin()

## Decide sign
def get_sign(v,e):
   """
     Aux workaround to find the closest element with sign
   """
   ind_p = find_nearest(v, e)
   ind_m = find_nearest(v,-e)
   if abs(v[ind_p]-e) < abs(v[ind_m]+e): return e
   else: return -e


def mean_val(v,M):
   """ Returns the expected value: <v|M|v> """
   return np.sum(np.conj(v)*M*v)

## Get ExchangeS ########################################
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


def get_ingap1(E,V,Ba,pos,vacs,cond,vale):
   #Check if in-gap
   v = E[E>vale]
   v= v[v<cond]
   if v.shape[0] == len(vacs):
      inds = np.where(E==v)[0]
      return inds
   else:
      print('WARNING, In-gap out of the gap!!!')
      inds = []
      for r0 in pos[vacs]:
         ls = []
         for i in range(V.shape[0]):
            v = V[i,:]
            l = get_lc(Ba,pos,v,r0)
            ls.append(l)
         inds.append( np.argmin(ls) )
      return inds

#def get_ingap(E,Nv=2,vb=False):
#   """
#     Returns the indices of the in-gap states as well as the conduction and
#     valence energies
#   """
#   absE = np.abs(E)
#   ind_E = np.argsort(absE)
#   ingap = [get_sign(E,e) for e in absE[ind_E[0:Nv]]]
#   EE = E
#   for ein in ingap:
#      EE = EE[EE!=ein]
#   cond = np.max(EE[EE<0])
#   vale = np.min(EE[EE>0])
#   ## indices
#   inds = [find_nearest(E,ie) for ie in ingap]
#   return inds,cond,vale

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

def get_exchanges(vL,vR):
   """
     Compute all the effective exchange couplings
   """
   JF, UL, UR, D, tLR, tRL = 0., 0., 0., 0., 0., 0.
   for iv in range(len(vL)):
      JF += ( np.conj(vL[iv])*vL[iv] ) * ( np.conj(vR[iv])*vR[iv] )
      UL += ( np.conj(vL[iv])*vL[iv] ) * ( np.conj(vL[iv])*vL[iv] )
      UR += ( np.conj(vR[iv])*vR[iv] ) * ( np.conj(vR[iv])*vR[iv] )
      D  += ( np.conj(vL[iv]) * vR[iv] )**2
      tLR += np.conj(vL[iv])*vL[iv] * np.conj(vL[iv])*vR[iv]
      tRL += np.conj(vR[iv])*vR[iv] * np.conj(vR[iv])*vL[iv]
   return UL, UR, 2*JF, D, tLR, tRL

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


## Basis
class basis(object):
   """
     Minimum dummy class to store indices of atom and orbital
   """
   def __init__(self,fname='dfct.basis'):
      self.ind,self.n = np.loadtxt(fname,usecols=(0,1),unpack=True,dtype=int)


class Spectrum(object):
   def __init__(self,fname,nv=None):
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

      # Basis for multiorbital systems
      self.Ba = basis(dfct_basis)

      # Electric Field
      self.elec = float(fname.split('/')[-2][1:]) # TODO Fix this!!!

      #Spectrum
      M = np.load(pris)
      self.Ep = M[:,0].real
      M = np.load(dfct)
      self.E = M[:,0].real
      self.V = M[:,1:]
      self.condP,self.valeP, self.cond,self.vale = get_cond_vale(self.E,self.Ep)
      self.gapP = self.condP - self.valeP
      self.gap = self.cond - self.vale

      if nv == None:
         self.Nv, self.dist, self.alpha, self.vacs = get_geo_stuff(fname)
      else:
         _, self.dist, self.alpha, self.vacs = get_geo_stuff(fname)
         self.Nv = nv
           ### Get conduction, valence, and in-gap states and indices
           #M = np.load(pris)
           #self.Ep = M[:,0].real
           #_,self.Pcond,self.Pvale = get_ingap(self.Ep, Nv=0)
           #self.gap_pris = self.Pvale-self.Pcond
           #M = np.load(dfct)
           #self.E = M[:,0].real
           #self.V = M[:,1:]
      # Position of atoms    #XXX why not use IO.xyz??
      self.pos = np.loadtxt(dfct_pos,skiprows=2,usecols=(1,2,3))
      self.Nc = self.pos.shape[0]
      self.sub = np.loadtxt(dfct_pos,skiprows=2,usecols=(4,))
      ##position
      self.inds = get_ingap1(self.E,self.V,self.Ba,self.pos,self.vacs,
                                                          self.cond,self.vale)
      ##inds,cond,vale = get_ingap(self.E, Nv=self.Nv)
      #self.inds = inds
      #self.cond = cond
      #self.vale = vale
      #self.gap_dfct = vale-cond
      self.E_ingap = self.E[self.inds]
      self.V_ingap = self.V[self.inds]
      self.properties()
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
               msg +='   SP:%.4f  SL: %.4f\n'%(self.SP_eig[ic],self.SL_eig[ic])
               ic += 1
         else: msg += '    %+.7f\n'%(ed)
      return msg
   def properties(self):
      Rx = self.pos[:,0][self.Ba.n]   # Operator position
      SUB = self.sub[self.Ba.n]   # Operator position
      ## Geometry and localization
      if self.V_ingap.shape[0] == 1:
         r0 = self.pos[self.vacs[0]]
         #self.lc = [get_lc(self.V_ingap[0], self.pos, r0)]
         self.lc = [get_lc(self.Ba,self.pos,self.V_ingap[0],r0,lim=0.9)]
         self.r0 = [r0]
      elif self.V_ingap.shape[0] == 2:
         ## In-gap splitting
         self.splitting = np.abs(self.E_ingap[0]-self.E_ingap[1])
         ## Define and sort left and right states
         v1 = (self.V_ingap[0] + self.V_ingap[1])/np.sqrt(2)
         v2 = (self.V_ingap[0] - self.V_ingap[1])/np.sqrt(2)
         A = mean_val(v1,Rx)
         B = mean_val(v2,Rx)
         ## Order of LR = vL, vR
         if A > B: self.LR = [v2,v1]      # Left, Right
         else: self.LR = [v1,v2]          # Left, Right
         try:  #XXX check if this is harmless
            ## position of left and right vacancies
            A,B = self.pos[self.vacs][:,0]
            if A > B:
               rL = self.pos[self.vacs[1]]
               rR = self.pos[self.vacs[0]]
            else:
               rL = self.pos[self.vacs[0]]
               rR = self.pos[self.vacs[1]]
            lcL = get_lc(self.Ba,self.pos,self.LR[0],rL)
            lcR = get_lc(self.Ba,self.pos,self.LR[1],rR)
            self.lc = [lcL,lcR]
            self.r0 = [rL,rR]
         except: pass
         UL, UR, JF, D, tLR, tRL = get_exchanges(self.LR[0], self.LR[1])
         self.J_f = JF
         self.UL = UL
         self.UR = UR
         self.D = D
         self.tLR = tLR
         self.tRL = tRL
         self.J_af = (self.E_ingap[0]-self.E_ingap[1])**2
      # Hexagon
      C = np.mean(self.pos,axis=0)
      pos = self.pos - C
      X = pos[:,0]
      Y = pos[:,1]
      Z = pos[:,2]
      ## Requires hexagon with 1 side parallel to X or Y
      lx = np.max(X) - np.min(X)   # diameter
      ly = np.max(Y) - np.min(Y)
      self.a = min((lx/2.,ly/2.))
      self.ap = self.a * np.sqrt(3)/2.
      ## Sublattice Polarization
      self.SP_eig = [mean_val(v,SUB) for v in self.V_ingap]
      if self.Nv == 2:
         self.SL_LR = [mean_val(v,SUB) for v in [self.LR[0],self.LR[1]]]
      ## Layer Polarization
      Z = self.pos[:,2]
      LA = np.where(Z>0,1,-1)  # Layer operator
      aux = []
      for i in range(len(LA)):
         for _ in self.Ba.ind[self.Ba.n==i]:
            aux.append(LA[i])
      LA = np.array(aux)
      self.SL_eig = [mean_val(v,LA) for v in self.V_ingap]
      if self.Nv == 2:
         self.SL_LR = [mean_val(v,LA) for v in [self.LR[0],self.LR[1]]]
   def plot_spectrum(self,inds=[]):   #,ax=None):
      def state_space(v,pos,S=100):
         """
         Plot real space distribution of a given wave-function
           v: vector to plot, with shape: (N,) 
           X,Y,Z: x,y,z coordinates of each of the atoms
         """
         X = pos[:,0]
         Y = pos[:,1]
         Z = np.where(pos[:,2]>0,1,-1)
         V = v * np.conj(v)
         V = np.bincount(self.Ba.n, weights=V)
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

      def onpick_wrapper(energies,vecs,base,plot):
         def onpick(event):
            try:
               ind = event.pickx[0,0]
               v = vecs[ind,:]
               print('Eigenvalue:',ind,'  Energy: %seV'%(energies[ind]))
               plot(v,base)
            except IndexError: pass   #print('no data selected')
         return onpick

      my_onpick = onpick_wrapper(self.E,self.V,self.pos,state_space)
      my_picker = picker_wrapper()

      fig, ax = plt.subplots()
      Xplt = range(len(self.E))
     # ax.scatter(Xplt,self.Ep,c='k', s=150,edgecolors='none',alpha=0.5,
     #                                                         label='Pristine')
      line, = ax.plot(range(len(self.E)),self.E, 'o', picker=my_picker,
                                                              label='Defected')
      ax.scatter(self.inds, self.E_ingap,c='r',s=100,edgecolors='none')
      ax.set_xlim([-1,len(self.E)+1])
      ax.grid()
      #ax.legend(loc=2)

      fig.canvas.mpl_connect('pick_event', my_onpick)
      fig.tight_layout()
      plt.show()

   def plot_state(self,inds=[]):   #,ax=None):
      if len(inds) == 0: inds = self.inds
      pos = self.pos
      vacs = self.vacs
      for ind in inds:
         v = self.V[ind]
         cc = np.conj(v) * v   # contains the square of the elements already
         cc = cc.real  # by construction cc is real
         X = pos[:,0]
         Y = pos[:,1]
         Z = pos[:,2]
         C = np.bincount(self.Ba.n, weights=cc)
         # Normalize the scalar
         C = C/np.max(C)
         zz = np.unique(Z)
         ## Plot Left and Right states
         fig, ax = plt.subplots()
         mks = ['s','o']
         cls = ['b','r']
         zord = 9
         for i in range(len(zz)):
            z = zz[i]
            try:
               ax.scatter(X[Z==z],Y[Z==z],c='k',s=10,marker=mks[i],
                                                        alpha=0.1,zorder=1)
               ax.scatter(X[Z==z],Y[Z==z],c=cls[i],s=500*C[Z==z],marker=mks[i],
                                     edgecolor='none',alpha=0.5,zorder=zord)
               zord += 1
            except: pass
         # Vacancies
         ax.scatter(X[vacs],Y[vacs],c='b',s=200,marker='>',zorder=21)
         ax.plot(X[vacs],Y[vacs],'k--',lw=2,zorder=22)
         # Confinement
         try:
            for lc,r0 in zip(self.lc, self.r0):
               circ = plt.Circle((r0[0],r0[1]), lc, color='g',lw=2,ls='--',\
                                  fill=False,zorder=23)
               ax.add_artist(circ)
         except AttributeError: pass
         # Figure settings
         ax.set_xlim([min(X),max(X)])
         ax.set_ylim([min(Y),max(Y)])
         ax.set_xticks([])
         ax.set_yticks([])
         if self.Nv == 1: ax.set_title('$\lambda_E = %s$'%(self.elec))
         else: ax.set_title('$E = %s$'%(self.E[ind]),fontsize=30)
         ax.set_aspect('equal')
         fig.tight_layout()
      plt.show()
   def ret_info(self):
      return self.Nc, self.Nv, self.elec, self.dist, self.alpha, self.gap_pris, self.gap_dfct, self.splitting, self.J_f, self.J_af, self.UL, self.UR, self.D
   def ret_LR_states(self):
      return self.LR,self.pos




if __name__ == '__main__':
   import sys
   fol = sys.argv[1]
   A = Spectrum(fol)
   print(A)
   v = np.conj(A.V_ingap[0]) * A.V_ingap[0]
   v = v.real   # a is real by construction
   print('hyper:',A.elec,'',v[-1]*1420)    # elec   hyperfine (MHz)
   A.plot_spectrum()

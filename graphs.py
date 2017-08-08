#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import gridspec
from itertools import product 
import numpy as np
#import util as ut
import logging
LG = logging.getLogger(__name__)


cdict={'red'  : ((0.,   0,   0),(0.6,0.0,0.0),(1, 1.0, 1.0)), 'green': ((0., 0.0, 0.0),(0.4,1.0,1.0),(0.6,1.0,1.0),(1, 0.0, 0.0)), 'blue' : ((0., 1.0, 1.0),(0.4,0.0,0.0),(1, 0.0, 0.0))}
my_cmap = LinearSegmentedColormap('my_colormap',cdict,256)


def crystal(pos,latt=[],fname=None,N=3,show=False):
   """
     Plots the Bravais lattice and the First Brillouin zone
   """
   if len(latt) == 0: LG.warning('No lattice vectors provided. FBZ will be a point')
   recip = geo.reciprocal(latt)
   fig = plt.figure(figsize=(20,10))
   gs = gridspec.GridSpec(1, 2)   
   fig.subplots_adjust(wspace=0.25,hspace=0.0)   
   ax1 = plt.subplot(gs[0,0])  # Original plot
   ax2 = plt.subplot(gs[0,1])  # Original plot
   ## Unit cell
   UCell(pos,latt,ax=ax1,show=False,tit='$Bravais$ $Lattice$')
   LG.debug('Plotted unit cell')
   ## FBZ
   FBZ(recip,N=N,ax=ax2,tit='$First$ $Brillouin$ $Zone$')
   LG.debug('Plotted FBZ')
   fig.tight_layout()
   if show: plt.show()
   if fname != None:
      LG.info('Bravais lattice saved in '+fname)
      fig.savefig(fname)



def FBZ(recip,N=3,ax=None,tit=None,show=False):
   """
     see the construction of the points to understand N
     could be the order of C symmetry????
   """
   if ax == None:
      fig = plt.figure() #figsize=(20,10))
      gs = gridspec.GridSpec(1, 1)
      fig.subplots_adjust(wspace=0.25,hspace=0.0)
      ax = plt.subplot(gs[0,0])  # Original plot

   # Empiric size of the arrow head
   if len(recip) > 0 :
      v_norm = np.nanmean([np.linalg.norm(v) for v in recip])
   else:
      v_norm = float('nan')
      LG.warning('No reciprocal lattice vectors received')
   hw = 0.025 * v_norm
   hl = 0.05 * v_norm
   for v in recip:
      vn = v/np.linalg.norm(v) # normalized vector
      vv = (np.linalg.norm(v)-hl)* vn #vector minus the length of the arrow
      ax.arrow(0,0,vv[0],vv[1],fc='b',ec='b',head_width=hw,head_length=hl)
   ## Special points
   LG.debug('Creating relevant points at a1/iN+a2/iN')
   lista = np.linspace(0,1,N+1)
   perms = [p for p in product(lista, repeat=len(recip)) ]
   points = []
   for p in perms:
      aux = np.array([0.,0.,0.])
      for i in range(len(p)):
         aux += p[i]*recip[i]
      points.append(aux)
   X = [p[0] for p in points]
   Y = [p[1] for p in points]
   ax.scatter(X,Y,c='r',s=100)

   if tit != None:
      if '$' in tit: ax.set_title(tit,fontsize=22)
      else: ax.set_title(tit)
   ax.axis('equal')
   ax.set_xlabel('$k_{x}$ $(2\pi / \AA)$',fontsize=17)
   ax.set_ylabel('$k_{y}$ $(2\pi / \AA)$',fontsize=17)
   ax.grid()
   if show: plt.show()


def vec_in_list(v,l,eps=0.000000001):
   """ Returns True if vector v is in the list of vectors l """
   for x in l:
      if np.linalg.norm(x-v) < eps: return True
   return False


def UCell(pos,latt=[],ax=None,tit=None,show=False):
   """
     Plots the unit cell, lattice vectors, and first neighbouring unit cells.
     ONLY VALID FOR 0D, 1D, 2D
   """
   if ax == None:
      fig = plt.figure() #figsize=(20,10))
      gs = gridspec.GridSpec(1, 1)
      fig.subplots_adjust(wspace=0.25,hspace=0.0)
      ax = plt.subplot(gs[0,0])  # Original plot

   ## Plot Unit cell
   X,Y = [],[]
   for r in pos:
      X.append(r[0])
      Y.append(r[1])
   sx = 0.075
   sy = 2* sx
   for i in range(len(pos)):
      r = pos[i]
      ax.text(r[0]+sx,r[1]+sy, str(i)) #, fontsize=20,bbox={'facecolor':'white', 'alpha':0.7, 'pad':5})
   ax.scatter(X,Y,c='k',s=100,edgecolors='none')
   LG.debug('Plotted %s atoms in the unitcell'%(len(pos)))
   ## Plot neighbouring cells
   if len(latt) > 0:
      fake_latt = []  # x, y, xmy, xy   XXX stores (vec,label)
      for i in range(len(latt)):
         fake_latt.append((latt[i],r'$\vec{a}_{%s}$'%(i+1)))
      aux = np.array([0.,0.,0.])
      ## XXX missing -a1+a2
      l = '$'
      for i in range(len(latt)):
         aux += latt[i]
         l += r'\vec{a}_{%s}+'%(i+1)
      l = l[0:-1] + '$'
      fake_latt.append((aux,l))
      cs = ['b','r','g','y','c','m']
      v_norm = np.mean([np.linalg.norm(v) for v in latt])
      ## Lattice vectors and cells
      hw = 0.04 * v_norm #v_norm * 0.17/3.46410161514
      hl = 0.05 * v_norm #hw * 0.3/0.2
      i = 0
      for vl in fake_latt:
         v,l = vl
         vn = v/np.linalg.norm(v) # normalized vector
         vv = (np.linalg.norm(v)-hl)* vn #vector minus the length of the arrow
         X,Y = [],[]
         for r in pos:
            w = r+v
            X.append(w[0])
            Y.append(w[1])
         ax.scatter(X,Y,c=cs[i],s=100,edgecolors='none')
         v_str = '(%.2f,%.2f,%.2f)'%(v[0],v[1],v[2])
         LG.debug('Plotted the unit cell shifted by: ' + v_str)
         ax.arrow(0,0,vv[0],vv[1],head_width=hw,head_length=hl,fc='b', ec='b')
         ax.text(v[0],v[1], l, fontsize=20,
                             bbox={'facecolor':'white', 'alpha':0.7, 'pad':5})
         i+=1

      ## Repited cells
      fake_latt = [v[0] for v in fake_latt]
      perms = [p for p in product(range(-1,2), repeat=len(latt)) ]
      perms = [p for p in perms if np.linalg.norm(p) != 0.]
      repited_vecs = []
      for p  in perms:
         aux = np.array([0.,0.,0.])
         for ip in range(len(p)):
            aux += p[ip] * latt[ip]
         if not vec_in_list(aux,fake_latt): repited_vecs.append(aux)
      for v in repited_vecs:
         X,Y = [],[]
         for r in pos:
            w = r+v
            X.append(w[0])
            Y.append(w[1])
         ax.scatter(X,Y,c=cs[i%(len(cs))],s=90,alpha=0.5,edgecolors='none')
         X = [0,v[0]]
         Y = [0,v[1]]
         ax.plot(X,Y,'b--',alpha=0.5)
         i+=1

   if tit != None:
      if '$' in tit: ax.set_title(tit,fontsize=22)
      else: ax.set_title(tit)
   ax.axis('equal')
   ax.set_xlabel('$x$ $(\AA)$',fontsize=17)
   ax.set_ylabel('$y$ $(\AA)$',fontsize=17)
   ax.grid()
   if show: plt.show()


def bands(X,Y,Z,show=False):
   LG.debug('Plotting %s points'%(len(X)))
   fig = plt.figure()
   gs = gridspec.GridSpec(1, 1)
   fig.subplots_adjust(wspace=0.,hspace=0.0)   
   ax = plt.subplot(gs[0,0])  # Original plot
   ax.scatter(X,Y,c=Z,s=20,cmap=my_cmap,edgecolors='none')
   ax.grid()
   ax.set_xlim([min(X),max(X)])
   ax.set_ylim([-10,10])
   if show:
      plt.tight_layout()
      plt.show()

from matplotlib import cm
def get_color(x,color_map='rainbow'):
   """ x needs to be normalized:  0<x<1 """
   cmap = cm.get_cmap(color_map)
   return cmap(x)  # rgba tuple

def spectrum(Es,Cs=[],ax=None,TOL=1e-5,Ef=0.0,y0=-10,y1=10,vb=False,show=False):
   """
    Plot the spectrum (as horizontal lines) with width proportional to the
    degeneracy of each state.
      TOL: energies that differ less than TOL are considered degenerate
      Ef: is used to plot a dashed line
      vb: verbose, to see states and degeneracy
      Cs: DEVELOPING to plot with color per state
   TODO: CHANGE!!!! better idea
   """
   if ax == None: fig, ax = plt.subplots()

   if isinstance(Es,list): Es = np.array(Es)
   if len(Cs) > 0: Cs = (np.array(Cs)-np.min(Cs))/abs(np.max(Cs)-np.min(Cs))
   else: pass
   ## Study degeneracy
   I = np.array([i for i in range(len(Es))])
   d = np.diff(Es)
   d = np.append(d,1000000)   # workaround for a fencepost error
   result = I[d>TOL]
   Es_set = Es[result]
   hist,bins = [],[]
   for e in Es_set:   # TODO Inefficient
      states = Es[abs(Es-e) < TOL]
      hist.append(e)
      bins.append(len(states))

   ## Plot
   x0,x1 = 0,1   # X limits (dummy)
   dx = abs(x0-x1)
   s=0.05*dx   # horizontal spacer between degenerated states
   dy = 0.1 * abs(y0 - y1)
   cont = 0    # if this works... it's magic
   for i in range(len(hist)):
      e = hist[i]
      w = bins[i]   # degeneracy
      LG.debug('%.5f / %s'%(e,int(w)))
      dxp = (dx - (w-1)*s)/w    # width of each of the degenerated states
      for j in range(w):
         if len(Cs) > 0: c = get_color( Cs[cont] )[0:-1]
         else: c = 'k'
         ax.plot([x0+j*(dxp+s),x0+(j+1)*dxp+j*s],[e,e],c=c,lw=2)
         cont += 1
   ax.axhline(Ef,color='k',ls='--',lw=1.5,alpha=0.2)
   ax.fill_between([x0-2*dx,x1+2*dx], Ef, -1000,facecolor='yellow',alpha=0.1)

   ax.set_xlim([x0-dx,x1+dx])
   ax.set_ylim([y0-dy,y1+dy])
   ax.set_ylim([-2,2])
   ax.set_xticklabels([])
   ax.set_ylabel('$E$ $(eV)$',fontsize=20)
   if show: plt.show()

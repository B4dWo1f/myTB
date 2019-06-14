#!/usr/bin/python3
# -*- coding: UTF-8 -*-


## Input file ##################################################################
## TODO maybe argparse this?
import sys
try: fini = sys.argv[1]
except IndexError: fini = 'SK.ini'

## Setup #######################################################################
#import setup
import load
##FP,HP,CP,SP,atoms,hoppings = setup.setup(fini)
FP,HP,CP,SP,atoms = load.setup(fini)
##SP.DOspin = False
##setup.compile_fortran('numeric.f95')


################################## LOGGING #####################################
import logging
import log_help
logging.basicConfig(level=logging.DEBUG,
                 format='%(asctime)s %(name)s:%(levelname)s - %(message)s',
                 datefmt='%Y/%m/%d-%H:%M:%S',
                 filename=FP.out+'main.log', filemode='w')
LG = logging.getLogger('main')
log_help.screen_handler(LG,lv='info')
################################################################################

print(FP)
print(HP)
print(CP)
print(SP)

print(' '*18,'-'*40,' '*20)

import IO
ats,pos,latt,sub = IO.read.xyz(SP.xyz_file)
if SP.force0D:
   LG.info('Forcing 0-dimensional system')
   latt =[]


## Base
# TODO make function in base.py
from basis import Base_Element,Base
elems = []
for i in range(len(pos)):
   a,r = ats[i],pos[i]
   elems.append( Base_Element(i,a,atoms,r) )

base = Base(elems,latt,atoms=atoms) #,cent=False)
base.get_neig(fol=FP.ham)
base.get_sublattice(sub)
base.get_layer()


LG.info('Duplicating basis')
base_pris = base.copy()
base_dfct = base.copy()
del base


LG.info('Defects for the basis')
if SP.vac.N > 0:
   IND_vac = base_dfct.vacancy(N=SP.vac.N,d=SP.vac.d,alpha=SP.vac.alpha) #,hollow=False)
   #IND_vac = base_dfct.vacancy(N=SP.vac.N,d=SP.vac.d,alpha=SP.vac.alpha) #,hollow=False)
else: IND_vac = []

if SP.ada.N >0:
   IND_ada = base_dfct.adatom(N=SP.ada.N,at='X', sp3=SP.ada.sp3)
   base_pris.adatom(N=SP.ada.N,at='X', dummy=True, sp3=SP.ada.sp3)


### Save basis
base_pris.save(FP.out+'pris.basis',FP.out+'pris.xyz')
base_dfct.save(FP.out+'dfct.basis',FP.out+'dfct.xyz')
LG.info('Base created')


import hamiltonian as ham
H_pris = ham.build_ham(base_pris,HP,'pris')
H_dfct = ham.build_ham(base_dfct,HP,'dfct')


LG.info('Check for spin')
if SP.DOspin:
   LG.warning('SPIN DOUBLING?')
   base_pris.dospin()
   base_dfct.dospin()
   # TODO testing
   base_pris.save(FP.out+'pris.basis.spin',FP.out+'pris.spin.xyz')
   base_dfct.save(FP.out+'dfct.basis.spin',FP.out+'dfct.spin.xyz')
   LG.info('Base with spin saved')


import operators as OP
import geometry as geo
from random import uniform, choice
#op = OP.orbital(base_dfct,'s')
import algebra as alg
#op = OP.spin(base_pris)
#op = OP.orbital(base_pris,'pz')
#op = OP.sublattice(base_pris)

if CP.spectrum:
   def aux(H,v0=None):
      """
        Dummy auxiliary function to calculate the spectrum in parallel
      """
      LG.info('Spectrum: %s'%(H.tag))
      es,v = H.get_N_states(Op=True,folder=FP.out,n=n_es,shw=Shw,pbc=SP.pbc,
                                                                         v0=v0)
      print('  ---- %s ----'%(H.tag))
      for e in es:
         print(e)
      return es,v
   Shw = False
   n_es = min([int(H_pris.dim//2),11])
   parallel = True
   import numpy as np
   if parallel:
      import multiprocessing as sub
      pool = sub.Pool(2)
      foo = pool.map(aux,[H_pris, H_dfct])
   else:
      v0 = None
      for h in [H_pris, H_dfct]:
         es,v = aux(h,v0)
         v0 = np.mean(v,axis=0)


if CP.bands:
   if len(latt) != 0:
      Shw = False
      full = False
      eigvec = True
      LG.info('Calculating bands')
      points = geo.get_points(base_pris.recip)
      G  = points[0]
      K  = points[6]
      Kp = points[9]
      M = (K+Kp)/2
      points = [G,K,Kp,G]
      #points = [K,G,M,Kp]
      #points = [0.9*K, K , 1.1*K]
      path = geo.recorrido(points,CP.nk)
      LG.debug('Calculating bands for Pristine')
      I,E,Z = H_pris.get_bands(path,V=eigvec,full=full,folder=FP.out,show=Shw)
      LG.info('Bands Pristine done')
      LG.debug('Calculating bands for Defected')
      v0 = Z[0]
      H_dfct.get_bands(path,V=eigvec,full=full,folder=FP.out,show=Shw,v0=v0)
      LG.info('Bands Defected done')
   else: LG.critical('No lattice vectors ==> No bands')


if CP.dos in ['full','window'] and len(latt) != 0:
   import numpy as np
   LG.info('Calculation DOS')
   Shw = False
   if CP.dos == 'full': full = True
   else: full = False
   if CP.nddos == None: Nddos = int(0.4*H_pris.dim)
   else: Nddos = CP.nddos
   if CP.nkdos == None: nk = 100
   else: nk = CP.nkdos
   path = []
   for ix in np.linspace(0,1,nk):
      for iy in np.linspace(0,1,nk):
         k = ix*base_pris.recip[0] + iy*base_pris.recip[1]
         path.append(k)
   H_pris.get_bands(path,V=CP.local,full=full,sigma=1e-6,k=Nddos,
                                              folder=FP.out,ext='dos',show=Shw)
   H_dfct.get_bands(path,V=CP.local,full=full,sigma=1e-6,k=Nddos,
                                              folder=FP.out,ext='dos',show=Shw)
else: LG.debug('No DOS calculations')

LG.info('All done. Bye!')
print('All done. Bye!')

print('='*80)
print('='*80)
exit()

#from calculations import get_DOS
##E = E[(E>-50) & (E<50)]
##mE,ME = min(E),max(E)
#mE,ME = -20,20
#nE = int((ME-mE)/0.1)
#E,Dp,Dd = get_DOS(mE,ME, H_dfct.intra,H_pris, path_slf=FP.slf,nE=nE,fol=FP.out,delta=0.001)
#
#import matplotlib.pyplot as plt
#fig, ax = plt.subplots()
#ax.plot(E,Dp,label='Pristine')
#ax.plot(E,Dd,label='Defected')
#ax.grid()
#ax.legend()
#plt.show()

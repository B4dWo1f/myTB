#!/usr/bin/python3
# -*- coding: UTF-8 -*-


## Input file ##################################################################
## TODO maybe argparse this?
import sys
try: fini = sys.argv[1]
except IndexError: fini = 'SK1.ini'

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
log_help.screen_handler(LG)
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

from time import time
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
   IND_vac = base_dfct.vacancy(N=SP.vac.N,d=None,alpha=SP.vac.alpha) #,hollow=False)
   #IND_vac = base_dfct.vacancy(N=SP.vac.N,d=SP.vac.d,alpha=SP.vac.alpha) #,hollow=False)
else: IND_vac = []

if SP.ada.N >0:
   IND_ada = base_dfct.adatom(N=SP.ada.N,at='H1')
   base_pris.adatom(N=SP.ada.N,at='H1', dummy=True)


### Save basis
told = time()
base_pris.save(FP.out+'pris.basis',FP.out+'base_pris.xyz')
base_dfct.save(FP.out+'dfct.basis',FP.out+'base_dfct.xyz')
LG.info('Base created')
print('            *** Base:',time()-told)


told = time()
import hamiltonian as ham
H_pris = ham.build_ham(base_pris,HP,'pris')
H_dfct = ham.build_ham(base_dfct,HP,'dfct')
print('     *** Hamiltonian:',time()-told)


LG.info('Check for spin')
if SP.DOspin:
   LG.warning('SPIN DOUBLING?')
   base_pris.dospin()
   base_dfct.dospin()


told = time()
import operators as OP
import geometry as geo
from random import uniform, choice
#op = OP.orbital(base_dfct,'s')
import algebra as alg
#op = OP.spin(base_pris)
#op = OP.orbital(base_pris,'pz')
#op = OP.sublattice(base_pris)
if CP.bands:
   if len(latt) != 0:
      Shw = False
      LG.info('Calculating bands')
      points = geo.get_points(base_pris.recip)
      points = [points[0],points[6],points[9], points[0]]
      path = geo.recorrido(points,CP.nk)
      LG.debug('Bands Pristine')
      I,E,Z = H_pris.get_bands(path,folder=FP.out,show=Shw)
      LG.info('Bands Pristine done')
      print('    ** Pris bands',time()-told)
      LG.debug('Bands Defected')
      H_dfct.get_bands(path,folder=FP.out,show=Shw)
      LG.info('Bands Defected done')
      print('    ** Dfct bands',time()-told)
   else: LG.critical('No lattice vectors ==> No bands')

if CP.spectrum:
   Shw = False
   LG.info('Spectrum: Pristine')
   import numpy as np
   n_es = min([int(H_pris.dim//2),8])
   #es,v = H_pris.get_N_states(Op=True,folder=FP.out,border=False,n=n_es)
   es,v = H_pris.get_N_states(Op=True,folder=FP.out,n=n_es,shw=Shw)
   #es,v = H_pris.get_spectrum(Op=True,folder=FP.out,border=False,shw=Shw)
   print('  ---- Pristine ----')
   for e in es:
      print(e)
   LG.info('Spectrum: Defected')
   n_es = min([int(H_dfct.dim//2),8])
   #es,v = H_dfct.get_N_states(Op=True,folder=FP.out,border=False,n=n_es)
   es,v = H_dfct.get_N_states(Op=True,folder=FP.out,n=n_es,shw=Shw)
   #es,v = H_dfct.get_spectrum(Op=True,folder=FP.out,border=False,shw=Shw)
   print('  ---- Defected ----')
   for e in es:
      print(e)
print('        *** Spectrum:',time()-told)

LG.info('All done. Bye!')
print('All done. Bye!')

#exit()
#print('='*80)
#print('='*80)
#
#
#from calculations import get_DOS
#E = E[(E>-50) & (E<50)]
#mE,ME = min(E),max(E)
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

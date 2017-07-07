#!/usr/bin/python3
# -*- coding: UTF-8 -*-


import sys
try: fini = sys.argv[1]
except IndexError: fini = 'SK1.ini'
import setup
FP,HP,CP,SP,atoms,hoppings = setup.setup(fini)

############################### LOGGING #####################################
import logging
import log_help
logging.basicConfig(level=logging.DEBUG,
                 format='%(asctime)s %(name)s:%(levelname)s - %(message)s',
                 datefmt='%Y/%m/%d-%H:%M:%S',
                 filename=FP.out+'main.log', filemode='w')
LG = logging.getLogger('main')
log_help.screen_handler(LG)
#############################################################################


print(FP)
print(HP)
print(CP)
print(SP)


## TODO thake this function somewhere else!!
import os
def compile_fortran(fname):
   def doit(fname):
      LG.debug('Backup file (.%s) not found'%(fname))
      LG.info('Compilando fortran con f2py')
      os.system('f2py -c -m %s %s'%(root_fname,fname))
      LG.info('   ...Compilado fortran con f2py')
      os.system('cp %s .%s'%(fname,fname))
      LG.warning('Hidden copy to avoid re-compiling')
   root_fname = '.'.join(fname.split('.')[0:-1])
   if not os.path.exists('.%s'%(fname)): doit(fname)
   else:
      LG.debug('Backup file (.%s) is present'%(fname))
      diff_for = os.popen('diff %s .%s'%(fname,fname)).read()
      diff_for = diff_for.lstrip().rstrip()
      diff_for.splitlines()
      so = os.popen('ls %s.*so 2> /dev/null'%(root_fname)).read()
      if len(diff_for) > 1 or len(so) == 0: doit(fname)
      else: LG.info('%s is already compiled'%(fname))

compile_fortran('numeric.f95')

print(' '*18,'-'*40,' '*20)

import IO
ats,pos,latt,sub = IO.xyz(SP.xyz_file)
latt = []


## Base
from basis import Base_Element,Base
elems = []
for i in range(len(pos)):
   a,r = ats[i],pos[i]
   elems.append( Base_Element(i,a,atoms[a],r) )

from time import time
base = Base(elems,latt)
base.get_neig(fol=FP.ham)
base.get_sublattice(sub)
base.get_layer()

LG.info('Duplicating basis')
base_pris = base.copy()
base_dfct = base.copy()
del base


if SP.vac.N > 0:
   IND_vac = base_dfct.vacancy(N=SP.vac.N,d=SP.vac.d,alpha=SP.vac.alpha)


## Save basis
told = time()
base_pris.save(FP.out+'pris.basis',FP.out+'base_pris.xyz')
base_dfct.save(FP.out+'dfct.basis',FP.out+'base_dfct.xyz')
LG.info('Base created')
print('            *** Base:',time()-told)



told = time()
## Create the Hamiltonians
import hamiltonian as ham
LG.info('Creating Pristine Hamiltonian')
# Pristine
LG.debug('Starting kinetic terms')
Htot = ham.kinetic(base_pris,hoppings)
if HP.lelec != 0.0:
   LG.info('Electric field: %s'%(HP.lelec))
   Htot.append( ham.electric(base_pris,HP.lelec) )
LG.info('Hamiltonian ready')
H_pris = ham.Hamiltonian(Htot,tag='pris')
H_pris.names()
del Htot
# Defected
LG.info('Creating Defected Hamiltonian')
LG.debug('Starting kinetic terms')
Htot = ham.kinetic(base_dfct,hoppings)
if HP.lelec != 0.0:
   Htot.append( ham.electric(base_dfct,HP.lelec) )
LG.info('Hamiltonian ready')
H_dfct = ham.Hamiltonian(Htot,tag='dfct')
H_dfct.names()
H_dfct.save_matrix(FP.ham)
LG.info('Hamiltonians done')
print('     *** Hamiltonian:',time()-told)


told = time()
import operators as OP
import geometry as geo
#from random import uniform, choice
#op = OP.orbital(base_dfct,'s')
if CP.bands:
   Shw = False
   LG.info('Calculating bands')
   points = geo.get_points(base_pris.recip)
   points = [points[0],points[6],points[9], points[0]]
   path = geo.recorrido(points,CP.nk)
   LG.debug('Bands Pristine')
   H_pris.get_bands(path,folder=FP.out,show=Shw)
   LG.info('Bands Pristine done')
   LG.debug('Bands Defected')
   #print(IND_vac)
   op = OP.atom(base_dfct,[IND_vac[0]])
   H_dfct.get_bands(path,Op=op,folder=FP.out,show=Shw)
   LG.info('Bands Defected done')

if CP.spectrum:
   #Shw = False
   LG.info('Spectrum: Pristine')
   #es,_ = H_pris.get_spectrum(Op=True,folder=FP.out)   #,show=Shw)
   es,v = H_pris.get_N_states(Op=True,folder=FP.out)
   LG.info('Spectrum: Defected')
   #op = OP.atom(base_dfct,[IND_vac[0]])
   #es,_ = H_dfct.get_spectrum(Op=True,folder=FP.out)   #,show=Shw)
   es,v = H_dfct.get_N_states(Op=True,folder=FP.out)
print('        *** Spectrum:',time()-told)

LG.info('All done. Bye!')
exit()
print('bla bla')

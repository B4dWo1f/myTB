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


## Base
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


if SP.vac.N > 0:
   IND_vac = base_dfct.vacancy(N=SP.vac.N,d=SP.vac.d,alpha=SP.vac.alpha)
else: IND_vac = []

if SP.ada.N >0:
   IND_ada = base_dfct.adatom(N=SP.ada.N)
   base_pris.adatom(N=SP.ada.N, dummy=True)


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
   Htot.append( ham.pseudo_rashba(base_pris,HP.lelec) )
   #Htot.append( ham.electric(base_pris,HP.lelec) )
LG.info('Hamiltonian ready')
H_pris = ham.Hamiltonian(Htot,tag='pris')
H_pris.names()
del Htot
# Defected
LG.info('Creating Defected Hamiltonian')
LG.debug('Starting kinetic terms')
Htot = ham.kinetic(base_dfct,hoppings)
if HP.lelec != 0.0:
   Htot.append( ham.pseudo_rashba(base_dfct,HP.lelec) )
   #Htot.append( ham.electric(base_dfct,HP.lelec) )
LG.info('Hamiltonian ready')
H_dfct = ham.Hamiltonian(Htot,tag='dfct')
H_dfct.names()
H_dfct.save_matrix(FP.ham)
LG.info('Hamiltonians done')
print('     *** Hamiltonian:',time()-told)
del Htot


told = time()
import operators as OP
import geometry as geo
#from random import uniform, choice
#op = OP.orbital(base_dfct,'s')
if CP.bands:
   if len(latt) != 0:
      Shw = True
      LG.info('Calculating bands')
      points = geo.get_points(base_pris.recip)
      points = [points[0],points[6],points[9], points[0]]
      path = geo.recorrido(points,CP.nk)
      LG.debug('Bands Pristine')
      I,E,Z = H_pris.get_bands(path,folder=FP.out,show=Shw)
      LG.info('Bands Pristine done')
      LG.debug('Bands Defected')
      H_dfct.get_bands(path,folder=FP.out,show=Shw)
      LG.info('Bands Defected done')
   else: LG.critical('No lattice vectors ==> No bands')

if CP.spectrum:
   Shw = False
   LG.info('Spectrum: Pristine')
   import numpy as np
   n_es = min([H_pris.dim,10])
   #es,v = H_pris.get_N_states(Op=True,folder=FP.out,border=False,n=n_es)
   es,v = H_pris.get_N_states(Op=True,folder=FP.out,n=n_es,shw=Shw)
   #es,v = H_pris.get_spectrum(Op=True,folder=FP.out,border=False,shw=Shw)
   print('  ---- Pristine ----')
   for e in es:
      print(e)
   LG.info('Spectrum: Defected')
   n_es = min([H_dfct.dim,10])
   #es,v = H_dfct.get_N_states(Op=True,folder=FP.out,border=False,n=n_es)
   es,v = H_dfct.get_N_states(Op=True,folder=FP.out,n=n_es,shw=Shw)
   #es,v = H_dfct.get_spectrum(Op=True,folder=FP.out,border=False,shw=Shw)
   print('  ---- Defected ----')
   for e in es:
      print(e)
print('        *** Spectrum:',time()-told)

LG.info('All done. Bye!')

print('='*80)
print('='*80)


#import mygreen_tools as green
#import numpy as np
#import os
#
#def get_DOS(Emin,Emax,vintra,h,path_slf='/tmp',nE=101,use_all=True,fol='./'):
#   de = 0.1* (Emax - Emin)
#   E = np.linspace(int(Emin-de),int(Emax+de),nE)
#   if use_all:
#      files = os.popen('ls -1 %s*.npy'%(path_slf)).read().splitlines()
#      es = [x.split('/')[-1].replace('.npy','') for x in files]
#      es = [x.split('/')[-1].replace('self','') for x in es]
#      es = np.array(sorted([float(x) for x in es]))
#      E = np.sort(np.append(E,es))
#   DOSd,DOSp = [],[]
#   f = open(fol+'dfct.dos','w')
#   g = open(fol+'pris.dos','w')
#   for e in E:
#      Gd = green.green_function(e,H_dfct.intra,H_pris,path_selfes=path_slf)
#      Gp = green.green_function(e,H_pris.intra,H_pris,path_selfes=path_slf)
#      dd = -Gd.trace()[0,0].imag/np.pi
#      dp = -Gp.trace()[0,0].imag/np.pi
#      f.write(str(e)+'   ')
#      g.write(str(e)+'   ')
#      vd = -Gd.diagonal().imag/np.pi
#      for ix in range(max(vd.shape)):
#         f.write(str(vd[0,ix])+'   ')
#      f.write('\n')
#      vp = -Gp.diagonal().imag/np.pi
#      for ix in range(max(vp.shape)):
#         g.write(str(vp[0,ix])+'   ')
#      g.write('\n')
#      DOSd.append(dd)
#      DOSp.append(dp)
#   f.close()
#   g.close()
#   return E,DOSp,DOSd

from calculations import get_DOS
E = E[(E>-50) & (E<50)]
mE,ME = min(E),max(E)
nE = 2*int((ME-mE)/0.1)
E,_,DOS = get_DOS(mE,ME, H_dfct.intra,H_pris, path_slf=FP.slf,nE=nE,fol=FP.out)
import matplotlib.pyplot as plt
fig, ax = plt.subplots()
ax.plot(E,DOS)
ax.grid()
plt.show()

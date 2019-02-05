#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import pp
import numpy as np
#from scipy.linalg import eig_banded
from scipy.sparse.linalg import eigsh
import algebra as alg
from tqdm import tqdm
import logging
import log_help
LG = logging.getLogger(__name__)

eps = 1e-5

def diagon(Hm,K,Op,sigma=0,n=0):
   """ Full diagonalization of the Hamiltonian for a given k """
   kx,ky,kz = K
   #H = Hamil(Htot,[kx,ky,kz])
   H = Hm.get_k(K).todense()
   if Op : return np.linalg.eigh(H)
   else: return np.linalg.eigvalsh(H)

def diagon_window(Hm,K,Op,sigma=0,n=5,v0=None):
   """ Diagonalize the hamiltonian for a given k in a window of energy
   Hk is a function H(k)
   """
   kx,ky,kz = K
   #H = Hamil(Htot,[kx,ky,kz])
   H = Hm.get_k(K)
   n = min([H.shape[0]-2,n])  # Protection for not enough eigvals
   if Op:
      return eigsh(H, k=n+1, sigma=sigma, which='LM',
                   return_eigenvectors=True,v0=v0)
   else:
      return eigsh(H, k=n+1, sigma=sigma, which='LM',
                   return_eigenvectors=False,v0=v0) #,maxiters=10000)


def bands(RECORRIDO,H,V=False,sigma=0,n=5,full=False,v0=None):
   """
      H is a Hamitlonian object
   """
   X, Y, Z = [], [], []
   cont=0
   #if n<=H.dim-1: diag = 
   print(n,H.dim,'',V)
   if n >= H.dim -1 or full:
      print('Full Diagonalization')
      diag = diagon
   else:
      print('Partial Diagonalization')
      diag = diagon_window
   for k in tqdm(RECORRIDO, unit='K-points'):
      if V:
        es, vs = diag(H,k,V,sigma,n,v0=v0)
        vs = vs.transpose()
        v0 = np.mean(vs,axis=0)
        for e,v in zip(es,vs):
           if e.imag > eps: sys.exit('Hamiltoniano no hermitico')
           X.append(cont)
           Y.append(e.real)
           Z.append(v)
      else:
        es = diag(H,k,V,sigma,n)
        #eig = diagon(H,k,Op)
        for e in es:
           if e.imag > eps: sys.exit('Hamiltoniano no hermitico')
           X.append(cont)
           Y.append(e.real)
           Z.append(0)
      cont += 1
   return X, Y, Z




#@log_help.log2screen(LG)
#def bandsPP(RECORRIDO,Htot,Op=False,sigma=None,n=None,ncpus=None,eps=0.00001):
#   """
#     Calculates the band in parallel along a certain path
#     Es = [E0, E1] determines the window in which to diagonalize
#   """
#   cont = 0
#   if sigma != None and n != None:
#      func = diagon_window
#      LG.info('Diagonalize the %s closest eigenvalues to %s'%(n,sigma))
#   else:
#      func = diagon
#      LG.info('Full diagonalization')
#   inputs = []
#   for k in RECORRIDO:
#      inputs.append([cont,k])
#      cont += 1
#   LG.debug('Diagonalize in %s k points'%(len(inputs)))
#
#   # tuple of all parallel python servers to connect with
#   ppservers = ()
#   if not ncpus:
#     try: ncpus = int(open('NCPUS','r').readlines()[0].rstrip('\n'))
#     except IOError:
#        ncpus = 1
#        LG.info('Not Parallel computation of the bands')
#        X, Y, V = bands(RECORRIDO,Htot,Op)
#        LG.info('   ...bands done')
#        return X, Y, V
#
#   js = pp.Server(ncpus, ppservers=ppservers)
#   msg='Creating PP server with %s CPUs in %s servers'%(ncpus,len(ppservers)+1)
#   LG.debug(msg)
#   jobs = [(i, js.submit(func, (Htot,k,Op,sigma,n), (),\
#           ("import numpy as np","from hamiltonian import Hamil",
#            "import algebra as alg",
#            "from scipy.sparse.linalg import eigsh"))) for i,k in inputs]
#   #job_server.print_stats()
#   X1, Y1 = [], []
#   for aux,job in jobs:
#      X1.append(aux)
#      Y1.append(job())
#
#   X, Y, V = [], [], []
#   if Op:
#     for x,y in zip(X1,Y1):
#        R = y[1].transpose()  # R tiene como filas los autovectores de H(k)
#        for E,vec in zip(y[0],R):
#           X.append(x)
#           Y.append(E)
#           V.append(np.matrix(vec))
#   else:
#     for x,y in zip(X1,Y1):
#        for E in y:
#           X.append(x)
#           if E.imag > eps: LG.critical('Complex eigenvalue!!')
#           Y.append(E.real)
#           V.append(0)
#   return X, Y, V
#
#
#def bandsPP2D(RECORRIDO,Htot,ind=0,ncpus=None,eps=0.00001):
#   """
#     Calculates the bands in 2D, returning Kx, Ky, ind-Eigval 
#   """
#   cont = 0
#   inputs = []
#   for k in RECORRIDO:
#      inputs.append([cont,k])
#      cont += 1
#
#   # tuple of all parallel python servers to connect with
#   ppservers = ()
#   if ncpus == None:  ncpus = 1
#   job_server = pp.Server(ncpus, ppservers=ppservers)
#   jobs = [(i, job_server.submit(diagon, (Htot,k,False,), (),\
#           ("import numpy as np","from hamiltonian import Hamil",))) for i,k in inputs]
#   #job_server.print_stats()
#   X1, Y1 = [], []
#   for aux,job in jobs:
#      X1.append(aux)
#      Y1.append(job())
#
#   X, Y, V = [], [], []
#   for i in range(len(Y1)):  #x,y in zip(X1,Y1):
#      k = RECORRIDO[X1[i]]
#      X.append(k[0])
#      Y.append(k[1])
#      V.append(Y1[i][ind].real)
#   return X, Y, V

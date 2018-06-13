#! /usr/bin/env python
# -*- coding: UTF-8 -*-

import pp
import numpy as np
#from scipy.linalg import eig_banded
from scipy.sparse.linalg import eigsh
import algebra as alg
import logging
import log_help
LG = logging.getLogger(__name__)

eps = 1e-5
####################################
#def Hamil(Hlist,k,chk=True):
#   """
#     Hlist: list of HTerms (class) which contains the name of the element,
#            the coupling, the matrix and the exponential
#     k: np.array k point in which we want to evaluate the hamiltonian
#   *** Assumes only intra,x,y,xy terms are provided
#   """
#   Hamiltoniano = np.matrix(np.zeros(Hlist[0].mat.shape,dtype=complex))
#   for Hterm in Hlist:
#      l = Hterm.coup
#      M = Hterm.mat #.todense()   # XXX
#      v = Hterm.exp
#      if np.linalg.norm(v) != 0.:
#         Hamiltoniano += l * M * np.exp(-1.j*np.dot(k,v))
#         Hamiltoniano += l * M.H * np.exp(-1.j*np.dot(k,-v))
#      else: Hamiltoniano += l * M * np.exp(-1.j*np.dot(k,v))
#   if chk:      # Check for Hermiticity
#      A = Hamiltoniano-Hamiltoniano.H
#      if np.allclose(A,np.zeros(A.shape,dtype=complex)): pass
#      else:
#         msg = 'Hamiltonian is not hermitian H(%.2f,%.2f,%.2f)'%(k[0],k[1],k[2])
#         LG.critical(msg)
#         sys.exit(1)
#   return Hamiltoniano
####################################

#def diagon(Htot,K,Op,sigma=0,n=0):
#   """ Full diagonalization of the Hamiltonian for a given k """
#   kx,ky,kz = K
#   H = Hamil(Htot,[kx,ky,kz])
#   if Op : return np.linalg.eigh(H)
#   else: return np.linalg.eigvalsh(H)

def diagon_window(Hm,K,Op,sigma=0,n=5):
   """ Diagonalize the hamiltonian for a given k in a window of energy
   Hk is a function H(k)
   """
   kx,ky,kz = K
   #H = Hamil(Htot,[kx,ky,kz])
   H = Hm.get_k(K)
   n = min([H.shape[0]-2,n])  # Protection for not enough eigvals
   if Op: return eigsh(H,k=n+1,sigma=sigma,which='LM',return_eigenvectors=True)
   else: return eigsh(H,k=n+1,sigma=sigma,which='LM',return_eigenvectors=False)


@log_help.log2screen(LG)
def bandsPP(RECORRIDO,Htot,Op=False,sigma=None,n=None,ncpus=None,eps=0.00001):
   """
     Calculates the band in parallel along a certain path
     Es = [E0, E1] determines the window in which to diagonalize
   """
   cont = 0
   if sigma != None and n != None:
      func = diagon_window
      LG.info('Diagonalize the %s closest eigenvalues to %s'%(n,sigma))
   else:
      func = diagon
      LG.info('Full diagonalization')
   inputs = []
   for k in RECORRIDO:
      inputs.append([cont,k])
      cont += 1
   LG.debug('Diagonalize in %s k points'%(len(inputs)))

   # tuple of all parallel python servers to connect with
   ppservers = ()
   if not ncpus:
     try: ncpus = int(open('NCPUS','r').readlines()[0].rstrip('\n'))
     except IOError:
        ncpus = 1
        LG.info('Not Parallel computation of the bands')
        X, Y, V = bands(RECORRIDO,Htot,Op)
        LG.info('   ...bands done')
        return X, Y, V

   js = pp.Server(ncpus, ppservers=ppservers)
   msg='Creating PP server with %s CPUs in %s servers'%(ncpus,len(ppservers)+1)
   LG.debug(msg)
   jobs = [(i, js.submit(func, (Htot,k,Op,sigma,n), (),\
           ("import numpy as np","from hamiltonian import Hamil",
            "import algebra as alg",
            "from scipy.sparse.linalg import eigsh"))) for i,k in inputs]
   #job_server.print_stats()
   X1, Y1 = [], []
   for aux,job in jobs:
      X1.append(aux)
      Y1.append(job())

   X, Y, V = [], [], []
   if Op:
     for x,y in zip(X1,Y1):
        R = y[1].transpose()  # R tiene como filas los autovectores de H(k)
        for E,vec in zip(y[0],R):
           X.append(x)
           Y.append(E)
           V.append(np.matrix(vec))
   else:
     for x,y in zip(X1,Y1):
        for E in y:
           X.append(x)
           if E.imag > eps: LG.critical('Complex eigenvalue!!')
           Y.append(E.real)
           V.append(0)
   return X, Y, V


def bandsPP2D(RECORRIDO,Htot,ind=0,ncpus=None,eps=0.00001):
   """
     Calculates the bands in 2D, returning Kx, Ky, ind-Eigval 
   """
   cont = 0
   inputs = []
   for k in RECORRIDO:
      inputs.append([cont,k])
      cont += 1

   # tuple of all parallel python servers to connect with
   ppservers = ()
   if ncpus == None:  ncpus = 1
   job_server = pp.Server(ncpus, ppservers=ppservers)
   jobs = [(i, job_server.submit(diagon, (Htot,k,False,), (),\
           ("import numpy as np","from hamiltonian import Hamil",))) for i,k in inputs]
   #job_server.print_stats()
   X1, Y1 = [], []
   for aux,job in jobs:
      X1.append(aux)
      Y1.append(job())

   X, Y, V = [], [], []
   for i in range(len(Y1)):  #x,y in zip(X1,Y1):
      k = RECORRIDO[X1[i]]
      X.append(k[0])
      Y.append(k[1])
      V.append(Y1[i][ind].real)
   return X, Y, V



def bands(RECORRIDO,H,Op=False):
   """
      H is a Hamitlonian object
   """
   X, Y, Z = [], [], []
   cont=0
   for k in RECORRIDO:
      if Op:
        eig, eigvec = diagon(H,k,Op)
        eigvec = eigvec.transpose()
        for eigval,eigvec in zip(eig,eigvec):
           v1 = eigvec
           if eigval.imag > eps: sys.exit('Hamiltoniano no hermitico')
           X.append(cont)
           Y.append(eigval.real)
           Z.append(v1)
      else:
        eig = diagon_window(H,k,Op)
        for eigval in eig:
           X.append(cont)
           Y.append(eigval.real)
           Z.append(0)
      cont += 1
   return X, Y, Z

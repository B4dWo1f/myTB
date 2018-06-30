#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np
from numpy import matrix, diag, zeros
from scipy.sparse import coo_matrix, bmat,csc_matrix
import bands
import graphs
import IO
import sys
import algebra as alg
from time import time
import logging
import log_help
LG = logging.getLogger(__name__)   # Logger for this module

l_orb = {'s':0,'p':1,'d':2}
orb_group = {'s':'s','px':'p','py':'p','pz':'p',
             'dxy':'d','dyz':'d','dzx':'d','dx2y2':'d','d3z2r2':'d'}

class HTerm(object):
   """ Elements of the Hamiltonian """
   def __init__(self, matrix, exp,coupling=1.,name=''):
      self.name = name
      self.coup = coupling
      self.mat = matrix
      self.exp = exp
   def __str__(self):
      if self.name != '': msg = 'Term: '+self.name+'\n'
      else: msg = 'Term: \n'
      msg += 'coupling: %s\n'%(self.coup)
      v = self.exp
      msg += '  vector: (%5.3f, %5.3f, %5.3f)\n'%(v[0],v[1],v[2])
      msg += str(self.mat)
      return msg

class Hamiltonian(object):
   @log_help.log2screen(LG)
   def __init__(self, Hlist,dospin=False,tag=''):
      ## Check dimensionality
      LG.info('Creating Hamiltonian: %s    (spin:%s)'%(tag,dospin))
      dims = [x.mat.shape[0] for x in Hlist]
      if len(set(dims)) == 1 and not dospin:
         pass  # same dimension for every term
      elif len(set(dims)) == 2 or dospin:
         LG.info('Spin doubling')
         md, Md = min(dims), max(max(dims),2*min(dims))
         for i in range(len(Hlist)):
            x = Hlist[i]
            if x.mat.shape[0] == md:
               LG.info('spin doubling term: %s'%(x.name))
               x.mat = alg.m2spin(x.mat)
               if x.mat.shape[0] != Md:
                  LG.critical('Error doubling spin in term %s'%(x.name))
      else:
         LG.critical('Different dimensions not coherent with spin')
         exit()
      ## Old stuff
      self.lista = Hlist
      try: self.dim = int(list(set([x.mat.shape[0] for x in Hlist]))[0])
      except TypeError:
         LG.critical('Different dimensions not coherent with spin')
         exit()
      self.tag = tag
   def __iter__(self): return (x for x in self.lista)
   def save_matrix(self,folder='./'):
      for i in range(len(self.lista)):
         H = self.lista[i]
         if len(H.name) > 0: fname = folder+H.name+'_%s.dat'%(self.tag)
         else: fname = folder+'H_%s%s.dat'%(self.tag,i)
         M = H.mat
         v = H.exp
         LG.info('Saving matrix %s'%(H.name))
         IO.write.mat(fname,M,v)
         LG.debug('  ...done')
   def get_hk_gen(self):
      """ Generate kdependent hamiltonian"""
      return hk_gen(self)
   def get_k(self,k): return Hamil(self.lista,k)
   def get_N_states(self,Op=False,border=True,n=7,sigma=0,
                                                        folder='./',shw=False):
      from scipy.sparse.linalg import eigsh
      try: self.intra
      except: self.names()
      LG.info('In get_N_states')
      if border: H = csc_matrix(self.intra)       # Island
      else: H = csc_matrix( self.get_k(np.array((0,0,0))) )  # Peri bound cond
      LG.info('H acquired')
      if Op:
         LG.info('Start Diagonalization')
         kn = n+1
         if not 1 < kn < H.shape[0]-1:
         #if n+1 >= H.shape[0]:   #XXX check
            print('FULL MATRIX')
            print('-=-=-=-=-=-=-=-=-=-=-=-=-=-')
            LG.warning('Diagonalizing full matrix')
            es,v = np.linalg.eigh(H.todense())
         else:
            es,v = eigsh(H,k=kn,sigma=sigma+0.000001,which='LM',
                                                    return_eigenvectors=True)
         v = v.transpose()
         ind_ord = np.argsort(es)   #XXX check that this works as expected
         es=es[ind_ord]
         v=v[ind_ord]
      else:
         LG.info('Start Diagonalization')
         es = eigsh(H,k=n+1,sigma=sigma,which='LM',return_eigenvectors=False)
         v = [0 for _ in es]
      bname = folder+'%s_spectrum'%(self.tag)
      LG.debug('Writing spectrum to: '+bname)
      np.save(bname,np.column_stack((es,v)))
      if shw: graphs.spectrum(es,show=True)
      return es,v
   def get_spectrum(self,Op=False,border=True,folder='./',shw=False):
      """
       border: if True ---> open bound condition
               if False ---> closed bound condition
       Op: save eigenvalues or not
      """
      try: self.intra
      except: self.names()
      LG.info('In get_spectrum')
      if border: H = self.intra                 # Island
      else: H = self.get_k(np.array((0,0,0)))   # Periodic boundary conditions
      H = H.todense()
      LG.info('H acquired')
      if not Op:
         LG.info('Diagonalize without eigenvectors')
         es = np.linalg.eigvalsh(H)
         v = [0 for _ in es]
      else:
         LG.info('Diagonalize with eigenvectors')
         es,v = np.linalg.eigh(H)
         v = v.transpose()
         ind_ord = np.argsort(es)   #XXX check that this works as expected
         es=es[ind_ord]
         v=v[ind_ord]
         #Cs = [(v * Op * v.H)[0,0].real for v in Z]
      bname = folder+'%s_spectrum'%(self.tag)
      LG.debug('Writing spectrum to: '+bname)
      np.save(bname,np.column_stack((es,v)))
      if shw:
         LG.info('Plotting')
         graphs.spectrum(es,show=True)
      return es,v
   def get_bands(self,path,Op=False,sigma=None,k=None,show=False,ncpus=4,
                                                                  folder='./'):
      if not isinstance(Op,bool):
         LG.info('Bands with eigevectors')
         Opp=True
      else:
         LG.info('Bands without eigevectors')
         Opp = Op
      #X,Y,Z = bands.bandsPP(path,self.lista,Op=Opp,sigma=sigma,n=k,ncpus=ncpus)
      X,Y,Z = bands.bands(path,self,Op=Opp) #,sigma=sigma,n=k,ncpus=ncpus)
      if Opp: Z = [(v * Op * v.H)[0,0].real for v in Z]
      bname = folder+'%s.bands'%(self.tag)
      LG.debug('Writing bands to: '+bname)
      f = open(bname,'w')
      for x,y,z in zip(X,Y,Z):
         f.write('%s   %s   %s\n'%(x,y,z))
      f.close()
      LG.info('Bands saved to: '+bname)
      if show: graphs.bands(X,Y,Z,show=True)
      return np.array(X),np.array(Y),np.array(Z)  #TODO check type compatib
   #DEPRECATED
   #def dospin(self):
   #   Hs = self.lista
   #   for h in Hs:
   #      h.mat = coo_matrix(alg.m2spin(h.mat))
   #   self.names()
   #   self.dim = self.intra.shape[0]
   #def disconnect(self,indices=[],inf=100000,hop=False):
   #   """
   #     This method will disconnect a given atom leaving inf for the onsite
   #     energy and putting 0 to every hopping.
   #     the provided indices are the rows/cols in the matrix to be put to
   #     zero or infinity
   #   """
   #   for ih in range(len(self.lista)):
   #      h = self.lista[ih]
   #      if h.name in ['intra']: #,'x','y','xy','xmy']:
   #         M = h.mat
   #         for i in indices:
   #            if hop:
   #               M[i,:] = 0.0   # Switch off hoppings as well
   #               M[:,i] = 0.0
   #            M[i,i] = inf
   #         self.lista[ih].mat = coo_matrix(M)
   #         self.intra = M   #XXX call names here?
   #def pick_orbitals(self,orbitals): #XXX DEPRECATED
   #   """
   #     DEPRECATED
   #     This method picks by index the elements of the Hamiltonian
   #   """
   #   for Ht in self.lista:
   #      Ht.mat = coo_matrix( Ht.mat.todense()[orbitals,:][:,orbitals] )
   #   self.dim = self.lista[0].mat.shape[0]
   def names(self,d=2):
      """ Assumes the correct order 0,a1,a2,a1+a2,a1-a2 """
      ## XXX WRONG!!!
      dd = 0
      zero = self.lista[0].mat * 0.0
      self.intra = zero  # initialize to zero
      for h in self.lista:
         if np.linalg.norm(h.exp) == 0.:
            LG.info('Adding %s term to intra with coupling %s'%(h.name,h.coup))
            self.intra += h.coup*h.mat
      try:
         self.tx = self.lista[1].mat
         dd += 1
      except IndexError: self.tx = zero
      try:
         self.ty = self.lista[2].mat
         dd += 1
      except IndexError: self.ty = zero
      try: self.txy = self.lista[3].mat
      except IndexError: self.txy = zero
      try: self.txmy = self.lista[4].mat
      except IndexError: self.txmy = zero
      LG.debug('Expected dimensionality: %s'%(dd))
      self.dimensionality = d


def hk_gen(h):
   """ Returns a function that generates a k dependent hamiltonian"""
   if h.dimensionality == 0: return None
   if h.dimensionality == 1:
      def hk(k):
         """k dependent hamiltonian, k goes from 0 to 1"""
         tk = h.inter * np.exp(1j*np.pi*2.*k)
         ho = h.intra + tk + tk.H
         return ho
      return hk  # return the function
   if h.dimensionality == 2:
      def hk(k):
         """k dependent hamiltonian, k goes from (0,0) to (1,1)"""
         k = np.array(k)
         ux = np.array([1.,0.])
         uy = np.array([0.,1.])
         ptk = [[h.tx,ux],[h.ty,uy],[h.txy,ux+uy],[h.txmy,ux-uy]]
         ho = (h.intra).copy() # intraterm
         for p in ptk: # loop over hoppings
            tk = p[0]*np.exp(1j*np.pi*2.*(p[1].dot(k)))  # add bloch hopping
            ho += tk + tk.H  # add bloch hopping
         return ho
      return hk

def Hamil(Hlist,k,chk=True):
   """
      Hlist: list of HTerms (class) which contains the name of the element,
        the coupling, the matrix and the exponential
      k: np.array k point in which we want to evaluate the hamiltonian
      *** Assumes only intra,x,y,xy terms are provided
   """
   Hamiltoniano = Hlist[0].mat * 0.0
   for Hterm in Hlist:
      l = Hterm.coup
      M = Hterm.mat #.todense()   # XXX
      v = Hterm.exp
      if np.linalg.norm(v) != 0.:
         Hamiltoniano += l * M * np.exp(-1.j*np.dot(k,v))
         Hamiltoniano += l * M.H * np.exp(-1.j*np.dot(k,-v))
      else: Hamiltoniano += l * M * np.exp(-1.j*np.dot(k,v))
   #if chk:      # Check for Hermiticity
   #   A = Hamiltoniano-Hamiltoniano.H
   #   if np.allclose(A,np.zeros(A.shape,dtype=complex)): pass
   #   else:
   #      msg = 'Hamiltonian is not hermitian H(%.2f,%.2f,%.2f)'%(k[0],k[1],k[2])
   #      LG.critical(msg)
   #      sys.exit(1)
   return Hamiltoniano

@log_help.disable2(LG)
@log_help.log2screen(LG)
def build_ham(base,hp,tag,dospin=False):
   LG.info('Creating %s Hamiltonian (spin:%s)'%(tag,dospin))
   LG.debug('Starting kinetic terms')
   Htot = kinetic(base,hp.hoppings)
   if np.linalg.norm(hp.lzee) != 0.0:
      LG.info('Zeeman Term: (%s,%s,%s)'%(hp.lzee[0],hp.lzee[1],hp.lzee[2]))
      Htot.append( zeeman(base,hp.lzee) )
   if hp.lmass != 0.0:
      LG.info('Mass Term: %s'%(hp.lmass))
      Htot.append( mass(base,hp.lmass) )
   if hp.lSO != 0.0:
      LG.info('Spin-Orbit coupling: %s'%(hp.lSO))
      Htot.append( soc(base,hp.lSO) )
   if hp.lelec != 0.0: ########################## layer on-site
      LG.info('Electric field: %s'%(hp.lelec))
      Htot.append( electric(base,hp.lelec) )
   if  hp.lrashba != 0.0: ####################### Rashba
      LG.info('Rashba coupling')
      Htot.append( pseudo_rashba(base,hp.lrashba) )
   LG.info('Hamiltonian %s ready'%(tag))
   H_pris = Hamiltonian(Htot,tag=tag,dospin=dospin)
   H_pris.names()
   del Htot
   return H_pris



#def dic2vec(d):
#   """
#     Given a dictionary with some Slater-Koster parameters:
#       {'Vpps': 7.48, 'Vsss': -7.76, 'Vsps': 8.16, 'Vppp': -3.59}
#     Returns a suitable vector for the Slater_Koster function:
#      [Vsss, Vsps, Vpps, Vppp, Vsds, Vpds, Vpdp, Vdds, Vddp, Vddd]
#   """
#   nam = ['Vsss', 'Vsps', 'Vpps', 'Vppp', 'Vsds', 'Vpds', 'Vpdp', 'Vdds',
#          'Vddp', 'Vddd']
#   vec = [0.0 for _ in nam]
#   for i in range(len(nam)):
#      n = nam[i]
#      try: vec[i] = d[n]
#      except KeyError: pass  #XXX check
#   return vec



#@log_help.log2screen(LG)
def kinetic(base,hoppings,func=None,coup=1):
   """
   TODO: dont run over the basis elements, instead run over the bonds
     Converts a matrix of neighbours in a mtrix of hoppings. NOTICE that
     the input and output size may differ
   """
   import newSK as SK
   base.get_indices()
   ndim = len(base.INDS)
   diag_onsite = []
   for i in range(len(base.INDS)):
      e = base.INDS[i], base.ATS[base.AUX_INDS[i]], base.ORBS[i]
      diag_onsite.append(base[e[0]].onsite[e[2]])
   LG.debug('Expecting to create a %sx%s Hamiltonian'%(ndim,ndim))
   try: bonds = base.bonds
   except AttributeError: bonds = base.get_neig()
   if np.linalg.norm(bonds[0][1]) != 0:
      LG.warning('Incorrect order of vectors (0, a1, a2, a1+a2..)')
   #names = [b[-1] for b in bonds]
   iname = 0
   Htot = []
   for ib in range(len(bonds)):
      #nam = names[ib]
      M,v,nam = bonds[ib]
      LG.info('Doing matrix: %s'%(nam))
      #II,JJ,DD = M.row,M.col,M.data
      if np.linalg.norm(v) == 0:  ## Add on-site energies
         II = list(range(len(diag_onsite)))
         JJ = list(range(len(diag_onsite)))
         DD = [e for e in diag_onsite]
         #II = np.append(II,[i for i in range(len(diag_onsite))])
         #JJ = np.append(JJ,[i for i in range(len(diag_onsite))])
         #DD = np.append(DD,[e for e in diag_onsite])
      else: II,JJ,DD = [],[],[]   # H_aux = np.zeros((ndim,ndim))
      auxII,auxJJ,auxDD = [],[],[]
      for i,j in zip(M.row,M.col):  # i,j label connected atoms
         at1 = base[i]  # we build here the Hamiltonian termS between these 2
         at2 = base[j]  # atoms
         r = at1.position - (at2.position+v)
         for ih in at1.indices:     # actually we only need half of these
            #ei = base.INDS[ih], base.AUX_INDS[ih], base.ORBS[ih]
            for jh in at2.indices:  # elements. It has to be symmetric
               #ej = base.INDS[jh], base.AUX_INDS[jh], base.ORBS[jh]
               bra = base.INDS[ih], base.ATS[base.AUX_INDS[ih]], base.ORBS[ih]
               ket = base.INDS[jh], base.ATS[base.AUX_INDS[jh]], base.ORBS[jh]
               #bra = base.basis[ih]
               #ket = base.basis[jh]
               ab = sorted([bra[1],ket[1]])
               hop_ato = '%s-%s'%(ab[0],ab[1])
               hop_orb = 't_%s_%s'%(bra[2],ket[2])
               try: SKp = hoppings[hop_ato][1]
               except KeyError:
                  LG.debug('No %s hopping'%(hop_ato))
                  continue
               ii = bra[0]  #base.basis[ih][0]
               jj = ket[0]  #base.basis[jh][0]
               #ii = base.basis[ih][0]
               #jj = base.basis[jh][0]
               t = SK.hoppings[hop_orb]
               if base[jj].layer != base[ii].layer:
                  # TODO generalize condition, or modify input
                  if base[jj].element == 'C' and base[ii].element == 'H':
                     f = 0.
                  elif base[jj].element == 'H' and base[ii].element == 'C':
                     f = 0.
                  else: f = hoppings['Interlayer']
               else: f = 1.
               auxII.append(ih)
               auxJJ.append(jh)
               auxDD.append( f*t(r,SKp) )
               #if hop_ato == 'C-H' or hop_ato == 'H-C':
               #   print(hop_ato,hop_orb,f,r,f*t(r,SKp))
               #print(hop_ato,hop_orb,f,r,f*t(r,SKp))
      II = np.append(II,auxII)
      JJ = np.append(JJ,auxJJ)
      DD = np.append(DD,auxDD)
      LG.info('  ...added %s term'%(nam))
      H_aux = csc_matrix( (DD, (II, JJ)), shape=(ndim,ndim) )
      H_aux.eliminate_zeros()
      Htot.append( HTerm(H_aux,v,coup,name=nam) )
   return Htot


def mass(base,lmass):
   """ Assumes sublattice attribute already calculated """
   v = np.array([0.,0.,0.])
   sub_dic = {'A':1,'B':-1}
   aux = [[None for _ in base.elements] for _ in base.elements]
   for i in range(len(base.elements)):
      E = base.elements[i]
      f = E.sublattice
      aux[i][i] = coo_matrix(f * np.identity(len(E.onsite)))
   return HTerm(csc_matrix(bmat(aux)),v,lmass,name='mass')


def soc(base,lso):
   """ Returns the SOC Hterm """
   v = np.array([0.,0.,0.])
   from SOC import soc_l
   #L_s = soc_l(0)
   #L_p = soc_l(1)
   #L_d = soc_l(2)
   #orb_mat = {'s':L_s,'p':L_p,'d':L_d}
   #s_ord = {'s':0}
   #p_ord = {'px':0,'py':1,'pz':2}
   #d_ord = {'dxy':0,'dyz':1,'dzx':2,'dx2y2':3,'d3r2y2':4}
   #ord_orb = {'s':s_ord, 'p':p_ord, 'd':d_ord}
   #print(L_s)
   #print(L_p)
   #aux = [] # will contain the blocks of the SOC matrix
   #for ei in base.ORBS:
   #   print(ei)
   #   Mi = orb_mat[ei[0]]
   #   dici = ord_orb[ei[0]]
   #   print(Mi)
   #   indi = dici[ei]
   #   print(M[ind])
   #exit()
   aux = [[None for _ in base.elements] for _ in base.elements] # Size=num atoms
   base.DOspin = True
   #SOC_s = np.matrix([[0,0],[0,0]])
   #SOC_p = 0.5*np.matrix([[0, 0,-1j, 0, 0,-1],
   #[0, 0, 0, 1j, 1, 0],
   #[1j, 0, 0, 0, 0, 1j],
   #[0, -1j, 0, 0, 1j, 0],
   #[0, 1, 0, -1j, 0, 0],
   #[-1, 0, -1j, 0, 0, 0]])
   #l_orb = {'s':SOC_s,'p':SOC_p}
   l_orb = {'s':0,'p':1,'d':2}
   for i in range(len(base.elements)):
      E = base.elements[i]
      orbs = []
      for x in E.orbitals:
         if x[0] not in orbs: orbs.append(x[0])
      auxx = [[None for _ in orbs] for _ in orbs]
      for j in range(len(auxx)):
         auxx[j][j] = soc_l(l_orb[orbs[j]])
      #auxx = [[None for _ in E.orbitals] for _ in E.orbitals]
      #for j in range(len(E.orbitals)):
      #   o = E.orbitals[j]
      #   auxx[j][j] = coo_matrix(soc_l(l_orb[o]))
      aux[i][i] = bmat(auxx)
   return HTerm(csc_matrix(bmat(aux)),v,lso,name='soc')


def sparse_pauli_matrix(n,tp=complex):
   """
     Sparse definition of N-range Pauli matrices
   """
   c = [2*i+1 for i in range(n)]
   r = [2*i for i in range(n)]
   row = np.array(r+c)
   col = np.array(c+r)
   data = np.array([1 for _ in row])
   Sx = coo_matrix((data, (row, col)), shape=(2*n, 2*n), dtype=tp)

   c = [2*i+1 for i in range(n)]
   r = [2*i for i in range(n)]
   row = np.array(r+c)
   col = np.array(c+r)
   data = np.array([1j*(-1)**(i+1) for i in row])
   Sy = coo_matrix((data, (row, col)), shape=(2*n, 2*n), dtype=tp)

   c = [i for i in range(2*n)]
   r = [i for i in range(2*n)]
   row = np.array(r)
   col = np.array(c)
   data = np.array([(-1)**(i) for i in row])
   Sz = coo_matrix((data, (row, col)), shape=(2*n, 2*n), dtype=tp)
   return Sx,Sy,Sz


def zeeman(base,lzee):
   """ Returns the Zeeman term """
   def pauli_matrix(n):
      """
        Builds 3 nxn matrices with the each of the Pauli Matrices in the
        diagonal
      """
      Sx = np.matrix([[0,1],[1,0]],dtype=complex)
      Sy = np.matrix([[0,-1j],[1j,0]],dtype=complex)
      Sz = np.matrix([[1,0],[0,-1]],dtype=complex)
      auxX = [[None for _ in range(n)] for _ in range(n)]
      auxY = [[None for _ in range(n)] for _ in range(n)]
      auxZ = [[None for _ in range(n)] for _ in range(n)]
      for i in range(n):
         auxX[i][i] = Sx
         auxY[i][i] = Sy
         auxZ[i][i] = Sz
      return bmat(auxX), bmat(auxY), bmat(auxZ)
   base.DOspin = True
   v = np.array([0.,0.,0.])
   coup = 1.0
   N = len(base.ORBS)
   sig = pauli_matrix(N)
   M = lzee[0]*sig[0] + lzee[1]*sig[1] + lzee[2]*sig[2]
   return HTerm(csc_matrix(M),v,coup,name='zeeman')


@log_help.log2screen(LG)
def electric(base,lElec):
   """ Returns the Electric Field term """
   #XXX may fail for multiorbital
   #XXX will use the layer number as a factor for the electric field.
   #    This may result in unexpected results for multilayers CHECK
   LG.info('Doing matrix for electric field. lelec=%s'%(lElec))
   v = np.array([0.,0.,0.])
   ndim = len(base.INDS)
   II = list(range(len(base.LAYS)))
   JJ = list(range(len(base.LAYS)))
   H_aux = csc_matrix( (base.LAYS, (II, JJ)), shape=(ndim,ndim) )
   H_aux.eliminate_zeros()
   LG.info('... added electric field')
   return HTerm(H_aux,v,lElec,name='electric')
   #return HTerm(csc_matrix(bmat(aux)),v,lElec,name='electric')


@log_help.log2screen(LG)
def pseudo_rashba(base,lElec):
   ## TODO use array slices to select the indices rather than running through
   # all the orbitals
   LG.info('Doing matrix for Rashba')
   v = np.array([0.,0.,0.])
   ndim = len(base.INDS)
   #M = np.matrix(np.zeros((ndim,ndim)))
   inds = list(range(len(base.ats)))
   II,JJ = [],[]
   for i_at in inds:   # Rashba is only intra-atomic
      v_orbs = base.ORBS[base.INDS==i_at]
      v_inds = base.AUX_INDS[base.INDS==i_at]
      s = v_orbs=='s'
      pz = v_orbs=='pz'
      inds_s = v_inds[s]
      inds_pz = v_inds[pz]
      if len(inds_s) > 0 and len(inds_pz) > 0:
         for i in inds_s:
            for j in inds_pz:
               II.append(i)
               JJ.append(j)
               II.append(j)
               JJ.append(i)
   data = [1 for _ in II]
   H_aux = csc_matrix( (data,(II,JJ)), shape=(ndim,ndim))
   #II,JJ = [],[]
   #for i in range(len(base.ORBS)):
   #   it = base.ORBS[i]
   #   for j in range(len(base.ORBS)):
   #      jt = base.ORBS[j]
   #      if (it=='s' and jt=='pz') or (it=='pz' and jt=='s'):
   #         #M[i,j] = 1
   #         II.append(i)
   #         JJ.append(j)
   ##   print(i,'/',len(base.ORBS),'  (%s)'%(len(II)))
   #data = [1 for _ in II]
   #H_aux = csc_matrix( (data,(II,JJ)), shape=(ndim,ndim))
   #LG.info('... added Rashba term')
   ##print('--------------')
   return HTerm(H_aux,v,lElec,name='rashba')


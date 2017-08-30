#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import logging
import log_help
LG = logging.getLogger(__name__)

import mygreen
import numpy as np
import algebra as alg


def green_function(e,v,h,delta=0.01,path_selfes='.',force=False,l=1.0):
   """
     Returns the Green's function of a given hamiltonian (v) with self-energy
     (s) and broadening (delta)
   ** Note that v is a matrix (on-site) but h is a hamiltonian class
   """
   e = round(e,8)
   LG = logging.getLogger('mygreen_tools')
   emat = np.matrix(np.identity(v.shape[0]),dtype=complex) * (e + delta*1j)
   if l == 0.:
      LG.debug('Self-Energy coupling = 0')
      #print('Self-Energy coupling = 0')
      return (emat - v ).I
   else:
      sname = path_selfes + 'self%s.npy'%(e)
      try:
         if force: raise FileNotFoundError
         s = np.load( sname )
         if v.shape != s.shape: s = alg.m2spin(s)
         LG.info('loaded from file '+sname)
         #print('loaded from file '+sname)
      except FileNotFoundError:
         LG.debug('Calculating for E=%s'%(e))
         print('Calculating for E=%s'%(e))
         #_,s = calc_selfe(e,h,delta=delta,path=path_selfes)
         _,s = mygreen.bloch_selfenergy(h,energy=e,delta=delta,\
                                       nk=100,error=10**(-5),\
                                       #mode="full")
                                       #mode="renormalization")
                                       mode="adaptative")
         if path_selfes != None:
            sname = path_selfes + 'self%s.npy'%(e)
            np.save(sname,s)
            LG.debug('Saved: %s'%(sname))
      return (emat - v - l*s ).I


#def calc_selfe(e,h,delta=0.01,path=None):
#   g,selfe = mygreen.bloch_selfenergy(h,energy=e,delta=delta,\
#                                    nk=100,error=10**(-5),\
#                                    mode="adaptative")
#                                    #mode="full")
#                                    #mode="renormalization")
#   if path != None:
#      sname = path+'self%s.npy'%(e)
#      np.save(sname,selfe)
#      #mygreen.Selfenergy(e,g,selfe).write(path)
#      LG.info('Writing: '+sname)
#   return g,selfe




#def integ_brute(vintra,selfes,delta=0.01):
#   def phase(C):
#      return cmath.phase(C.energy+15)
#   # get Epath (ordered)
#   recorrido = sorted(selfes,key=phase,reverse=True)
#   integ = 0.0*1j
#   for i in range(1,len(recorrido)):
#      e1 = recorrido[i-1].energy
#      e = recorrido[i].energy
#      s = recorrido[i].selfe
#      integ += green_function(e,vintra,s,delta=delta)*(e-e1)
#   #integ = integ/len(selfes)
#   return 0,integ
#
#
#def integr_simpson(vintra,selfes,delta=0.01):
#   """
#      Calculates the integral over the path selfes.
#    Selfes must be a list of objects with the attributes energy, green, selfe
#   """
#   #def f(e,v,s):
#   #   emat = np.matrix(np.identity(len(v))) * (e + delta*1j)
#   #   return (emat - v - s ).I
#   def phase(C):
#      return cmath.phase(C+15)
#   # get Epath (ordered)
#   recorrido = []
#   for SLF in selfes:
#      recorrido.append(SLF.energy)
#   recorrido = list(set(recorrido))
#   recorrido = sorted(recorrido,key=phase,reverse=True)
#
#   # Defines the points where the func is evaluated.
#   eutil,selfeutil = [],[]
#   for i in range(0,len(recorrido),2):
#      Eaux = recorrido[i]
#      for SLF in selfes:
#         if SLF.energy == Eaux and Eaux not in eutil:
#            selfeutil.append(SLF.selfe)
#            eutil.append(Eaux)
#            continue   #XXX break ??
#   #eutil = [ recorrido[i] for i in range(0,len(recorrido),2) ]
#   #selfeutil = [ selfes[i] for i in range(0,len(selfes),2) ]
#
#   # aux points are middle pts in the simpson rule
#   eaux,selfeaux = [],[]
#   for i in range(1,len(recorrido)-1,2):
#      Eaux = recorrido[i]
#      for SLF in selfes:
#         if SLF.energy == Eaux and Eaux not in eutil:
#            selfeaux.append(SLF.selfe)
#            eaux.append(Eaux)
#            continue
#   #eaux = [ recorrido[i] for i in range(1,len(recorrido)-1,2) ]
#   #selfeaux = [ selfes[i] for i in range(1,len(selfes)-1,2) ]
#
#   gvs = []
#   integ = 0.0*1j
#   # TODO This loop can be parallel ?
#   for ix in range(len(eutil)-1):
#      a = eutil[ix]
#      b = eutil[ix+1]
#      c = eaux[ix]
#      selfea = selfeutil[ix]
#      selfeb = selfeutil[ix+1]
#      selfec = selfeaux[ix]
#
#      fa = green_function(a,vintra,selfea)  # XXX this should be green func
#      fb = green_function(b,vintra,selfeb)  # XXX this should be green func
#      fc = green_function(c,vintra,selfec)  # XXX this should be green func
#      integ += (fa+4*fc+fb) * (b-a)/6.  #(f(a) + 4*f(c) + f(b)) * (b-a)/6.
#      gvs.append(fa)
#   gvs.append(fb)
#   return gvs,integ
#
#
#def calc_selfes(h,es,delta=0.01,ncpus=4,path=None):
#   import pp
#   """ Parallel calculation of the self-energies """
#   # tuple of all parallel python servers to connect with
#   ppservers = ()
#   #ppservers = ("192.168.8.58:1234","192.168.8.49:1234",)
#   server = pp.Server(ncpus, ppservers=ppservers)
#   jobs =  []
#   for i in es:
#      jobs.append((i,
#                   server.submit(blafunc, (i,h,delta,path), (),
#                   ('import mygreen','numpy as np',))))
#   #res = [ (e,job()) for e, job in jobs]
#   def readncpus(fname='NCPUS'):
#      try: return int(open(fname,'r').read())
#      except IOError: return ncpus
#   res = []
#   for e, job in jobs:
#      res.append( (e,job()) )
#      if server.get_ncpus() != readncpus(): server.set_ncpus(readncpus())
#   selfesspinless, greensspinless = [],[]
#   resultado = []
#   for r in res:
#      E = r[0]
#      G = r[1][0]
#      S = r[1][1]
#      greensspinless.append(G) #r[1][0])
#      selfesspinless.append(S) #r[1][1])
#      resultado.append(Selfenergy(E,G,S))
#   server.print_stats()
#   return resultado
#   #return selfesspinless, greensspinless
#
#
#def get_selfes(energies,path_selfes,h,delta=0.001,ncpus=4,force=False):
#   """
#      Read or calculate the Self energies.
#      force forces the calculation of all the self-energies
#   """
#   check = check_file   # local rename of function
#   ## check if already calculated
#   calculated,missing = [],[]
#   for e in energies:
#      gname = 'green%s.npy'%(e)
#      sname = 'self%s.npy'%(e)
#      if check(path_selfes+gname) and check(path_selfes+sname):
#         G = np.load(path_selfes+gname)
#         S = np.load(path_selfes+sname)
#         A = Selfenergy(e,G,S)  # object with E,self,green atributes
#         if not force: calculated.append(A)
#         else: missing.append(e)
#      else:
#         missing.append(e)
#   ## calculate (and write) missing matrices
#   if len(missing) > 0:
#      print('==>',len(missing))
#      print('Calculating %s missing energies'%(len(missing)))
#      rest = calc_selfes(h,missing,delta=delta,ncpus=ncpus,path=path_selfes)
#   else: rest = []
#   #for OBJ in rest:
#   #   OBJ.write(path_selfes)
#   #   print 'Writing:',path_selfes+'green%s.npy'%(OBJ.energy)
#   return calculated + rest
#
#
#def calc_DOS(base,hamil,list_selfes,fname_DOS='OUT/',delta=0.001,ncpus=4):
#   """
#     Returns the spin up and down DOS as well as the energies in which they
#     were calculated
#   """
#   es,dos = [],[]
#   dos_up,dos_down = [],[]
#   for Selfe in list_selfes:
#      E = Selfe.energy
#      SLF = Selfe.selfe
#      g = green_function(E,hamil,SLF,delta)
#      dia = np.diagonal(g)
#      ## Total
#      auxUP,auxDOWN = 0.0, 0.0
#      #TODO use indexing instead of loop
#      for i in range(g.shape[0]//2): # 0,len(g)-1,2):
#         auxUP += g[2*i,2*i]
#         auxDOWN += g[2*i+1,2*i+1]
#      es.append(E.real)
#      aux = -g.trace()[0,0].imag/np.pi
#      if abs(-(auxUP + auxDOWN).imag/np.pi - aux) > 0.000001 :
#         exit('-*-'*10)
#      dos.append(-g.trace()[0,0].imag/np.pi)
#      dos_up.append(-auxUP.imag/np.pi)
#      dos_down.append(-auxDOWN.imag/np.pi)
#   f = open(fname_DOS,'w')
#   for i in range(len(es)):
#      f.write(str(es[i])+'   '+str(dos_up[i])+'   '+str(dos_down[i])+'\n')
#   f.close()
#   return es,dos_up,dos_down,dos
#
#def check_file(path):
#   """" Checks if a given file exists or not """
#   if isinstance(path,str):
#      if os.path.isfile(path):
#         return True
#      else:
#         return False
#
#def check_folder(fol):
#   """ Creates a given folder if it didn't exist """
#   if isinstance(fol,str):
#      if os.path.isdir(fol): pass
#      else: os.system('mkdir -p %s'%(fol))
#
#def read_hamil(root='./',norm=True):
#   """
#     Read the Hamiltonian terms intra,x,y,xy,xmy.
#     The norm options makes ALL hoppings equal to 1.0
#   """
#   files = ['intra.dat','x.dat','y.dat','xy.dat','xmy.dat']
#   ## check dimension
#   existing_files = []
#   for fi in files:
#      if check_file(root+fi): existing_files.append(fi)
#   dim_aux = 0
#   for i in range(len(existing_files)):
#      fi = existing_files[i]
#      f = open(root+fi,'r')
#      ndim = f.readline().lstrip().rstrip().split()[1].lstrip().rstrip()
#      ndim = int(ndim.replace('(','').replace(',',''))
#      f.close()
#      if i > 0:
#         if ndim != dim_aux:
#            sys.exit('%s and %s have different dimension!'%(existing_files[i],\
#                                                          existing_files[i-1]))
#      dim_aux = ndim
#   if len(files) != len(existing_files):
#      print('WARNING: Missing hoppings. Using 0-matrices for them.')
#
#   ## Initialize all matrices to 0
#   mats = [np.matrix(np.zeros((ndim,ndim),dtype=complex)) for _ in files]
#   for im in range(len(files)):
#      fi = files[im]
#      if fi in existing_files:
#         I,J,D = [],[],[] #row,col,data
#         f = open(root+fi,'r')
#         lines = f.readlines()
#         for l in lines:
#            line = l.lstrip().rstrip()
#            if line[0] != '#':
#               li = line.rstrip('\n').split()
#               I.append(int(li[0]))
#               J.append(int(li[1]))
#               D.append(float(li[2]))
#         f.close()
#         aux = mats[im]
#         for i,j,d in zip(I,J,D):
#            if norm:
#               aux[i,j] = d/d
#            else:
#               aux[i,j] = d
#   return mats
#
##def get_distances(f,es,path,h,delta=0.001,ncpus=4):
##   """
##    Returns de list of distances between all the points in the mesh.
##    The distance is defined as the difference of a given function in each point
##        dists[i] = |e[i+1]-e[i]| * sum( |f[i] - f[i+1| ) / len(f)
##    The element dists[i] correndond to the difference of the
##    selfenergies between es[i] and es[i+1]
##   """
##   selfenergies = f(es,path,h,delta=delta,ncpus=ncpus,force=False)
##   dists = [0 for _ in range(len(es)-1)]
##   ## compruebo si hacen falta más
##   for ie in range(len(es)-1):
##      e = es[ie]
##      e1 = es[ie+1]
##      selfe = selfenergies[ie].selfe
##      selfe1 = selfenergies[ie+1].selfe
##      dists[ie] = np.abs(e1-e)*np.sum(np.abs(selfe-selfe1))/len(selfe)
##   return dists
##
##def adjust_mesh(es,fu,path,h,thres=0.00001,delta=0.001,ncpus=4):
##   """
##     Refine the mesh based only in the function. NOT THE INTEGRAL.
##   """
##   cont = 0
##   out = False
##   while not out:
##      #LG.debug('Adjust mesh: %s'%(cont))
##      print('Adjust mesh: %s'%(cont))
##      es_new = []
##      dists = get_distances(fu,es,path,h,delta=delta,ncpus=ncpus)
##      for ie in range(len(es)-1):
##         if dists[ie] < thres:
##            es_new.append(es[ie])
##         else:
##            es_new.append(es[ie])
##            es_new.append((es[ie]+es[ie+1])/2.)
##      es_new.append(es[ie+1]) #because last point would never be included
##
##      if len(es) == len(es_new): out = True   # No new points
##      else: out = False  # New points
##      es = deepcopy(es_new)  # just in case
##      cont += 1
##   return es
#
#
##def energy_path(func,path,h,a=-30,b=0,c=0,d=30,limDOS=3,Ne=100,
##                                           delta=0.001,thres=0.00001,ncpus=4):
##   """
##     Given a path  it returns another one in wich the function to be evaluated
##   does not vary more than thres between neighbouring points.
##     The returning mesh ensures that the Riemman sumation for the integral
##   should be enough (yet we will use Simpson)
##     func: function to adjust for integration
##     path: path to folder containing the self-energies
##   """
##   Ne_start = 10
##   es1 = np.linspace(a+c*1j, a+d*1j,Ne_start*2)  # Because the simpson rule
##   es2 = np.linspace(a+d*1j, b+d*1j,Ne_start*2)  # uses f(a) + f(a+b/2) + f(b)
##   es3 = np.linspace(b+d*1j, b+c*1j,Ne_start)
##   es = np.append(es1,es2)
##   es = np.append(es,es3)
##   ## remove duplicates
##   es_bien = []
##   for c in es:
##      if c not in es_bien: es_bien.append(c)
##   es = es_bien
##   es_new = adjust_mesh(es,func,path,h,thres,delta=delta,ncpus=ncpus)
##   ## remove duplicates
##   es = [ c for c in es_new if c not in es ]
##   #es = []
##   #for c in es_new:
##   #   if c not in es: es.append(c)
##   try: ini,fin = complex(limDOS[0]), complex(limDOS[1])
##   except TypeError: ini,fin = complex(-limDOS), complex(limDOS)
##   # Include an extra point trying to include 0.0
##   esR = list(np.linspace(complex(ini), complex(fin), Ne+1-Ne%2))
##   return es,esR
#
#
#def check_selfe(es,path):
#   onlyfiles = [ f for f in listdir(path) if isfile(join(path,f)) ]
#   ES = []
#   for f in onlyfiles:
#      if '.npy' in f:
#         aux = f.replace('.npy','').replace('green','').replace('self','')
#         aux = complex(aux)
#      else: pass
#      if aux not in ES: ES.append(aux)
#   ESset = set(ES)
#   esset = set(es)
#   print('total of %s Energies'%(len(ES)))
#   print('Getting Selfes')
#   exit()
#   return True
#   
#
#def read_selfes(self_folder):
#   fg = open(self_folder+'/greens.txt')
#   fs = open(self_folder+'/selfs.txt')
#   greens,selfs = [],[]
#   for g,s in zip(fg.readlines(),fs.readlines()):
#      greens.append(g.rstrip()+'.npy')
#      selfs.append(s.rstrip()+'.npy')
#   fg.close()
#   fs.close()
#   es,selfesspinless,greensspinless = [],[],[]
#   for g,s in zip(greens,selfs):
#      Eg = complex(g.split('/')[-1].replace('green','').replace('.npy',''))
#      Es = complex(s.split('/')[-1].replace('self','').replace('.npy','') )
#      if Eg != Es: sys.exit('Problem reading self/green files')
#      es.append(Eg)
#      G = np.load(g)
#      S = np.load(s)
#      greensspinless.append(np.matrix(G))
#      selfesspinless.append(np.matrix(S))
#   return es,selfesspinless,greensspinless
#
#def write_mag(integ,base,it=0,fname=None):
#   pos = []
#   for P in base.elems:
#      pos.append(P.position)
#
#   mags = []
#   for i in range(len(base.elems)):
#      nup = -integ[i,i].imag/np.pi
#      ndown = -integ[i+1,i+1].imag/np.pi
#      mags.append(np.array([0,0,(nup-ndown)/2.]))  # Only Sz magnetization!!!
#
#   if fname:
#      print('Final: %s'%(fname))
#      f = open('%s'%(fname),'w')
#   else:
#      print('Final: magnetization%s.txt'%(it))
#      f = open('magnetization%s.txt'%(it),'w')
#   for P,M in zip(pos,mags):
#      if M[0].imag > 0.01: print('M0 imag!!!')
#      if M[1].imag > 0.01: print('M1 imag!!!')
#      if M[2].imag > 0.01: print('M2 imag!!!')
#      f.write(str(P[0])+'   '+str(P[1])+'   '+str(P[2])+'   '+
#              str(M[0].real)+'   '+str(M[1].real)+'   '+str(M[2].real)+'\n')
#   f.close()
#
#
#def write_selfes(es,selfesspinless,greensspinless,self_folder):
#   print('Salvar en folder:',self_folder,'\n')
#   fg = open('%s/greens.txt'%(self_folder),'w')
#   fs = open('%s/selfs.txt'%(self_folder),'w')
#   for E,S,G in zip(es,selfesspinless,greensspinless):
#      ## self
#      fname = self_folder+'/self'+str(E)
#      fs.write(fname+'\n')
#      np.save(fname,S)
#      ## green
#      fname = self_folder+'/green'+str(E)
#      fg.write(fname+'\n')
#      np.save(fname,G)
#   fg.close()
#   fs.close()
#
##def mail(msg,to='noel.alberto.garcia@gmail.com',subject='DONE'):
##   os.system('echo %s | mutt %s -s %s'%(msg,to,subject))
##
##
##def send_mail(to,body,subject='',attach=None):
##   command = 'echo "%s" | mutt %s '%(body,to)
##   if len(subject) > 0:
##      command += '-s %s '%(subject)
##   else: pass
##   if attach != None:
##      if isinstance(attach,list):
##         for fil in attach:
##            command += '-a %s '%(fil)
##      else:
##         command += '-a %s '%(attach)
##   os.system(command)

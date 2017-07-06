#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np

"""
   This library provides all the Slater-Koster hoppings as detailed in the paper
 (Phys. Rev. 94, 1498 (1954)).
 ** NOTE: the d orbitals have not been tested (21/12/2016)
"""

def vec(r):
   rn = np.linalg.norm(r)
   l = r[0]/rn
   m = r[1]/rn
   n = r[2]/rn
   return l,m,n

def dic2vec(d):
   """
     Given a dictionary with some Slater-Koster parameters:
         {'Vpps': 7.48, 'Vsss': -7.76, 'Vsps': 8.16, 'Vppp': -3.59}
     Returns a suitable vector for the Slater_Koster function:
        [Vsss, Vsps, Vpps, Vppp, Vsds, Vpds, Vpdp, Vdds, Vddp, Vddd]
   """
   nam = ['Vsss', 'Vsps', 'Vpps', 'Vppp', 'Vsds', 'Vpds', 'Vpdp', 'Vdds',
                                                                'Vddp', 'Vddd']
   vec = [0.0 for _ in nam]
   for i in range(len(nam)):
      n = nam[i]
      try: vec[i] = d[n]
      except KeyError: pass  #XXX check
   return vec


# s-s
def t_s_s(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return Vsss

# p-p
def t_px_px(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (l**2)*Vpps+(1-l**2)*Vppp

def t_py_py(r,SKp):
   l,m,n = vec(r)
   return t_px_px((m,-l,n),SKp)

def t_pz_pz(r,SKp):
   l,m,n = vec(r)
   return t_px_px((n,m,-l),SKp)

def t_px_py(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return l*m*(Vpps-Vppp)
def t_py_px(r,SKp): return t_px_py(r,SKp)

def t_px_pz(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return l*n*(Vpps-Vppp)
def t_pz_px(r,SKp): return t_px_pz(r,SKp)
def t_pz_py(r,SKp): return t_py_pz(r,SKp)

def t_py_pz(r,SKp):
   l,m,n = vec(r)
   return t_px_py((m,n,l),SKp)

# d-d  # order: xy, yz, zx, x2y2, 3z2r2
def t_dxy_dxy(r,SKp):  # 0,0
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return 3.*(l**2)*(m**2)*Vdds + (l**2+m**2-4.*(l**2)*(m**2))*Vddp +\
                                        (n**2+((l**2)*(m**2)))*Vddd
def t_dyz_dyz(r,SKp):  # 1,1
   l,m,n = vec(r)
   return t_dxy_dxy((n,m,-l),SKp)
def t_dzx_dzx(r,SKp):  # 2,2
   l,m,n = vec(r)
   return t_dxy_dxy((l,n,-m),SKp)

def t_dx2y2_dx2y2(r,SKp):  # 3,3
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (3./4.)*((l**2-m**2)**2)*Vdds + (l**2+m**2-(l**2-m**2)**2)*Vddp +\
                                       (n**2+(1./4.)*(l**2-m**2)**2)*Vddd

def t_d3z2r2_d3z2r2(r,SKp):  # 4,4
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return ((n**2-(1./2.)*(l**2+m**2))**2)*Vdds + 3.*(n**2)*(l**2+m**2)*Vddp +\
                                              (3./4.)*((l**2+m**2)**2)*Vddd
def t_dxy_dyz(r,SKp):  # 0,1
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return 3.*l*(m**2)*n*Vdds+l*n*(1.-4.*m**2)*Vddp+l*n*(m**2-1.)*Vddd

def t_dxy_dzx(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return 3.*(l**2)*m*n*Vdds+m*n*(1.-4.*l**2)*Vddp+m*n*(l**2-1.)*Vddd

def t_dxy_dx2y2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (3./2.)*l*m*(l**2-m**2)*Vdds + 2*l*m*(m**2-l**2)*Vddp +\
                                   (1./2.)*l*m*(l**2-m**2)*Vddd

def t_dxy_d3z2r2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return np.sqrt(3)*l*m*(n**2-(1./2.)*(l**2+m**2))*Vdds - \
          2*np.sqrt(3)*l*m*(n**2)*Vddp+(1./2.)*np.sqrt(3)*l*m*(1+n**2)*Vddd
def t_dyz_dzx(r,SKp):
   l,m,n = vec(r)
   return t_dxy_dzx((-n,-m,-l),SKp)

def t_dyz_dx2y2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (3./2.)*m*n*(l**2-m**2)*Vdds - m*n*(1.+2.*(l**2-m**2))*Vddp +\
                                    m*n*(1.+(1./2.)*(l**2-m**2))*Vddd

def t_dyz_d3z2r2(r,SKp):   # 1,4
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return np.sqrt(3)*m*n*(n**2-(1./2.)*(l**2+m**2))*Vdds +\
          np.sqrt(3)*m*n*(l**2+m**2-n**2)*Vddp -\
          (1./2.)*np.sqrt(3)*m*n*(l**2+m**2)*Vddd

def t_dzx_dx2y2(r,SKp):   # 2,3
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (3./2.)*n*l*(l**2-m**2)*Vdds+n*l*(1.-2.*(l**2-m**2))*Vddp -\
                                  n*l*(1.-(1./2.)*(l**2-m**2))*Vddd

def t_dzx_d3z2r2(r,SKp):   # 2,4
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return np.sqrt(3)*l*n*(n**2-(1./2.)*(l**2+m**2))*Vdds +\
          np.sqrt(3)*l*n*(l**2+m**2-n**2)*Vddp -\
          (1./2.)*np.sqrt(3)*l*n*(l**2+m**2)*Vddd

def t_dx2y2_d3z2r2(r,SKp):   # 3,4
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (1./2.)*np.sqrt(3)*(l**2-m**2)*(n**2-(1./2.)*(l**2+m**2))*Vdds +\
          np.sqrt(3)*(n**2)*(m**2-l**2)*Vddp +\
          (1./4.)*np.sqrt(3)*(1+n**2)*(l**2-m**2)*Vddd

def t_dyz_dxy(r,SKp): return t_dxy_dyz(r,SKp)  # 1,0
def t_dzx_dxy(r,SKp): return t_dxy_dzx(r,SKp)  # 2,0
def t_dzx_dyz(r,SKp): return t_dyz_dzx(r,SKp)  # 2,1
def t_dx2y2_dxy(r,SKp): return t_dxy_dx2y2(r,SKp)  # 3,0
def t_dx2y2_dyz(r,SKp): return t_dyz_dx2y2(r,SKp)  # 3,1
def t_dx2y2_dzx(r,SKp): return t_dzx_dx2y2(r,SKp)  # 3,2
def t_d3z2r2_dxy(r,SKp): return t_dxy_d3z2r2(r,SKp)  # 4,0
def t_d3z2r2_dyz(r,SKp): return t_dyz_d3z2r2(r,SKp)  # 4,1
def t_d3z2r2_dzx(r,SKp): return t_dzx_d3z2r2(r,SKp)  # 4,2
def t_d3z2r2_dx2y2(r,SKp): return t_dx2y2_d3z2r2(r,SKp)  # 4,3

# s-p // p-s
def t_s_px(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return l*Vsps

def t_s_py(r,SKp):
   l,m,n = vec(r)
   return t_s_px((m,-l,n),SKp)
def t_s_pz(r,SKp):
   l,m,n = vec(r)
   return t_s_px((n,m,-l),SKp)

def t_px_s(r,SKp): return -t_s_px(r,SKp)   # CHECK
def t_py_s(r,SKp): return -t_s_py(r,SKp)   # CHECK
def t_pz_s(r,SKp): return -t_s_pz(r,SKp)   # CHECK

# s-d // d-s
def t_s_dxy(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return np.sqrt(3)*l*m*Vsds
def t_s_dyz(r,SKp):
   l,m,n = vec(r)
   return t_s_dxy((n,m,-l),SKp)
def t_s_dzx(r,SKp):
   l,m,n = vec(r)
   return t_s_dxy((l,n,-m),SKp)
def t_s_dx2y2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (1./2.)*np.sqrt(3)*(l**2-m**2)*Vsds
def t_s_d3z2r2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (n**2-(1./2.)*(l**2+m**2))*Vsds
def t_dxy_s(r,SKp): return t_s_dxy(r,SKp)   # CHECK
def t_dyz_s(r,SKp): return t_s_dyz(r,SKp)   # CHECK
def t_dzx_s(r,SKp): return t_s_dzx(r,SKp)   # CHECK
def t_dx2y2_s(r,SKp): return t_s_dx2y2(r,SKp)   # CHECK
def t_d3z2r2_s(r,SKp): return t_s_d3z2r2(r,SKp)   # CHECK

# p-d
def t_px_dxy(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return np.sqrt(3)*(l**2)*m*Vpds+m*(1-2*l**2)*Vpdp
def t_px_dyz(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return np.sqrt(3)*l*m*n*Vpds-2*l*m*n*Vpdp
def t_px_dzx(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return np.sqrt(3)*(l**2)*n*Vpds+n*(1-2*l**2)*Vpdp
def t_px_dx2y2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (1./2.)*np.sqrt(3)*l*(l**2-m**2)*Vpds+l*(1-l**2+m**2)*Vpdp
def t_px_d3z2r2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return l*(n**2-(1./2.)*(l**2+m**2))*Vpds-np.sqrt(3)*l*(n**2)*Vpdp
def t_py_dxy(r,SKp):
   l,m,n = vec(r)
   return t_px_dxy((m,l,-n),SKp)
def t_py_dyz(r,SKp):
   l,m,n = vec(r)
   return t_px_dzx((m,-l,n),SKp)
def t_py_dzx(r,SKp):
   l,m,n = vec(r)
   return t_px_dyz((m,n,l),SKp)
def t_py_dx2y2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (1./2.)*np.sqrt(3)*m*(l**2-m**2)*Vpds-m*(1+l**2-m**2)*Vpdp
def t_py_d3z2r2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return m*(n**2-(1./2.)*(l**2+m**2))*Vpds-np.sqrt(3)*m*(n**2)*Vpdp
def t_pz_dxy(r,SKp):
   l,m,n = vec(r)
   return t_px_dyz((n,-l,-m),SKp)
def t_pz_dyz(r,SKp):
   l,m,n = vec(r)
   return t_px_dxy((n,m,-l),SKp)
def t_pz_dzx(r,SKp):
   l,m,n = vec(r)
   return t_px_dzx((n,-m,l),SKp)
def t_pz_dx2y2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return (1./2.)*np.sqrt(3)*n*(l**2-m**2)*Vpds-n*(l**2-m**2)*Vpdp
def t_pz_d3z2r2(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return n*(n**2-(1./2.)*(l**2+m**2))*Vpds+np.sqrt(3)*n*(l**2+m**2)*Vpdp

def t_dxy_px(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -(np.sqrt(3)*(l**2)*m*Vpds+m*(1-2*l**2)*Vpdp)
def t_dyz_px(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -(np.sqrt(3)*l*m*n*Vpds-2*l*m*n*Vpdp)
def t_dzx_px(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -(np.sqrt(3)*(l**2)*n*Vpds+n*(1-2*l**2)*Vpdp)
def t_dx2y2_px(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -((1./2.)*np.sqrt(3)*l*(l**2-m**2)*Vpds+l*(1-l**2+m**2)*Vpdp)
def t_d3z2r2_px(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -(l*(n**2-(1./2.)*(l**2+m**2))*Vpds-np.sqrt(3)*l*(n**2)*Vpdp)
def t_dxy_py(r,SKp):
   l,m,n = vec(r)
   return t_px_dxy((m,l,-n),SKp)
def t_dyz_py(r,SKp):
   l,m,n = vec(r)
   return t_px_dzx((m,-l,n),SKp)
def t_dzx_py(r,SKp):
   l,m,n = vec(r)
   return t_px_dyz((m,n,l),SKp)
def t_dx2y2_py(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -((1./2.)*np.sqrt(3)*m*(l**2-m**2)*Vpds-m*(1+l**2-m**2)*Vpdp)
def t_d3z2r2_py(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -(m*(n**2-(1./2.)*(l**2+m**2))*Vpds-np.sqrt(3)*m*(n**2)*Vpdp)
def t_dxy_pz(r,SKp):
   l,m,n = vec(r)
   return t_px_dyz((n,-l,-m),SKp)
def t_dyz_pz(r,SKp):
   l,m,n = vec(r)
   return t_px_dxy((n,m,-l),SKp)
def t_dzx_pz(r,SKp):
   l,m,n = vec(r)
   return t_px_dzx((n,-m,l),SKp)
def t_dx2y2_pz(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -((1./2.)*np.sqrt(3)*n*(l**2-m**2)*Vpds-n*(l**2-m**2)*Vpdp)
def t_d3z2r2_pz(r,SKp):
   Vsss,Vsps,Vpps,Vppp,Vsds,Vpds,Vpdp,Vdds,Vddp,Vddd = dic2vec(SKp)
   l,m,n = vec(r)
   return -(n*(n**2-(1./2.)*(l**2+m**2))*Vpds+np.sqrt(3)*n*(l**2+m**2)*Vpdp)


hoppings = {'t_s_s':t_s_s,
't_px_px':t_px_px, 't_py_py':t_py_py, 't_pz_pz':t_pz_pz,
't_px_py':t_px_py, 't_py_px':t_py_px, 't_px_pz':t_px_pz, 't_pz_px':t_pz_px,
't_pz_py':t_pz_py, 't_py_pz':t_py_pz,
't_dxy_dxy':t_dxy_dxy, 't_dyz_dyz':t_dyz_dyz, 't_dzx_dzx':t_dzx_dzx,
't_dx2y2_dx2y2':t_dx2y2_dx2y2, 't_d3z2r2_d3z2r2':t_d3z2r2_d3z2r2,
't_dxy_dyz':t_dxy_dyz, 't_dxy_dzx':t_dxy_dzx, 't_dxy_dx2y2':t_dxy_dx2y2,
't_dxy_d3z2r2':t_dxy_d3z2r2, 't_dyz_dzx':t_dyz_dzx,
't_dyz_dx2y2':t_dyz_dx2y2, 't_dyz_d3z2r2':t_dyz_d3z2r2,
't_dzx_dx2y2':t_dzx_dx2y2, 't_dzx_d3z2r2':t_dzx_d3z2r2,
't_dx2y2_d3z2r2':t_dx2y2_d3z2r2, 't_dyz_dxy':t_dyz_dxy, 't_dzx_dxy':t_dzx_dxy,
't_dzx_dyz':t_dzx_dyz, 't_dx2y2_dxy':t_dx2y2_dxy, 't_dx2y2_dyz':t_dx2y2_dyz,
't_dx2y2_dzx':t_dx2y2_dzx, 't_d3z2r2_dxy':t_d3z2r2_dxy,
't_d3z2r2_dyz':t_d3z2r2_dyz, 't_d3z2r2_dzx':t_d3z2r2_dzx,
't_d3z2r2_dx2y2':t_d3z2r2_dx2y2,
't_s_px':t_s_px, 't_s_py':t_s_py, 't_s_pz':t_s_pz,
't_px_s':t_px_s, 't_py_s':t_py_s, 't_pz_s':t_pz_s,
't_s_dxy':t_s_dxy, 't_s_dyz':t_s_dyz, 't_s_dzx':t_s_dzx, 't_s_dx2y2':t_s_dx2y2,
't_s_d3z2r2':t_s_d3z2r2,
't_dxy_s':t_dxy_s, 't_dyz_s':t_dyz_s, 't_dzx_s':t_dzx_s, 't_dx2y2_s':t_dx2y2_s,
't_d3z2r2_s':t_d3z2r2_s,
't_px_dxy':t_px_dxy, 't_px_dyz':t_px_dyz, 't_px_dzx':t_px_dzx,
't_px_dx2y2':t_px_dx2y2, 't_px_d3z2r2':t_px_d3z2r2,
't_py_dxy':t_py_dxy, 't_py_dyz':t_py_dyz, 't_py_dzx':t_py_dzx,
't_py_dx2y2':t_py_dx2y2, 't_py_d3z2r2':t_py_d3z2r2,
't_pz_dxy':t_pz_dxy, 't_pz_dyz':t_pz_dyz, 't_pz_dzx':t_pz_dzx,
't_pz_dx2y2':t_pz_dx2y2, 't_pz_d3z2r2':t_pz_d3z2r2,
't_dxy_px':t_dxy_px, 't_dyz_px':t_dyz_px, 't_dzx_px':t_dzx_px,
't_dx2y2_px':t_dx2y2_px, 't_d3z2r2_px':t_d3z2r2_px, 't_dxy_py':t_dxy_py,
't_dyz_py':t_dyz_py, 't_dzx_py':t_dzx_py, 't_dx2y2_py':t_dx2y2_py,
't_d3z2r2_py':t_d3z2r2_py,
't_dxy_pz':t_dxy_pz, 't_dyz_pz':t_dyz_pz, 't_dzx_pz':t_dzx_pz,
't_dx2y2_pz':t_dx2y2_pz, 't_d3z2r2_pz':t_d3z2r2_pz}


if __name__ == '__main__':
   pos = [np.array([-0.7,0.,0.]), np.array([0.7,0.,0.])]
   #orbs = ['s','px','py','pz', 'dxy', 'dyz', 'dzx', 'dx2y2', 'd3z2r2']
   orbs = ['s','px','py','pz']
   from itertools import product
   combis = [p for p in product(orbs, repeat=2)]
   SK = {'Vpps': 7.48, 'Vsss': -7.76, 'Vsps': 8.16, 'Vppp': -3.59}
   SK = {'Vpps': 1, 'Vsss': 2, 'Vsps': 3, 'Vppp': 4}
   ## Complicao
   rs = [np.array([1.4,0.,0.]), np.array([0,1.4,0.]), np.array([0.,0.,1.4]),
         np.array([1.4,1.4,0.]),np.array([0.,1.4,1.4]),np.array([1.4,0.,1.4]),
         np.array([1.4,1.4,1.4])]
   rs = [np.array([1.4,0.,0.]), np.array([0,1.4,0.]), np.array([0.,0.,1.4])]
   for r in rs:
      print('===>',r)
      for c in combis:
         name = 't_'
         name += '_'.join(c)
         t = locals()[name](r,SK)
         if abs(t) > 0.0: print(name,'-->',t)
      print('')

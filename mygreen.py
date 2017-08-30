#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import numpy as np


def green_renormalization(intra,inter,energy=0.0,nite=None,
                                       error=0.0001,info=False,delta=0.001):
   """
     Calculates bulk and surface Green function by a renormalization
     algorithm, as described in I. Phys. F: Met. Phys. 15 (1985) 851-858
   """
   e = np.matrix(np.identity(intra.shape[0])) * (energy + 1j*delta)
   ite = 0
   alpha = inter
   beta = inter.H
   epsilon = intra
   epsilon_s = intra + inter * (e-intra).I * inter.H
   epsilon_s = intra
   while True: # implementation of Eq 11
      einv = (e - epsilon).I # inverse
      epsilon_s = epsilon_s + alpha * einv * beta
      epsilon = epsilon + alpha * einv * beta + beta * einv * alpha
      alpha = alpha * einv * alpha  # new alpha
      beta = beta * einv * beta  # new beta
      ite += 1
      # stop conditions
      if not nite == None:
         if ite > nite:  break
      else:
         if np.max(np.abs(alpha))<error and np.max(np.abs(beta))<error: break
   if info: print("Converged in ",ite,"iterations")
   g_surf = (e - epsilon_s).I # surface green function
   g_bulk = (e - epsilon).I  # bulk green function 
   return g_bulk,g_surf


def green_kchain(h,k=0.,energy=0.,delta=0.01,only_bulk=True,error=0.0001):
   """ Calculates the green function of a kdependent chain for a 2d system """
   def gr(ons,hop):
      """ Calculates G by renormalization"""
      gf,sf = green_renormalization(ons,hop,energy=energy,nite=None,
                                    error=error,info=False,delta=delta)
      if only_bulk:  return gf
      else:  return gf,sf
   tky = h.ty*np.exp(1j*np.pi*2.*k)
   #tkx = h.ty*np.exp(1j*np.pi*2.*k)
   tkxy = h.txy*np.exp(1j*np.pi*2.*k)
   tkxmy = h.txmy*np.exp(-1j*np.pi*2.*k)  # notice the minus sign !!!!
   # chain in the x direction
   ons = h.intra + tky + tky.H  # intra of k dependent chain
   hop = h.tx + tkxy + tkxmy  # hopping of k-dependent chain
   return gr(ons,hop)  # return green function


def bloch_selfenergy(h,nk=100,energy=0.0,delta=0.01,mode="full",error=0.00001):
   """
     Calculates the selfenergy of a cell coupled to a pristine environement.
     input is a hamiltonian class
   """
   def gr(ons,hop):
      """ Calculates G by renormalization"""
      gf,sf = green_renormalization(ons,hop,energy=energy,nite=None,
                             error=error,info=False,delta=delta)
      return gf,sf
   d = h.dimensionality # dimensionality of the system
   g = h.intra *0.0j # initialize green function
   e = np.matrix(np.identity(h.dim))*(energy + delta*1j) # complex energy
   if mode=="full":  # full integration
      if d==1: # one dimensional
         ks = np.linspace(0.,1.,nk,endpoint=False)
      elif d==2: # two dimensional
         ks = []
         kk = np.linspace(0.,1.,nk,endpoint=False)  # interval 0,1
         for ikx in kk:
            for iky in kk:
               ks.append([ikx,iky])
         ks = np.array(ks)  # all the kpoints
      else: raise # raise error
      hk_gen = h.get_hk_gen() #()  # generator of k dependent hamiltonian
      for k in ks:  # loop in BZ
         g += (e - hk_gen(k)).I  # add green function  
      g = g/len(ks)  # normalize
   #####################################################
   #####################################################
   elif mode=="renormalization":
      if d==1: # full renormalization
         g,s = gr(h.intra,h.inter)  # perform renormalization
      elif d==2: # two dimensional, loop over k's
         ks = np.linspace(0.,1.,nk,endpoint=False)
         for k in ks:  # loop over k in y direction
            #  add contribution to green function
            m = green_kchain(h,k=k,energy=energy,delta=delta,error=error)
            g += m
         g = g/len(ks)
   ######################################################
   ######################################################
   elif mode=="adaptative":
      if d==1: # full renormalization
         g,s = gr(h.intra,h.inter)  # perform renormalization
      elif d==2: # two dimensional, loop over k's
         #ks = np.linspace(0.,1.,nk,endpoint=False)
         import integration
         def fint(k):
            """ Function to integrate """
            return green_kchain(h,k=k,energy=energy,delta=delta,error=error)
         # eps is error, might work....
         g = integration.integrate_matrix(fint,xlim=[0.,1.],eps=error)
         ##XXX g = integration.integrate_matrix(fint,xlim=[0.,1.],eps=error)
          # chain in the y direction
          #ons = h.intra + tkx + tkx.H  # intra of k dependent chain
          #hop = h.ty + tkxy + tkxmy.H  # hopping of k-dependent chain
          #g += gr(ons,hop)  # add contribution to green function
      else: raise
   else: raise
   # now calculate selfenergy
   selfenergy = e - h.intra - g.I
   return g,selfenergy

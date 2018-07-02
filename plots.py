#!/usr/bin/python3
# -*- coding: UTF-8 -*-

"""
 This script should be in charge of walk through the files and decide what to
 plot. It will use then the library graphs.py to plot the results
"""

import numpy as np
import matplotlib.pyplot as plt
import graphs


def bands(folder):
   try:
      pris = np.load(folder + 'pris.bands.npy')
      dfct = np.load(folder + 'dfct.bands.npy')
   except FileNotFoundError:
      pris = np.loadtxt(folder + 'pris.bands')   
      dfct = np.loadtxt(folder + 'dfct.bands')

   Xp = pris[:,0]
   Yp = pris[:,1]
   Zp = pris[:,2]
   
   Xd = dfct[:,0]
   Yd = dfct[:,1]
   Zd = dfct[:,2]
   
   fig, ax = plt.subplots()
   graphs.bands(Xp,Yp,ax=ax,alpha=0.2)
   graphs.bands(Xd,Yd,ax=ax)
   plt.show()


def spectrum(folder):
   try:
      pris = np.load(folder + 'pris.spectrum.npy')
      dfct = np.load(folder + 'dfct.spectrum.npy')
   except FileNotFoundError:
      print('Files pris.spectrum.npy or dfct.spectrum.npy not found')
      exit()
   Ep = pris[:,0]
   Ed = dfct[:,0]
   fig, ax = plt.subplots()
   graphs.spectrum(Ep,ax=ax)
   graphs.spectrum(Ed,shift=1.1,ax=ax)
   plt.show()


if __name__ == '__main__':
   import argparse

   parser = argparse.ArgumentParser(description='Show TB results')

   help_msg = 'Plot bands'
   parser.add_argument('-b',action='store_true', help=help_msg)
   help_msg = 'Plot spectrum'
   parser.add_argument('-s',action='store_true', help=help_msg)
   help_msg = 'Path to results folder'
   parser.add_argument('-i',nargs='*',type=str, help=help_msg)


   args = parser.parse_args()

   if args.b==False and args.s==False:
      parser.print_help()

   for fol in args.i:
      print('')
      print(fol)
      print('')
      if args.b: bands(fol)
      if args.s: spectrum(fol)

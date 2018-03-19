#!/usr/bin/python3
# -*- coding: UTF-8 -*-

import logging

LlG = logging.getLogger('perform')
fh = logging.FileHandler('/tmp/performance.log',mode='w')
fmt = logging.Formatter('%(asctime)s %(name)s:%(levelname)s - %(message)s')
fh.setFormatter(fmt)
fh.setLevel(logging.DEBUG)
LlG.addHandler(fh)



def screen_handler(lg=None,lv='debug',fmt='%(name)s -%(levelname)s- %(message)s'):
   """
     This function adds a screenHandler to a given logger. If no logger is
     specified, just return the handler
   """
   if lv == 'debug': lv = logging.DEBUG
   elif lv == 'info': lv = logging.INFO
   elif lv == 'warning': lv = logging.WARNING
   elif lv == 'error': lv = logging.ERROR
   elif lv == 'critical': lv = logging.CRITICAL
   sh = logging.StreamHandler()
   sh.setLevel(lv)
   fmt = logging.Formatter(fmt)
   sh.setFormatter(fmt)
   if lg != None: lg.addHandler(sh)
   return sh


def log2screen(lg):
   """ This decorator adds *temporarily* a screenHandler to a given logger """
   def do_it(wrapped):
      def inner(*args, **kwargs):
         sh = screen_handler(lg)  # add ScreenHandler
         ret = wrapped(*args, **kwargs)
         lg.removeHandler(sh)     # remove ScreenHandler
         return ret
      return inner
   return do_it


def disable(lg):
   """
    Temporarily raise the log level to CRITICAL to avoid over logging
   """
   def do_it(wrapped):
      def inner(*args, **kwargs):
         lv = lg.getEffectiveLevel()
         lg.setLevel(logging.CRITICAL)
         ret = wrapped(*args, **kwargs)
         lg.setLevel(lv)
         return ret
      return inner
   return do_it

def disable2(lg):
   """
    Temporarily raise the log level to INFO to avoid over logging
   """
   def do_it(wrapped):
      def inner(*args, **kwargs):
         lv = lg.getEffectiveLevel()
         lg.setLevel(logging.INFO)
         ret = wrapped(*args, **kwargs)
         lg.setLevel(lv)
         return ret
      return inner
   return do_it

## Timer Decorator
from time import time
def timer(wrapped):
   """
   Logs the execution time of a certain funcion
   """
   def inner(*args, **kwargs):
      t = time()
      ret = wrapped(*args, **kwargs)
      LG.info(wrapped.__name__+' in %s'%(time()-t))
      return ret
   return inner


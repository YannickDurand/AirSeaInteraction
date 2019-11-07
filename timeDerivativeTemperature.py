#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 11:17:49 2019

@author: ydurand
"""

#%%
import numpy as np

from stockDiagKEs import *

#%%

def timeDerivativeTemperature(dataNetPE,t,areaGrid,adv=-1):
     dtdt, advrt, advxt, advyt= dataDiagPETemp(dataNetPE,t)
     
     dtdt = dtdt/86400.
     
     # t0 = np.where(t0 == 0,np.NaN,t0)
     # t1 = np.where(t1 == 0,np.NaN,t1)
     if adv==-1:
          
          advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
          
          adv_at = (np.diff(advrt,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
          adv_at =  adv_at/areaGrid

          tot=dtdt+adv_at
          
     else:
          tot=dtdt

     return tot

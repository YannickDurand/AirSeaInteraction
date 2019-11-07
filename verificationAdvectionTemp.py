#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 28 18:26:30 2019

@author: ydurand
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#import MITgcmutils as mit
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

from stockDiagKEs import *
from netCDF4 import Dataset
from fonctionH import *
from plotKES import *
#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=6*3600.                #en s
maxT=0.
alphat = 2e-4
rho0=999.8
tStart=5
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================
#
###############################################################################
#
# =============================================================================
# ###Without Wind
PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_wh_Wqt_Cyclone'
# =============================================================================
#
###############################################################################
#
# =============================================================================
###Different Wind
#PATH='/media/durand/Yannick/19_OAC/Control_Case/fplane_diff_Wind/conf1/Case_wh_qT_withW_ACyclone_V7'
PATH='/dsk1/ydurand/T03R01N3_2'
# =============================================================================


# =============================================================================
# ###Without Vortex
#PATH='/media/durand/Yannick/19_OAC/Control_Case/Case_withoutWqT_NoVortex'
# =============================================================================

# =============================================================================
# ###Without Atmosphere
#PATH='/media/durand/Yannick/19_OAC/Control_Case/Case_withoutAtm_ACyclone'
# =============================================================================


fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)

#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter2,time2=dimRes(dataNetCDFS)
nbIter,time=getTime(dataNetCDFK)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

#%%
# =============================================================================
# ---------------------------Evolution Potential Energy--------------------------
# =============================================================================
#%%
evoPE=np.array([])
evoPE2=np.array([])
evoWb2=np.array([])
evoDiffT=np.array([])
ix = np.int(si_x/2)

def comp_b (temp,salt):
  return alphat*temp #+ betas*salt

for i in range(tStart,nbIter):
     t=i
     dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys= dataDiagPE(dataNetPE,t)
     
     dtdt = dtdt/86400
     dsdt = dsdt/86400
     
     # t0 = np.where(t0 == 0,np.NaN,t0)
     # t1 = np.where(t1 == 0,np.NaN,t1)
          
     advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
     #advrs = np.append(advrs,advrs[None,0,:,:],axis=0)
     
     adv_at = (np.diff(advrt,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_at =  adv_at/areaGrid
     
     #adv_as = -( np.diff(advrs,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_as=0
     #adv_as = adv_as/(areaGrid)*0     

     ####################################
     
     dbdt = comp_b (dtdt,dsdt)
     advb = comp_b(adv_at,adv_as)
     tot=dbdt+advb
     dpedt = -areaGrid*9.81*999.8*(tot)*(zzc[:,np.newaxis,np.newaxis])
     
     evoPE=np.append(evoPE,(np.sum(dpedt)))
          
     wb2, wb = getWb(dataNetPE,i)
     wb2=wb2*9.81*areaGrid
     evoWb2=np.append(evoWb2,(np.sum(wb2)))
     print(i)
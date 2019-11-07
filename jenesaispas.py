#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  5 11:44:06 2019

@author: ydurand


"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#import MITgcmutils as mit
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

from stockDiagKEs import *
from stockData import *
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
f=1e-4
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
#PATH='/media/durand/Yannick/19_OAC/Control_Case/T03R01N3H60'
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
#define the pdt
pdt=6*3600
period=18*3600*2.                                                              #Two coriolis period
start = 100
end   = start + int(period/pdt)                                                #Compute the end in term of pdt (value of i)
www,rhoGW=0.,0.

for i in range(start,end):
    t=i
    www=www+getDoubleV(dataNetCDFS,t)
    rhoGW=rhoGW+getUPRESS(dataNetCDFK,t)
    ttt=getUPRESS(dataNetCDFK,t)
    
    plt.figure(i)
    plt.contourf(xxx*1e-3,zzc,www[:,165,:])

    plt.figure(i+100)
    plt.contourf(xxx*1e-3,zzc,ttt[:,165,:])
    
    
www=www/(end-start)
rhoGW=rhoGW/(end-start)
#%%
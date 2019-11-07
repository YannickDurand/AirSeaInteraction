#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 17:12:17 2019

@author: ydurand
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#import MITgcmutils as mit
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import scipy.stats

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
tStart=0.
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
PATH='/dsk1/ydurand/T03R01N3'
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
temp=getTemp(dataNetCDFS,tStart)
ujk,vik=getSpeed2(dataNetCDFS,tStart)

vik=vik[:,165,:]
temp2=temp[:,165,:]

rho=Rho(rho0,maxT,temp)
tempGrad=np.diff(temp2,axis=1)/1200/(1e-4*rho[:,:-1,165])*rho0
#%%
plt.contourf(xxx,zzc,vik,level=5)
plt.contour(xxx[:-1],zzc,tempGrad,level=5)
plt.colorbar()
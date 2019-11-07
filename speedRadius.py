#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 09:47:35 2019

@author: durand

A fonction to find the speed radius witch represent the radius of the eddy
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from stockData import *
from fonctionRadius import *

#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=21600.             #en s
maxT=0.
alphat = 2e-4
rho0=999.8
tStart=1
dx=1200.
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

PATH='/home/durand/OAC/Control_Case/T03R01N3'


fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
filePhy='/phiHyd.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)
dataNetPhy=(PATH+filePhy)


#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter,time=dimRes(dataNetCDFS)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

#%%
evoSpeedRadius=np.zeros(nbIter-tStart)
j=0

for i in range(tStart, nbIter):
   
    evoSpeedRadius[j]=speedRadius(dataNetCDFS,i,dx)
    print(j)
    j=j+1
    
#%%
t=158

verticalSpeedRadius=speedRadiusLevel(dataNetPhy,dataNetCDFS,t,dx)

evoVerticalSpeedRadius=np.zeros((nbIter-tStart,30))
j=0

for i in range(tStart, nbIter-1):
    t=i
    evoVerticalSpeedRadius[j]=speedRadiusLevel(dataNetPhy,dataNetCDFS,t,dx)
    print(j)
    j=j+1
#%%































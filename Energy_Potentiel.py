#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:15:04 2019

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
from timeDerivativeTemperature import *
#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=10800             #en s
maxT=0.
alphat = 2e-4
rho0=998.
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qt_withW_ACyclone_npS'
PATH='/media/ydurand/Yannick/19_OAC/Control_Case/fplane/conf1/Case_withoutWqT_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_withoutWqT_ACyclone'
#PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'

#PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_withoutWqT_ACyclone'
#PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_wh_qT_withW_ACyclone'



fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)


file5=PATH+fileKEs

#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter,time=dimRes(dataNetCDFS)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]


#%%
# =============================================================================
# ---------------------------Initial condition---------------------------------
# =============================================================================

#%%
# =============================================================================
# ---------------------------Evolution Potential Energy--------------------------
# =============================================================================

##----> timeDerivativeTemperature(dataNetPE,t,Grid,adv=-1):
#%%

totalEP,totalRHOGW  =  np.array([]),np.array([])
evoMeanTemp=np.array([])
evoRhoGW=np.array([])
evoDbDt=np.array([])

for i in range(0,nbIter-1):
     t=i
     temp=getTemp(dataNetCDFS,t)
     rho=Rho(rho0,maxT,temp)
     
     EPGrid=rho*9.81*areaGrid*((zzc[:,np.newaxis,np.newaxis]))
     totalEP=np.append(totalEP,np.sum(EPGrid))     
     
     wij=getDoubleV(dataNetCDFS,t)
     rhoGW=wij*9.81*rho*areaGrid
     evoRhoGW=np.append(evoRhoGW,(np.sum(rhoGW)))
     
     dTdt=timeDerivativeTemperature(dataNetPE,t,areaGrid,1)
     dbdt=alphat*dTdt*rho0
     dbdtVol=9.81*((zzc[:,np.newaxis,np.newaxis]))*dbdt*areaGrid
     evoDbDt=np.append(evoDbDt,np.sum(dbdtVol))
     
     print(i)
     
sumTotalEPGrid=(totalEP-totalEP[0])
derEPgrid=np.diff(sumTotalEPGrid/10800.)
sumEvoRhoGW=np.cumsum(evoRhoGW)
sumEvoDbDt=np.cumsum(evoDbDt)
#%%
# =============================================================================
# ##########--------------Difference -----#######################"# 
# =============================================================================
sumEvoTot=sumEvoRhoGW[:]-sumEvoDbDt
diffSumEvoTotSumEPGrid=sumTotalEPGrid[:]-sumEvoTot


#%%
# =============================================================================
# ------------------------Plot Plot Plots-------------------------------------
# =============================================================================
# =============================================================================
# plotComp(i,time,A,B,C,nameA,nameB,nameC)
# =============================================================================
zeros=np.zeros(240)
plotComp(900,time[:-1],evoRhoGW[:-1],evoDbDt[:-1],derEPgrid,"A","B","C")
#%%
plotComp(1001,time[:],sumEvoRhoGW[:],sumEvoDbDt,sumEvoTot,"A","B","C")
#%%
plotComp(2001,time[:],sumEvoTot,sumTotalEPGrid[:],diffSumEvoTotSumEPGrid,"A","B","C")
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:15:04 2019

@author: ydurand
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm


#import MITgcmutils as mit
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import scipy.stats

from stockDiagKEsMean import *
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
tStart=159-6*3*5
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================
#
###############################################################################
#PATH3='/dsk1/ydurand/T03R01N3'
PATH3='/home/durand/OAC/Control_Case/T03R01N3'
PATH2='/media/durand/Yannick/19_OAC/Control_Case/T03R01N2'
PATH1='/media/durand/Yannick/19_OAC/Control_Case/T03R01N1'

fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
fileTave='/diag_ocnTave.nc'


dataNetCDFK=(PATH1+fileKEs)
dataNetCDFS=(PATH1+fileS)
dataGrid=(PATH3+grid)
dataNetPE=(PATH1+filePE)

dataNetCDFK2=(PATH2+fileKEs)
dataNetCDFS2=(PATH2+fileS)
dataNetPE2=(PATH2+filePE)

dataNetCDFK3=(PATH3+fileKEs)
dataNetCDFS3=(PATH3+fileS)
dataNetPE3=(PATH3+filePE)
dataNetTave3=(PATH3+fileTave)

#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter2,time2=dimRes(dataNetCDFS3)
nbIter,time=getTime(dataNetCDFK3)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]


#%%

KEx3=getKEonX(dataNetCDFK3)
meanKEx3=np.mean(KEx3,axis=0)
meanKEx3=meanKEx3/86400.
advKE=getAdvonX(dataNetCDFK3)
meanAdvKE=np.mean(advKE,axis=0)

dissVx3=getDissVonY(dataNetCDFK3)
dissVx3=dissVx3[(160-6*3*5-28):]
meanDissVx3=np.mean(dissVx3,axis=0)

dissHx=getDissHonX(dataNetCDFK3)
meanDissHx3=np.mean(dissHx,axis=0)

########################for Test

termBarot=getBarotroponX(dataNetCDFK3)
meanTermBarot=np.mean(termBarot,axis=0)
termBaroc=getBarocliniconX(dataNetCDFK3)
meanTermBaroc=np.mean(termBaroc,axis=0)

#%%
wwx=getWonY(dataNetCDFS3)
meanWWX=np.mean(wwx,axis=0)
wbx=getWbonX(dataNetPE3)
meanWBX=np.mean(wbx,axis=0)
#%%
"""
# =============================================================================
# =============================================================================
# # plot
# =============================================================================
# =============================================================================
middle=np.zeros((30,330))
middle[:,int((330/2))]=1
#%%
plt.figure(1000)
plt.contourf(xxx*1e-3,zzc,meanKEx3, cmap=cm.jet, levels=100)
plt.title("Mean of the kinetic energy variation")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.legend()
plt.grid()
plt.show()

plt.figure(10001)
plt.contourf(xxx*1e-3,zzc[:10],meanAdvKE[:10,:], cmap=cm.jet, levels=100)
plt.title("Mean of the kinetic energy advection")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.legend()
plt.grid()
plt.show()

plt.figure(10)
plt.contourf(xxx[90:-90]*1e-3,zzc[:15],dissVx3[80,:15,90:-90], cmap=cm.jet, levels=100)
plt.title("Mean of the vertical dissipation variation")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.legend()
plt.grid()
plt.show()

plt.figure(1002)
plt.contourf(xxx*1e-3,zzc[:],meanDissHx3[:,:], cmap=cm.jet, levels=100)
plt.title("Mean of the horizontal dissipation variation")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.legend()
plt.grid()
plt.show()
#%%
#Closure
plt.figure(1004)
plt.contourf(xxx*1e-3,zzc[:],meanKEx3-meanAdvKE-meanDissHx3-meanDissVx3-meanTermBarot-meanTermBaroc, cmap=cm.jet, levels=100)
plt.title("Mean of the horizontal dissipation variation")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.grid()
plt.legend()
plt.title("Mean of the difference variation")
plt.show()
#%%
#plot for w
plt.figure(2004)
plt.contourf(xxx*1e-3,zzc[:],meanTermBaroc, cmap=cm.jet, levels=100)
plt.title("Mean on 30 days of the vertical speed")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.legend()
plt.show()

plt.figure(2005)
plt.contourf(xxx[120:200]*1e-3,zzc[:14],meanWWX[:14,120:200], cmap=cm.jet, levels=100,vmin=0.)
plt.title("Mean on 30 days of the vertical density speed")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.legend()
plt.show()

#%%
plt.figure(2006)
plt.contourf(xxx*1e-3,zzc[:],meanWWX, cmap=cm.jet, levels=100)
plt.title("Mean of the vertical dissipation variation")
plt.colorbar(label='Joules')
plt.xlabel("km")
plt.ylabel("m")
plt.legend()
plt.title("Mean of the vertical dissipation variation")
plt.show()
#%%
ujk,vik=getUVVEL(dataNetTave3)
#%%
vmax=0.5
vcont = np.linspace(-np.abs(vmax),np.abs(vmax),6)
CS=plt.contour(xxx[90:-90]*1e-3,zzc[:15],ujk[159,:15,90:-90,164],colors='b',levels=6,linewidths=1)
plt.clabel(CS,inline=1, fontsize=10)

#%%
"""
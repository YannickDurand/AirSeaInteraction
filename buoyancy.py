#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 14:16:37 2019

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

from stockDataBuoyancy import *
from stockData import *
from netCDF4 import Dataset
from Operator import *
from fonctionH import Rho
from fonctionH import vort_rel
from plotKES import *

#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=6*3600.                #en s
maxT=0.
alphat = 2e-4
rho0=999.8
tStart=7*4
dx=1200.
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================
#
###############################################################################
#PATH3='/dsk1/ydurand/T03R01N3'
PATH='/home/durand/OAC/Control_Case/T03R01N3'

fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
fileTave='/diag_ocnTave.nc'
filePhy='/phiHyd.nc'

dataGrid=(PATH+grid)
dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataNetPE=(PATH+filePE)
dataNetTave=(PATH+fileTave)
dataNetPhi=(PATH+filePhy)

#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter2,time2=dimRes(dataNetCDFS)
nbIter,time=getTime(dataNetCDFK)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

dz=np.ones((30,330,330))
dz=dz*zzc[:,np.newaxis,np.newaxis]

#%%
t=159
temp=getTemp(dataNetCDFS,t+1)
rho=Rho(rho0,maxT,temp)

phy=getPhi(dataNetPhi,t)
buoyancy=9.81*phy+9.81/rho0*rho*dz


stockMeanB=np.zeros(nbIter)
for t in range(nbIter):
     
     temp=getTemp(dataNetCDFS,t)
     rho=Rho(rho0,maxT,temp)
     drhodz=partwz02(rho,dz,30,330,330)
     #drhodz2=drhodz*9.81/rho0
     stockMeanB[t]=np.mean(drhodz)
     print(t)

phy=getPhi(dataNetPhi,t)
buoyancyEnd=9.81*phy+9.81/rho0*rho*dz


#%%
"""
#Temperature comparaison
t=0
temp=getTemp(dataNetCDFS,t)
ujk=getU(dataNetTave,t)
#ujk=ucc(ujk)
vik=getV(dataNetTave,t)
#vik=vcc(vik)
dy=dx
qij0=vort_rel(dx, dy, ujk, vik)


t=159
temp2=getTemp(dataNetCDFS,t)
ujk=getU(dataNetTave,t)
#ujk=ucc(ujk)
vik=getV(dataNetTave,t)
#vik=vcc(vik)
dy=dx
qij1=vort_rel(dx, dy, ujk, vik)



diffT=temp-temp2
diffQij=qij0-qij1

vmin=np.min(qij)#*1e-5
vmin=vmin*1e5
vmax=np.max(qij)
vcont = np.linspace(vmin,-0.1,6)

plt.figure()
plt.contourf(xxx[90:-91]*1e-3,xxx[90:-91]*1e-3,(temp2[0,90:-91,90:-91]-temp[0,90:-91,90:-91])/(40*24*3600))
plt.colorbar()
CS=plt.contour(xxx[90:-91]*1e-3,xxx[90:-91]*1e-3,qij1[0,90:-90,90:-90]*1e5,vcont,linewidths=1)
plt.clabel(CS,inline=1, fontsize=10)
plt.show()
"""
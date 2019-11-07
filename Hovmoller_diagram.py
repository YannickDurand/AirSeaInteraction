#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  7 15:36:50 2019

@author: ydurand
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from stockData import *
from stockDiagKEs import *
from Operator import *
from fonctionH import *
from fonctionRadius import *


#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=21600.             #en s
maxT=0.
alphat = 2e-4
rho0=999.8
rhoAtm=1.
f=1e-4
tStart=100
dx=1200.
UatmS=3.
kappa=0.4
DragCoeff=2.37e-3
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'
PATH='/dsk1/ydurand/T03R01N3_2'

fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
filePhy='/phiHyd.nc'
fileKPP='/kpp_state.nc'
fileAML='/cheapAML.nc'
fileTAU='/TAU.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)
dataNetPhy=(PATH+filePhy)
dataNetKPP=(PATH+fileKPP)
dataNetCheapAML=(PATH+fileAML)
dataNetTAU=(PATH+fileTAU)


#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter,time=dimRes(dataNetCDFS)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

#%%

evoSpeedRadius=np.zeros(nbIter-tStart)
evoCenterEddy=np.zeros((nbIter-tStart,2))
j=0

for i in range(tStart+1, nbIter):
    t=i
    evoSpeedRadius[j]=speedRadius(dataNetCDFS,i,dx)
    aaa=list(getCenterEta(dataNetCDFS,t))
    evoCenterEddy[j,0]=aaa[0]
    evoCenterEddy[j,1]=aaa[1]
    j=j+1
    
del aaa
#%%
HovmollerVert=np.zeros(((nbIter-tStart),30))
HovmollerHori=np.zeros(((nbIter-tStart),330))

j=0
for i in range(tStart,nbIter-10):
     t=i
     www=vitesseW(dataNetCDFS,t)
     centerGridy=evoCenterEddy[j,0]
     centerGridx=evoCenterEddy[j,1]
     nbGridPoint=int(evoSpeedRadius[j]/dx)
     boundaryLeft=int(centerGridx-nbGridPoint)-2

     www0=www[:,boundaryLeft-50,165]
     HovmollerVert[j,:]=www0
     j=j+1

     www1=www[20,165,:]
     HovmollerHori[j,:]=www1






















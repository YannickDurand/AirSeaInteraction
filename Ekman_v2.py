#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 11:04:08 2019

@author: durand

Script to compute the Ekman Velocity follow Stern :
     - The mean is done on 2 periode of Coriolis Parameter

"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from stockData import *
from stockDiagKEs import *
from Operator import *
from fonctionH import *


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
tStart=5
dx=1200.
UatmS=3.
kappa=0.4
DragCoeff=2.37e-3
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

PATH='/media/ydurand/Yannick/19_OAC/Control_Case/T03R01N3_Barotrope'

fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
filePhy='/phiHyd.nc'
fileKPP='/kpp_state.nc'
fileAML='/cheapAML.nc'
fileTAU='/TAU.nc'
fileTave='/diag_ocnTave.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)
dataNetPhy=(PATH+filePhy)
dataNetKPP=(PATH+fileKPP)
dataNetCheapAML=(PATH+fileAML)
dataNetTAU=(PATH+fileTAU)
dataNetTave=(PATH+fileTave)


#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter,time=dimRes(dataNetCDFS)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

#%%
#define the pdt
pdt=6*3600
period=18*3600*30.                                                              #Two coriolis period
start = 33
end   = 160-start                                                #Compute the end in term of pdt (value of i)
ekmanVelocity=0.

for i in range(start,end):
     t=i
     
     ujk,vik=getSpeed2(dataNetCDFS,t)
     qij=vort_rel(dx, dx, ujk, vik)
     qijend=qij[:,1,:]
     qijstart=qij[:,-1,:]
     qij=np.append(qij,qijend[:,np.newaxis,:],axis=1)
     qij=np.insert(qij,0,qijstart,axis=1)
     qijend=qij[:,:,1]
     qijstart=qij[:,:,-1]
     qij=np.append(qij,qijend[:,:,np.newaxis],axis=2)
     qij=np.insert(qij,0,qijstart,axis=2)
     
     qijCenterY=vcc(qij)
     qijCenterX=ucc(qij)
     
     windStressSurfaceX, windStressSurfaceY=getTAU(dataNetTAU,t)
     
     
     termX=windStressSurfaceX[:,:]/(rho0*(f+qijCenterY[1,:,:]))
     termY=windStressSurfaceY[:,:]/(rho0*(f+qijCenterX[1,:,:]))
     
     gradxWSSY=np.diff(termY,axis=1)/dx
     gradyWSSX=np.diff(termX,axis=0)/dx
     ekmanVelocity=-(gradyWSSX[:,1:-1]-gradxWSSY[1:-1,:])+ekmanVelocity
     

ekmanVelocity=ekmanVelocity/(end-start)
#%%
#for i in range
www=0
#del www2
for i in range(start,end):
     www2=0
     www2=vitesseWTave(dataNetTave,t)
     www2=www2[8]
     www=www+www2
     print(i)
www=www/(end-start)
#%%
diffEkmanVelocityWWW=ekmanVelocity-www[:-1,:-1]


#%%

plt.figure(11)
plt.contourf(ekmanVelocity*1e5)
plt.colorbar(label="$10^{-5}.m.s^{-1}$")
plt.title("Ekman Velocity from Stern")
plt.show()

plt.figure(12)
plt.contourf(www*1e5)
plt.colorbar(label="$10^{-5}.m.s^{-1}$")
plt.title("www at z="+str(int(zzc[8])))
plt.show()

plt.figure(13)
plt.contourf(diffEkmanVelocityWWW*1e5)
plt.colorbar(label="$10^{-5}.m.s^{-1}$")
plt.title("diff")
plt.show()


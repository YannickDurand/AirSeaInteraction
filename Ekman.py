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

PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'
PATH='/dsk1/ydurand/T03R01N3_2'
PATH='/home/durand/OAC/Control_Case/T03R01N3'

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
#define the pdt
pdt=6*3600
period=18*3600*5.                                                              #Two coriolis period
start = 60
end   = start + int(period/pdt)                                                #Compute the end in term of pdt (value of i)
ekmanVelocity=0.

for i in range(start,end):
     t=i
     
     ujk,vik=getSpeed2(dataNetCDFS,t)
     qij=vort_rel(dx, dx, ujk, vik)

     windStressSurfaceX, windStressSurfaceY=getTAU(dataNetTAU,t)
     
     termX=windStressSurfaceX[1:,1:-1]
     termY=windStressSurfaceY[1:-1,1:]
     
     gradxWSSY=np.gradient(termY,dx,axis=1)
     gradyWSSX=np.gradient(termX,dx,axis=0)
     ekmanVelocity=(gradyWSSX[:,:]-gradxWSSY[:,:])/(rho0*(f+qij[0,:,:]))+ekmanVelocity
     

ekmanVelocity=ekmanVelocity/(end-start)
#%%
#for i in range
www=0
#del www2
for i in range(start,end):
     www2=0
     www2=vitesseW(dataNetCDFS,t)
     #www2=www2[4]
     www=www+www2
     print(i)
www=www/(end-start)
#%%
diffEkmanVelocityWWW=ekmanVelocity-www[1:,1:]


#%%

plt.figure(11)
plt.contourf(ekmanVelocity*1e5)
plt.colorbar(label="$10^{-5}.m.s^{-1}$")
plt.title("Ekman Velocity from Stern")
plt.show()

plt.figure(12)
plt.contourf(www*1e5)
plt.colorbar(label="$10^{-5}.m.s^{-1}$")
plt.title("www at z="+str(int(zzc[4])))
plt.show()

plt.figure(13)
plt.contourf(diffEkmanVelocityWWW*1e5)
plt.colorbar(label="$10^{-5}.m.s^{-1}$")
plt.title("diff")
plt.show()


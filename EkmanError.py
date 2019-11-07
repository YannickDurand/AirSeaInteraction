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
PATH='/media/durand/Yannick/19_OAC/Control_Case/T03R01N3'

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
#Distributiuon of the wind stress
t=150



#define the pdt
pdt=6*3600
period=18*3600*7.                                                              #Two coriolis period
start = 100
end   = start + int(period/pdt)                                                #Compute the end in term of pdt (value of i)
divWind,ekmanVelocity=0.,0.

for i in range(start,end):
     t=i
     
     ujk,vik=getSpeed2(dataNetCDFS,t)
     qij=vort_rel(dx, dx, ujk, vik)
     
     divUJK=np.diff(ujk[3],axis=1)/1200.
     divVIK=np.diff(vik[3],axis=0)/1200.
     divWind=divUJK+divVIK+divWind
     
     ujk,vik=getSpeedSurface(dataNetCDFS,t)
     eta=getEta(dataNetCDFS,t)
     Uatm=getUatm(dataNetCheapAML,t)
     Vatm=np.zeros((331,330))
     
     ucc=ucc1L(ujk)
     vcc=vcc1L(vik)
     
     temp=getTemp(dataNetCDFS,t+1)
     rho=Rho(rho0,maxT,temp)
     rho=rho[0]
     
     
     Udiffx=Uatm[:,:]-ujk
     UdiffxC=ucc1L(Udiffx)
     Udiffy=Vatm[:,:]-vik
     UdiffyC=vcc1L(Udiffy)
     
     Udiff=Udiffx[:,:-1]+Udiffy[:-1,:]
     normUdiff=(UdiffxC**2+UdiffyC**2)**(0.5)
     windStressSurfaceX=rhoAtm*DragCoeff*normUdiff*UdiffxC
     windStressSurfaceY=rhoAtm*DragCoeff*normUdiff*UdiffyC
     
     termX=windStressSurfaceX[:-1,:-1]
     termY=windStressSurfaceY[:-1,:-1]
     
     gradxWSSY=np.gradient(termY,dx,axis=1)
     gradyWSSX=np.gradient(termX,dx,axis=0)
     ekmanVelocity=-(gradyWSSX[:,:]-gradxWSSY[:,:])/(rho0*(f+qij[0,:,:]))+ekmanVelocity
     

ekmanVelocity=ekmanVelocity/(end-start)
divWind=divWind/(end-start)
#%%
#for i in range
www=0

for i in range(start,end):
     www2=0
     www2=vitesseW(dataNetCDFS,t)
     www2=www2[3]
     www=www+www2
     print(i)
www=www/(end-start)
#%%
diffEkmanVelocityWWW=ekmanVelocity-www[:-1,:-1]
diffEQ2=diffEkmanVelocityWWW+divWind[:-1,:-1]
#%%
plt.figure(11)
plt.contourf(ekmanVelocity*1e5)
plt.colorbar(label="$m.s^{-1}$")
plt.title("Wrong term")
plt.show()

plt.figure(12)
plt.contourf(www*1e5)
plt.colorbar(label="$m.s^{-1}$")
plt.title("www at z=21 m")
plt.show()

plt.figure(13)
plt.contourf(diffEkmanVelocityWWW*1e5)
plt.colorbar(label="$m.s^{-1}$")
plt.title("diff")
plt.show()

plt.figure(14)
plt.contourf(diffEQ2*1e5)
plt.colorbar(label="$m.s^{-1}$")
plt.title("diff")
plt.show()
#%%
#Compute of the divergence of the wind

#%%                                             #Compute the end in term of pdt (value of i)
ekmanVelocity2=0.
"""
for i in range(start,end):
     t=i
     
     ujk,vik=getSpeed2(dataNetCDFS,t)
     qij=vort_rel(dx, dx, ujk, vik)

     
     windStressSurfaceX, windStressSurfaceY=getTAU(dataNetTAU,t)
     
     termX=windStressSurfaceX[1:,1:-1]/(rho0*(f+qij[0,:,:]))
     termY=windStressSurfaceY[1:-1,1:]/(rho0*(f+qij[0,:,:]))
     
     gradxWSSY=np.gradient(termY,dx,axis=1)
     gradyWSSX=np.gradient(termX,dx,axis=0)
     ekmanVelocity2=-(gradyWSSX[:,:]-gradxWSSY[:,:])+ekmanVelocity2
     

ekmanVelocity2=ekmanVelocity2/(end-start)

diffdeladiff=diffEkmanVelocityWWW+abs(ekmanVelocity2)

#%%
plt.figure(14)
plt.contourf(diffdeladiff)
plt.colorbar()
plt.show()
"""
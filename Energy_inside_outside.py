#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 25 11:07:00 2019

@author: durand

Fonction to vertify the kinetic energy budget inside and outside the eddy

The radius of the eddy is define by the speed radius.

The center of the eddy is define by the maximum of eta at the surface
or the maximum of the phy pertubation into each layer.
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from fonctionRadius import *

from stockDiagKEs import *
from netCDF4 import Dataset
from fonctionH import *
from plotKES import *
import matplotlib.cm as cm


#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=21600.             #en s
maxT=0.
alphat = 2e-4
rho0=999.8
tStart=5
dx=1200.
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

PATH='/dsk1/ydurand/T03R01N3'



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
# =============================================================================
#   #####-----  Stockage of the speed radius and the center  -----#####
# =============================================================================
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
"""
evoVerticalSpeedRadius=np.zeros((nbIter-tStart,30))
j=0

for i in range(tStart, nbIter-1):
    t=i
    evoVerticalSpeedRadius[j]=speedRadiusLevel(dataNetPhy,dataNetCDFS,t,dx)
    print(j)
    j=j+1
"""
#%%
# =============================================================================
#          #####-----     Stockage of Depht     -----#####
# =============================================================================
evoDepth=np.zeros(nbIter-tStart)
j=0
for i in range(tStart, nbIter-1):
     t=i
     evoDepth[j]=depthLevel2(dataNetCDFS,t,dx)
     j=j+1
#%%
# =============================================================================
# ---------------------------Evolution Kinetic Energy--------------------------
# =============================================================================
    
#%%
evoEK,totFEXT,totUDISS,totUPRESS=np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart)
totUDISSh,totUDISSv=np.zeros(nbIter-tStart),np.zeros(nbIter-tStart)

evoInsideEK,evoInsideFEXT,evoInsideUDISS,evoInsideUPRESS=np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart)
evoInsideUDISSh,evoInsideUDISSv=np.zeros(nbIter-tStart),np.zeros(nbIter-tStart)

evoInsideDISSv=np.zeros(nbIter-tStart)
j=0
for i in range(tStart,nbIter-1):
     t=i
     utot, uadv, upress, udissh, udissv, uext = meshcons2(dataNetCDFK,i)
     utot=utot/86400.
     
     
     temp=getTemp(dataNetCDFS,t+1)
     rho=Rho(rho0,maxT,temp)
     pro=int(evoDepth[j])
     
     
     #Calculation on the volume
     utot, uadv, upress, udissh, udissv, uext =utot*areaGrid, uadv*areaGrid, upress*areaGrid, udissh*areaGrid, udissv*areaGrid, uext*areaGrid
         
     #add rho
     t=i     
     #multiply by rho    add today
     utot=utot*rho0+uadv*rho0
     
     #total kinetic Energy
     evoEK[j]=(np.sum(utot))
     
     #external Force
     uext=uext*rho0
     totFEXT[j]=np.sum(uext)
     
     #dissipation
     dissip=udissh*rho0+udissv*rho0
     totUDISSh[j]=np.sum(udissh*rho0)
     totUDISSv[j]=np.sum(udissv*rho0)
     totUDISS[j]=(np.sum(udissh*rho0)+np.sum(udissv*rho0))
     
     #force de pression
     upress=upress*rho0
     totUPRESS[j]=np.sum(upress)  
     
     #####-----    Create grid     -----#####
     centerGridy=evoCenterEddy[j,0]
     centerGridx=evoCenterEddy[j,1]
     nbGridPoint=int(evoSpeedRadius[j]/dx)
     
     boundaryLeft=int(centerGridx-nbGridPoint)-2
     boundaryRight=int(centerGridx+nbGridPoint)+2
     boundaryTop=int(centerGridy+nbGridPoint)+2
     boundaryBottom=int(centerGridy-nbGridPoint)-2
     
     centerXY=(boundaryRight-boundaryLeft)/2
     #utotInside=utot[:,boundaryLeft:boundaryRight,boundaryBottom:boundaryTop]
     
     ####----     divide inside and outside     -----######
     for m in range(boundaryLeft,boundaryRight):
          for k in range(boundaryBottom,boundaryTop):
               if (((k-centerGridx)*dx)**2+((m-centerGridy)*dx)**2)**(1/2)<evoSpeedRadius[j]:
                    evoInsideEK[j]=evoInsideEK[j]+np.sum(utot[pro:,k,m])
                    evoInsideUPRESS[j]=evoInsideUPRESS[j]+np.sum(upress[pro:,k,m])
                    evoInsideFEXT[j]=evoInsideFEXT[j]+np.sum(uext[pro:,k,m])
                    evoInsideUDISS[j]=evoInsideUDISS[j]+np.sum(dissip[pro:,k,m])
                    evoInsideUDISSv[j]=evoInsideUDISSv[j]+np.sum(udissv[pro:,k,m]*rho0)
                    evoInsideUDISSh[j]=evoInsideUDISSh[j]+np.sum(udissh[pro:,k,m]*rho0)
     
     print(i)
     j=j+1
del utot, uadv, upress, udissh, udissv, uext, temp, rho

#%%
evoOutsideEK=evoEK-evoInsideEK
evoOutsideUPRESS=totUPRESS-evoInsideUPRESS 
evoOutsideFEXT=totUPRESS-evoInsideFEXT
evoOutsideUDISS=totUDISS-evoInsideUDISS   
evoOutsideUDISSv=totUDISSv-evoInsideUDISSv   
evoOutsideUDISSh=totUDISSh-evoInsideUDISSh  


sumTotInsideEK=np.cumsum(evoInsideEK*deltaT)
sumTotInsideUPRESS=np.cumsum(evoInsideUPRESS*deltaT)
sumTotInsideFEXT=np.cumsum(evoInsideFEXT*deltaT)
sumTotInsideUDISS=np.cumsum(evoInsideUDISS*deltaT)
sumTotInsideUDISSv=np.cumsum(evoInsideUDISSv*deltaT)
sumTotInsideUDISSh=np.cumsum(evoInsideUDISSh*deltaT)

sumTotOutsideEK=np.cumsum(evoOutsideEK*deltaT)
sumTotOutsideUPRESS=np.cumsum(evoOutsideUPRESS*deltaT)
sumTotOutsideFEXT=np.cumsum(evoOutsideFEXT)
sumTotOutsideUDISS=np.cumsum(evoOutsideUDISS*deltaT)
sumTotOutsideUDISSv=np.cumsum(evoOutsideUDISSv*deltaT)
sumTotOutsideUDISSh=np.cumsum(evoOutsideUDISSh*deltaT)


sumTotEK=np.cumsum(evoEK*deltaT)
sumTotDISS=np.cumsum(totUDISS*deltaT)
sumTotFEXT=np.cumsum(totFEXT*deltaT)
sumTotUPRESS=np.cumsum(totUPRESS*deltaT)
#%%
# =============================================================================
# #######------      Divide inside and outside for rhogw     -----#####
# =============================================================================

#%%
# =============================================================================
# #####-----      Plot plot plots    -----#####
# plotCompT(i,time,A,B,C,nameA,nameB,nameC,title):
# =============================================================================

plotComp2T(100,time[tStart-1:],sumTotInsideEK,sumTotOutsideEK,'KE inside','KE outside', 'Evolution of Kinetic enrgy outside and inside')
plotComp2T(101,time[tStart-1:],sumTotInsideUPRESS,sumTotOutsideUPRESS,'$rogw$ inside','$rhogw$ outside', 'Evolution of the transfer therm outside and inside')
plotComp2T(102,time[tStart-1:],sumTotInsideUDISS,sumTotOutsideUDISS,'$diss$ inside x','$diss$ outside x', 'Evolution of the total dissipation outside and inside')
plotComp2T(103,time[tStart-1:],sumTotInsideUDISSv,sumTotOutsideUDISSv,'$dissV$ inside x','$dissV$ outside x', 'Evolution of the vertical dissipation outside and inside')
plotComp2T(104,time[tStart-1:],sumTotInsideFEXT,sumTotOutsideFEXT,'$F_{ext}$ inside x','$F_{ext}$ outside x', 'Evolution of the external forcing outside and inside')
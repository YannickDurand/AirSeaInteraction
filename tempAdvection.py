#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 14:06:20 2019

@author: durand
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage as scp

from stockDiagKEs import *
from fonctionH import *
from Sauvegarde import *
from Operator import *

#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=21600.           #en s
maxT=0.
dh=1200.
alphaT = 2e-4
rho0=999.8

#####-----     diferrent for each case     -----#####
# 

# 1.8e-5*6  for A-Cyclone without wind
# 3.0e-5*6  for Cyclone   without wind

# 2.5e-4*6  for No-Vortex with wind     10m/s                                  #sum inverse
# 1.1e-4*6  for No-Vortex without wind

# 5.0e-4*6  for A-Cyclone without Atm
coefDiff=1e-5*6
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================


PATH='/dsk1/ydurand/T03R01N3_2'


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

dz=np.ones((30,330,330))
dzzz=np.ones((30,330,330))

dz=dz*zzc[:,np.newaxis,np.newaxis]
dzzz=dzzz*zzc[:,np.newaxis,np.newaxis]
dz=np.diff(dz,axis=0)


#%%
evoAdvXYZ=np.array([])
evoPEAdvXYZ=np.array([])

j=0
"""

advTot=np.zeros(((int(nbIter/8)),19,329,329))
advXt=np.zeros(((int(nbIter/8)),19,329,329))
advYt=np.zeros(((int(nbIter/8)),19,329,329))
advZt=np.zeros(((int(nbIter/8)),19,329,329))

for i in range(39,nbIter-1):
     t=i
     hhd,ujk,vik,wij,temp,eta=vitesse(dataNetCDFS,zzt,HHH=-1,tstart=i)
     rho=Rho(rho0,maxT,temp)
     rho=rho[:-1,:-1,:-1]
     
     gradTempX=np.gradient(temp,dh,axis=2)
     gradTempY=np.gradient(temp,dh,axis=1)
     
     dTdx=partui(temp,dh)
     dTdy=partvj(temp,dh)
     
     dTdz=-partwz(temp,dz)
     
     uuu=ucc(ujk)
     vvv=vcc(vik)
     
     advX=uuu[:,:,:-1]*dTdx
     advY=vvv[:,:-1,:]*dTdy
     advZ=wij[:-1,:,:]*dTdz
     advXYZ=advX[:-1,:-1,:]+advY[:-1,:,:-1]+advZ[:,:-1,:-1]
     
     #advXt=np.append(advXt,advX)
     #advYt=np.append(advYt,advY)
     #advZt=np.append(advZt,advZ)

     if i%300==1:
         
         advTot[j,:,:,:]=advXYZ[:,:,:]
         advXt[j,:8.5,:,:]=advX[:-1,:-1,:]
         advYt[j,:,:,:]=advY[:-1,:,:-1]
         advZt[j,:,:,:]=advZ[:,:-1,:-1]
         j=j+1
         print(j)

     PEAdvXYZ=advXYZ*rho*9.81*alphaT*zzc[:-1,np.newaxis,np.newaxis]*areaGrid[:-1,:-1,:-1]
     evoAdvXYZ=np.append(evoAdvXYZ,np.sum(advXYZ))
     evoPEAdvXYZ=np.append(evoPEAdvXYZ,np.sum(PEAdvXYZ))
     
sumEvoAdvXYZ=np.cumsum(evoAdvXYZ*deltaT)   
sumEvoPEAdvXYZ=np.cumsum(evoPEAdvXYZ*deltaT)   

"""  
#%%
evoDiffTherm=np.array([])

###Here we compute the Laplacian
for i in range(5,nbIter-1):
     t=i
     temp=getTemp(dataNetCDFS,i)
     rho=Rho(rho0,maxT,temp)

     dTdx=np.gradient(temp,dh,axis=2)
     dTdy=np.gradient(temp,dh,axis=1)
     
     dT2dx=np.gradient(dTdx,dh,axis=2)
     dT2dy=np.gradient(dTdy,dh,axis=1)

     
     dTdz=-partwz02(temp,dzzz,30,330,330)
     dT2dz=-partwz02(dTdz,dzzz,30,330,330)
     
     dT2dv=dT2dx+dT2dy+dT2dz
     
     PEdiffTherm=-rho*9.81*alphaT*zzc[:,np.newaxis,np.newaxis]
     PEdiffTherm=PEdiffTherm*areaGrid[:,:,:]*dT2dv*coefDiff
     evoDiffTherm=np.append(evoDiffTherm,np.sum(PEdiffTherm))
     
     print(i)
     
sumEvoDiffTerm=np.cumsum(abs(evoDiffTherm*deltaT))
print(sumEvoDiffTerm[154])


#%%
"""
# =============================================================================
# ##########--------------Difference -----#######################"# 
# =============================================================================
sumDiffPEEvo=sumEvoPE-sumTotalEPGrid                                           # Here sumEvoPE is only the variation of d(rho)/dt

sumDiffPEEvo2=sumDiffPEEvo-sumEvoWb                                            # Here we add (rhogw)

sumDiffPEEvo2a=sumDiffPEEvo2+sumEvoRhoGW

sumDiffPEEvo3=sumDiffPEEvo2-sumEvoDiffTerm
####    We have almost the same thing but one parameter diffe
####    Differ by 1.3e10 at the end
####    Where this Oscillation come from?    

####     ----     For Total Energy     ---- ####
sumEvoTot=sumEvoPE+sumTotEK
sumEvoTot2=sumEvoTot-sumTotDISS-sumTotDISSModel
sumEvoTot3=sumEvoTot2-sumTotFEXT
sumEvoTot4=sumEvoTot3-sumEvoDiffTerm
sumEvoTot5=sumEvoTot4-sumTotUETA

plotComp(919,time[39:],sumEvoTot3,sumEvoDiffTerm,sumEvoTot4,"A","B","C")
"""
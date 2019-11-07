#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 15:16:06 2019

@author: ydurand

Fonction to check the depth of the isotherm of 12 °C
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from fonctionRadius import *

from stockData import *
from stockDiagKEs import *


#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=21600.             #en s
maxT=0.
alphat = 2e-4
rho0=999.8
tStart=40
dx=1200.
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

PATH='/dsk1/ydurand/T03R01N3_2'

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
"""
evoSpeedRadius=np.zeros(nbIter-tStart)
evoCenterEddy=np.zeros((nbIter-tStart,2))
j=0

for i in range(tStart, nbIter):
    t=i
    evoSpeedRadius[j]=speedRadius(dataNetCDFS,i,dx)
    aaa=list(getCenterEta(dataNetCDFS,t))
    evoCenterEddy[j,0]=aaa[0]
    evoCenterEddy[j,1]=aaa[1]
    
    #print(j)
    j=j+1
    
"""
#%%


chooseTemp=15.
isoThermInit=np.zeros((330,330))
isoThermInit2=np.zeros((330,330))

isoTherm=np.zeros((330,330))
#%%
t=0
tempInit=getTemp(dataNetCDFS,t)
plt.contourf(xxx,zzc,tempInit[:,165,:])
plt.colorbar()
plt.show()
for i in range(330):
     for j in range(330):
          k=0
          while temp[k,j,i]>chooseTemp:
               k=k+1
               
          k=k                                       #level number of the layer over 15°C
          
          #compute the coefficient of T between the two layer close to 15°C
          alphaT=tempInit[k+1,j,i]-tempInit[k,j,i]
          alphaT=alphaT/(zzc[k+1]-zzc[k])
          
          zzz=(chooseTemp-tempInit[k,j,i])/alphaT
          zzz2=zzz+zzc[k]
          #print(zzc[k])
          isoThermInit[j,i]=zzc[k]
          isoThermInit2[j,i]=zzz2
     print(i)


#%%
t=130
temp=getTemp(dataNetCDFS,t)
for i in range(330):
     for j in range(330):
          k=0
          while temp[k,j,i]>chooseTemp:
               k=k+1
          
          k=k                                       #level number of the layer over 15°C
          
          #compute the coefficient of T between the two layer close to 15°C
          alphaT=temp[k+1,j,i]-temp[k,j,i]
          alphaT=alphaT/(zzc[k+1]-zzc[k])
          
          zzz=(chooseTemp-temp[k,j,i])/alphaT+zzc[k]
          isoTherm[j,i]=zzz
     print(i)

#%%
diffDepthIsoTherm=isoTherm-isoThermInit
diffTempVertical=tempInit[:,165,:]-temp[:,165,:]
#%%
plt.figure(100)
plt.contourf(xxx*1e-3,xxx*1e-3,isoThermInit2)
plt.title("Depth of the isotherme 10°C at t= 5 days")
plt.colorbar(label='m')
plt.show()

plt.figure(101)
plt.contourf(xxx*1e-3,xxx*1e-3,isoTherm)
plt.title("Depth of the isotherme 10°C at t="+str(int(time[150])))
plt.colorbar(label='m')
plt.show()

plt.figure(102)
plt.contourf(xxx*1e-3,xxx*1e-3,diffDepthIsoTherm)
plt.title("Difference of the isotherme 10°C between t=0 and t="+str(int(time[150])))
plt.colorbar(label='m')
plt.show()

#%%
plt.contourf(xxx,zzc,diffTempVertical)
plt.colorbar()
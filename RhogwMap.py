#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 08:47:18 2019

@author: durand

Map of Rho gw
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from stockData import *
from stockDiagKEs import *
from stockKPP import *
from netCDF4 import Dataset
from fonctionH import *
from plotKES import *

from fonctionRadius import *


import matplotlib.cm as cm

plt.close('all')
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

PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'
PATH='/dsk1/ydurand/T03R01N3_2'


fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
filePhy='/phiHyd.nc'
fileKPP='/kpp_state.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)
dataNetPhy=(PATH+filePhy)
dataNetKPP=(PATH+fileKPP)


#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter,time=dimRes(dataNetCDFS)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

#%%
evoMeanBoundaryDepth,evoMeanBoundaryDepthLevel=np.array([]),np.array([])
t=100

for i in range(5,nbIter-1):
    t=i
    boundaryDepth=getDepthBoundary(dataNetKPP,t)
    boundaryDepthMean=np.mean(boundaryDepth)
    j=0
    while boundaryDepthMean>abs(zzc[j]):
        layerLevel=j
        j=j+1
    evoMeanBoundaryDepthLevel=np.append(evoMeanBoundaryDepthLevel,int(layerLevel))
    evoMeanBoundaryDepth=np.append(evoMeanBoundaryDepth,boundaryDepthMean)
        
    
plt.figure(i)
plt.contourf(xxx[:]*1e-3,xxx[:]*1e-3,boundaryDepth[:,:],cmap=cm.jet, levels=100)
plt.title('Boundary Depth at t='+str(time[t])+' days')
plt.colorbar()
plt.show()
#%%
period=18*3600*2. 
pdt=6.*3600.
"""
for i in range(28,nbIter,28):
    t=i
    start = i
    end   = start + int(period/pdt)
    rhogw=np.zeros((30,330,330))
    for j in range(start,end):
         rhogw=rhogw+getUPRESS(dataNetCDFK,t)*1e8
    rhogw=rhogw/(end-start)    
    
    zzzSurface=1
    zzzBellow=int(evoMeanBoundaryDepthLevel[i]-1)
    
    #Plot horizontal
    plt.figure((i*10)+2000)
    plt.xlabel('x (km)')
    plt.ylabel("y (km)")
    CS=plt.contour(xxx[:]*1e-3,xxx[:]*1e-3,rhogw[zzzSurface,:,:],cmap='gray', levels=5)
    plt.clabel(CS,fmt='%2.1f',fontsize=9, inline=1,)
    plt.contourf(xxx[:]*1e-3,xxx[:]*1e-3,rhogw[zzzBellow,:,:],cmap=cm.jet, levels=100)
    plt.colorbar(label='$kg.m^{-2}.s^{-3}x10^{-8}$')
    plt.title('Rho gw at the surface(contours) z= '+str(int(evoMeanBoundaryDepth[i]))+'m (colors), t='+str(time[i])+' days')
    plt.legend()
    plt.show()


    plt.figure((i*10)+1000)
    CS=plt.contour(xxx[100:230]*1e-3,xxx[100:230]*1e-3,rhogw[zzzSurface,100:230,100:230],cmap='gray', levels=5)
    plt.clabel(CS,fmt='%2.1f',fontsize=9, inline=1,)
    plt.contourf(xxx[100:230]*1e-3,xxx[100:230]*1e-3,rhogw[zzzBellow,100:230,100:230],cmap=cm.jet, levels=100)
    plt.colorbar(label='$kg.m^{-2}.s^{-3}x10^{-8}$')
    plt.xlabel('x (km)')
    plt.ylabel("y (km)")
    plt.title('Rho gw at the surface(contours) z= '+str(int(evoMeanBoundaryDepth[i]))+'m (colors), t='+str(time[i])+' days')
    plt.show()

     #plot vertical
     #Can we find a relation between the gain of Kinetic Energy and rho gw,
    plt.figure((i*10)+5000)
    plt.contourf(xxx[:]*1e-3,zzc,rhogw[:,165,:],cmap=cm.jet, levels=100)
    plt.colorbar(label='$kg.m^{-2}.s^{-3}x10^{-8}$')
    plt.xlabel('x (km)')
    plt.ylabel("z (m)")
    plt.legend()
    plt.title('Rho gw vertical in relation to x at, t='+str(time[i])+' days')
    plt.show()
     
    plt.figure((i*10)+5100)
    plt.contourf(xxx[:]*1e-3,zzc,rhogw[:,:,165],cmap=cm.jet, levels=100)
    plt.colorbar(label='$kg.m^{-2}.s^{-3}x10^{-8}$')
    plt.title('Rho gw vertical in relation to y at, t='+str(time[i])+' days')
    plt.xlabel('y (km)')
    plt.ylabel("z (m)")
    plt.show()



#%%
# =============================================================================
# =============================================================================
                         # # Vertical Plot # # 
# =============================================================================
# =============================================================================
#%%



#%%
###############################################################################
# =============================================================================
# # Here we separate the plot between inside and outside!
# =============================================================================
###############################################################################
#%%
"""
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
# =============================================================================
#          #####-----     Stockage of Depht     -----#####
# =============================================================================

evoDepth=np.zeros(nbIter-tStart)
j=0
for i in range(tStart, nbIter-1):
     t=i
     evoDepth[j]=depthLevel(dataNetPhy,dataNetCDFS,t,dx)
     j=j+1
     

#%%
#define the pdt
pdt=6*3600
period=18*3600*2.                                                              #Two coriolis period

"""
#####-----    Create grid     -----#####
centerGridy=evoCenterEddy[start,0]
centerGridx=evoCenterEddy[start,1]
nbGridPoint=int(evoSpeedRadius[start]/dx)

boundaryLeft=int(centerGridx-nbGridPoint)-2     
boundaryRight=int(centerGridx+nbGridPoint)+2
boundaryTop=int(centerGridy+nbGridPoint)+2
boundaryBottom=int(centerGridy-nbGridPoint)-2

"""

for j in range(23,nbIter,28):
     
     
     
     start = j

     end   = start + int(period/pdt)                                                #Compute the end in term of pdt (value of i)
     www,rhogwInside=0.,0.
     pro=int(evoDepth[j])
     
     centerGridy=evoCenterEddy[start,0]
     centerGridx=evoCenterEddy[start,1]
     nbGridPoint=int(evoSpeedRadius[start]/dx)
     
     boundaryLeft=int(centerGridx-nbGridPoint)-2     
     boundaryRight=int(centerGridx+nbGridPoint)+2
     boundaryTop=int(centerGridy+nbGridPoint)+2
     boundaryBottom=int(centerGridy-nbGridPoint)-2
     print('boudnary lesft',boundaryLeft,'   boundary right',boundaryRight)
     
     midmid=int((boundaryRight-boundaryLeft)/2.)-1
     
     ujk,vik=getSpeed2(dataNetCDFS,start)
     ujk,vik=ujk[:pro,boundaryLeft:boundaryRight,boundaryBottom:boundaryTop],vik[:pro,boundaryLeft:boundaryRight,boundaryBottom:boundaryTop]
     
     for i in range(start,end):
          t=i     
          rhoGW=getUPRESS(dataNetCDFK,t)
          rhogwInside=rhogwInside+rhoGW[:pro,boundaryLeft:boundaryRight,boundaryBottom:boundaryTop]

     rhogwInside=rhogwInside/(end-start)
    
     plt.figure((i*11)+100)
     plt.contourf(xxx[boundaryLeft:boundaryRight]*1e-3,zzc[:pro],rhogwInside[:,midmid,:],cmap=cm.jet, levels=100)
     plt.colorbar(label='$kg.m^{-2}.s^{-3}x10^{-8}$')
     plt.title('Rho gw vertical in relation to x at, t='+str(time[i])+' days')
     plt.xlabel('x (km)')
     plt.ylabel("z (m)")
     plt.show()

     plt.figure((i*13)+100)
     plt.contourf(xxx[boundaryLeft:boundaryRight]*1e-3,zzc[:pro],rhogwInside[:,:,midmid],cmap=cm.jet, levels=100)
     plt.colorbar(label='$kg.m^{-2}.s^{-3}x10^{-8}$')
     plt.title('Rho gw vertical in relation to y at, t='+str(time[i])+' days')
     plt.xlabel('y (km)')
     plt.ylabel("z (m)")
     plt.show()
     
     plt.figure((i*13)+110)
     plt.contourf(xxx[boundaryLeft:boundaryRight]*1e-3,zzc[:pro],ujk[:,:,midmid],cmap=cm.jet, levels=100)
     plt.colorbar(label='$m.s^{-1}$')
     plt.title('velocity v vertical in relation to y at, t='+str(time[i])+' days')
     plt.xlabel('y (km)')
     plt.ylabel("z (m)")
     plt.show()
    
     
     """
     zzzBellow=int(evoMeanBoundaryDepthLevel[j]-1)
    
     plt.figure((i*17)+100)
     plt.xlabel('x (km)')
     plt.ylabel("y (km)")
     CS=plt.contour(xxx[boundaryLeft:boundaryRight]*1e-3,xxx[boundaryLeft:boundaryRight]*1e-3,rhogwInside[zzzSurface,:,:],cmap='gray', levels=5)
     plt.clabel(CS,fmt='%2.1f',fontsize=9, inline=1,)
     plt.contourf(xxx[boundaryLeft:boundaryRight]*1e-3,xxx[boundaryLeft:boundaryRight]*1e-3,rhogwInside[zzzBellow,:,:],cmap=cm.jet, levels=100)
     plt.colorbar(label='$kg.m^{-2}.s^{-3}x10^{-8}$')
     plt.title('Rho gw at the surface(contours) z= '+str(int(evoMeanBoundaryDepth[i]))+'m (colors), t='+str(time[i])+' days')
     plt.legend()
     plt.show()
     """






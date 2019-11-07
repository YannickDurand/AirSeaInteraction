#!/usr/bin/env python
# coding: utf-8

#%%


#!/usr/bin/env python3

"""
Created on Fri Mar 15 2019

@author: ydurand
"""


#%%

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
tStart=7*4
HHH=0
dx=1200.
dy=dx
f=1e-4
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================
#
###############################################################################
#PATH3='/dsk1/ydurand/T03R01N3'
PATH3='/media/ydurand/Yannick/19_OAC/Control_Case/TNaNR01N3'
PATH2='/media/durand/Yannick/19_OAC/Control_Case/T03R01N2'
PATH1='/media/ydurand/Yannick/19_OAC/Control_Case/T03R01N1'

fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'


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

#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter2,time2=dimRes(dataNetCDFS3)
nbIter,time=getTime(dataNetCDFK3)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]


#%%

#%%    
###########################################################################
#-------------------Plots--------------------------------------------------
###########################################################################
planetaryVorticity=np.ones((30,329,329))*1e-4
   #print(ujk[i].shape)
#hik=diffH[:,:,j]
t=0
ujk,vik=getSpeed2(dataNetCDFS,t)
qij=vort_rel(dx, dy, ujk, vik)

vortAbs=qij#+planetaryVorticity
   #Calculation beforehand
#diver=conv(dx, dy, ujk[:,:,HHH,i], vik[:,:,HHH,i])*1E7
   #Ek=energyK(uck[:,:,i], vck[:,:,i], hik)
#enPot=energyP(dx, dy, zzd)
   #potEns=potentialEnstrophy(qij[:,:,j,i],hik[1:-2,1:-2])

#%%
#les plots

   #####-----Eta-----#####
   plt.figure(figsize=(12.5,10*Ly/Lx))
   plt.contour(xxx,yyy, hik, 20)
   plt.colorbar()
   plt.title('Eta' + ' ; t='+str(t))
   plt.show()
   
#####-----Voritcity------#####
plt.figure(100)
plt.contourf(xxx[:-1]*1e-3,xxx[:-1]*1e-3, vortAbs[0]/1E-4, cmap=cm.jet, vmax=0, levels=100)
plt.colorbar(label='$10^{-4}s^{-1}$')
plt.title('vorticite absolue' + ' ; t='+str(t))
plt.xlabel('x km')
plt.ylabel('y km')
plt.legend()
plt.show()
   
   #####-----divergence-----#####
   plt.figure(figsize=(12.5,10*Ly/Lx))
   plt.contour(xxx,yyy,diver,20)
   plt.colorbar()
   plt.title('Convergence' + ' ; t='+str(t))
   plt.show()
   #plt.close(901)
   
   #####-----Kinetic Energy-----#####
   plt.figure(figsize=(12.5,10*Ly/Lx))
   Ek=energyK(uck[i], vck[i], hik)
   plt.contour(xxx,yyy,Ek,20)
   plt.colorbar()
   plt.title('Kinetic Energy' + ' ; t='+str(t))
   plt.show()
   
   #####-----Potential Energy-----#####
   
   plt.figure(figsize=(12.5,10*Ly/Lx))
   plt.contour(xxx,yyy,energyP(dx, dy, hik),20)
   plt.colorbar()
   plt.title('Potential Energy' + ' ; t='+str(t))
   plt.show()
   
   #####-----Potential Enstrophy-----#####
   PE=potentialEnstrophy(qij, hik[1:-1,1:-1])
   plt.figure(figsize=(12.5,10*Ly/Lx))
   plt.contour(PE,20)
   plt.colorbar()
   plt.title('Potential Enstrophy ' + ' ; t='+str(t))
   plt.show()  
"""
   
   
"""
   #####-----Voritcity------#####
   qij=vort_pot(mesh, dx, dy, f, ujk[i], vik[i])
   if t==0 and zzt[j]>-0.1:
  qij=np.random.randn((res-1), (res-1))*1.E-2;
   figure1=plt.figure(1, figsize=(12.5,10*Ly/Lx))
   plt.contour(x,y, qij/1E-4, 20)
   plt.colorbar()
   plt.title('vorticite potentielle' + ' ; t='+str(t))
   figure1.savefig(str(vorPath)+'vorticity_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
   plt.close(1)

   #####-----divergence-----#####
   diver=conv(dx, dy, ujk[i], vik[i])
   if t==0 and zzt[j]>-0.1:
  diver=np.random.randn((250), (250))*1.E-3;
   figure2=plt.figure(2, figsize=(12.5,10*Ly/Lx))
   plt.contour(dx*mesh.xik,dy*mesh.yik,diver,20)
   plt.colorbar()
   plt.title('Divergence' + ' ; t='+str(t))
   figure2.savefig(str(divPath)+'divergence_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
   plt.close(2)

   #####-----Kinetic Energy-----#####
   figure3=plt.figure(3, figsize=(12.5,10*Ly/Lx))
   Ek=energyK(uck[i], vck[i], hik)
   if t==0 and zzt[j]>-0.1:
  Ek=np.random.randn((250), (250))*1.E-2;
   plt.contour(dx*mesh.xik,dy*mesh.yik,Ek,20)
   plt.colorbar()
   plt.title('Kinetic Energy' + ' ; t='+str(t))
   figure3.savefig(str(ekPath)+'EK_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
   plt.close(3)

   #####-----Potential Energy-----#####
   
   figure4=plt.figure(4, figsize=(12.5,10*Ly/Lx))
   plt.contour(dx*mesh.xik,dy*mesh.yik,energyP(dx, dy, hik),20)
   plt.colorbar()
   plt.title('Potential Energy' + ' ; t='+str(t))
   figure4.savefig(str(epPath)+'EP_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
   plt.close(4)
   
   #####-----Potential Enstrophy-----#####
   PE=potentialEnstrophy(qij, hik[1:-1,1:-1])
   figure5=plt.figure(5, figsize=(12.5,10*Ly/Lx))
   plt.contour(PE,20)
   plt.colorbar()
   plt.title('Potential Enstrophy ' + ' ; t='+str(t))
   figure5.savefig(str(ensPath)+'Ens_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
   plt.close(5)  
   """

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 16:14:25 2019

@author: ydurand
"""



import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.cm as cm


#import MITgcmutils as mit
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

from stockDiagKEsMean import *
from stockData import *
from netCDF4 import Dataset
from fonctionH import *
from plotKES import *

from fonctionH import *
from Operator import *
from fonctionRadius import *

#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=6*3600.                #en s
maxT=0.
alphat = 2e-4
rho0=999.8
tStart=5
f=1e-4
dx=1200.
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

# ###Without Wind
PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_wh_Wqt_Cyclone'
# =============================================================================
#
###############################################################################
#
# =============================================================================
###Different Wind

PATH='/home/durand/OAC/Control_Case/T03R01N3'
#PATH='/media/durand/Yannick/19_OAC/Control_Case/T03R01N3H60'
# =============================================================================



fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
fileTave='/diag_ocnTave.nc'

dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)
dataNetTave=(PATH+fileTave)
#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter2,time2=dimRes(dataNetCDFS)
nbIter,time=getTime(dataNetCDFK)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]
dz=np.ones((30,330,331))
dz2=np.ones((30,331,330))

dz=dz*zzc[:,np.newaxis,np.newaxis]
dz2=dz2*zzc[:,np.newaxis,np.newaxis]

#dz=np.diff(dz,axis=0)
#%%
evoSpeedRadius=np.zeros(nbIter)
evoCenterEddy=np.zeros((nbIter,2))
j=0

for i in range(0, nbIter):
    t=i
    evoSpeedRadius[j]=speedRadius(dataNetCDFS,i,dx)
    aaa=list(getCenterEta(dataNetCDFS,t))
    evoCenterEddy[j,0]=aaa[0]
    evoCenterEddy[j,1]=aaa[1]
    j=j+1
    
del aaa
diffy=int(evoCenterEddy[159,0]-evoCenterEddy[0,0])
diffx=int(evoCenterEddy[159,1]-evoCenterEddy[0,1])
#%%
#normeU=np.zeros((160,30,330,330))
start=0
ujk,vik=getUVVEL(dataNetTave)
wij=getWVEL(dataNetTave)
#%%
for i in range(160):
    print('coucou')
    #normeU[i]=(ucc(ujk[i])**2+vcc(vik[i])**2)**(1./2.)
#%%
ujk2=np.mean(ujk,axis=0)
dUdzFinal=partwz(ujk[159],dz)#
dUdzFinal2=partwz02(ujk2,dz**2,30,330,331)
dUdzInit=partwz(ujk[0],dz)#,30,330,331)

dVdzFinal=partwz(vik[159],dz2)#,30,331,330)
#dVdzInit=partwz(vik[0],dz2)#,30,331,330)
wik2=np.mean(wij,axis=0)

dwdz=partwz(wik2,dz[:,:,:-1])
dwdx=np.gradient(wik2,1200,axis=2)
dwdy=np.gradient(wik2,1200,axis=1)


termTest=(dwdx[:-1]+dwdy[:-1])#+dwdz)
#%%
start=100
end=160
wijPlot=np.mean(wij[start:end],axis=0)


#%%
vmax=0.5
vcont = np.linspace(-np.abs(vmax),np.abs(vmax),10)
plt.contourf(xxx,zzc[:],vik[0,:,int(evoCenterEddy[0,1]),:])
plt.colorbar()
plt.contour(xxx*1e-3,zzc,dVdzInit[:,int(evoCenterEddy[0,1]),:],levels=10,color='k',linewidths=1)
#%%
plt.figure(3)
plt.contourf(xxx,zzc[1:15],(dUdzFinal[1:15,:,int(evoCenterEddy[159,1])]/5*1e-5+termTest[1:15,:,int(evoCenterEddy[159,1])]))
plt.grid()
plt.colorbar()
plt.show()
plt.plot((dUdzFinal[:,int(evoCenterEddy[159,0])+int(evoSpeedRadius[159]/1200.),int(evoCenterEddy[159,1])]),zzc[:-1],label='t=40 days')

plt.plot([0,-0.003],[zzl[2],zzl[2]], linewidth=1, linestyle=':',color='k', label='OBL t=0 days')
plt.plot([0,-0.003],[zzl[6],zzl[6]], linewidth=1, linestyle='-.',color='k', label='OBL t=40 days')

plt.xlabel('$m.s^{-2}$')
plt.title('Vertical derivative of the horizontal speed')
plt.xlabel('$m.s^{-2}$')
plt.ylabel('depht (m)')
plt.legend()
plt.show()
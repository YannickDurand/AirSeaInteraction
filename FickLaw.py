#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 10:37:41 2019

@author: ydurand
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import scipy.ndimage as scp

import matplotlib.animation as animation


from stockDiagKEs import *
from fonctionH import *
from Operator import *
from Sauvegarde import *
import matplotlib.cm as cm


#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=10800.           #en s
maxT=0.
dh=1200.
alphaT = 2e-4
rho0=998.
coefDiff= 0.003
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

#PATH='/media/ydurand/Yannick/19_OAC/Control_Case/Case_withoutWqT_ACyclone'
#PATH='/home/ydurand/OAC/Case_wh_qT_withW_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qt_withW_ACyclone_npS'
PATH='/media/ydurand/Yannick/19_OAC/Control_Case/fplane/conf1/Case_withoutWqT_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_withoutWqT_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'

PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_withoutWqT_ACyclone'
PATH='/home/ydurand/OAC/Control_Case/fplane/conf1/Case_wh_qT_withW_ACyclone'
#PATH='/media/ydurand/Yannick/19_OAC/Control_Case/Case_withoutAtm_ACyclone'

PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'


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
xxx=(xxx/1000)

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

dz=np.ones((20,330,330))
dz=dz*zzc[:,np.newaxis,np.newaxis]
#dz=np.diff(dz,axis=0)
#%%
lapSnap=np.ones((4,20,330,330))
j=0
#Laplacian of the temperature field
for i in range(40,nbIter-1,50):
     t=i
     temp=getTemp(dataNetCDFS,i)

     dTdx=np.gradient(temp,dh,axis=2)
     dTdy=np.gradient(temp,dh,axis=1)
     
     #Horizontal Laplacian
     dT2dx=np.gradient(dTdx,dh,axis=2)
     dT2dy=np.gradient(dTdy,dh,axis=1)

     #Vertical Laplacian
     dTdz=partwz02(temp,dz,20,330,330)
     dT2dz=partwz02(dTdz,dz,20,330,330)

     
     dT2dv=dT2dx+dT2dy+dT2dz
     lapSnap[j,:,:,:]=dT2dv
     j=j+1
     print(j)
#%%
advSnap=np.ones((4,20,330,330))
j=0
#Compute of advection
for i in range(40,nbIter-1,50):
     t=i
     hhd,ujk,vik,wij,temp,eta=vitesse(dataNetCDFS,zzt,HHH=-1,tstart=i)

     dTdx=np.gradient(temp,dh,axis=2)
     dTdy=np.gradient(temp,dh,axis=1)
     
     dTdz=partwz02(temp,dz,20,330,330)
     
     uuu=ucc(ujk)
     vvv=vcc(vik)
     
     #Horizontal Advection
     advX=uuu[:,:,:]*dTdx
     advY=vvv[:,:,:]*dTdy
     
     #Vertical Advection
     advZ=wij[:,:,:]*dTdz

     
     advXYZ=advX+advY+advZ
     
     advSnap[j,:,:,:]=advXYZ
     j=j+1
     print(j)
#%%
dTdtSnap=np.ones((4,20,330,330))
j=0
#Compute of euleurienne time derivative
for i in range(40,nbIter-1,60):
     t=i
     temp=getTemp(dataNetCDFS,t)
     temp2=getTemp(dataNetCDFS,(t+1))
     
     dTdt=(temp2-temp)/(deltaT)
     dTdtSnap[j,:,:,:]=dTdt
     j=j+1
     print(j)
#%%
#Flick Law, find the coefficient. 
#Is there a constant coefficient???
DTDt=(dTdtSnap+advSnap)*1e7
diffTherm=lapSnap*coefDiff*1e7/3
diffTherm=np.true_divide(diffTherm, 4)
diff=DTDt-diffTherm  

#%%
###plot
plt.figure(100)
plt.contourf(xxx[:],zzc[:],diffTherm[1,:,165,:],cmap=cm.jet, levels=100)   
plt.colorbar()
plt.show()
#%%
plt.figure(101)
plt.contourf(xxx[:],xxx[:],lapSnap[3,2,:,:],cmap=cm.jet, levels=100)   
plt.colorbar()
plt.show()
#%%
plt.figure(103)
plt.contourf(xxx[:],zzc[:],diff[1,:,165,:],cmap=cm.jet, levels=100)   
plt.colorbar()
plt.show()
#%%

dt = deltaT

fig = plt.figure()
i=0
im = plt.imshow(DTDt[0,:,165,:], animated=True)
def updatefig(*args):
    global i
    if (i<99):
        i += 1
    else:
        i=0
    im.set_array(DTDt[i,:,165,:])
    return im,
ani = animation.FuncAnimation(fig, updatefig,  blit=True)
plt.show()


fig, ax = plt.subplots()
#ax.set(xlim=(xxx[0], xxx[-2]), ylim=(zzt[0], zzt[-2]))
cax = ax.contourf(xxx[:-2],zzc[:-2],diffTherm[0,:,165,:], cmap='jet',levels=100)
fig.colorbar(cax)
 
def animate(i):
    ax.collections = []
    ax.contourf(xxx[:-2],zzc[:-2],diffTherm[i,:,165,:], cmap='jet',levels=100)

    
anim = animation.FuncAnimation(fig, animate, frames=100)
 
plt.draw()
plt.show()
#%%
fig, ax = plt.subplots()
#ax.set(xlim=(xxx[0], xxx[-2]), ylim=(zzt[0], zzt[-2]))
cax = ax.contourf(xxx[:],zzc[:],DTDt[0,:,165,:]-diffTherm[0,:,165,:], cmap='jet',levels=100)
fig.colorbar(cax)
 
def animate(i):
    ax.collections = []
    ax.contourf(xxx[:],zzc[:],DTDt[i,:,165,:]-diffTherm[i,:,165,:], cmap='jet',levels=100)
    ax.set_title('Time (day): ' + str(i*(10800*2)/(3600*24)+5))
    
anim = animation.FuncAnimation(fig, animate, frames=50)

plt.draw()
plt.show()
#%%
"""
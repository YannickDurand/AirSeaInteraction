#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 16 14:02:19 2019

@author: ydurand
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 15:15:04 2019

@author: ydurand
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#import MITgcmutils as mit
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

from stockDiagKEs import *
from netCDF4 import Dataset
from fonctionH import *
from plotKES import *

import os
#%%

# =============================================================================
# ------------------Create Folder to stock-------------------------------------
# =============================================================================
PATH=                                                                          ###--- Choose of the path to stock
NAME=

os.makedirs((PATH+NAME),exist_ok=True)
Folder=(PATH+NAME)
#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=10800.             #en s
maxT=0.
alphat = 2e-4
rho0=998.
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================
#
###############################################################################
#
# =============================================================================
# ###Without Wind
PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_wh_Wqt_Cyclone'
# =============================================================================
#
###############################################################################
#
# =============================================================================
###Different Wind
#PATH='/media/durand/Yannick/19_OAC/Control_Case/fplane_diff_Wind/conf1/Case_wh_qT_withW_ACyclone_V7'
PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'
# =============================================================================


# =============================================================================
# ###Without Vortex
#PATH='/media/durand/Yannick/19_OAC/Control_Case/Case_withoutWqT_NoVortex'
# =============================================================================

# =============================================================================
# ###Without Atmosphere
#PATH='/media/durand/Yannick/19_OAC/Control_Case/Case_withoutAtm_ACyclone'
# =============================================================================


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

#%%
# =============================================================================
# ---------------------------Initial condition---------------------------------
# =============================================================================

utot_0, uadv_0, upress_0, u_eta_0, udissh_0, udissv_0, uext_0, u_ab_0 = meshcons(dataNetCDFK,-1,0)
totEk0=np.sum(utot_0)

#%%
# =============================================================================
# ---------------------------Evolution Kinetic Energy--------------------------
# =============================================================================
#%%
sumTotEK,evoEK,totFEXT,totUDISS,totUPRESS,totUAB=np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40)
totUDISSh,totUDISSv,totUETA=np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40)

j=0
for i in range(39,nbIter-1):
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=meshcons(dataNetCDFK,-1,i)
     utot=utot/86400.
     
     temp=getTemp(dataNetCDFS,i)
     rho=Rho(rho0,maxT,temp)
     #Calculation on the volume
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=utot*areaGrid, uadv*areaGrid, upress*areaGrid, u_eta*areaGrid, udissh*areaGrid, udissv*areaGrid, uext*areaGrid, u_ab*areaGrid
         
     #add rho
     t=i     
          #multiply by rho    add today
     utot=utot*rho+uadv*rho
     
     evoEK[j]=(np.sum(utot))
     
     #external Force
     totFEXT[j]=np.sum(uext[:,:,:]*rho)
     
     #dissipation
     totUDISSh[j]=np.sum(udissh*rho)
     totUDISSv[j]=np.sum(udissv*rho)
     totUDISS[j]=(np.sum(udissh*rho)+np.sum(udissv*rho))
     
     #force de pression
     totUPRESS[j]=np.sum(upress*rho)
     totUETA[j]=np.sum(u_eta*rho)

     #error Model
     totUAB[j]=np.sum(u_ab*rho)
     
     print(j)
     j=j+1
     

sumTotEK=np.cumsum(evoEK*deltaT)
sumTotDISS=np.cumsum(totUDISS*deltaT)
sumTotDISSModel=np.cumsum(totUAB*deltaT)
sumTotFEXT=np.cumsum(totFEXT*deltaT)
sumTotUPRESS=np.cumsum(totUPRESS*deltaT)
sumTotUETA=np.cumsum(totUETA*deltaT)
#%%
plt.figure(1004)
plt.plot(time[40:],sumTotEK, label = "Dt(EK)")
plt.plot(time[40:],sumTotFEXT, label = "Fext")
#plt.plot(time[40:],totUDISSh, label = 'horizontal dissipation')
plt.plot(time[40:],sumTotDISS, label = 'Total disspation')
#plt.plot(time[40:],totUDISSv, label ='vertical dissipation')
plt.plot(time[40:],sumTotUPRESS, label = 'Hydrostatic Vertical Pressure grad')
#plt.plot(time[40:],sumTotUETA, label = 'Pression Force Vertical')
plt.plot(time[40:],sumTotDISSModel, label = 'from error')

tototot=sumTotEK-sumTotFEXT-sumTotDISS-sumTotUPRESS
plt.plot(time[40:],tototot, label = 'Total')

plt.title("Evolution ok Kinetic Energy")
plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-3}$ ]")
plt.grid(True)
plt.legend()
plt.show()

#     plt.savefig(str(Folder)+str(nameA)+'_'+str(nameB)+'_'+str(nameC)+'s.png')

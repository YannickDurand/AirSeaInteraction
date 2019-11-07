#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:57:05 2019

@author: ydurand
"""

#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from stockDiagKEs import *
from stockData import *
from netCDF4 import Dataset
from fonctionH import *
import matplotlib.cm as cm

#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

PATH='/dsk1/ydurand/T03R01N0'

fileS='/state.nc'
grid='/grid.nc'

#file for Potential Energy
fileOcnTave='/diag_ocnTave.nc'
fileOcnSnap='/diag_ocnSnap.nc'
filePE='/tracer_wb.nc'

#destination of the file
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)

dataOcnTave=(PATH+fileOcnTave)
dataOcnSnap=(PATH+fileOcnSnap)
dataNetPE=(PATH+filePE)

# =============================================================================
# tropEddy=Dataset(dataNetCDFK)
# utot=np.array(tropEddy.variables['TOTKETEN'][:,:,:,:])
# =============================================================================

#%%
#-------------------------load grid-------------------------------------------#
xxx,si_x,nbIter,ttt=dimRes(dataNetCDFS)
si_y=si_x
ttt=ttt/(3600*24)

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]
si_z=dimZZ

hFacC, DRF = verticalGrid(dataGrid)

alphat = 2e-4
betas = 7.4e-4
ix = np.int(si_x/2)


#%%
flag_grid = 0
def comp_b (temp,salt):
  return alphat*temp + betas*salt
#%%
# =============================================================================
# =============================================================================
                         # # START YANNICK # #
# =============================================================================
# =============================================================================
#%%
#Initial State
maxT=16

t=100
temp=getTemp(dataNetCDFS,t)
rho0=999.8
rho=Rho(rho0,maxT,temp)
EPGrid=rho*9.81*areaGrid
EP_0=EPGrid
sumEP_0=np.sum(EP_0)

#%%
# =============================================================================
# ####################-----Time Evolution-----########################
# =============================================================================

totalEP,totalRHOGW  =  np.array([]),np.array([])
evoMeanTemp=np.array([])
     
#%%
def comp_EP(dataNetPE,t,ix):
     t=i
     dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys= dataDiagPE(dataNetPE,t)
     
     dtdt = dtdt/86400
     dsdt = dsdt/86400
     
     # t0 = np.where(t0 == 0,np.NaN,t0)
     # t1 = np.where(t1 == 0,np.NaN,t1)
          
     advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
     #advrs = np.append(advrs,advrs[None,0,:,:],axis=0)
     
     adv_at = (np.diff(advrt,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_at =  adv_at/areaGrid*0
     
     #adv_as = -( np.diff(advrs,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_as=0
     #adv_as = adv_as/(areaGrid)*0     

     ####################################
     
     dbdt = comp_b(dtdt,dsdt)
     advb = comp_b(adv_at,adv_as)
     tot=dbdt+advb
     dpedt = -areaGrid*9.81*999.8*(tot)*(zzc[:,np.newaxis,np.newaxis])


     return dpedt[:,:,ix]


#%%
# =============================================================================
# ########################Mean over five 1/f###################################
# =============================================================================
xx = np.arange(si_x)
tab=(len(zzt),len(xx))


filterf=5*18*3600
pdt=6*3600
t=80



pdtFilter=int(filterf/(pdt*2.))
wbFilter=np.zeros(tab)
epFilter=np.zeros(tab)

for i in range(t-pdtFilter,t+pdtFilter):
     wb2, wb = getWb(dataNetPE,i,ix)
     wb2, wb=wb2*9.81, wb*9.81
     wbFilter=np.add(wbFilter,wb)
     
     epFilter=comp_EP(dataNetPE,i,ix)

     
wbFilter=1/(pdtFilter*2.)*wbFilter
epFilter=1/(pdtFilter*2.)*epFilter

#%%
# =============================================================================
# #####-----Comparaison with W_b-----#####
# =============================================================================
vmin=-2e-5
vmax=2e-5
wbFilter = np.where(wbFilter<vmin,vmin,wbFilter)
wbFilter = np.where(wbFilter>vmax,vmax,wbFilter)
  

plt.figure(803)
xx = np.arange(si_x)
plt.contourf(xx,zzt[1:],wbFilter[1:,:],cmap=cm.jet, levels=100)
plt.colorbar(format='%.0e')
CS = plt.contour(xx,zzt[1:],epFilter[1:,:],colors='k')
plt.clabel(CS, fontsize=9, inline=1)
plt.show()











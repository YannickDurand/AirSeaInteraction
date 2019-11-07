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
from netCDF4 import Dataset
from fonctionH import *
import matplotlib.cm as cm

#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

#PATH='/media/ydurand/Yannick/19_OAC/Control_Case/Case_withoutWqT_ACyclone'
PATH='/home/ydurand/OAC/Case_wh_qT_withW_ACyclone'
#PATH='/media/ydurand/Yannick/19_OAC/Control_Case/fplane/conf1/Case_withoutWqT_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_ACyclone'
#PATH='/dsk1/ydurand/fplane/CONF1/Case_withoutWqT_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_withoutWqT_Cyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_Cyclone'


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

RAC,zzt,zzl,zzd,RC,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]
si_z=dimZZ

hFacC, DRF = verticalGrid(dataGrid)

alphat = 2e-4
betas = 7.4e-4
ix = np.int(si_x/2)



#%%
t=150
dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys= dataDiagPE(dataNetPE,t)
dtdt = dtdt/86400
dsdt = dsdt/86400

# t0 = np.where(t0 == 0,np.NaN,t0)
# t1 = np.where(t1 == 0,np.NaN,t1)



advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
advrs = np.append(advrs,advrs[None,0,:,:],axis=0)


adv_at = -( np.diff(advrt,axis=0) \
  + np.diff(advyt,axis=1)      \
  + np.diff(advxt,axis=2)     \
  /(areaGrid))


adv_as = -( np.diff(advrs,axis=0) \
  + np.diff(advyt,axis=1)      \
  + np.diff(advxt,axis=2))     \
  /(areaGrid)
#%%
flag_grid = 0
def comp_b (temp,salt):
  return alphat*temp + betas*salt

# b0 = comp_b(t0,s0)
# b1 = comp_b(t1,s1)

# pe0 = RC*b0
# pe1 = RC*b1
#%%
dbdt = comp_b (dtdt,dsdt)
advb = comp_b(adv_at, adv_as)
dpedt = RC[:,np.newaxis,np.newaxis]*(dbdt)#-advb)


def yzplot(psi,*args, **kwargs):
  
  psi = np.where(np.isnan(psi),0.,psi)
  vmax = np.max(np.abs((psi)))
  vmax = kwargs.get('vmax', vmax)
  vmin = -vmax
  psi = np.where(psi<vmin,vmin,psi)
  psi = np.where(psi>vmax,vmax,psi)
  
  title = kwargs.get('title',None)

  fgrid = kwargs.get('fgrid', 0)

  if fgrid:
    xx = YC[:,ix]*1e-3
    yy = RC[:,0,0]
  else:
    si_y,si_x = psi.shape
    xx = np.arange(si_x)
    yy = np.arange(si_y)
    
  plt.figure()
  plt.contourf(xx,zzt,psi,100,cmap=plt.cm.seismic)#,vmin=vmin,vmax=vmax,extend='both')
  plt.colorbar(format='%.0e')
  plt.title(title)
  if fgrid:
    plt.xlabel('x (km)')
    plt.ylabel('z (m)')

#psi = (pe1-pe0)/4500
#psi = psi[:,:,ix]
# yzplot(psi,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-5)

psi2 = dpedt[:,:,ix]
yzplot(psi2,title=r"tottend (m\,s$^{-2}$)")#,vmax=1e-5)

#psi3 = (psi - psi2)/psi2
# yzplot(psi3,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax = 1e-3)

#%%
psi4 = adv_at[:,:,ix]
yzplot(psi4,title='value of advection of T',vmax = 1e-4)

psi5 = dtdt[:,:,ix]
yzplot(psi5,title='variation of T (°C.s$^{-1})',vmax = 1e-4)

psi6 = psi5 - psi4
yzplot(psi6,vmax = 1e-4)
#%%
# =============================================================================
# =============================================================================
                         # # START YANNICK # #
# =============================================================================
# =============================================================================
#%%
#Initial State
maxT=16

t=0
temp=getTemp(dataNetCDFS,t)
rho0=998.
rho=Rho(rho0,maxT,temp)
EPGrid=rho*9.81*areaGrid
EP_0=EPGrid
sumEP_0=np.sum(EP_0)
#%%
evoPE,evoPE2=np.array([]),np.array([])

for i in range(nbIter-1):
     t=i
     dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys= dataDiagPE(dataNetPE,t)
     
     dtdt = dtdt/86400
     dsdt = dsdt/86400
     
     # t0 = np.where(t0 == 0,np.NaN,t0)
     # t1 = np.where(t1 == 0,np.NaN,t1)
     
     ix = np.int(si_x/2)
     
     advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
     advrs = np.append(advrs,advrs[None,0,:,:],axis=0)
     
     adv_at = -( np.diff(advrt,axis=0) \
       + np.diff(advyt,axis=1)      \
       + np.diff(advxt,axis=2)     \
       /(areaGrid))
     
     
     adv_as = -( np.diff(advrs,axis=0) \
       + np.diff(advyt,axis=1)      \
       + np.diff(advxt,axis=2))     \
       /(areaGrid)     
       
     dbdt = comp_b (dtdt,dsdt)
     advb = comp_b(adv_at, adv_as)
     dpedt = areaGrid*9.81*(dbdt-advb)*(zzt[:,np.newaxis,np.newaxis]-zzl[:,np.newaxis,np.newaxis])

     evoPE=np.append(evoPE,(np.sum(dpedt)))
     
     
     wb2, wb = getWb(dataNetPE,i)
     temp=np.multiply(wb2,dpedt)
          
     evoPE2=np.append(evoPE,(np.sum(temp)))
     
     
#%%     
plt.figure(901)
plt.plot(evoPE2)
plt.show()


#%%
# =============================================================================
# ####################-----Time Evolution-----########################
# =============================================================================

totalEP,totalRHOGW  =  np.array([]),np.array([])
evoMeanTemp=np.array([])
#%%     
for i in range(nbIter-1):
     t=i
     temp=getTemp(dataNetCDFS,t)
     rho0=998.
     rho=Rho(rho0,maxT,temp)
     
     wb2, wb = getWb(dataNetPE,i)
     rhogW=rho*9.81*wb2*areaGrid
     print('wb',np.sum(rhogW))
     totalRHOGW=np.append(totalRHOGW,np.sum(rhogW))

#%%
plt.figure(902)
plt.plot(totalRHOGW,label='$rho$gw')
plt.plot(evoPE2, label = 'D/Dt(b)-adv(b)')
plt.legend()
plt.show()

#%%
#####----- For the total Potential Energy -----#####

for i in range(nbIter-1):
     t=i
     temp=getTemp(dataNetCDFS,t)
     rho0=998.
     rho=Rho(rho0,maxT,temp)
     EPGrid=rho*9.81*areaGrid
     totalEP=np.append(totalEP,np.sum(EPGrid))#/sumEP_0)
     evoMeanTemp=np.append(evoMeanTemp,np.mean(temp))
#%%
plt.figure(903)
plt.plot(ttt[1:],totalEP)
plt.title('Evolution of the Potential Energy')
plt.show()

plt.figure(904)
plt.plot(evoMeanTemp)
plt.show()     
     
     
#%%
def comp_EP(dataNetPE,t,ix):
     dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys= dataDiagPE(dataNetPE,t)
     dtdt = dtdt/86400
     dsdt = dsdt/86400
     
     # t0 = np.where(t0 == 0,np.NaN,t0)
     # t1 = np.where(t1 == 0,np.NaN,t1)
     
     
     
     advrt = np.append(advrt,advrt[None,0,:],axis=0)
     advrs = np.append(advrs,advrs[None,0,:],axis=0)
     
     
     adv_at = -( np.diff(advrt,axis=0) \
       + np.diff(advyt,axis=1)      \
       + np.diff(advxt,axis=2)     \
       /(areaGrid))
     
     
     adv_as = -( np.diff(advrs,axis=0) \
       + np.diff(advyt,axis=1)      \
       + np.diff(advxt,axis=2))     \
       /(areaGrid)

     dbdt = comp_b (dtdt,dsdt)
     advb = comp_b(adv_at, adv_as)
     dpedt = RC[:,np.newaxis,np.newaxis]*(dbdt)#-advb)

     return dpedt[:,:,ix]


#%%
# =============================================================================
# ########################Mean over five 1/f###################################
# =============================================================================
xx = np.arange(si_x)
tab=(len(zzt),len(xx))


filterf=5*18*3600
pdt=3*3600
t=80



pdtFilter=int(filterf/(pdt*2.))
wbFilter=np.zeros(tab)
epFilter=np.zeros(tab)

for i in range(t-pdtFilter,t+pdtFilter):
     wb2, wb = getWb(dataNetPE,i,ix)
     wbFilter=np.add(wbFilter,wb)
     
     epFilter=np.add(epFilter,comp_EP(dataNetPE,i,ix))

     
wbFilter=1/(pdtFilter*2.)*wbFilter
epFilter=1/(pdtFilter*2.)*epFilter*1e8

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











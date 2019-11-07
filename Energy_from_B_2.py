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
from stockDiagKEs import *
from netCDF4 import Dataset
from fonctionH import *


#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================

#PATH='/media/ydurand/Yannick/19_OAC/Control_Case/Case_withoutWqT_ACyclone'
#PATH='/home/ydurand/OAC/Case_wh_qT_withW_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qt_withW_ACyclone_npS'
PATH='/media/ydurand/Yannick/19_OAC/Control_Case/fplane/conf1/Case_withoutWqT_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_withoutWqT_ACyclone'
PATH='/dsk1/ydurand/fplane/CONF1/Case_wh_qT_withW_Cyclone'

PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_wh_qT_withW_Cyclone'



fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)


file5=PATH+fileKEs
# =============================================================================
# tropEddy=Dataset(dataNetCDFK)
# utot=np.array(tropEddy.variables['TOTKETEN'][:,:,:,:])
# =============================================================================


#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter,time=dimRes(dataNetCDFS)
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

#%%
# =============================================================================
# #-----------------------Cell for choice----------------------------------------
# =============================================================================
flag_uv = 5 # 1: u , 2: v, 3:KEu, 4:KEv, 5:KEE
flag_grid = 1

#%%
# =============================================================================
# #%==================== LOAD FIELDS ===================================
# =============================================================================
# load grid
# =============================================================================
# if flag_grid:
#   XC    = mit.rdmds(dir0+'XC*')
#   YC    = mit.rdmds(dir0+'YC*')
#   XG    = mit.rdmds(dir0+'XG*')
#   YG    = mit.rdmds(dir0+'YG*')
#   DXC   = mit.rdmds(dir0+'DXC*')
#   DYC   = mit.rdmds(dir0+'DYC*')
#   hFacC = mit.rdmds(dir0+'hFacC*')
#   hFacS = mit.rdmds(dir0+'hFacS*')
#   hFacW = mit.rdmds(dir0+'hFacW*')
#   RAS   = mit.rdmds(dir0+'RAS*')
#   RAW   = mit.rdmds(dir0+'RAW*')
#   RAC   = mit.rdmds(dir0+'RAC*')
#   RAZ   = mit.rdmds(dir0+'RAZ*')
#   RC    = mit.rdmds(dir0+'RC*')
#   RF    = mit.rdmds(dir0+'RF*')
#   DRC   = mit.rdmds(dir0+'DRC*')
#   DRF   = mit.rdmds(dir0+'DRF*')
#   Depth = mit.rdmds(dir0+'Depth*')
# =============================================================================


if flag_uv == 1:
  filer = file1
elif flag_uv == 2:
  filer = file2
elif flag_uv == 3:
  filer = file3
elif flag_uv == 4:
  filer = file4
elif flag_uv == 5:
  filer = file5

i = 7

#iters1 = mit.mds.scanforfiles(dir0 + filer)
utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab =meshcons(dataNetCDFK,-1,i)

#%%
plt.close("all")
utot = utot /86400.

#si_y,si_x = utot.shape
ix = np.int(si_x/2)


def yzplot(psi,*args, **kwargs):
  
  vmax = np.max(np.abs((psi)))
  vmax = kwargs.get('vmax', vmax)
  vmin = -vmax
  psi = np.where(psi<vmin,vmin,psi)
  psi = np.where(psi>vmax,vmax,psi)
  
  title = kwargs.get('title',None)

  fgrid = kwargs.get('fgrid', 0)

  if fgrid==0:
    xx = YC[:,ix]*1e-3
    yy = RC[:,0,0]
  else:
    si_y,si_x = psi.shape
    xx = np.arange(si_x)
    yy = np.arange(si_y)
    
  plt.figure()
  plt.contourf(xx,zzc,psi,100,cmap=plt.cm.seismic,levels=256)#vmin=vmin,vmax=vmax,extend='both')
  plt.colorbar(format='%.0e')
  plt.title(title)
  if fgrid:
    plt.xlabel('x (km)')
    plt.ylabel('z (m)')
# =============================================================================
# 
# # =============================================================================
# psi = utot[:,:,ix]
# yzplot(psi,title=r"tottend (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-6)
# # 
# # 
# psi2 = uadv + upress + udissv + udissh + u_eta + u_ab + uext
# psi2 = psi2[:,:,ix]
# yzplot(psi,title=r"sum (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-6)
# # 
# error = np.abs(uadv) + np.abs(upress) + np.abs(udissv) + np.abs(udissh) + np.abs(u_eta) + np.abs(u_ab) + np.abs(uext)
# error2 = error[:,:,ix]/np.abs(utot[:,:,ix])
# # 
# psi3 = (psi2 - utot[:,:,ix])/utot[:,:,ix]/error2
# yzplot(psi3,title=r"relative error (m\,s$^{-2}$)",fgrid=flag_grid,vmax=1e-7)
# # =============================================================================
# =============================================================================

#%%
# =============================================================================
# ---------------------------Initial condition---------------------------------
# =============================================================================

utot_0, uadv_0, upress_0, u_eta_0, udissh_0, udissv_0, uext_0, u_ab_0 = meshcons(dataNetCDFK,-1,0)
totEk0=np.sum(utot_0)

#%%
# =============================================================================
# ---------------------------Evolution-----------------------------------------
# =============================================================================

sumTotEK,evoEK,totUEXT,totUDISS,totUPRESS=np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
totUDISSh,totUDISSv,totUETA=np.array([]),np.array([]),np.array([])

sumPrec=0
for i in range(nbIter-1):
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=meshcons(dataNetCDFK,-1,i)
     
     #Calculation on the volume
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=utot*areaGrid, uadv*areaGrid, upress*areaGrid, u_eta*areaGrid, udissh*areaGrid, udissv*areaGrid, uext*areaGrid, u_ab*areaGrid
     
     evoEK=np.append(evoEK,np.sum(utot)/86400.)
     
     """
     totUEXT=np.append(totUEXT,np.sum(uext))
     #dissipation
     totUDISS=np.append(totUDISS,np.sum(udissh)+np.sum(udissv))
     totUDISSh=np.append(totUDISSh,np.sum(udissh))
     totUDISSv=np.append(totUDISSv,np.sum(udissv))
     totUPRESS=np.append(totUPRESS,np.sum(upress))
     totUETA=np.append(totUETA,np.sum(u_eta))
     """
     
     #for the 
     varEK=evoEK[i]*18400.
     #print("1 = ",varEK)
     sumTotEK=np.append(sumTotEK,sumPrec+varEK)
     #print("2 = ",np.sum(sumTotEK))
     sumPrec=sumTotEK[i]
     
     
#%%
maxT=20
"""


#%%
t=1
ix = np.int(si_x/2)
dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys = dataDiagPE(dataNetPE,t)

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
"""
#%%
evoPE=np.array([])
alphat = 2e-4
ix = np.int(si_x/2)
maxT=16.

def comp_b (temp,salt):
  return alphat*temp #+ betas*salt

for i in range(nbIter-1):
     t=i
     #dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys = dataDiagPE(dataNetPE,t)
     dtdt,dsdt=getTOT(dataNetPE,t)
     dtdt = dtdt/86400.
     dsdt = dsdt/86400.
     """
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
     """  
     ###---From the Buayancy
     dbdt = comp_b (dtdt,dsdt)
     #advb = comp_b(adv_at, adv_as)
     dpedt = areaGrid*(dbdt)*9.81#-advb)
     
     ###---From the surface
     temp=getTemp(dataNetCDFS,t)
     temp=temp[0,:,:]
     rho0=998.
     rho=Rho(rho0,maxT,temp)
     eta=getEta(dataNetCDFS,t)
     eta2=getEta(dataNetCDFS,(t+1))
     
     deltaEta=(eta2-eta)*rho*9.81

     sumDeltaEta=np.sum(dpedt)+np.sum(deltaEta)     
     evoPE=np.append(evoPE,sumDeltaEta)

#%%
plt.figure(10)
plt.plot(evoEK, label = "Dt(EK)")
plt.plot(evoPE, label = "Dt(EP)")
plt.title("Evolution of Dt(EP) and Dt(EK) in fonction of time")
plt.legend()
plt.show()
#%%
# =============================================================================
# --------------------------Try filtering--------------------------------------
# =============================================================================

# =============================================================================
# totEK_f = savgol_filter(totEK, 25 , 2)
# totUEXT_f = savgol_filter(totUEXT, 25 , 2)
# totUDISS_f = savgol_filter(totUDISS, 25 , 2)
# totUPRESS_f = savgol_filter(totUPRESS, 25 , 2)
# 
# =============================================================================
#%%
# =============================================================================
# ------------------------Plot Plot Plots-------------------------------------
# =============================================================================
# =============================================================================
#plt.close("all")
# 
# =============================================================================
plt.figure(2)
plt.plot(totEK, label = "Dt(EK)")
plt.plot(totUEXT, label = "Fext")
plt.plot(totUDISS, label = 'Dissipation')
plt.plot(totUPRESS, label = 'Pression Force')
plt.plot(totUETA, label = 'from the surface')
plt.legend()
plt.show()
#%%
# =============================================================================
# plt.figure(3)
# plt.plot(totEK_f, label = "EK")
# plt.plot(totUEXT_f, label = "Fext")
# plt.plot(totUDISS_f, label = 'Dissipation')
# plt.plot(totUPRESS_f, label = 'Pression Force')
# plt.legend()
# plt.show()
# =============================================================================

plt.figure(4)
plt.plot(totUDISS, label = 'Total disspation')
plt.plot(totUDISSh, label='horizontal dissipation')
plt.plot(totUDISSv, label ='vertical dissipation')
plt.legend()
plt.show()

plt.figure(5)
plt.plot(totUPRESS, label="pression force")
plt.legend()
plt.show()
#%%
plt.figure(1)
plt.plot(sumTotEK, label="from diag")
plt.title('sum EC variation from t=0 to t=x')
plt.show()
#%%
# =============================================================================
# ---------------------------Evolution total----------------------------------
# =============================================================================

#for i in range (40):
     
     







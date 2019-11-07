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

#PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_withoutWqT_ACyclone'
#PATH='/home/durand/OAC/Control_Case/fplane/conf1/Case_wh_qT_withW_ACyclone'



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
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

#%%
# =============================================================================
# #%==================== LOAD FIELDS ===================================
# =============================================================================


t = 7

#iters1 = mit.mds.scanforfiles(dir0 + filer)
utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab =meshcons(dataNetCDFK,-1,t)

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
# ---------------------------Evolution Kinetic Energy--------------------------
# =============================================================================
#%%

sumTotEK,evoEK,totUEXT,totUDISS,totUPRESS,totUAB=np.array([]),np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
totUDISSh,totUDISSv,totUETA=np.array([]),np.array([]),np.array([])

sumTotDISS,sumTotDISSModel=np.array([]),np.array([])
sumTotFEXT=np.array([])
sumTotUPRESS=np.array([])

sumPrec,sumPrecDiss,sumPrecDissModel,sumPrecTotFEXT,sumPrecTotUPRESS=0,0,0,0,0
maxT=0.

for i in range(nbIter-1):
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=meshcons(dataNetCDFK,-1,i)
     
     #Calculation on the volume
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=utot*areaGrid, uadv*areaGrid, upress*areaGrid, u_eta*areaGrid, udissh*areaGrid, udissv*areaGrid, uext*areaGrid, u_ab*areaGrid
         
     #add rho
     t=i
     temp=getTemp(dataNetCDFS,t+1)
     rho0=998.
     rho=Rho(rho0,maxT,temp)
     
     #multiply by rho
     utot=utot*rho
     
     evoEK=np.append(evoEK,np.sum(utot)/86400.)
    
     totUEXT=np.append(totUEXT,np.sum(uext[0,:,:]*rho[0,:,:]))
     #dissipation
     totUDISSh=np.append(totUDISSh,np.sum(udissh*rho))
     totUDISSv=np.append(totUDISSv,np.sum(udissv*rho))
     totUDISS=np.append(totUDISS,np.sum(udissh*rho)+np.sum(udissv*rho))
     
     #force de pression
     totUPRESS=np.append(totUPRESS,np.sum(upress*rho))
     totUETA=np.append(totUETA,np.sum(u_eta*rho))

     #error Model
     totUAB=np.append(totUAB,np.sum(u_ab*rho))
     
     #for the time evolution
     varEK=evoEK[i]*10800.
     sumTotEK=np.append(sumTotEK,sumPrec+varEK)
     sumPrec=sumTotEK[i]
     
     varDISS=totUDISS[i]*10800.
     sumTotDISS=np.append(sumTotDISS,sumPrecDiss+varDISS)
     sumPrecDiss=sumTotDISS[i]
     
     varDISSModel=totUAB[i]*10800. ############ATTTENTION
     sumTotDISSModel=np.append(sumTotDISSModel,sumPrecDissModel+varDISSModel)
     sumPrecDissModel=sumTotDISSModel[i]
     
     varFEXT=totUEXT[i]*10800.
     sumTotFEXT=np.append(sumTotFEXT,sumPrecTotFEXT+varFEXT)
     sumPrecTotFEXT=sumTotFEXT[i]
     
     varUPRESS=totUPRESS[i]*10800
     sumTotUPRESS=np.append(sumTotUPRESS,sumPrecTotUPRESS+varUPRESS)
     sumPrecTotUPRESS=sumTotUPRESS[i]
     print(i)
#%%
# =============================================================================
# ---------------------------Evolution Potential Energy--------------------------
# =============================================================================
#%%
evoPE=np.array([])
evoPE2=np.array([])
evoWb2=np.array([])
evoDiffT=np.array([])

sumEvoPE=np.array([])
sumEvoWb=np.array([])


alphat = 2e-4
ix = np.int(si_x/2)
maxT=16.
sumprec=0
sumprecWb,sumPrecDiffT=0,0
rho0=998.


def comp_b (temp,salt):
  return alphat*temp #+ betas*salt

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
     
    
     
     adv_at = (np.diff(advrt,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_at =  adv_at/areaGrid
     
     adv_as=np.zeros((20,330,330))
     
     """
     adv_as = -( np.diff(advrs,axis=0) \
       + np.diff(advyt,axis=1)      \
       + np.diff(advxt,axis=2))     \
       /(areaGrid)     
     """
     ####################################
     
     dbdt = comp_b (dtdt,dsdt)
     advb = comp_b(adv_at,adv_as)
     tot=dbdt+advb
     dpedt = -areaGrid*9.81*998.*(tot)*(zzc[:,np.newaxis,np.newaxis])
     
     evoPE=np.append(evoPE,(np.sum(dpedt)))
     
     total=sumprec+evoPE[i]*10800
     sumEvoPE=np.append(sumEvoPE,total)
     sumprec=sumEvoPE[i]
          
     wb2, wb = getWb(dataNetPE,i)
     wb2=wb2*9.81*areaGrid
     evoWb2=np.append(evoWb2,(np.sum(wb2)))
     
     total2=sumprecWb+np.sum(wb2)*10800.
     sumEvoWb=np.append(sumEvoWb,total2)
     sumprecWb=sumEvoWb[i]     
     print(i)


#%%

totalEP,totalRHOGW  =  np.array([]),np.array([])
evoMeanTemp=np.array([])
evoRhoGW=np.array([])
sumTotalEPGrid=np.array([])

for i in range(nbIter-1):
     t=i
     temp=getTemp(dataNetCDFS,t)
     rho=Rho(rho0,maxT,temp)
     
     EPGrid=rho*9.81*areaGrid*(abs(zzc[:,np.newaxis,np.newaxis]))
     totalEP=np.append(totalEP,np.sum(EPGrid))#/sumEP_0)
     sumTotalEPGrid=np.append(sumTotalEPGrid,(totalEP[i]-totalEP[0]))   
     
     evoMeanTemp=np.append(evoMeanTemp,np.mean(temp))
     
     tempDY=(temp[:,-1:,:]-2*temp[:,:-1,:]+temp[:,1:,:])/(1200**2)
     
     
     #laplacianTEMP=
     
     wij=getDoubleV(dataNetCDFS,t)
     rhoGW=wij*9.81*rho*areaGrid
     evoRhoGW=np.append(evoRhoGW,(np.sum(rhoGW)))
#%%
# =============================================================================
# ##########--------------Difference -----#######################"# 
# =============================================================================
totalDEP=evoPE-evoWb2
sumTotalDEP=sumEvoWb+sumEvoPE

#Potential Energy
sumDiffEvoPE=sumEvoPE-sumTotalEPGrid
diffEvoPE=evoPE-totalEP

#total Energie
evoEnergieTotal=totalDEP+evoEK
sumEvoTot=sumTotEK+sumEvoPE#sumTotalEPGrid+sumEvoWb#

###disspation
evoEner=evoEnergieTotal+totUDISS
sumEvoEner=sumEvoTot-sumTotDISS

sumEvoEner2=sumEvoEner+sumTotDISSModel
evoEner2=evoEner+totUAB
evoEnerLiss=savgol_filter(evoEner2, 25 , 6)

sumEvoEner3=sumEvoEner2+sumEvoWb
evoEner3=evoEner2+evoWb2

sumEvoEner4=sumEvoEner3+sumTotFEXT
evoEner4=evoEner2+totUEXT

sumEvoEner5=sumTotUPRESS
#%%
# =============================================================================
# ------------------------Plot Plot Plots-------------------------------------
# =============================================================================
# =============================================================================
# plotComp(i,time,A,B,C,nameA,nameB,nameC)
# =============================================================================
zeros=np.zeros(240)
#%%
plotComp(900,time[10:],evoEK[10:],evoPE[10:],evoEnergieTotal[10:],"$D_t(KE)$","$D_t(PE)$","$D_t(TE)$")
#%%
plotComp(901,time[10:],evoEnergieTotal[10:],totUDISS[10:],evoEner[10:],"$D_t(TE)$","$Dissipation$","$D_t(ET)=D_t(TE)-D_t(Diss)$")

#%%
plotComp(901,time,evoEner2,evoWb2,evoEner3,"$D_t(ET)$","$D_t(rhogw)$","$D_t(ET)-D_t(rhogw)$")

#%%
plotComp(914,time,evoPE,evoWb2,totalDEP,"$D_t(EP)$","$RhoWG$","$D_t(EP)-RhoWG$")


#%%
#%%
plotComp(911,time,sumEvoPE,totalEP-totalEP[0],sumTotalDEP,"$\sum{D_tEP_i x dt}$","$EP(t)-EP(t=0)$","$\sum{D_tEP_i x dt}-rhogw$")

# =============================================================================
# TEST TO HAVE ZERO
# =============================================================================
######
## Here we try to find something near zero for total evolution
######
#%%
plotCompT(912,time[:],sumEvoPE,sumTotalEPGrid,sumEvoWb,"$PE(t)-PE(t=0)$","$PEG(t)-PEG(t=0)$","$rhogw(t)-rhogw(t=0)")
plotCompT(911,time[:],sumTotEK,sumEvoPE,sumEvoTot,"$KE(t)-KE(t=0)$","$PE(t)-PE(t=0)$","$ME(t)-ME(t=0)$")

#%%
plotCompT(913,time,sumEvoTot,sumTotDISS,sumEvoEner,"$ME(t)-ME(t=0)$","$KEdiss(t)-KEdiss(t=0)$","$TO(t)-TO(t=0)$")

#%%
plotCompT(915,time,sumEvoEner,sumTotDISSModel,sumEvoEner2,"$TO(t)-TO(t=0)$","$KEdiss(t)-KEdiss(t=0)$","$TO_M(t)-TO_M(t=0)$")

#%%
plotCompT(916,time,sumEvoEner2,sumEvoWb,sumEvoEner3,"$TO_M(t)-TO_M(t=0)$","$rhogw$","$diff$")

#%%
plotCompT(917,time,sumEvoEner3,sumTotFEXT,sumEvoEner4,"$TO_M(t)-TO_M(t=0)$","$Fext(t)-Fext(t=0)$","$diff$")

#%%
plotCompT(917,time,sumEvoEner4,sumTotUPRESS,sumEvoEner5,"$T-Fext","$div(up)$","$diff$")

#%%
# =============================================================================
# #############---------------Transfer Analysis-----------###########"
# =============================================================================
#%%
transEvo=evoWb2+totUPRESS
plotComp(2001,time,evoWb2,totUPRESS,transEvo,'Tg*wb2 from tracerWB','rho*UPRESS from diagKE','difference between')

#%%
plotComp(2002,time,evoEK,evoPE,evoWb2,'$D_t(KE)$','$D_t(PE)$','$g*rho*w$')

#%%
# =============================================================================
# #############---------------Kinetic Energy-----------###########"
# =============================================================================
#%%
plt.figure(1004)
plt.plot(time[1:],evoEK, label = "Dt(EK)")
plt.plot(time[1:],totUEXT, label = "Fext")
plt.plot(time[1:],totUDISSh, label = 'horizontal dissipation')
plt.plot(time[1:],totUDISS, label = 'Total disspation')
plt.plot(time[1:],totUDISSv, label ='vertical dissipation')
plt.plot(time[1:],totUPRESS, label = 'Hydrostatic Pressure grad')
plt.plot(time[1:],totUETA, label = 'Pression Force Vertical')
plt.plot(time[1:],totUAB, label = 'from error')

plt.title("Evolution ok Kinetic Energy")
plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-3}$ ]")
plt.legend()
plt.show()

#%%
plotComp(1006,time,totUDISS,totUDISSh,totUDISSv,'Total disspation','horizontal dissipation','vertical dissipation')

#%%
plotComp(1006,time,totUPRESS,totUETA,totUETA+totUPRESS,'Pression Force Horizontal','Pression Force Vertical','diff')
#%%
plt.figure(1007)
plt.plot(time[1:],totUEXT, label = "Fext")
plt.plot(time[1:],totUETA, label = 'UETA')

plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-3}$ ]")
plt.title("Loose of Kinetic Energy by the UETA")
plt.legend()
plt.show()

#%%
plt.figure(1008)
plt.plot(time[1:],totUEXT, label = "Fext")

plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-3}$ ]")
plt.title("Part of Kinetic Energy by External Forcing")
plt.legend()
plt.show()
#%%
###########################
###Try to separate by order
###############################################################################
#%%
###Operation
totEKWHT=evoEK-totUPRESS ###3004
totEKWHD=totEKWHT-totUDISS ###3005
totEKWSU=totEKWHD- totUETA ###3006
NoneEner=totEKWSU -  totUEXT###3007

###############################################################################
#%%
plotComp(3004,time,evoEK,totUPRESS,totEKWHT,'$D_t(EK)$','Hydrostatic Pressure grad','$D_t(EK)-Hydrostatic Pressure grad')

#%%
plotComp(3005,time,totEKWHT,totUDISS,totEKWHD,'$D_t(EK)-D_t(rhogw)$','Total disspation','$Dt_(EK)-D_t(rhogw)-D_t(diss)$')

#%%
plotComp(3006,time,totEKWHD,totUETA,totEKWSU,'$D_t(EK)$-Pression Force-Total disspation','from the surface','$D_t(EK)-rhogw$-Total disspation-from Surface')

#%%
plotComp(3008,time,totEKWSU,totUEXT,NoneEner,'$D_t(EK)$-Pression Force-Total disspation-from Surface','Fext','$D_t(EK)-rhogw$-Total disspation-Fext')

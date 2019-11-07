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
rho0=999.8
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
# ---------------------------Evolution Kinetic Energy--------------------------
# =============================================================================
#%%
sumTotEK,evoEK,totFEXT,totUDISS,totUPRESS,totUAB=np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40)
totUDISSh,totUDISSv,totUETA=np.zeros(nbIter-40),np.zeros(nbIter-40),np.zeros(nbIter-40)

j=0
for i in range(39,nbIter-1):
     t=i
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=meshcons(dataNetCDFK,i)
     utot=utot/86400.
     temp=getTemp(dataNetCDFS,t)
     rho=Rho(rho0,maxT,temp)
     #Calculation on the volume
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=utot*areaGrid, uadv*areaGrid, upress*areaGrid, u_eta*areaGrid, udissh*areaGrid, udissv*areaGrid, uext*areaGrid, u_ab*areaGrid
         
     #add rho
     t=i     
          #multiply by rho    add today
     utot=utot*rho0+uadv*rho0
     
     evoEK[j]=(np.sum(utot))
     
     #external Force
     totFEXT[j]=np.sum(uext[:,:,:]*rho0)
     
     #dissipation
     totUDISSh[j]=np.sum(udissh*rho0)
     totUDISSv[j]=np.sum(udissv*rho0)
     totUDISS[j]=(np.sum(udissh*rho0)+np.sum(udissv*rho0))
     
     #force de pression
     totUPRESS[j]=np.sum(upress*rho0)
     totUETA[j]=np.sum(u_eta*rho0)

     #error Model
     totUAB[j]=np.sum(u_ab*rho0)
     
     print(j)
     j=j+1
     

sumTotEK=np.cumsum(evoEK*deltaT)
sumTotDISS=np.cumsum(totUDISS*deltaT)
sumTotDISSModel=np.cumsum(totUAB*deltaT)
sumTotFEXT=np.cumsum(totFEXT*deltaT)
sumTotUPRESS=np.cumsum(totUPRESS*deltaT)
sumTotUETA=np.cumsum(totUETA*deltaT)

#%%
# =============================================================================
# ---------------------------Evolution Potential Energy--------------------------
# =============================================================================
#%%
evoPE=np.array([])
evoPE2=np.array([])
evoWb2=np.array([])
evoDiffT=np.array([])
ix = np.int(si_x/2)

def comp_b (temp,salt):
  return alphat*temp #+ betas*salt

for i in range(39,nbIter-1):
     t=i
     dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys= dataDiagPE(dataNetPE,t)
     
     dtdt = dtdt/86400
     dsdt = dsdt/86400
     
     # t0 = np.where(t0 == 0,np.NaN,t0)
     # t1 = np.where(t1 == 0,np.NaN,t1)
          
     advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
     advrs = np.append(advrs,advrs[None,0,:,:],axis=0)
     
     adv_at = (np.diff(advrt,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_at =  adv_at/areaGrid*0
     
     adv_as = -( np.diff(advrs,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_as = adv_as/(areaGrid)*0     

     ####################################
     
     dbdt = comp_b (dtdt,dsdt)
     advb = comp_b(adv_at,adv_as)
     tot=dbdt+advb
     dpedt = -areaGrid*9.81*998.*(tot)*(zzc[:,np.newaxis,np.newaxis])
     
     evoPE=np.append(evoPE,(np.sum(dpedt)))
          
     wb2, wb = getWb(dataNetPE,i)
     wb2=wb2*9.81*areaGrid
     evoWb2=np.append(evoWb2,(np.sum(wb2)))

sumEvoPE=np.cumsum(evoPE*deltaT)
sumEvoWb=np.cumsum(evoWb2*deltaT)

#%%

totalEP,totalRHOGW  =  np.array([]),np.array([])
evoMeanTemp=np.array([])
evoRhoGW=np.array([])
sumTotalEPGrid=np.array([])

for i in range(40,nbIter):
     t=i
     temp=getTemp(dataNetCDFS,t)
     rho=Rho(rho0,maxT,temp)
     
     EPGrid=rho*9.81*areaGrid*((zzc[:,np.newaxis,np.newaxis]))
     totalEP=np.append(totalEP,np.sum(EPGrid))     
     evoMeanTemp=np.append(evoMeanTemp,np.mean(temp))

     wij=getDoubleV(dataNetCDFS,t)
     rhoGW=wij*9.81*rho*areaGrid
     evoRhoGW=np.append(evoRhoGW,(np.sum(rhoGW)))

sumTotalEPGrid=(totalEP-totalEP[0])
sumEvoRhoGW=np.cumsum(evoRhoGW*deltaT)
#%%

# =============================================================================
# ##########--------------Difference -----#######################"# 
# =============================================================================
sumDiffPEEvo=sumEvoPE-sumTotalEPGrid                                           # Here sumEvoPE is only the variation of d(rho)/dt

sumDiffPEEvo2=sumDiffPEEvo-sumEvoWb                                            # Here we add (rhogw)

sumDiffPEEvo2a=sumDiffPEEvo2-sumEvoRhoGW

#sumDiffPEEvo3=sumDiffPEEvo2-sumEvoDiffTerm
####    We have almost the same thing but one parameter diffe
####    Differ by 1.3e10 at the end
####    Where this Oscillation come from?    

####     ----     For Total Energy     ---- ####
sumEvoTot=sumEvoPE+sumTotEK
sumEvoTot2=sumEvoTot-sumTotDISS-sumTotDISSModel
sumEvoTot3=sumEvoTot2-sumTotFEXT
#sumEvoTot4=sumEvoTot3-sumEvoDiffTerm
#sumEvoTot5=sumEvoTot4-sumTotUETA

#%%
# =============================================================================
# ------------------------Plot Plot Plots-------------------------------------
# =============================================================================
#to save ==> plt.savefig(str(ekPath)+'EK_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
# =============================================================================
# plotComp(i,time,A,B,C,nameA,nameB,nameC)
# =============================================================================

zeros=np.zeros(240)
plotComp(915,time[39:],sumEvoPE,sumTotalEPGrid,sumDiffPEEvo,"$ alpha D_t(T)$","PE(t=t)-PE(t=0)","Difference")
plotCompSauv(915,time[39:],sumEvoPE,sumTotalEPGrid,sumDiffPEEvo,"$ alpha D_t(T)$","PE(t=t)-PE(t=0)","Difference",Folder)



#%%
plotComp(917,time[39:],sumEvoPE,sumTotEK,sumEvoTot,"$PE(t=t)-PE(t=0)","KE(t=t)-KE(t=0)","Difference")


#%%8.5
# =============================================================================
# TEST TO HAVE ZERO
# =============================================================================
######
## Here we try to find something near zero for total evolution
######
#%%
plotCompT(911,time[39:],sumTotEK,sumEvoPE,sumEvoTot,"$KE(t)-KE(t=0)$","$PE(t)-PE(t=0)$","$ME(t)-ME(t=0)$")
plotComp(915,time[39:],sumEvoTot,sumTotDISS,sumEvoTot2,"$ME(t)-ME(t=0)$","Dissipation","new==>$ME(t)-ME(t=0)$")
plotComp(917,time[39:],sumEvoTot2,sumTotFEXT,sumEvoTot3,"$ME(t)-ME(t=0)$","External Forcing","new==>$ME(t)-ME(t=0)$")
#plotComp(919,time[39:],sumEvoTot3,sumEvoDiffTerm,sumEvoTot4,"$ME(t)-ME(t=0)$","Therm Diffusion","new==>$ME(t)-ME(t=0)$")



# =============================================================================
# #############---------------Transfer Analysis-----------###########"
# =============================================================================
#%%
transEvo=evoRhoGW+totUPRESS
sumTransEvo=sumTotUPRESS+sumEvoRhoGW
plotComp(2001,time[39:],evoWb2,totUPRESS,transEvo,'rhog*w from PE','rho*g*w from KE','difference between')
plotCompT(2003,time[39:],sumEvoWb,sumTotUPRESS,sumTransEvo,'rhog*w from PE','rho*g*w from KE','difference between')


#%%
plotComp(2002,time,evoEK,evoPE,evoWb2,'$D_t(KE)$','$D_t(PE)$','$g*rho*w$')

#%%
# =============================================================================
# #############---------------Kinetic Energy-----------###########"
# =============================================================================
#%%
plt.figure(1004)
plt.plot(time[40:],sumTotEK, label = "Dt(EK)")
plt.plot(time[40:],sumTotFEXT, label = "Fext")
#plt.plot(time[40:],totUDISSh, label = 'horizontal dissipation')
plt.plot(time[40:],sumTotDISS, label = 'Total disspation')
#plt.plot(time[40:],totUDISSv, label ='vertical dissipation')
plt.plot(time[40:],sumTotUPRESS, label = 'Hydrostatic Pressure grad')
#plt.plot(time[40:],totUETA, label = 'Pression Force Vertical')
#plt.plot(time[40:],totUAB, label = 'from error')

tototot=sumTotEK+sumTotFEXT+sumTotDISS+sumTotUPRESS+sumTotUAB
plt.plot(time[40:],tototot, label = 'Hydrostatic Pressure grad')

plt.title("Evolution ok Kinetic Energy")
plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-3}$ ]")
plt.legend()
plt.show()

#%%
plotComp(1006,time,totUDISS,totUDISSh,totUDISSv,'Total disspation','horizontal dissipation','vertical dissipation')

#%%
#%%
plt.figure(1007)
plt.plot(time[1:],totUEXT, label = "Fext")
plt.plot(time[1:],totUDISS, label = 'UDISS')

plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-3}$ ]")
plt.title("External Forcing and Dissipation")
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

WhereItGo=totUDISS+totUEXT
###############################################################################
#%%
plotComp(3004,time[39:],evoEK,totUPRESS,totEKWHT,'$D_t(EK)$','Hydrostatic Pressure grad','$D_t(EK)-Hydrostatic Pressure grad')

#%%
plotComp(3005,time[39:],totEKWHT,totUDISS,totEKWHD,'$D_t(EK)-D_t(rhogw)$','Total disspation','$Dt_(EK)-D_t(rhogw)-D_t(diss)$')

#%%
plotComp(3006,time[39:],totEKWHD,totUETA,totEKWSU,'$D_t(EK)$-Pression Force-Total disspation','from the surface','$D_t(EK)-rhogw$-Total disspation-from Surface')

#%%
plotComp(3008,time[39:],totEKWSU,totUEXT,NoneEner,'$D_t(EK)$-Pression Force-Total disspation-from Surface','Fext','$D_t(EK)-rhogw$-Total disspation-Fext')

#%%
plotComp(5000,time[39:],sumEvoRhoGW,sumTotEK,sumEvoRhoGW-sumTotEK,'sumEvoRhoGW','sum Evo KE','difference')

#%%
"""
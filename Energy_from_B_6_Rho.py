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
import scipy.stats

from stockDiagKEs import *
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
tStart=0
#%%
# =============================================================================
# #-------------------------Cell for NetCDF------------------------------------
# =============================================================================
#
###############################################################################

PATH='/home/durand/OAC/Control_Case/T03R01N3_60km'

fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'
fileTave='/diag_ocnTave.nc'
filePhy='/phiHyd.nc'

dataGrid=(PATH+grid)
dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataNetPE=(PATH+filePE)
dataNetTave=(PATH+fileTave)
dataNetPhi=(PATH+filePhy)

#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter2,time2=dimRes(dataNetCDFS)
nbIter,time=getTime(dataNetCDFK)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]

"""
#%%
# =============================================================================
# ---------------------------Evolution Kinetic Energy--------------------------
# =============================================================================
#%%
sumTotEK,evoEK,totFEXT,totUDISS,totUPRESS,totUAB=np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart)
totUDISSh,totUDISSv,totUETA=np.zeros(nbIter-tStart),np.zeros(nbIter-tStart),np.zeros(nbIter-tStart)

j=0
for i in range(tStart,nbIter):
     t=i
     utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab=meshcons(dataNetCDFK,i)
     utot=utot/86400.
     temp=getTemp(dataNetCDFS,t+1)
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
     
     print(i)
     j=j+1
del utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab
     

sumTotEK=np.cumsum(evoEK*deltaT)
sumTotDISS=np.cumsum(totUDISS*deltaT)
sumTotDISSh=np.cumsum(totUDISSh*deltaT)
sumTotDISSv=np.cumsum(totUDISSv*deltaT)
sumTotDISSModel=np.cumsum(totUAB*deltaT)
sumTotFEXT=np.cumsum(totFEXT*deltaT)
sumTotUPRESS=np.cumsum(totUPRESS*deltaT)
sumTotUETA=np.cumsum(totUETA*deltaT)
#%%

plot(9001,time[tStart:],sumTotEK,"Kinetic Energy N3","Kinetic Energy")
plot(9002,time[tStart:],sumTotUPRESS,"rhogw N3","Rho gw")
plot(9003,time[tStart:],sumTotDISS,"dissipation N3","Dissipation")


#%%
# =============================================================================
# #Regression linear
# =============================================================================

absComp=np.arange(tStart,nbIter)-tStart/4

#Kinetic Energy
lr = scipy.stats.linregress(time[tStart:],sumTotEK)
coeffN1=lr[0]
trendN1=coeffN1*(absComp/4.)
plot(9001,time[tStart:],trendN1,"Trend Kinetic Energy N3","Kinetic Energy")

#rhogw
lr2 = scipy.stats.linregress(time[tStart:],sumTotUPRESS)
coeffN1=lr2[0]
trendN1=coeffN1*(absComp/4.)
plot(9002,time[tStart:],trendN1,"Trend rhogw N3","Rho gw")

#dissipation
lr3 = scipy.stats.linregress(time[tStart:],sumTotDISS)
coeffN1=lr3[0]
trendN1=coeffN1*(absComp/4.)
plot(9003,time[tStart:],trendN1,"Trend dissipation N3","Dissipation")

#external Forcing
lr4=scipy.stats.linregress(time[tStart:],sumTotFEXT)
"""
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

for i in range(tStart,nbIter):
     t=i
     dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys= dataDiagPE(dataNetPE,t)
     
     dtdt = dtdt/86400
     dsdt = dsdt/86400
     
     # t0 = np.where(t0 == 0,np.NaN,t0)
     # t1 = np.where(t1 == 0,np.NaN,t1)
          
     advrt = np.append(advrt,advrt[None,0,:,:],axis=0)
     #advrs = np.append(advrs,advrs[None,0,:,:],axis=0)
     
     adv_at = (np.diff(advrt,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_at =  adv_at/areaGrid
     
     #adv_as = -( np.diff(advrs,axis=0)+ np.diff(advyt,axis=1)+ np.diff(advxt,axis=2))
     adv_as=0
     #adv_as = adv_as/(areaGrid)*0     

     ####################################
     
     dbdt = comp_b(dtdt,dsdt)
     advb = comp_b(adv_at,adv_as)
     #tot=dbdt-advb
     dpedt = -areaGrid*9.81*999.8*(dbdt)*(zzc[:,np.newaxis,np.newaxis])
     dpedt2 = -areaGrid*9.81*999.8*(advb)*(zzc[:,np.newaxis,np.newaxis])
     
     evoPE=np.append(evoPE,(np.sum(dpedt)))
     evoPE2=np.append(evoPE2,(np.sum(dpedt2)))
          
     wb2, wb = getWb(dataNetPE,i)
     wb2=wb2*9.81*areaGrid
     evoWb2=np.append(evoWb2,(np.sum(wb2)))
     print(i)
del dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys,dbdt, advb, dpedt, wb, wb2

sumEvoPE=np.cumsum(evoPE[100:]*deltaT)
sumEvoPE2=np.cumsum(evoPE2[100:]*deltaT)

sumEvoWb=np.cumsum(evoWb2[100:]*deltaT)

#%%

totalEP,totalEPMean,totalEP2,totalEP3,totalRHOGW  =  np.array([]),np.array([]),np.array([]),np.array([]),np.array([])
totalRHOAnoma=np.array([])
totalRHOAnoma2=np.array([])
totalRHOAnoma3=np.array([])

evoEta=np.array([])
evoRhoGW=np.array([])
sumTotalEPGrid=np.array([])

rhoMean=np.zeros((30,330,330))

for i in range(tStart+1,nbIter+1):
     t=i
     temp=getTemp(dataNetTave,t-1)
     rho=Rho(rho0,maxT,temp)
     eta=getEta(dataNetCDFS,t)
     eta=eta*rho[0]*9.81*(1200**2)
     evoEta=np.append(evoEta,np.sum(eta))
     
     EPGrid=rho*9.81*areaGrid*((zzc[:,np.newaxis,np.newaxis]))
     tempo=(np.sum(EPGrid)+evoEta[i-1])
     totalEP=np.append(totalEP,tempo)
     
     RHOAnoma=getRHOAnoma(dataNetTave,t-1)
     RHOAnoma2=np.where(RHOAnoma<0,RHOAnoma,0)
     RHOAnoma3=np.where(RHOAnoma<=0.0001,0,EPGrid)
     RHOAnoma3=np.where(RHOAnoma>=0.0001,0,EPGrid)
          
     RHOAnoma01=RHOAnoma*9.81*areaGrid*((zzc[:,np.newaxis,np.newaxis]))
     RHOAnoma21=RHOAnoma2*9.81*areaGrid*((zzc[:,np.newaxis,np.newaxis]))
     
     totalRHOAnoma=np.append(totalRHOAnoma,np.sum(RHOAnoma01))
     totalRHOAnoma2=np.append(totalRHOAnoma2,np.sum(RHOAnoma21))
     totalRHOAnoma3=np.append(totalRHOAnoma3,np.sum(EPGrid))
     
     wij=getDoubleV(dataNetCDFS,t)
     rhoGW=wij*9.81*rho*areaGrid
     evoRhoGW=np.append(evoRhoGW,(np.sum(rhoGW)))
     print(i)

sumTotalEPGrid=(totalEP[100:]-totalEP[100])
sumTotalEPGridMean=(totalEPMean[100:]-totalEPMean[100])

sumEvoRhoGW=np.cumsum(evoRhoGW[100:]*deltaT)
sumRHOAnoma=(totalRHOAnoma[100:]-totalRHOAnoma[100])
sumRHOAnoma2=(totalRHOAnoma2[100:]-totalRHOAnoma2[100])
sumRHOAnoma3=(totalRHOAnoma3[100:]-totalRHOAnoma3[100])

#%%
"""
# =============================================================================
# ##########--------------Difference -----#######################"# 
# =============================================================================
sumDiffPEEvo=sumEvoPE+sumEvoWb-sumTotalEPGrid                                  # Here sumEvoPE is only the variation of d(rho)/dt

sumDiffPEEvo2=sumDiffPEEvo-sumEvoWb                                            # Here we add (rhogw)

sumDiffPEEvo2a=sumDiffPEEvo2-sumEvoRhoGW

#sumDiffPEEvo3=sumDiffPEEvo2-sumEvoDiffTerm
####    We have almost the same thing but one parameter diffe
####    Differ by 1.3e10 at the end
####    Where this Oscillation come from?    

####     ----     For Total Energy     ---- ####
sumEvoTot=sumEvoPE+sumTotEK
sumEvoTot2=sumEvoTot-sumTotDISS-sumTotDISSModel
sumEvoTot3=sumEvoPE+sumTotEK-sumTotFEXT-sumTotDISS+sumTotDISSModel
#sumEvoTot4=sumEvoTot3-sumEvoDiffTerm
#sumEvoTot5=sumEvoTot4-sumTotUETA
#%%
# =============================================================================
# ------------------------Plot Plot Plots-------------------------------------
# =============================================================================
# =============================================================================
# plotComp(i,time,A,B,C,nameA,nameB,nameC)
# =============================================================================
zeros=np.zeros(240)
plotCompT(915,time[4:],sumEvoPE,sumTotalEPGrid,sumEvoWb,"$ alpha D_t(T)$","PE(t=t)-PE(t=0)","rhogw",'Evolution of potential energy')
plotCompT(916,time[4:],sumEvoPE+sumEvoWb,sumTotalEPGrid,sumDiffPEEvo,"$ alpha D_t(T) + rhogw $","PE(t=t)-PE(t=0)","sumDiffPEEvo",'Evolution of potential energy')



#%%
plotCompT(917,time[4:],sumEvoPE,sumTotEK[:],sumEvoTot,"$PE(t=t)-PE(t=0)$","KE(t=t)-KE(t=0)","$ME(t)-ME(t=0)$","Evolution of the Kinetic, Potential and Mecanical Energy")

#%%8.5
# =============================================================================
# TEST TO HAVE ZERO
# =============================================================================
######
## Here we try to find something near zero for total evolution
######
#%%
plotCompT(911,time[4:],sumTotEK[:],sumEvoPE,sumEvoTot,"$KE(t)-KE(t=0)$","$PE(t)-PE(t=0)$","$ME(t)-ME(t=0)$",)
plotCompT(915,time[4:],sumEvoTot,sumTotDISS[:],sumEvoTot2,"$ME(t)-ME(t=0)$","Dissipation","new==>$ME(t)-ME(t=0)$")
plotCompT(917,time[4:],sumEvoTot2,sumTotFEXT[:],sumEvoTot3,"$ME(t)-ME(t=0)$","External Forcing","new==>$ME(t)-ME(t=0)$")
#plotComp(921,time[5:],sumEvoTot3,sumEvoDiffTerm,sumEvoTot4,"$ME(t)-ME(t=0)$","Therm Diffusion","new==>$ME(t)-ME(t=0)$")
plotComp4T(919,time[4:],sumEvoTot,sumTotDISS[:],sumTotFEXT[:],sumEvoTot3,"$ME(t)-ME(t=0)$","Dissipation","External Forcing","residu",'Total energy evolution')


# =============================================================================
# #############---------------Transfer Analysis-----------###########"
# =============================================================================
#%%
transEvo=evoRhoGW+totUPRESS[:]
plotComp(2001,time[4:],evoWb2,totUPRESS[:],transEvo,'rhog*w from PE','rho*g*w from KE','difference between')


#%%
plotComp(2002,time[4:],sumTotEK[:],sumTotUPRESS,sumTotEK[:]+sumTotUPRESS,'$KE(t)-KE(t=0)$','$g*rho*w$','difference')

#%%
# =============================================================================
# #############---------------Kinetic Energy-----------###########"
# =============================================================================
#%%
plt.figure(1004)
plt.plot(time[tStart:],sumTotUETA[:], label = "$\sum D_t(EK)$")
plt.plot(time[tStart:],sumTotUPRESS[:], label = "Fext")
#plt.plot(time[40:],totUDISSh, label = 'horizontal dissipation')
#plt.plot(time[tStart:],sumTotDISS[:], label = 'Total disspation')
#plt.plot(time[40:],totUDISSv, label ='vertical dissipation')
plt.plot(time[tStart:],sumEvoRhoGW[:], label = 'Hydrostatic Pressure grad')
#plt.plot(time[40:],totUETA, label = 'Pression Force Vertical')
#plt.plot(time[40:],totUAB, label = 'from error')

tototot=sumEvoRhoGW+sumTotUPRESS+sumTotUETA
plt.plot(time[tStart:],tototot[:], label = 'KE-Fext-Dissip')

plt.title("Evolution of Kinetic Energy")
plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-2}$ ]")
plt.legend()
plt.grid(True)
plt.show()
#%%
plt.figure(1005)
plt.plot(time[5:],sumTotUPRESS[:], label = 'Hydrostatic Pressure grad')

plt.title("Evolution of Transfer term from kinetic energy to potential Energy")
plt.xlabel("t [ day ]")
plt.ylabel("Energy [ kg.m$^{-2}$.s$^{-2}$ ]")
plt.legend()
plt.grid(True)
plt.show()
#%%
plotCompT(1006,time[4:],sumTotDISS[:],sumTotDISSh[:],sumTotDISSv[:],'Total disspation','horizontal dissipation','vertical dissipation','Evolution of dissipation')

#%%
plotCompT(911,time[4:],sumTotEK[:],sumEvoPE,sumEvoTot,"$KE(t)-KE(t=0)$","$PE(t)-PE(t=0)$","$ME(t)-ME(t=0)$")

#%%
plt.figure(1007)
plt.plot(time[5:],totFEXT, label = "Fext")
plt.plot(time[5:],totUDISS, label = 'UDISS')

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
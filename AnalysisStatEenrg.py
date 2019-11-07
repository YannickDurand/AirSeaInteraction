#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 11:25:33 2019

@author: ydurand

Data Analysis for Energy
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.stats

from stockData import *
from stockDiagKEs import *

#%%
# =============================================================================
# ------------------------Selection Variable------------------------
# =============================================================================
deltaT=6*3600.                #en s
maxT=0.
alphat = 2e-4
rho0=999.8
tStart=7*4
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
PATH='/dsk1/ydurand/T03R01N3'
#PATH='/media/durand/Yannick/19_OAC/Control_Case/T03R01N3H60'
# =============================================================================



fileS='/state.nc'
fileKEs='/diagKEs.nc'
grid='/grid.nc'
filePE='/tracer_wb.nc'


dataNetCDFK=(PATH+fileKEs)
dataNetCDFS=(PATH+fileS)
dataGrid=(PATH+grid)
dataNetPE=(PATH+filePE)

#-------------------------------load grid-------------------------------------#
xxx,si_x,nbIter2,time2=dimRes(dataNetCDFS)
nbIter,time=getTime(dataNetCDFK)
time=np.divide(time[:],(24*3600))
si_y=si_x

RAC,zzt,zzl,zzd,zzc,dimZZ=dimGrid(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]


#%%
KEplusDiss=np.array([])

stratValue        =np.array(([3.e-5], [2.e-5], [1.e-5]))

kineticEnergyStrat=np.array(([3.e-5,-1.23e13], [2.e-5,-1.25e13], [1.e-5,-1.29e13]))
rhogwStrat        =np.array(([3.e-5, 7.27e12], [2.e-5, 9.64e12], [1.e-5, 1.44e13]))
dissStrat         =np.array(([3.e-5,-4.78e13], [2.e-5,-5.14e13], [1.e-5,-5.75e13]))
fExtStrat         =np.array(([3.e-5, 2.80e13], [2.e-5, 2.94e13], [1.e-5, 3.03e13]))


dissVStrat        =np.array(([3.e-5,-3.74e13], [2.e-5,-40.4e13], [1.e-5,-4.67e13]))
dissHStrat        =np.array(([3.e-5,-1.10e13], [2.e-5,-1.09e13], [1.e-5,-1.08e13]))


Residue           =np.array(([3.e-5, 5.89e13], [2.e-5, 6.21e13], [1.e-5, 6.27e13]))

KEplusDisstemp=dissStrat[:,1]+fExtStrat[:,1]
KEplusDiss=np.insert(stratValue,1,KEplusDisstemp,axis=1)
KEplusDiss=KEplusDiss[:-1,:]
#%%
# =============================================================================
# Linear regression 
# =============================================================================
absComp=np.arange(1,4)*1e-5


#Kinetic Energy
lr = scipy.stats.linregress((kineticEnergyStrat))
coeffN1=lr[0]
trendKE=coeffN1*(absComp)+lr[1]

#Rhogw
lr2 = scipy.stats.linregress((rhogwStrat))
coeffN1=lr2[0]
trendRhogw=coeffN1*(absComp)+lr2[1]

#diss
     #total
lr3 = scipy.stats.linregress((dissStrat))
coeffN1=lr3[0]
trendDiss=coeffN1*(absComp)+lr3[1]
     #vertical
lr4 = scipy.stats.linregress((dissVStrat))
coeffN1=lr4[0]
trendDissV=coeffN1*(absComp)

#ExternalForcing
lr5 = scipy.stats.linregress((fExtStrat))
coeffN1=lr5[0]
trendFext=coeffN1*(absComp)+lr5[1]

#Reidue
lr6 = scipy.stats.linregress((Residue))
coeffN1=lr6[0]
trendResid=coeffN1*(absComp)+lr6[1]

#TryTotry
lr7=scipy.stats.linregress((KEplusDiss))
#%%
# =============================================================================
# # plot
# =============================================================================
plt.figure(1001)
plt.scatter(kineticEnergyStrat[:,0]*1e5,kineticEnergyStrat[:,1])
plt.plot(absComp*1e5,trendKE,'r')
plt.title("Evolution of Kinetic Energy as a function of stratification")
plt.xlim(0,4e-5*1e5)
plt.ylabel("Energy variation after 40 days [ kg.m$^{-2}$.s$^{-2}$ ]")
plt.xlabel("Brunt-Vaïsaïla square [$10^{-5}$]")
plt.legend()
plt.grid(True)
plt.text(2.5,-1.225e13,"trend curve = "+str(int(lr[0]/1e15))+"$x10^{15}N²$ + ("+str(int(lr[1]/1e11))+")x$10^{11}$",style = 'italic', fontsize = 8,)
plt.text(2.5,-1.228e13,"Correlation = "+str(lr[2]),style = 'italic', fontsize = 8,)
plt.show()
#%%
plt.figure(1003)
plt.scatter(kineticEnergyStrat[:,0]*1e5,rhogwStrat[:,1])
plt.plot(absComp*1e5,trendRhogw,'r')
plt.title("Evolution of $rho gw$ as a function of stratification")
plt.xlim(0,4e-5*1e5)
plt.ylabel("Energy variation after 40 days [ kg.m$^{-2}$.s$^{-2}$ ]")
plt.xlabel("Brunt-Vaïsaïla square [$10^{-5}.s^{-2}$]")
plt.legend()
plt.grid(True)
plt.text(2.5,9.6e12,"trend curve = "+str(int(lr2[0]/1e15))+"$x10^{15}N²$ + ("+str(int(lr2[1]/1e11))+")x$10^{11}$",style = 'italic', fontsize = 8,)
plt.text(2.5,9.3e12,"Correlation = "+str(lr2[2]),style = 'italic', fontsize = 8,)
plt.show()

#%%
plt.figure(1005)
plt.scatter(kineticEnergyStrat[:,0]*1e5,dissStrat[:,1])
plt.plot(absComp*1e5,trendDiss,'r')
plt.title("Evolution of dissipation as a function of stratification")
plt.xlim(0,4e-5*1e5)
plt.ylabel("Energy variation after 40 days [ kg.m$^{-2}$.s$^{-2}$ ]")
plt.xlabel("Brunt-Vaïsaïla square [$10^{-5}.s^{-2}$]")
plt.legend()
plt.grid(True)
plt.text(2.5,-5.1e13,"trend curve = "+str(int(lr3[0]/1e15))+"$x10^{15}N²$ + ("+str(int(lr3[1]/1e11))+")x$10^{11}$",style = 'italic', fontsize = 8,)
plt.text(2.5,-5.14e13,"Correlation = "+str(lr3[2]),style = 'italic', fontsize = 8,)
plt.show()

#%%
plt.figure(1007)
plt.scatter(kineticEnergyStrat[:,0]*1e5,fExtStrat[:,1])
plt.plot(absComp*1e5,trendFext,'r')
plt.title("Evolution of external forcing as a function of stratification")
plt.xlim(0,4e-5*1e5)
plt.ylabel("Energy variation after 40 days [ kg.m$^{-2}$.s$^{-2}$ ]")
plt.xlabel("Brunt-Vaïsaïla square [$10^{-5}.s^{-2}$]")
plt.legend()
plt.grid(True)
plt.text(2.5,2.88e13,"trend curve = "+str(int(lr5[0]/1e15))+"$x10^{15}N²$ + ("+str(int(lr5[1]/1e11))+")x$10^{11}$",style = 'italic', fontsize = 8,)
plt.text(2.5,2.87e13,"Correlation = "+str(lr5[2]),style = 'italic', fontsize = 8,)
plt.show()

#%%
plt.figure(1009)
plt.scatter(kineticEnergyStrat[:,0]*1e5,Residue[:,1])
plt.plot(absComp*1e5,trendResid,'r')
plt.title("Evolution of the Residue as a function of stratification")
plt.xlim(0,4e-5*1e5)
plt.ylabel("Energy variation after 40 days [ kg.m$^{-2}$.s$^{-2}$ ]")
plt.xlabel("Brunt-Vaïsaïla square [$10^{-5}.s^{-2}$]")
plt.legend()
plt.grid(True)
plt.text(2.5,6.135e13,"trend curve = "+str(int(lr6[0]/1e15))+"$x10^{15}N²$ + ("+str(int(lr6[1]/1e11))+")x$10^{11}$",style = 'italic', fontsize = 8,)
plt.text(2.5,6.12e13,"Correlation = "+str(lr6[2]),style = 'italic', fontsize = 8,)
plt.show()
#%%
# =============================================================================
# Pourcent with the case N3 
# =============================================================================
fracStratification=kineticEnergyStrat[:,0]/kineticEnergyStrat[0,0]

kineticEnergyPourcent=(kineticEnergyStrat[:,1]-kineticEnergyStrat[2,1])/(kineticEnergyStrat[2,1])

rhogwPourcent=(rhogwStrat[:,1])/rhogwStrat[2,1]

dissPourcent=dissStrat[:,1]/dissStrat[2,1]

fExtPourcent=fExtStrat[:,1]/fExtStrat[2,1]

dissHPourcent=dissHStrat[:,1]/dissHStrat[2,1]
dissVPourcent=dissVStrat[:,1]/dissVStrat[2,1]

"""

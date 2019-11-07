#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 12:06:36 2019

@author: ydurand
"""
import numpy as np

#from fonctionH import *
#from stockData import *

import matplotlib.pyplot as plt
import matplotlib.cm as cm


from fonctionH import *
from Sauvegarde import *
from stockDiagKEs import *



# =============================================================================
#%%
###############################################################################
#-------------------------Folder for NetCDF------------------------------------
###############################################################################

PATH='/dsk1/ydurand/T03R01N3'
PATH2='/dsk1/ydurand/T03R01N1'



fileS='/state.nc'
fileKEs='/diagKEs.nc'
filePE='/tracer_wb.nc'
grid='/grid.nc'

dataNetCDFK=(PATH+fileKEs)
dataGrid=(PATH+grid)
dataNetCDFS=(PATH+fileS)
dataNetCDFS2=(PATH2+fileS)

dataNetPE=(PATH+filePE)

#%%
###############################################################################
#------------------------Cell to stock data------------------------------------
###############################################################################

xxx,xxx1,yyy,zzt,dimZZT,zzc,time,nbIter,res=meshconsState(dataNetCDFS)
time=time/(24*3600)


RAC,zzt,zzl,zzd,dimZZ=dimGridState(dataGrid)
areaGrid=RAC[np.newaxis,:,:]*zzd[:,np.newaxis,np.newaxis]
#%%
#%%
###############################################################################
#-----------------------Initial state------------------------------------------
###############################################################################

hhd0,ujk0,vik0,wij0,temp0,eta0=vitesse(dataNetCDFS,zzt,HHH=-1,tstart=0)
hhd0,ujk0,vik0,wij0,temp0,eta2=vitesse(dataNetCDFS2,zzt,HHH=-1,tstart=0)

eta1_0_max=np.max(eta0)
eta2_0_max=np.max(eta2)
#%%
###############################################################################
#------------------------Calculation-------------------------------------------
#---Fonction use : 
#---vitesse2(dataNetCDF,HHH=0,tstart=0,tend=0)
##########################################################  #####################

sumEk=[]
sumEk2=[]
eta1=np.zeros(nbIter)
eta2=np.zeros(nbIter)

for i in range(nbIter):
     hhd,ujk,vik,wij,temp,eta=vitesse(dataNetCDFS,zzt,HHH=-1,tstart=i)
     hhd0,ujk0,vik0,wij0,temp0,eta22=vitesse(dataNetCDFS2,zzt,HHH=-1,tstart=i)
     
     Ek=energyK(ujk,vik,hhd,i)*999.8
     Ek=np.multiply(Ek,areaGrid)#*999.8
# =============================================================================
     Ek2=energyK(ujk0,vik0,hhd0,i)
     Ek2=np.multiply(Ek2,areaGrid)
# =============================================================================
     
     sumEk.append(np.sum(Ek))
# =============================================================================
     sumEk2.append(np.sum(Ek2))
# =============================================================================

     #for elevation of the height
     eta1[i]=(np.max(eta)/eta1_0_max)
# =============================================================================
#      eta22=eta((fileNetCdf2),tstart=i,tend=i+1)
     eta2[i]=(np.max(eta22)/eta2_0_max)
# =============================================================================
     print(i)
     
#%%  
plt.figure(10)     
plt.plot(time,sumEk,label='N3')
plt.xlabel("t [ day ]")
plt.ylabel("Energy variation [ m$^{-1}$.s$^{-2}$ ]")
plt.grid(True)
plt.legend()
plt.title("Total Kinetic Energy")
# =============================================================================
# plt.plot(sumEk2)
# =============================================================================
plt.show()
#%%
plt.figure(11)
plt.plot(time,eta1,label='N3')
plt.ylabel("h(t)/h(t=0)")
plt.xlabel("t [ day ]")
plt.title("Evolution of the height pertubation")
plt.grid(True)
plt.legend()
# =============================================================================
# plt.plot(eta2)
# =============================================================================
plt.show()
#%%
# =============================================================================
# --------------------------Time Evolution-------------------------------------
# =============================================================================

#%%
# =============================================================================
# ##############----------Temp evolution-------------------############
# =============================================================================
"""

######
#Take a look on the evolution of the temperature of each layer
######
tempMean=np.zeros((nbIter-1,20))
for i in range(nbIter-1):
     t=i
     temp=getTemp(dataNetCDFS,t)
     tempMean[i,:]=[np.mean(temp[j,:,:]) for j in range(20)]
     print(i)

#%%
plt.figure(700)
plt.plot(time[1:],tempMean[:,0],label="level "+str(int(zzt[0]))+" m")
plt.plot(time[1:],tempMean[:,1],label="level "+str(int(zzt[1]))+" m")
plt.plot(time[1:],tempMean[:,2],label="level "+str(int(zzt[2]))+" m")
plt.plot(time[1:],tempMean[:,3],label="level "+str(int(zzt[3]))+" m")
plt.plot(time[1:],tempMean[:,4],label="level "+str(int(zzt[4]))+" m")
plt.plot(time[1:],tempMean[:,5],label="level "+str(int(zzt[5]))+" m")


#plt.plot(time[1:],tempMean[:,7],label="level"+str(int(zzt[7])))
#plt.plot(time[1:],tempMean[:,10],label="level"+str(int(zzt[10])))
#plt.plot(time[1:],tempMean[:,15],label="level"+str(int(zzt[15])))

plt.ylabel("T [°C]")
plt.xlabel("t [ day ]")
plt.title("Evolution in time of temperature at specific level")
plt.legend()

plt.show()


#%%
# =============================================================================
# ##############----------Enstrophy evolution-------------------############
# =============================================================================

dx,dy=1200.,1200.
t=49
evoTotEns=np.array([])
#initial state
ujk,vik=getVarVor(dataNetCDFS,t)
vort_0=vort_pot(dx, dy, ujk, vik)
enst_0=potentialEnstrophy(vort_0, zzd)*dx**2
sumEnst_0=np.sum(enst_0)

for i in range(50,nbIter-1):
     t=i
     ujk,vik=getVarVor(dataNetCDFS,t)
     vort=vort_pot(dx, dy, ujk, vik)
     enst=potentialEnstrophy(vort, zzd)*dx**2
     evoTotEns=np.append(evoTotEns,(np.sum(enst)/sumEnst_0))
     
     if i%50==0:
          plt.figure(i)
          plt.contourf(enst[1,:,:],cmap=cm.jet, levels=100)
          plt.show()

#%%
plt.figure(701)
plt.plot(evoTotEns)
plt.show()


"""


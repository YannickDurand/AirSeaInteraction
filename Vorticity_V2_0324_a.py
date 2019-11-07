#!/usr/bin/env python
# coding: utf-8

#%%


#!/usr/bin/env python3

"""
Created on Fri Mar 15 2019

@author: ydurand
"""


#%%


import numpy as np
import matplotlib.pyplot as plt
import math as math
import os
#import scipy.io.netcdf as netcdf
from netCDF4 import Dataset


from Operator import *
from fonctionH import *
from Sauvegarde import *
from lesPlots import *





#%%
###########################################################################
#-------------Cell for Images Stockage-----------------------------------------
###########################################################################

stock=0
if stock==1:
    PATH='/home/durand/OAC/5_Analysis'
    NAME='/Case_Barotrope_AC/contour/' 
    #shutil.rmtree(PATH+NAME)
    os.makedirs((PATH+NAME),exist_ok=True)
    
    folderV='/vorticity/'
    folderD='/divergence/'
    folderEk='/Kinetic_Energy/'
    folderEp='/Potential_Energy/'
    folderEns='/Potential_Enstrophy/'




#%%
##############################################################################
#--------------------------Stock Variable----------------------------------------
###############################################################################

#dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Modon_Bruno/Case_withoutAtm/Case_1/Global'
#dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Dipole_Yannick/Case_Dipole_withAtm/Case_2_onlyV/Global_Vn'
#dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Modon_Bruno/Case_withoutAtm/Case_1/Global'
#dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Modon_Bruno/Case_withAtm/Case_2_U/Global'
#dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Baroclinic_AC/WithAtm/Case_1_UandV'
#dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Baroclinic_AC/withoutAtm'
#dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Baroclinic_AC/betaPlane/Global'
dir0='/media/ydurand/Yannick/19_OAC/4_Test/Case_Baroclinic_AC/betaPlane/withU/uP'

file1 = '/state.nc'
#file1= '/state.0000000000.t001.nc'
file2 = '/state.0000004980.t001.nc'
file3 = '/state.0000009960.t001.nc'
    
#get variable Grid
xxx,xxx1,yyy,zzt,dimZZT,zzc,time,nbIter,res=meshcons((dir0 + file1))
#get variable parameter
HHH=1
j=HHH
print("H = ",j, "la hauteur est ", zzt[HHH])


hhd,ujk,vik,wij,temp,eta=vitesse((dir0 + file1),zzt, HHH)

#%%
stockIter=np.array([50,80])

#%%
##########################################################################
#------------------Other Data need----------------------------------------
##########################################################################

#information about the model
f=1E-4
dx,dy=1500.,1500.
Nx,Ny=res,res
Lx,Ly=Nx*dx,Ny*dy


x,y=np.ones(res-1),np.ones(res-1)
for i in range((res-1)):
    x[i]=i*dx
    y[i]=i*dy

#Calculation beforehand


tempo=np.ones(((res+1),(res+1),19),dtype=float)
diffH=np.multiply(tempo,hhd)


#%%

#Quelle profondeur?


#%%
#####----create folder to stock-----#####
if stock == 1:
    hauteur=(str(int(zzt[HHH])))+("/")
    os.makedirs(PATH+NAME+folderV+hauteur,exist_ok=True)
    vorPath=PATH+NAME+folderV+hauteur
    os.makedirs(PATH+NAME+folderD+hauteur,exist_ok=True)
    divPath=PATH+NAME+folderD+hauteur
    os.makedirs(PATH+NAME+folderEk+hauteur,exist_ok=True)
    ekPath=PATH+NAME+folderEk+hauteur
    os.makedirs(PATH+NAME+folderEp+hauteur,exist_ok=True)
    epPath=PATH+NAME+folderEp+hauteur
    os.makedirs(PATH+NAME+folderEns+hauteur,exist_ok=True)
    ensPath=PATH+NAME+folderEns+hauteur

#%%    
###########################################################################
#-------------------Plots--------------------------------------------------
###########################################################################
for j in range(len(stockIter)):
        i=stockIter[j]
        #print(ujk[i].shape)
        hik=diffH[:,:,j]
        t=time[i]
        qij=vort_pot(dx, dy, f, ujk[:,:,HHH], vik[:,:,HHH])
        #Calculation beforehand
        diver=conv(dx, dy, ujk[:,:,HHH,i], vik[:,:,HHH,i])*1E7
        #Ek=energyK(uck[:,:,i], vck[:,:,i], hik)
        enPot=energyP(dx, dy, hik)
        #potEns=potentialEnstrophy(qij[:,:,j,i],hik[1:-2,1:-2])
        
        plotEta(xxx,yyy,eta[:,:,i],t) 
        #plotVort(x,y,qij[:,:,i]*(1e5),t)
        #plotDiv(xxx,yyy,diver,t)
        """
        #####-----Eta-----#####
        plt.figure(figsize=(12.5,10*Ly/Lx))
        plt.contour(xxx,yyy, hik, 20)
        plt.colorbar()
        plt.title('Eta' + ' ; t='+str(t))
        plt.show()
        
        #####-----Voritcity------#####
        plt.figure(figsize=(12.5,10*Ly/Lx))
        plt.contour(x,y, qij/1E-4, 20)
        plt.colorbar()
        plt.title('vorticite potentielle' + ' ; t='+str(t))
        plt.show()
        plt.close()
        
        #####-----divergence-----#####
        plt.figure(figsize=(12.5,10*Ly/Lx))
        plt.contour(xxx,yyy,diver,20)
        plt.colorbar()
        plt.title('Convergence' + ' ; t='+str(t))
        plt.show()
        #plt.close(901)
        
        #####-----Kinetic Energy-----#####
        plt.figure(figsize=(12.5,10*Ly/Lx))
        Ek=energyK(uck[i], vck[i], hik)
        plt.contour(xxx,yyy,Ek,20)
        plt.colorbar()
        plt.title('Kinetic Energy' + ' ; t='+str(t))
        plt.show()
        
        #####-----Potential Energy-----#####
        
        plt.figure(figsize=(12.5,10*Ly/Lx))
        plt.contour(xxx,yyy,energyP(dx, dy, hik),20)
        plt.colorbar()
        plt.title('Potential Energy' + ' ; t='+str(t))
        plt.show()
        
        #####-----Potential Enstrophy-----#####
        PE=potentialEnstrophy(qij, hik[1:-1,1:-1])
        plt.figure(figsize=(12.5,10*Ly/Lx))
        plt.contour(PE,20)
        plt.colorbar()
        plt.title('Potential Enstrophy ' + ' ; t='+str(t))
        plt.show()  
        """
        
        
        """
        #####-----Voritcity------#####
        qij=vort_pot(mesh, dx, dy, f, ujk[i], vik[i])
        if t==0 and zzt[j]>-0.1:
            qij=np.random.randn((res-1), (res-1))*1.E-2;
        figure1=plt.figure(1, figsize=(12.5,10*Ly/Lx))
        plt.contour(x,y, qij/1E-4, 20)
        plt.colorbar()
        plt.title('vorticite potentielle' + ' ; t='+str(t))
        figure1.savefig(str(vorPath)+'vorticity_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
        plt.close(1)

        #####-----divergence-----#####
        diver=conv(dx, dy, ujk[i], vik[i])
        if t==0 and zzt[j]>-0.1:
            diver=np.random.randn((250), (250))*1.E-3;
        figure2=plt.figure(2, figsize=(12.5,10*Ly/Lx))
        plt.contour(dx*mesh.xik,dy*mesh.yik,diver,20)
        plt.colorbar()
        plt.title('Divergence' + ' ; t='+str(t))
        figure2.savefig(str(divPath)+'divergence_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
        plt.close(2)

        #####-----Kinetic Energy-----#####
        figure3=plt.figure(3, figsize=(12.5,10*Ly/Lx))
        Ek=energyK(uck[i], vck[i], hik)
        if t==0 and zzt[j]>-0.1:
            Ek=np.random.randn((250), (250))*1.E-2;
        plt.contour(dx*mesh.xik,dy*mesh.yik,Ek,20)
        plt.colorbar()
        plt.title('Kinetic Energy' + ' ; t='+str(t))
        figure3.savefig(str(ekPath)+'EK_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
        plt.close(3)

        #####-----Potential Energy-----#####
        
        figure4=plt.figure(4, figsize=(12.5,10*Ly/Lx))
        plt.contour(dx*mesh.xik,dy*mesh.yik,energyP(dx, dy, hik),20)
        plt.colorbar()
        plt.title('Potential Energy' + ' ; t='+str(t))
        figure4.savefig(str(epPath)+'EP_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
        plt.close(4)
        
        #####-----Potential Enstrophy-----#####
        PE=potentialEnstrophy(qij, hik[1:-1,1:-1])
        figure5=plt.figure(5, figsize=(12.5,10*Ly/Lx))
        plt.contour(PE,20)
        plt.colorbar()
        plt.title('Potential Enstrophy ' + ' ; t='+str(t))
        figure5.savefig(str(ensPath)+'Ens_at_'+str(int(zzt[HHH]))+'m_time_'+str(t)+'s.png')
        plt.close(5)  
        """

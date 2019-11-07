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
import shutil
#import scipy.io.netcdf as netcdf
from netCDF4 import Dataset


from grille import *
from Operator import *
from fonctionH import *





#%%
###########################################################################
#-------------Cell for Images Stockage-----------------------------------------
###########################################################################
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


###########################################################################
#-------------Cell for import Data-----------------------------------------
###########################################################################
dir0 = '/home/durand/OAC/4_Test/Case_Barotrope_AC'

#dir0 = '/dsk1/ydurand/OAC/Test_Barotrop'
#dir0='/home/durand/OneDrive/1_Mes_Documents/1_Physique/1_Meteorology/19_OAC/4_Test/testBarotrope'
file1 = '/state.0000000000.t001.nc'
file2 = '/state.0000004980.t001.nc'
file3 = '/state.0000009960.t001.nc'


tropEddy=Dataset(dir0 + file3)

print ("Dimension of the data")
print (tropEddy.dimensions.keys()) 
print (tropEddy.variables['Zl']) 

print ("interrogate variable")
print (tropEddy.variables.keys()) 
print("\n Pour la variable u")
print (tropEddy.variables['U'])
print("\n Pour la variable v")
print (tropEddy.variables['V'])
print("\n Pour la variable Eta ")
print (tropEddy.variables['Eta'])
print("\n Pour la variable Z ")
print (tropEddy.variables['Z'])
print("\n Pour la variable Zl ")
print (tropEddy.variables['Zl'])
print("\n Pour la variable Yp1 ")
print (tropEddy.variables['Yp1'])


#%%


##########################################################################
#------------------Other Data need----------------------------------------
##########################################################################

#information about the model
res=250
f=1E-4
dx,dy=1200.,1200.
Nx,Ny=res,res
Lx,Ly=Nx*dx,Ny*dy
x,y=np.ones(249),np.ones(249)
for i in range(249):
    x[i]=i*dx
    y[i]=i*dy

xmesh = Periodic_FD(Nx)
ymesh = Periodic_FD(Ny)
mesh = Rectangle(xmesh,ymesh)
#%%


############################################################################
#--------------Cell to get the variable from Netcdf File--------------------
############################################################################
############################################################################
#--------------Cell to get the variable from Netcdf File--------------------
############################################################################

#get the grid
xxx=tropEddy.variables['X'][:]
yyy=tropEddy.variables['Y'][:]
zzt=tropEddy.variables['Zl'][:]
zzc=tropEddy.variables['Z'][:]
eta=tropEddy.variables['Eta'][:]
time=tropEddy.variables['T'][:]
nbIter=len(tropEddy.dimensions['T'])
nbZ=len(tropEddy.dimensions['Z'])
nbIter=int(nbIter)
print("Le nombre d'itération est ",nbIter)




#A remplacer par des niveaux de densité
#print("\n", zzz[1])
#Quelle profondeur?
print("les hauteurs sont :")
for j in range (nbZ):
    print("H = ",j, "la hauteur est ", zzt[j])
    HHH=j
    #print(eta[1,120:140,120:140])
    ujk=[]
    uck,vck=[],[]
    tt=[]
    vik=[]
    hhh=[]
    #get variable
    for i in range(nbIter):
        ujk.append(tropEddy.variables['U'][i,HHH,:,:])
        vik.append(tropEddy.variables['V'][i,HHH,:,:])
        if HHH<=0:
            hhh.append(tropEddy.variables['Eta'][i,:,:]-zzt[HHH])
        else:
            Ht=HHH-1
            hhh.append(tropEddy.variables['Eta'][i,:,:]-zzt[HHH]+zzt[Ht])
        ujk[i]=np.transpose(ujk[i])
        vik[i]=np.transpose(vik[i])
        hhh[i]=np.transpose(hhh[i])
        uck.append(ucc(ujk[i],res))  #create center variable
        vck.append(vcc(vik[i],res))
    print("La hauteur de la couche est en moyenne de :", np.mean(hhh), " m")
    #tropEddy.close()       #to create space in memory we close the netcdf file
    
        #####----create folder to stock-----#####
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
    
    ###########################################################################
    #-------------------Plots--------------------------------------------------
    ###########################################################################
    for i in range(1,nbIter):
            
            #print(ujk[i].shape)
            hik=hhh[i]
            t=time[i]  
            """
            #####-----Eta-----#####
            plt.figure(figsize=(12.5,10*Ly/Lx))
            plt.contour(xxx,yyy, hik, 20)
            plt.colorbar()
            plt.title('Eta' + ' ; t='+str(t))
            plt.show()
            
            #####-----Voritcity------#####
            qij=vort_pot(mesh, dx, dy, f, ujk[i], vik[i])
            plt.figure(figsize=(12.5,10*Ly/Lx))
            plt.contour(x,y, qij/1E-4, 20)
            plt.colorbar()
            plt.title('vorticite potentielle' + ' ; t='+str(t))
            plt.show()
            
            #####-----divergence-----#####
            diver=conv(dx, dy, ujk[i], vik[i])
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
            temp=hik[1:-1,1:-1]
            print(temp.shape)
            PE=potentialEnstrophy(qij, hik[1:-1,1:-1])
            plt.figure(figsize=(12.5,10*Ly/Lx))
            plt.contour(PE,20)
            plt.colorbar()
            plt.title('Potential Enstrophy ' + ' ; t='+str(t))
            plt.show()  
            
            
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

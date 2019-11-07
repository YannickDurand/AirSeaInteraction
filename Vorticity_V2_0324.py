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


# In[3]:


#########################################################################
#----------------Grid and operator---------------------------------------
#########################################################################
def partui(ujk,res,dx):
    diffu=np.diff(ujk, axis=0)
    return diffu

def partvj(vik,res,dy):
    diffv=np.diff(vik, axis=1)
    return diffv

def partuj(ujk,res,dy):
    ujk=ujk[1:res,:]
    difuj=np.zeros((res-1,res-1))
    for i in range(1,res-1):
        difuj[:,i]=(ujk[:,i+1]-ujk[:,i])/dy
    return difuj

def partvi(vik,res,dx):
    vik=vik[:,1:res]
    difvi=np.zeros((res-1,res-1))
    for i in range(1,res-1):
        difvi[i,:]=(vik[i+1,:]-vik[i,:])/dx
    return difvi
        

def ucc(ujk,res):
    uci=np.ones((res,res))
    for i in range(res):
        left=ujk[i,:]
        right=ujk[i+1,:]
        uci[i,:]=(left+right)/2
    return uci

def vcc(vik,res):
    vcj=np.ones((res,res))
    for i in range(res):
        left=vik[:,i]
        right=vik[:,i+1]
        vcj[:,i]=(left+right)/2
    return vcj


# In[4]:


##########################################################################
#-------------Cell for subroutine-----------------------------------------
##########################################################################

#####-----Vorticty-----#####
def vort_pot(mesh, dx, dy, f, ujk, vik):
    qij=-partuj(ujk,res,dy)+partvi(vik,res,dx)
    return qij

#####-----Divergence-----#####
def conv(mesh, dx, dy, ujk, vik):
    conv=partui(ujk,res,dx)+partvj(vik,res,dy)
    return conv



#####-----Kinetic Energy-----#####
def energyK(mesh, ujk, vil, hik):
    Kx = mesh.mj(ujk*ujk)*hik
    Ky = mesh.ml(vil*vil)*hik
    return Kx+Ky

#####-----Potential Energy-----#####
def energyP(mesh, dx, dy, hik):
    Pe=(1./2.)*(hik)*(hik)
    return Pe

#####-----Mecanic Energy-----#####
def sumEnergyK(mesh, ujk, vil, hik):
    KxT = (mesh.mi(ujk*ujk)*hik).sum()
    KyT = (mesh.mk(vil*vil)*hik).sum()
    PeT=((1./2.)*(hik)*(hik)).sum()
    ekT = KxT + KyT + PeT
    return ekT

def sumEnergyP(hik):
    PeT=((1./2.)*(hik)*(hik)).sum()
    return PeT

#####-----Potential Enstrophy-----#####
def potentialEnstrophy(mesh, qik, hik):
    Z=mesh.mk(qik*qik)*hik*(1./2.)
    return Z

def sumPE(mesh,qik,hik):
    Zt=(-mesh.mi(qik*qik)*hik*(1./2.)).sum()
    return Zt

#%%
###########################################################################
#-------------Cell for Images Stockage-----------------------------------------
###########################################################################
stock=0
if stock: 
    PATH='/home/durand/OneDrive/1_Mes_Documents/1_Physique/1_Meteorology/19_OAC/5_Analysis'
    NAME='/test_Barotropic' 
    shutil.rmtree(PATH+NAME)
    os.makedirs(PATH+NAME)#,exist_ok=True)
    
    folderV='/vorticity/'
    os.mkdir(PATH+NAME+folderV)
    vorPath=PATH+NAME+folderV
    
    folderD='/divergence/'
    os.mkdir(PATH+NAME+folderD)
    divPath=PATH+NAME+folderD
    
    folderEk='/Kinetic_Energy/'
    os.mkdir(PATH+NAME+folderEk)
    ekPath=PATH+NAME+folderEk
    
    folderEp='/Potential_Energy/'
    os.mkdir(PATH+NAME+folderEp)
    epPath=PATH+NAME+folderEp
    
    folderEns='/Potential_Enstrophy/'
    os.mkdir(PATH+NAME+folderEns)
    ensPath=PATH+NAME+folderEns

print(ensPath)

#%%


###########################################################################
#-------------Cell for import Data-----------------------------------------
###########################################################################
#dir0 = '/media/durand/6ABF2172220D421E/OneDrive/1_Mes_Documents/1_Physique/1_Meteorology/19_OAC/4_Test/1_test/mnc_test_0002/'
#dir0 = '/dsk1/ydurand/OAC/Test_Barotrop'
#dir0='/home/durand/OneDrive/1_Mes_Documents/1_Physique/1_Meteorology/19_OAC/4_Test/testBarotrope'
dir0='/home/durand/OneDrive/1_Mes_Documents/1_Physique/1_Meteorology/19_OAC/4_Test/Case_2Eddies'

file1 = '/state.0000000000.t001.nc'
file2 = '/state.0000002490.t001.nc'
file3 = '/state.0000004980.t001.nc'


tropEddy=Dataset(dir0 + file1)

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

xmesh = Periodic_FD(Nx)
ymesh = Periodic_FD(Ny)
mesh = Rectangle(xmesh,ymesh)


#%%


############################################################################
#--------------Cell to get the variable from Netcdf File--------------------
############################################################################

#get the grid
xxx=tropEddy.variables['X'][:]
yyy=tropEddy.variables['Y'][:]
zzt=tropEddy.variables['Zl'][:]
zzc=tropEddy.variables['Z'][:]
eta=tropEddy.variables['Eta'][:]
time,nbIter=tropEddy.variables['T'][:],len(tropEddy.variables['T'])

ujk=[]
uck,vck=[],[]
tt=[]
vik=[]
hhh=[]

#A remplacer par des niveaux de densité
#print("\n", zzz[1])
#Quelle profondeur?
print("les hauteurs sont :")
for i in range (20):
    print("H = ",i, "la hauteur est ", zzt[i])
HHH=7

#print(eta[1,120:140,120:140])

#get variable
for i in range(nbIter):
    ujk.append(tropEddy.variables['U'][i,HHH,:,:])
    vik.append(tropEddy.variables['V'][i,HHH,:,:])
    hhh.append(tropEddy.variables['Eta'][i,:,:]-zzt[HHH])
    ujk[i]=np.transpose(ujk[i])
    vik[i]=np.transpose(vik[i])
    uck.append(ucc(ujk[i],res))  #create center variable
    vck.append(vcc(vik[i],res))
print("La hauteur de la couche est en moyenne de :", np.mean(hhh), " m")
#tropEddy.close()       #to create space in memory we close the netcdf file


# %%


###########################################################################
#---------------Information about the Data---------------------------------
###########################################################################


# %%


###########################################################################
#-------------------Plots--------------------------------------------------
###########################################################################
x,y=np.ones(249),np.ones(249)
for i in range(249):
    x[i]=i*dx
    y[i]=i*dy
for i in range(nbIter):
        print(i)
        #print(ujk[i].shape)
        hik=hhh[i]
        t=time[i]

        
        #####-----Voritcity------#####
        qij=vort_pot(mesh, dx, dy, f, ujk[i], vik[i])
        figure1=plt.figure(1, figsize=(12.5,10*Ly/Lx))
        plt.contour(x,y, qij/1E-4, 20)
        plt.colorbar()
        plt.title('vorticite potentielle' + ' ; t='+str(t))
        figure1.savefig(str(vorPath)+'vorticity_at_time_'+str(t)+'s.png')
        plt.close(1)
        
        #####-----divergence-----#####
        diver=conv(mesh, dx, dy, ujk[i], vik[i])
        figure2=plt.figure(2, figsize=(14,10*Ly/Lx))
        plt.contourf(dx*mesh.xik,dy*mesh.yik,diver,100)
        plt.colorbar()
        plt.title('Convergence' + ' ; t='+str(t))
        figure2.savefig(str(divPath)+'divergence_at_time_'+str(t)+'s.png')
        plt.close(2)
       
        #####-----Kinetic Energy-----#####
        figure3=plt.figure(3, figsize=(12.5,10*Ly/Lx))
        Ek=energyK(mesh, uck[i], vck[i], hik)
        plt.contour(dx*mesh.xik,dy*mesh.yik,Ek,20)
        plt.colorbar()
        plt.title('Kinetic Energy' + ' ; t='+str(t))
        figure3.savefig(str(ekPath)+'EK_at_time_'+str(t)+'s.png')
        plt.close(3)
        
        #####-----Potential Energy-----#####
        
        figure4=plt.figure(4, figsize=(12.5,10*Ly/Lx))
        plt.contour(dx*mesh.xik,dy*mesh.yik,energyP(mesh, dx, dy, hik),20)
        plt.colorbar()
        plt.title('Potential Energy' + ' ; t='+str(t))
        figure4.savefig(str(epPath)+'EP_at_time_'+str(t)+'s.png')
        plt.close(4)
        """
        #####-----Potential Enstrophy-----#####
        plt.figure(figsize=(12.5,10*Ly/Lx))
        plt.contour(dx*mesh.xik,dy*mesh.yik,potentialEnstrophy(mesh, qij, hik),20)
        plt.colorbar()
        plt.title('Potential Enstrophy ' + ' ; t='+str(t))
        plt.show()  
        """


# In[ ]:





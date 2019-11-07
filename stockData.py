#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 10:27:13 2019

@author: ydurand
"""

from netCDF4 import Dataset

import numpy as np

def vitesse2(dataNetCDF,H=-1,tstart=-1,tend=0):
    tropEddy=Dataset(dataNetCDF)
    hhd=np.array(np.diff(tropEddy.variables['Zl'][:]))
    
    if tstart == -1 and H == -1 :
        ujk=np.array(tropEddy.variables['U'][:,:,:,:])
        vik=np.array(tropEddy.variables['V'][:,:,:,:])

    elif H == -1 and tstart != -1 :
        ujk=np.array(tropEddy.variables['U'][tstart:tend,:,:,:])
        vik=np.array(tropEddy.variables['V'][tstart:tend,:,:,:])
        
    elif tstart == -1 and H != -1 :
        ujk=np.array(tropEddy.variables['U'][:,H,:,:])
        vik=np.array(tropEddy.variables['V'][:,H,:,:])

    ujk=np.transpose(ujk)
    vik=np.transpose(vik)

    tropEddy.close()
    return hhd,ujk,vik

def eta(dataNetCDF,tstart=-1,tend=0):
     tropEddy=Dataset(dataNetCDF)
     if tstart == -1 :
        eta=np.array(tropEddy.variables['Eta'][:,:,:])

     else :
        eta=np.array(tropEddy.variables['Eta'][tstart:tend,:,:])        
        
     tropEddy.close()   
     return eta

def vitesseW(dataNetCDF,tstart=-1):
    tropEddy=Dataset(dataNetCDF)
    
    if tstart == -1 :
        wij=np.array(tropEddy.variables['W'][:,:,:,:])

    else :
        wij=np.array(tropEddy.variables['W'][tstart,:,:,:])


    tropEddy.close()
    return wij

def vitesseWTave(dataNetCDF,tstart=-1):
    tropEddy=Dataset(dataNetCDF)
    
    if tstart == -1 :
        wij=np.array(tropEddy.variables['W'][:,:,:,:])

    else :
        wij=np.array(tropEddy.variables['WVEL'][tstart,:,:,:])


    tropEddy.close()
    return wij

def getEta(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)
     
     eta=np.array(tropEddy.variables['Eta'][tstart,:,:])        
        
     tropEddy.close()   
     return eta
 
    
def dimRes(dataNetCDF):
    tropEddy=Dataset(dataNetCDF) 

    xxx,res=tropEddy.variables['X'][:],len(tropEddy.dimensions['X'])
    
    
    tday,nbIter=np.array(tropEddy.variables['T']),len(tropEddy.dimensions['T'])
    
    tropEddy.close()
    return xxx,res,nbIter,tday

def dimGrid(dataNetCDF):
     tropEddy=Dataset(dataNetCDF) 
     
     rAz=np.array(tropEddy.variables['rA'][:,:])
     zzt,dimZZT=np.array(tropEddy.variables['Zl'][:]),len(tropEddy.dimensions['Zl'])
     zzl=np.array(tropEddy.variables['Zu'][:])
     zzc=np.array(tropEddy.variables['RC'][:])
     
     diffZZ=abs(np.diff(np.array(tropEddy.variables['Zp1'][:])))
     
     tropEddy.close()
     return rAz,zzt,zzl,diffZZ,zzc,dimZZT
    
    
# =============================================================================
# 
# =============================================================================
def getSpeed(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)

    ujk=np.array(tropEddy.variables['U'][t,1,:,:])
    vik=np.array(tropEddy.variables['V'][t,1,:,:])

    tropEddy.close()
    return ujk,vik    

def getSpeedSurface(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)

    ujk=np.array(tropEddy.variables['U'][t,0,:,:])
    vik=np.array(tropEddy.variables['V'][t,0,:,:])

    tropEddy.close()
    return ujk,vik    
    
def getPhy(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
    
    phy=np.array(tropEddy.variables['phiHyd'][t,:,:,:])
    
    tropEddy.close()
    return phy
    
def getSpeed2(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)

    ujk=np.array(tropEddy.variables['U'][t,:,:,:])
    vik=np.array(tropEddy.variables['V'][t,:,:,:])

    tropEddy.close()
    return ujk,vik

# =============================================================================
# 
# =============================================================================
def getUatm(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
    
    Uatm=np.array(tropEddy.variables['CH_Uwind'][t,0,:,:])
    
    tropEddy.close()
    return Uatm

def getTAU(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
    
    TAUX=np.array(tropEddy.variables['oceTAUX'][t,0,:,:])
    TAUY=np.array(tropEddy.variables['oceTAUY'][t,0,:,:])
    
    tropEddy.close()
    return TAUX, TAUY


















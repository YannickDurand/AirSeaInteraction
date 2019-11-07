#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 12:01:53 2019

@author: ydurand
"""


import numpy as np
from netCDF4 import Dataset

def getTime(dataNetCDF):
    tropEddy=Dataset(dataNetCDF)
    tday,nbIter=np.array(tropEddy.variables['T']),len(tropEddy.dimensions['T'])
    return nbIter,tday

def dimRes(dataNetCDF):
    tropEddy=Dataset(dataNetCDF) 

    print("\n Pour la variable Z ")
    print (tropEddy.variables['Z'])
    print("\n Pour la variable Zl ")
    print (tropEddy.variables['Zl'])
    print("\n Pour la variable T ")
    print (tropEddy.variables['T'])
    
    xxx,res=tropEddy.variables['X'][:],len(tropEddy.dimensions['X'])
    
    
    tday,nbIter=np.array(tropEddy.variables['T']),len(tropEddy.dimensions['T'])
    
    tropEddy.close()
    return xxx,res,nbIter,tday
# =============================================================================
# 
# =============================================================================
# =============================================================================
# 
# =============================================================================
def getU(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
  
    ujk=np.array(tropEddy.variables['UVEL'][t,:,:,:])


    tropEddy.close()
    return ujk

def getV(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
  
    vik=np.array(tropEddy.variables['VVEL'][t,:,:,:])


    tropEddy.close()
    return vik

def getWVEL(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
  
    wij=np.array(tropEddy.variables['WVEL'][t,:,:,:])


    tropEddy.close()
    return wij
 
def getWonX(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
  
    wij=np.array(tropEddy.variables['W'][t,:,:,:])


    tropEddy.close()
    return wij

def getPhi(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
  
    wij=np.array(tropEddy.variables['phiHyd'][t,:,:,:])


    tropEddy.close()
    return wij

def getWbonX(dataNetCDF,t):
     tropEddy=Dataset(dataNetCDF)

     wb2=np.array(tropEddy.variables['WRHOMAS2'][t,:,:,:]) 
             
     tropEddy.close()
     return wb2
# =============================================================================
# 
# =============================================================================

def getTemp(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)
     
     temp=np.array((tropEddy.variables['Temp'][tstart,:,:,:]))
     
     tropEddy.close()
     return temp



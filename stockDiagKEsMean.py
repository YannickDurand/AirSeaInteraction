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
def getDissVonY(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     udissv=np.array(tropEddy.variables['KEm_Disv'][28:,:,:,165])
   
     tropEddy.close()  
     return udissv
 
def getDissHonY(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     udissv=np.array(tropEddy.variables['KEm_Diss'][28:,:,:,165])
   
     tropEddy.close()  
     return udissv

def getKEonY(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     utot=np.array(tropEddy.variables['TOTKETEN'][28:,:,:,165])
   
     tropEddy.close()  
     return utot

def getBarocliniconY(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)

     upress=np.array(tropEddy.variables['KEm_dPg'][28:,:,:,165])

   
     tropEddy.close()  
     return upress
 
def getBarotroponY(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)

     upress=np.array(tropEddy.variables['KEm_dEg'][28:,:,:,165])

   
     tropEddy.close()  
     return upress
 
def getAdvonY(dataNetCDF):
    
    tropEddy=Dataset(dataNetCDF)

    data=np.array(tropEddy.variables['KEm_Adve'][28:,:,:,165])

   
    tropEddy.close()  
    return data
# =============================================================================
# 
# =============================================================================
def getDissVonX(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     udissv=np.array(tropEddy.variables['KEm_Disv'][28:,:,165,:])
   
     tropEddy.close()  
     return udissv
 
def getDissHonX(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     udissv=np.array(tropEddy.variables['KEm_Diss'][28:,:,165,:])
   
     tropEddy.close()  
     return udissv

def getKEonX(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     utot=np.array(tropEddy.variables['TOTKETEN'][28:,:,165,:])
   
     tropEddy.close()  
     return utot

def getBarocliniconX(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)

     upress=np.array(tropEddy.variables['KEm_dPg'][28:,:,165,:])

   
     tropEddy.close()  
     return upress

def getBarotroponX(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)

     upress=np.array(tropEddy.variables['KEm_dEg'][28:,:,165,:])

   
     tropEddy.close()  
     return upress
 
def getAdvonX(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
    
     data=np.array(tropEddy.variables['KEm_Adve'][28:,:,165,:])
    
       
     tropEddy.close()  
     return data
 
def getWonX(dataNetCDF):
    tropEddy=Dataset(dataNetCDF)
  
    wij=np.array(tropEddy.variables['W'][28:,:,165,:])


    tropEddy.close()
    return wij

def getWonY(dataNetCDF):
    tropEddy=Dataset(dataNetCDF)
  
    wij=np.array(tropEddy.variables['W'][28:,:,:,165])


    tropEddy.close()
    return wij


def getWbonX(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)

     wb2=np.array(tropEddy.variables['WRHOMAS2'][28:,:,165,:]) 
             
     tropEddy.close()
     return wb2
# =============================================================================
# 
# =============================================================================
def getDissVonZ(dataNetCDF,zzz):
     tropEddy=Dataset(dataNetCDF)
     
     udissv=np.array(tropEddy.variables['KEm_Disv'][28:,zzz,:,:])
   
     tropEddy.close()  
     return udissv
 
def getDissHonZ(dataNetCDF,zzz):
     tropEddy=Dataset(dataNetCDF)
     
     udissv=np.array(tropEddy.variables['KEm_Diss'][28:,zzz,:,:])
   
     tropEddy.close()  
     return udissv

def getKEonZ(dataNetCDF,z):
     tropEddy=Dataset(dataNetCDF)
     
     utot=np.array(tropEddy.variables['TOTKETEN'][28:,zzz,:,:])
   
     tropEddy.close()  
     return utot

def getRhoGWonZ(dataNetCDF,z):
     tropEddy=Dataset(dataNetCDF)

     upress=np.array(tropEddy.variables['KEm_dPg'][28:,zzz,:,:])
   
     tropEddy.close()  
     return upress
# =============================================================================
# 
# =============================================================================
def getFext(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)

     uext=np.array(tropEddy.variables['KEm_Ext'][28:,0,:,:])
   
     tropEddy.close()  
     return uext
 
def getWVEL(dataNetCDF):
    tropEddy=Dataset(dataNetCDF)
  
    wij=np.array(tropEddy.variables['WVEL'][:,:,:,:])


    tropEddy.close()
    return wij

def getUVVEL(dataNetCDF):
    tropEddy=Dataset(dataNetCDF)
  
    ujk=np.array(tropEddy.variables['UVEL'][:,:,:,:])
    vik=np.array(tropEddy.variables['VVEL'][:,:,:,:])


    tropEddy.close()
    return ujk, vik

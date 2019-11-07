#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:59:44 2019

@author: ydurand
"""


import numpy as np
from netCDF4 import Dataset

def meshcons(dataNetCDF,tstart=-1):
     tropEddy=Dataset(dataNetCDF)
     
     utot=np.array(tropEddy.variables['TOTKETEN'][tstart,:,:,:])
     uadv=np.array(tropEddy.variables['KEm_Adve'][tstart,:,:,:])
     upress=np.array(tropEddy.variables['KEm_dPg'][tstart,:,:,:])
     u_eta=np.array((tropEddy.variables['KEm_dEg'][tstart,:,:,:]))
     udissh=np.array(tropEddy.variables['KEm_Diss'][tstart,:,:,:])
     udissv=np.array(tropEddy.variables['KEm_Disv'][tstart,:,:,:])
     uext=np.array(tropEddy.variables['KEm_Ext'][tstart,:,:,:])
     u_ab=np.array(tropEddy.variables['KEAB_g'][tstart,:,:,:])
   
     tropEddy.close()  
     return utot, uadv, upress, u_eta, udissh, udissv, uext, u_ab

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

def getUPRESS(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)

     upress=np.array(tropEddy.variables['KEm_dPg'][tstart,:,:,:])
   
     tropEddy.close()  
     return upress


def verticalGrid(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     """
     print ("Dimension of the data")
     print (tropEddy.dimensions.keys())  
    
     print ("interrogate variable")
     print (tropEddy.variables.keys())
     print (tropEddy.variables['HFacC'])
     """
     hFacC=np.array(tropEddy.variables['HFacC'][:,:,:])
     drF=np.array(tropEddy.variables['drF'][:])
     
     tropEddy.close()
     return hFacC,drF
     

def diagOcnTave(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     tropEddy.close()
     return

def diagOcnSnap(dataNetCDF):
     tropEddy=Dataset(dataNetCDF)
     
     tropEddy.close()
     return

def dataDiagPE(dataNetCDF,t=-1):
     tropEddy=Dataset(dataNetCDF)
     # print("\n Pour la variable ADVx_TH ")
     # print (tropEddy.variables['ADVx_TH'])

     dtdt=np.array(tropEddy.variables['TOTTTEND'][t,:,:,:])
     #dsdt=np.array(tropEddy.variables['TOTSTEND'][t,:,78:-78,78:-78])
     dsdt=0
     
     advrt=np.array(tropEddy.variables['ADVr_TH'][t,:,:,:])
     advxt=np.array(tropEddy.variables['ADVx_TH'][t,:,:,:])
     advyt=np.array(tropEddy.variables['ADVy_TH'][t,:,:,:])
     
     #advrs=np.array(tropEddy.variables['ADVr_SLT'][t,:,78:-78,78:-78])
     #advxs=np.array(tropEddy.variables['ADVx_SLT'][t,:,78:-78,78:-78])
     #advys=np.array(tropEddy.variables['ADVy_SLT'][t,:,78:-78,78:-78]) 
     advrs,advxs,advys=0,0,0

     tropEddy.close()
     return dtdt, dsdt, advrt, advxt, advyt, advrs, advxs, advys

def getTemp(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)
     
     temp=np.array((tropEddy.variables['THETA'][tstart,:,:,:]))
     
     tropEddy.close()
     return temp

def getWb(dataNetCDF,tstart,ix=-1):
     tropEddy=Dataset(dataNetCDF)
               
     if ix!=-1:
          wb2=np.array(tropEddy.variables['WRHOMAS2'][tstart,:,:,ix]) 
          wb=np.array(tropEddy.variables['WRHOMASS'][tstart,:,:,ix]) 
     else :
          wb2=np.array(tropEddy.variables['WRHOMAS2'][tstart,:,:,:]) 
          wb=np.array(tropEddy.variables['WRHOMASS'][tstart,:,:,:])           
          
          
     tropEddy.close()
     return wb2, wb

def getTOT(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)

     
     dtdt=np.array(tropEddy.variables['TOTTTEND'][tstart,:,:,:])
     dsdt=np.array(tropEddy.variables['TOTSTEND'][tstart,:,:,:])

     tropEddy.close()
     return dtdt, dsdt

def getEta(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)
     
     eta=np.array(tropEddy.variables['Eta'][tstart,:,:])
     
     tropEddy.close()
     return eta

def getDoubleV(dataNetCDF,tstart):
    tropEddy=Dataset(dataNetCDF)
    
    wij=np.array(tropEddy.variables['W'][tstart,:,:,:])
    
    return wij



def dataDiagPETemp(dataNetCDF,t=-1,x=-1):
     tropEddy=Dataset(dataNetCDF)
     # print("\n Pour la variable ADVx_TH ")
     # print (tropEddy.variables['ADVx_TH'])
     
     if t==-1 and x==-1:
          dtdt=np.array(tropEddy.variables['TOTTTEND'][:,:,:,:])
          
          advrt=np.array(tropEddy.variables['ADVr_TH'][:,:,:,:])
          advxt=np.array(tropEddy.variables['ADVx_TH'][:,:,:,:])
          advyt=np.array(tropEddy.variables['ADVy_TH'][:,:,:,:])
          


     elif t!=-1 and x==-1:
          dtdt=np.array(tropEddy.variables['TOTTTEND'][t,:,:,:])
          
          advrt=np.array(tropEddy.variables['ADVr_TH'][t,:,:,:])
          advxt=np.array(tropEddy.variables['ADVx_TH'][t,:,:,:])
          advyt=np.array(tropEddy.variables['ADVy_TH'][t,:,:,:])
          
          
     elif t!=-1 and x!=-1:
          dtdt=np.array(tropEddy.variables['TOTTTEND'][t,:,:,x])
          
          advrt=np.array(tropEddy.variables['ADVr_TH'][t,:,:,x])
          advxt=np.array(tropEddy.variables['ADVx_TH'][t,:,:,x])
          advyt=np.array(tropEddy.variables['ADVy_TH'][t,:,:,x])
          


     tropEddy.close()
     return dtdt, advrt, advxt, advyt


def dataDiagPE3(dataNetCDF,t):
     tropEddy=Dataset(dataNetCDF)

     dtdt=np.array(tropEddy.variables['TOTTTEND'][t,:,:,:])
     dsdt=np.array(tropEddy.variables['TOTSTEND'][t,:,:,:])

          
     tropEddy.close()
     return dtdt, dsdt

def meshcons2(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)

     utot=np.array(tropEddy.variables['TOTKETEN'][tstart,:,:,:])
     uadv=np.array(tropEddy.variables['KEm_Adve'][tstart,:,:,:])
     upress=np.array(tropEddy.variables['KEm_dPg'][tstart,:,:,:])
     udissh=np.array(tropEddy.variables['KEm_Diss'][tstart,:,:,:])
     udissv=np.array(tropEddy.variables['KEm_Disv'][tstart,:,:,:])
     uext=np.array(tropEddy.variables['KEm_Ext'][tstart,:,:,:])

     return utot, uadv, upress, udissh, udissv, uext

def getTime(dataNetCDF):
    tropEddy=Dataset(dataNetCDF)
    tday,nbIter=np.array(tropEddy.variables['T']),len(tropEddy.dimensions['T'])
    return nbIter,tday


def getRHOAnoma(dataNetCDF,t):
    tropEddy=Dataset(dataNetCDF)
    RHOAnoma=np.array(tropEddy.variables['RHOAnoma'][t,:,:,:])
    return RHOAnoma
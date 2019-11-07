#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 08:42:14 2019

@author: durand
"""

import numpy as np
from netCDF4 import Dataset


def meshconsState(dataNetCDF,tstart=0,tend=0):
    tropEddy=Dataset(dataNetCDF)
    '''
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
    '''

    
    xxx,res=tropEddy.variables['X'][:],len(tropEddy.dimensions['X'])
    xxx1=tropEddy.variables['Xp1'][:]
    yyy=tropEddy.variables['Y'][:]
    zzt,dimZZT=tropEddy.variables['Zl'][:],len(tropEddy.dimensions['Zl'])
    zzc=np.array(tropEddy.variables['Z'][:])
    time,nbIter=tropEddy.variables['T'][:],len(tropEddy.dimensions['T'])
    
    
    if tstart ==0 :
        time,nbIter=np.array(tropEddy.variables['T'][:]),len(tropEddy.dimensions['T'])
    else :
        time,nbIter=tropEddy.variables['T'][tstart:tend],len(tropEddy.dimensions['T'])
        
        
    tropEddy.close()        
    return xxx,xxx1,yyy,zzt,dimZZT,zzc,time,nbIter,res


def vitesse(dataNetCDF,zzt,HHH=-1,tstart=-1,tend=-1):
    tropEddy=Dataset(dataNetCDF)
    zzc=np.array(tropEddy.variables['Z'][:])

    H=HHH
    if tstart == -1 and H == -1 :
        ujk=np.array(tropEddy.variables['U'][:,:,:,:])
        vik=np.array(tropEddy.variables['V'][:,:,:,:])
        wij=np.array(tropEddy.variables['W'][:,:,:,:])
        temp=np.array((tropEddy.variables['Temp'][:,:,:,:]))
        eta=np.array(tropEddy.variables['Eta'][:,:,:])
    elif H== -1 and tstart != -1 and tend == -1 :
        ujk=np.array(tropEddy.variables['U'][tstart,:,:,:])
        vik=np.array(tropEddy.variables['V'][tstart,:,:,:])
        wij=np.array(tropEddy.variables['W'][tstart,:,:,:])
        temp=np.array((tropEddy.variables['Temp'][tstart,:,:,:]))
        eta=np.array(tropEddy.variables['Eta'][tstart,:,:])
    elif tstart == -1 and H != -1 :
        ujk=np.array(tropEddy.variables['U'][:,H,:,:])
        vik=np.array(tropEddy.variables['V'][:,H,:,:])
        wij=np.array(tropEddy.variables['W'][:,H,:,:])
        temp=np.array((tropEddy.variables['Temp'][:,H,:,:]))
        eta=np.array(tropEddy.variables['Eta'][:,:,:])

    hhd=np.array(np.diff(zzt))
    hhd=np.append(hhd, [(zzt[19]-zzc[19])*2],axis=0)
    tropEddy.close()
    return hhd,ujk,vik,wij,temp,eta
    
def dimGridState(dataNetCDF):
     tropEddy=Dataset(dataNetCDF) 
     
     rAz=np.array(tropEddy.variables['rA'][:,:])
     zzt,dimZZT=np.array(tropEddy.variables['Zl'][:]),len(tropEddy.dimensions['Zl'])
     zzl=np.array(tropEddy.variables['Zu'][:])
     
     diffZZ=abs(np.diff(np.array(tropEddy.variables['Zp1'][:])))
     
     return rAz,zzt,zzl,diffZZ,dimZZT

def getTemp(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)
     
     temp=np.array((tropEddy.variables['Temp'][tstart,:,:,:]))
     
     tropEddy.close()
     return temp

def getVarVor(dataNetCDFS,tstart):
     tropEddy=Dataset(dataNetCDFS)
     ujk=np.array(tropEddy.variables['U'][tstart,:,:,:])
     vik=np.array(tropEddy.variables['V'][tstart,:,:,:])
     
     tropEddy.close()
     return ujk, vik
     
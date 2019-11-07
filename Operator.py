#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 2019

@author: durand
"""
import numpy as np
from grille import *


#%%

#########################################################################
#----------------Grid and operator---------------------------------------
#########################################################################
res=250

def partui(ujk,dx):
    diffu=np.diff(ujk,axis=2)
    diffu=diffu[:,:]/dx
    return diffu

def partui1L(ujk,dx):
    diffu=np.diff(ujk,axis=1)
    diffu=diffu[:,:]/dx
    return diffu

def partvj(vik,dy):
    diffv=np.diff(vik,axis=1)
    diffv=diffv[:,:]/dy

    return diffv

def partvj1L(vik,dy):
    diffv=np.diff(vik,axis=0)
    diffv=diffv[:,:]/dy

    return diffv

def partwz(wij,dz):
    diffz=np.diff(wij,axis=0)
    diffz=diffz/dz[:-1]
    
    return diffz

def partwz02(wij,dz,zzz,yyy,xxx):
    diffz=np.zeros((zzz,yyy,xxx))
    for i in range(1,zzz-2):
        diffzBack=(wij[i,:,:]-wij[i-1,:,:])/(dz[(i),:,:]-dz[i-1,:,:])
        diffzForw=(-wij[i,:,:]+wij[i+1,:,:])/(dz[(i+1),:,:]-dz[i,:,:])
        diffz[i,:,:]=(diffzBack+diffzForw)/2.
    diffz[0,:,:]=-(wij[0,:,:]-wij[1,:,:])/(dz[0,:,:]-dz[1,:,:])
    diffz[zzz-1,:,:]=(wij[zzz-2,:,:]-wij[zzz-1,:,:])/(dz[zzz-2,:,:]-dz[zzz-1,:,:])
    
    return diffz

def partuj(ujk,dy):
    difuj=np.diff(ujk,axis=1)
    difuj=difuj[::,::,1:-1]
    return difuj

def partvi(vik,dx):
    difvi=np.diff(vik,axis=2)
    difvi=difvi[::,1:-1:,::]
    return difvi
        

def ucc(ujk):
    uci=(ujk[:,:,:-1]+ujk[:,:,1:])/2.
    return uci

def vcc(vik):
    vcj=(vik[:,:-1,:]+vik[:,1:,:])/2.
         
    return vcj

def ucc1L(ujk):
    uci=(ujk[:,:-1]+ujk[:,1:])/2.
    return uci

def vcc1L(vik):
    vcj=(vik[:-1,:]+vik[1:,:])/2.
         
    return vcj
#def lapla()
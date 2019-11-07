#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:54:05 2019

@author: durand
"""


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from stockData import *

#%%
def speedRadius(dataNetCDFS,t,dx):
    eta=getEta(dataNetCDFS,t)
    ujk,vik=getSpeed(dataNetCDFS,t)
    ujk=ujk[:,:]
    vik=vik[:,:]
    eta=eta[:,:]
    
    
    ####---- Find the maximum of Eta ----####
    maxEta=np.unravel_index(np.argmax(eta, axis=None), eta.shape)
    
    ####---- Find the maximum and the minimum of the speed ----####
    maxUUU=np.unravel_index(np.argmax(ujk, axis=None), ujk.shape)
    maxVVV=np.unravel_index(np.argmax(vik, axis=None), vik.shape)
    minUUU=np.unravel_index(np.argmin(ujk, axis=None), ujk.shape)
    minVVV=np.unravel_index(np.argmin(vik, axis=None), vik.shape)
    
    maxUUUy,maxUUUx=maxUUU[0],maxUUU[1]
    maxVVVy,maxVVVx=maxVVV[0],maxVVV[1]
    minUUUy,minUUUx=minUUU[0],minUUU[1]
    minVVVy,minVVVx=minVVV[0],minVVV[1]
    
    maxEtay,maxEtax=maxEta[0],maxEta[1]
    ####---- Compute the distance by the mean----####
    speedRadius1=(((maxUUUx-maxEtax)*dx)**2+((maxUUUy-maxEtay)*dx)**2)**(1/2)
    speedRadius2=(((minUUUx-maxEtax)*dx)**2+((minUUUy-maxEtay)*dx)**2)**(1/2)
    speedRadius3=(((maxVVVx-maxEtax)*dx)**2+((maxVVVy-maxEtay)*dx)**2)**(1/2)
    speedRadius4=(((minVVVx-maxEtax)*dx)**2+((minVVVy-maxEtay)*dx)**2)**(1/2)
    
    speedRadius=(1/4)*(speedRadius1+speedRadius2+speedRadius3+speedRadius4)

    return speedRadius

def speedRadiusLevel(dataNetPhy,dataNetCDFS,t,dx):
    phy2=getPhy(dataNetPhy,t)
    ujk2,vik2=getSpeed2(dataNetCDFS,t)

    vertMaxPhy=np.arange(60).reshape(30,2)
    maxPhyArray=np.array([])
    maxPhyValue=np.array([])
    verticalSpeedRadius=np.array([])
    
    
    for j in range(0,30):
        ujk=ujk2[j,40:-40,40:-40]
        vik=vik2[j,40:-40,40:-40]
        phy=phy2[j,40:-40,40:-40]
        
        ####---- Find the maximum of Phy ----####
        maxPhy=np.unravel_index(np.argmax(phy, axis=None), phy.shape)
        maxPhyValue=np.append(maxPhyValue,np.max(phy[j]))
        maxPhyy,maxPhyx=maxPhy[0],maxPhy[1]
        vertMaxPhy[j,0],vertMaxPhy[j,1]=maxPhyy,maxPhyx
        #,maxPhyx
        
        ####---- Find the maximum and the minimum of the speed ----####
        maxUUU=np.unravel_index(np.argmax(ujk, axis=None), ujk.shape)
        maxVVV=np.unravel_index(np.argmax(vik, axis=None), vik.shape)
        minUUU=np.unravel_index(np.argmin(ujk, axis=None), ujk.shape)
        minVVV=np.unravel_index(np.argmin(vik, axis=None), vik.shape)
        
            
        maxUUUy,maxUUUx=maxUUU[0],maxUUU[1]
        maxVVVy,maxVVVx=maxVVV[0],maxVVV[1]
        minUUUy,minUUUx=minUUU[0],minUUU[1]
        minVVVy,minVVVx=minVVV[0],minVVV[1]
        
        ####---- Compute the speed Radius at each level ----####
        speedRadius1=(((maxUUUx-maxPhyx)*dx)**2+((maxUUUy-maxPhyy)*dx)**2)**(1/2)
        speedRadius2=(((minUUUx-maxPhyx)*dx)**2+((minUUUy-maxPhyy)*dx)**2)**(1/2)
        speedRadius3=(((maxVVVx-maxPhyx)*dx)**2+((maxVVVy-maxPhyy)*dx)**2)**(1/2)
        speedRadius4=(((minVVVx-maxPhyx)*dx)**2+((minVVVy-maxPhyy)*dx)**2)**(1/2)
        
        speedRadius=(1/4)*(speedRadius1+speedRadius2+speedRadius3+speedRadius4)
        
        verticalSpeedRadius=np.append(verticalSpeedRadius,speedRadius)

    return verticalSpeedRadius

def getCenterEta(dataNetCdf,t):
     eta=getEta(dataNetCdf,t)
     maxEta=np.arange(1).reshape(1,1)
                
     maxEta2=np.unravel_index(np.argmax(eta, axis=None), eta.shape)
     #maxEta[0,0]=np.array([maxEta2[0],maxEta2[1]])
     return maxEta2


def depthLevel(dataNetPhy,dataNetCDFS,t,dx):
    phy2=getPhy(dataNetPhy,t)
    ujk2,vik2=getSpeed2(dataNetCDFS,t)

    vertMaxPhy=np.arange(60).reshape(30,2)
    maxPhyArray=np.array([])
    maxPhyValue=np.array([])
    verticalSpeedRadius=np.array([])
    
    
    for j in range(0,30):
        ujk=ujk2[j,40:-40,40:-40]
        vik=vik2[j,40:-40,40:-40]
        phy=phy2[j,40:-40,40:-40]
        
        ####---- Find the maximum of Phy ----####
        maxPhy=np.unravel_index(np.argmax(phy, axis=None), phy.shape)
        maxPhyValue=np.append(maxPhyValue,np.max(phy[j]))
        maxPhyy,maxPhyx=maxPhy[0],maxPhy[1]
        vertMaxPhy[j,0],vertMaxPhy[j,1]=maxPhyy,maxPhyx
        #,maxPhyx
        
        ####---- Find the maximum and the minimum of the speed ----####
        maxUUU=np.unravel_index(np.argmax(ujk, axis=None), ujk.shape)
        maxVVV=np.unravel_index(np.argmax(vik, axis=None), vik.shape)
        minUUU=np.unravel_index(np.argmin(ujk, axis=None), ujk.shape)
        minVVV=np.unravel_index(np.argmin(vik, axis=None), vik.shape)
        
            
        maxUUUy,maxUUUx=maxUUU[0],maxUUU[1]
        maxVVVy,maxVVVx=maxVVV[0],maxVVV[1]
        minUUUy,minUUUx=minUUU[0],minUUU[1]
        minVVVy,minVVVx=minVVV[0],minVVV[1]
        
        ####---- Compute the speed Radius at each level ----####
        speedRadius1=(((maxUUUx-maxPhyx)*dx)**2+((maxUUUy-maxPhyy)*dx)**2)**(1/2)
        speedRadius2=(((minUUUx-maxPhyx)*dx)**2+((minUUUy-maxPhyy)*dx)**2)**(1/2)
        speedRadius3=(((maxVVVx-maxPhyx)*dx)**2+((maxVVVy-maxPhyy)*dx)**2)**(1/2)
        speedRadius4=(((minVVVx-maxPhyx)*dx)**2+((minVVVy-maxPhyy)*dx)**2)**(1/2)
        
        speedRadius=(1/4)*(speedRadius1+speedRadius2+speedRadius3+speedRadius4)
        verticalSpeedRadius=np.append(verticalSpeedRadius,speedRadius)
        
    k=1
    while verticalSpeedRadius[k]<verticalSpeedRadius[k-1]+13e3:
         depth=k-1
         k=k+1
          
    return depth


def depthLevel2(dataNetCDFS,t,dx):
     ujk,vik=getSpeed2(dataNetCDFS,t)
     www=vitesseW(dataNetCDFS,t)
     speedMax=np.max(ujk)
     
     maxUUU=np.unravel_index(np.argmax(ujk, axis=None), ujk.shape)
     maxUUUy,maxUUUx=maxUUU[1],maxUUU[2]
     
     k=0
     depth2=0
     while abs(ujk[k,maxUUUy,maxUUUx])>0.37*speedMax:
          depth2=k
          k=k+1
     
     return depth2
     
     
     
     
     
     
     
     
     
     
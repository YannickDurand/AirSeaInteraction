#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 24 19:28:20 2019

@author: durand
"""

from Operator import *
import numpy as np
#%%


##########################################################################
#-------------Cell for subroutine-----------------------------------------
##########################################################################

#####-----Vorticty-----#####
def vort_rel(dx, dy, ujk, vik):
    qij=-partuj(ujk,dy)/dy
    qij=qij+partvi(vik,dx)/dx
    return qij

#####-----Divergence-----#####
def conv(dx, dy, ujk, vik):
    conv=partui(ujk,res,dx)+partvj(vik,res,dy)
    return conv


#####-----Kinetic Energy-----#####
def energyK(ujk, vik, hik, tstart,tend=-1):
    uck,vck=ucc(ujk),vcc(vik)
    #print(np.max(uck))
    if tstart==-1:
         Kx = (uck[:,:,:,:]*uck[:,:,:,:])*abs(hik[np.newaxis,:,np.newaxis,np.newaxis])
         Ky = (vck[:,:,:,:]*vck[:,:,:,:])*abs(hik[np.newaxis,:,np.newaxis,np.newaxis])
    elif tstart !=-1 and tend==-1 :
         Kx = (uck[:,:,:]*uck[:,:,:])#*abs(hik[:,np.newaxis,np.newaxis])
         Ky = (vck[:,:,:]*vck[:,:,:])#*abs(hik[:,np.newaxis,np.newaxis])
    else :
         Kx = (uck[tstart:tend,:,:,:]*uck[tstart:tend,:,:,:])*abs(hik[np.newaxis,:,np.newaxis,np.newaxis])
         Ky = (vck[tstart:tend,:,:,:]*vck[tstart:tend,:,:,:])*abs(hik[np.newaxis,:,np.newaxis,np.newaxis])
    return (Kx+Ky)*0.5

#####-----Potential Energy-----#####
def energyP(eta, hik):
    Pe=9.81*(1./2.)*(hik[np.newaxis,:,np.newaxis,np.newaxis]-eta[:,np.newaxis,:,:])**(2)
    return Pe

#####-----Mecanic Energy-----#####
def sumEnergyK(ujk, vil, hik):
    KxT = (ujk*ujk)*hik.sum()
    KyT = (vil*vil)*hik.sum()
    PeT=((1./2.)*(hik)*(hik)).sum()
    ekT = KxT + KyT + PeT
    return ekT

def sumEnergyP(hik):
    PeT=((1./2.)*(hik)*(hik)).sum()
    return PeT

#####-----Potential Enstrophy-----#####
def potentialEnstrophy(qik, hik):
    qikt=np.roll(qik,-1,axis=1)
    qikt=np.roll(qikt,-1,axis=2)
    potEns2=qik+np.roll(qik,-1,axis=1)+np.roll(qik,-1,axis=2)+qikt
    potEns2=potEns2[::,:-1:,:-1:]
    return potEns2/4.*hik[:,np.newaxis,np.newaxis]

def sumPE(mesh,qik,hik):
    Zt=(-mesh.mi(qik*qik)*hik*(1./2.)).sum()
    return Zt

def Rho(Rho0,maxT,Temp, Salinity=0):
    coeffAlpha=2.E-4
    
    diffTemp=Temp-maxT
    rho=Rho0-Rho0*coeffAlpha*diffTemp
    return rho

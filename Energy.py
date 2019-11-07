#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:33:59 2019

@author: ydurand
"""

import numpy as np

from Sauvegarde import *
from fonctionH import *
from lesPlots import *

#%%
###############################################################################
#-------------------------Cell for NetCDF------------------------------------
###############################################################################

PATH='/media/ydurand/Yannick/19_OAC/Control_Case/Case_withoutWqT_ACyclone'
PATH='/media/ydurand/Yannick/19_OAC/Control_Case/Case_wh_qT_withW_ACyclone'


fileS='/state.nc'
fileKEs='/diagKEs.nc'
dataNetCDFU=PATH+fileS

#%%
###############################################################################
#------------------------Cell to stock data------------------------------------
###############################################################################

xxx,xxx1,yyy,zzt,dimZZT,zzc,time,nbIter,res=meshcons(dataNetCDFU,tstart=0,tend=10)
hhd,ujk,vik,wij,temp,eta=vitesse(dataNetCDFU,zzt,HHH=-1,tstart=0,tend=10)

#%%
###############################################################################
#--------------------------Calculation-----------------------------------------
###############################################################################

EKs=energyK(ujk, vik, hhd)
EPs=energyP(eta,hhd)

EMs=EKs+EPs

partEKs=(np.roll(EKs,-1,axis=0)-np.roll(EKs,1,axis=0))/(86400.)


#%%
###############################################################################
#----------------------------Plots---------------------------------------------
###############################################################################
HHH=1
t=7

plotEk(xxx,yyy,partEKs[t,HHH,:,:],t)
plotEp(xxx,yyy,EPs[t,HHH,:,:],t)
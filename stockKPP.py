#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 09:01:39 2019

@author: durand
"""


import numpy as np
from netCDF4 import Dataset

def getDepthBoundary(dataNetCDF,tstart):
     tropEddy=Dataset(dataNetCDF)

     upress=np.array(tropEddy.variables['KPPHBL'][tstart,:,:])
   
     tropEddy.close()  
     return upress
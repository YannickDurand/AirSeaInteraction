#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 10:12:04 2019

@author: durand
"""
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

 
    
    
def plotEta(xxx,yyy,eta,t):            
    #####-----Eta-----#####
    plt.figure(figsize=(12.5,10))
    plt.contourf(xxx,yyy, eta, cmap=cm.jet, levels=100)
    plt.colorbar()
    plt.title('Eta' + ' ; t='+str(t))
    plt.show()

def plotVort(x,y,qij,t):
    #####-----Voritcity------#####
    plt.figure(figsize=(12.5,10))
    plt.contourf(x,y, qij, cmap=cm.jet, levels=100)
    plt.colorbar()
    plt.title('vorticite potentielle' + ' ; t='+str(t))
    plt.show()
    #plt.close()
    
def plotDiv(xxx,yyy,diver,t):
    #####-----divergence-----#####
    plt.figure(figsize=(12.5,10))
    plt.contourf(xxx,yyy,diver,cmap=cm.jet, levels=100)
    plt.colorbar()
    plt.title('Convergence' + ' ; t='+str(t))
    plt.show()

def plotEk(xxx,yyy,Ek,t):
    #####-----Kinetic Energy-----#####
    plt.figure(figsize=(12.5,10))
    plt.contourf(xxx,yyy,Ek,cmap=cm.jet, levels=100)
    plt.colorbar()
    plt.title('Kinetic Energy' + ' ; t='+str(t))
    plt.show()
    
def plotEp(xxx,yyy,enPot,t):
    #####-----Potential Energy-----#####
    plt.figure(figsize=(12.5,10))
    plt.contourf(xxx,yyy,enPot,cmap=cm.jet, levels=100)
    plt.colorbar()
    plt.title('Potential Energy' + ' ; t='+str(t))
    plt.show()
    
def plotEns(x,y,potEns,t):
    #####-----Potential Enstrophy-----#####
    plt.figure(figsize=(12.5,10))
    plt.contour(potEns,20)
    plt.colorbar()
    plt.title('Potential Enstrophy ' + ' ; t='+str(t))
    plt.show()  

    
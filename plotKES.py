#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 10:16:01 2019

@author: ydurand
"""

import matplotlib.pyplot as plt

def plot(i,time,A,nameA,title):
     plt.figure(i)
     plt.plot(time[:],A, label = str(nameA))
     
     plt.title(str(title))
     plt.xlabel("t [ day ]")
     plt.ylabel("Energy variation [ kg.m$^{-2}$.s$^{-2}$ ]")
     plt.grid(True)
     plt.legend()
     plt.show()


def plotComp(i,time,A,B,C,nameA,nameB,nameC):
     plt.figure(i)
     plt.plot(time[1:],A, label = str(nameA))
     plt.plot(time[1:],B, label = str(nameB))
     plt.plot(time[1:],C, label = str(nameC))
     
     plt.title("Evolution in time")
     plt.xlabel("t [ day ]")
     plt.ylabel("Energy variation [ kg.m$^{-2}$.s$^{-3}$ ]")
     plt.grid(True)
     plt.legend()
     plt.show()
     
     

def plotCompT(i,time,A,B,C,nameA,nameB,nameC,title):
     plt.figure(i)
     plt.plot(time[1:],A, label = str(nameA))
     plt.plot(time[1:],B, label = str(nameB))
     plt.plot(time[1:],C, label = str(nameC))
     
     plt.title(str(title))
     plt.xlabel("t [ day ]")
     plt.ylabel("Energy variation [ kg.m$^{-2}$.s$^{-2}$ ]")
     plt.grid(True)
     plt.legend()
     plt.show()
     
def plotComp2T(i,time,A,B,nameA,nameB,title):
     plt.figure(i)
     plt.plot(time[1:],A, label = str(nameA))
     plt.plot(time[1:],B, label = str(nameB))
     
     plt.title(str(title))
     plt.xlabel("t [ day ]")
     plt.ylabel("Energy variation [ kg.m$^{-2}$.s$^{-2}$ ]")
     plt.grid(True)
     plt.legend()
     plt.show()
     
def plotComp4T(i,time,A,B,C,D,nameA,nameB,nameC,nameD,title):
     plt.figure(i)
     plt.plot(time[1:],A, label = str(nameA))
     plt.plot(time[1:],B, label = str(nameB))
     plt.plot(time[1:],C, label = str(nameC))
     plt.plot(time[1:],D, label = str(nameD))
     
     plt.title(str(title))
     plt.xlabel("t [ day ]")
     plt.ylabel("Energy variation [ kg.m$^{-2}$.s$^{-2}$ ]")
     plt.grid(True)
     plt.legend()
     plt.show()
     
def plotCompSauv(i,time,A,B,C,nameA,nameB,nameC,Folder):
     plt.figure(i)
     plt.plot(time[1:],A, label = str(nameA))
     plt.plot(time[1:],B, label = str(nameB))
     plt.plot(time[1:],C, label = str(nameC))
     
     plt.title("Evolution in time")
     plt.xlabel("t [ day ]")
     plt.ylabel("Energy variation [ kg.m$^{-2}$.s$^{-3}$ ]")
     plt.grid(True)
     plt.legend()
     plt.savefig(str(Folder)+str(nameA)+'_'+str(nameB)+'_'+str(nameC)+'s.png')
     
     

def plotCompTSauv(i,time,A,B,C,nameA,nameB,nameC,Folder):
     plt.figure(i)
     plt.plot(time[1:],A, label = str(nameA))
     plt.plot(time[1:],B, label = str(nameB))
     plt.plot(time[1:],C, label = str(nameC))
     
     plt.title("Evolution in time")
     plt.xlabel("t [ day ]")
     plt.ylabel("Energy variation [ kg.m$^{-2}$.s$^{-2}$ ]")
     plt.grid(True)
     plt.legend()
     plt.savefig(str(Folder)+str(nameA)+'_'+str(nameB)+'_'+str(nameC)+'s.png')
     

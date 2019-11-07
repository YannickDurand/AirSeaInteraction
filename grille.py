#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Contient les outils de définition de la grille et opérateurs associés."""


import numpy as np
from scipy import sparse as sparse


#Permet de retourner un vecteur (pas 'reshape')
def col_vec(data):
    vec = np.zeros((data.size, 1))
    vec[:, 0] = data
    return vec
def row_vec(data):
    vec = np.zeros((1, data.size))
    vec[0, :] = data
    return vec


####################### Horizontal mesh #########################

class Periodic_FD:
    """Contient grille 1D et opérateurs 1D diff et moyennage (périodiques) sur grille 2D."""
    # staggering
    # m0  m1  m2  m3
    #   u0  u1  u2  u3
    def __init__(self, N):
        """ Construit les opérateurs de différentiation et de moyennage de la grille périodique.

        Définit les variables self.N, self.xi, self.onesi, self.zeroi, self.xj, self.onesj,
           self.onesk, self.onesl, self.xk, self.xl
        et les fonctions self.di, self.mi, self.dj, self.mj
        """

        self.N = N
        ones = np.ones(N)

        ### Définition des opérateurs matriciels de différentiations horizontale et de moyennages ###
        ### Et des vecteurs définissant les coordonnées de la grille     ###
        right = sparse.spdiags([ones], [0], N, N, format='csr') #matrice unité
        left = sparse.spdiags([ones], [-1], N, N, format='csr') #matrice sous-diagonale
        left[0, -1] = 1
        di = right-left
        mi = .5*(right+left)
        self.di = lambda f: di.dot(f)
        self.mi = lambda f: mi.dot(f)
        self.xi = col_vec(np.linspace(0, N-1, N))-N/2.
        self.onesi = np.ones((N, 1))
        self.zeroi = np.zeros((N, 1))

        right = sparse.spdiags([ones], [1], N, N, format='csr')
        right[-1, 0] = 1
        left = sparse.spdiags([ones], [0], N, N, format='csr')
        dj = right-left
        mj = .5*(right+left)
        self.dj = lambda f: dj.dot(f)
        self.mj = lambda f: mj.dot(f)
        self.xj = col_vec(np.linspace(0.5, N-.5, N))-N/2.
        self.onesj = np.ones((N, 1))

        #Définition grille verticale
        self.onesk = np.ones((1, N))
        self.onesl = np.ones((1, N))
        self.xk = row_vec(np.linspace(0, N-1, N))-N/2.
        self.xl = row_vec(np.linspace(0.5, N-.5, N))-N/2.

    ### Définition des opérateurs matriciels de différentiations verticaux et de moyennages ###
    #Outils
    def hstack(self, fields):
        return np.array(np.hstack(fields))
    def up(self, f):
        """Place la dernière collone au début"""
        N = self.N
        return self.hstack((f[:, N-1:N], f[:, :-1]))
    def down(self, f):
        """Place la première collone à la fin"""
        return self.hstack((f[:, 1:], f[:, 0:1]))
    #Définitions
    def dk(self, fl):
        return fl-self.up(fl)
    def mk(self, fl):
        return .5*(fl+self.up(fl))
    def dl(self, fk):
        return self.down(fk)-fk
    def ml(self, fk):
        return .5*(fk+self.down(fk))


############################ x-y rectangular mesh ##########################

class Rectangle:
    """Grille 2D"""
    def __init__(self, xmesh, ymesh):
        """Construit la grille 2D à partir des deux grilles 1D"""
        self.xmesh, self.ymesh = xmesh, ymesh
        self.Nx, self.Ny = xmesh.N, ymesh.N
        self.di, self.dj, self.mi, self.mj = xmesh.di, xmesh.dj, xmesh.mi, xmesh.mj
        self.dk, self.dl, self.mk, self.ml = ymesh.dk, ymesh.dl, ymesh.mk, ymesh.ml

        self.xik = xmesh.xi*ymesh.onesk
        self.xjk = xmesh.xj*ymesh.onesk
        self.xil = xmesh.xi*ymesh.onesl
        self.xjl = xmesh.xj*ymesh.onesl

        self.yik = xmesh.onesi*ymesh.xk
        self.yil = xmesh.onesi*ymesh.xl
        self.yjk = xmesh.onesj*ymesh.xk
        self.yjl = xmesh.onesj*ymesh.xl


####### Sanity checks ########
#Nx, Ny = 10 ,5
#xmesh = Periodic_FD(Nx)
#ymesh = Periodic_FD(Ny)
#mesh = Rectangle(xmesh, ymesh)

#print(xmesh.xi - xmesh.mi(xmesh.xj))
#print(xmesh.xj - xmesh.mj(xmesh.xi))
#print(xmesh.xk - xmesh.mk(xmesh.xl))
#print(xmesh.xl - xmesh.ml(xmesh.xk))
#print('xmesh.di(xmesh.xj)')
#print(xmesh.di(xmesh.xj))
#print('xmesh.dj(xmesh.xi)')
#print(xmesh.dj(xmesh.xi))
#print('xmesh.dk(mesh.yil)')
#print(mesh.dk(mesh.yil))
#print('xmesh.dl(mesh.yik)')
#print(mesh.dl(mesh.yik))

#print(mesh.xik - mesh.mi(mesh.xjk))
#print(mesh.xjk - mesh.mj(mesh.xik))
#print(mesh.yik - mesh.mk(mesh.yil))
#print(mesh.yil - mesh.ml(mesh.yik))
#print(mesh.di(mesh.xjk))
#print(mesh.dj(mesh.xik))
#print(mesh.dk(mesh.yjl))
#print(mesh.dl(mesh.yjk))




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 29 01:55:08 2019

@author: paramveer
"""
import numpy as np
import matplotlib.pyplot as plt
nu = .3
modu = 300000
g = modu/(2*(1+nu))
k = modu/(3*(1-2*nu))
een = np.zeros([3,3])
sn = np.zeros([3,3])
epn = np.zeros([3,3])
sigmay = 200
dsigmay = 20
epeqn = 0
sappend = np.array([])
eappend = np.array([])
ee11 = 0
def sy(ep):
    a = 200
    b = 100
    c = .9
    syield = a + b*(1-np.exp(-c*ep))
    dsyield = b*c*np.exp(-c*ep)
    return syield, dsyield
for i in range(100):

    if i==0:
        dde = 0
    else:
        dde = .00001
        
    ee11 = ee11 + dde
    dele = np.zeros([3,3])
    dele[0,0] = dde
    dele[1,1] = -nu*dde
    dele[2,2] = -nu*dde
    eetrial = een + dele
    eevtrial = np.trace(eetrial)
    eedtrial = eetrial - (1/3.0)*eevtrial
    p = k*eevtrial
    sdevtrial = 2*g*eedtrial
    qtrial = np.sqrt(3/2)*(sdevtrial[0,0]**2 + sdevtrial[1,1]**2 + sdevtrial[2,2]**2 \
                    + 2*(sdevtrial[0,1]**2 + sdevtrial[0,2]**2 + sdevtrial[1,2]**2))**(1/2)
    sdevtrialnorm = (sdevtrial[0,0]**2 + sdevtrial[1,1]**2 + sdevtrial[2,2]**2 \
                + 2*(sdevtrial[0,1]**2 + sdevtrial[0,2]**2 + sdevtrial[1,2]**2))**(1/2)
    sigmay, dsigmay = sy(epeqn)
    
    ftrial = qtrial - sigmay
    if ftrial <= 0:
        snn =  sdevtrial+ p
        epnn = epn
        eenn = eetrial
        epeqnn = epeqn
    else:
        tol = 1*10**-6
        delgama, ddgama = 0, 0#ftrial/(3*g + dsigmay)
        fntrial = qtrial - 3*g*delgama - sigmay
        while abs(fntrial/sigmay)>tol:
            #delgama = delgama +ddgama 
            #eqepnn = epeqn + delgama
            #sigma, dsigmay = sy(eqepnn)
            denom =  - 3*g - dsigmay
            ddgama = -fntrial/denom
            delgama = delgama + ddgama
            epeqnn = epeqn + delgama
            sigmay, dsigmay = sy(epeqnn)
            fntrial = qtrial - 3*g*delgama - sigmay
            print(fntrial)
        factor = (1- delgama*3*g/qtrial)
        sdevnn = factor*sdevtrial
        snn = sdevnn + p
        eenn = (1/(2*g))*sdevnn + (1/3)*eevtrial*np.eye(3)
        epnn =  epn + delgama * sdevtrial/sdevtrialnorm
    een = eenn
    epn = epnn
    epeqn = epeqnn
    sappend = np.append(sappend, snn[0,0])
    eappend = np.append(eappend, ee11)
    plt.plot(eappend, sappend)
    #print(een)
        
        
    
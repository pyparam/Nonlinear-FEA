#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  1 22:38:38 2020

@author: paramveer
"""
import numpy as np
import matplotlib.pyplot as plt
E = 300000
nu = .3
k = E/(3*(1-2*nu))
g = E/(2*(1+nu))

een = np.zeros([3,3])
epn = np.zeros([3,3])
sn = np.zeros([3,3])
dele = np.zeros([3,3])
eqepn = 0
tol = 1*10**(-8)
#sigmay = 200
#H = 10
sarray = np.array([])
earray = np.array([])
e11array = np.array([])
ee11 = 0
#delgama = 0
def sy(ep):
    a = 600
    b = 600
    c = 6
    syield = a + b*(1-np.exp(-c*ep))
    dsyield = b*c*np.exp(-c*ep)
    return syield, dsyield
for i in range(5000):
    
    if i==0:
        dde = 0
    else:
        dde= 0.0002
    ee11 = ee11+ dde
    dele[0,0] = dde
    dele[1,1]= -nu*dde
    dele[2,2]= -nu*dde
    
    eetrial = een + dele
    eevtrial = np.trace(eetrial)
    eedtrial = eetrial - (1/3)*eevtrial*np.eye(3)
    p = k*eevtrial
    sdtrial = 2*g*eedtrial
    sdtrialnorm = np.linalg.norm(sdtrial)
    qtrial = np.sqrt(3/2)*sdtrialnorm
    sigmay, dsigmay = sy(eqepn)
    ftrial = qtrial - sigmay
    if ftrial<tol:
        snn = sdtrial+p*np.eye(3)
        eenn = eetrial
        epnn = epn
        eqepnn = eqepn
        eqvmises = qtrial
    else:
        delgama = 0 #ftrial/(3*g + H)
        for ii in range(100):
            ddgama = ftrial/(3*g + dsigmay)
            delgama = delgama +  ddgama
            eqepnn = eqepn+ delgama
            sigmay, dsigmay = sy(eqepnn)
            ftrial = qtrial - 3*g*delgama - sigmay
            
            if abs(ftrial)<tol:
                break
            if ii>98:
                print('failed', abs(ftrial))
        #eqepnn = eqepn + delgama
        factor = (1- 3*g*delgama/qtrial)
        sdnn = factor*sdtrial
        snn = sdnn + p*np.eye(3)
        eenn = (1/(2*g))*sdnn+ (1/3)*eevtrial*np.eye(3)
        sdnorm = np.linalg.norm(sdnn)
        epnn = epn + np.sqrt(3/2)*delgama*sdnn/(sdnorm)
        eqvmises = np.sqrt(3/2)*sdnorm
        
    sn = snn
    een = eenn
    epn = epnn
    eqepn =  eqepnn
        
    earray= np.append(earray, epnn[0, 0]+ eenn[0,0])
    sarray= np.append(sarray, eqvmises)
plt.plot(earray, sarray, label =r'$\sigma_{eq}= \sqrt{\frac{3}{2} \bf{S:S}}$')
plt.ylabel(r'$\sigma \,  \, \, (MPa)$')
plt.xlabel(r'$\epsilon$')
plt.legend(loc=0) 
        
        
        
        
        
    
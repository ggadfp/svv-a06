# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 18:28:18 2020

@author: halms
"""

import numpy as np
import matplotlib.pyplot as plt

#%%


Nx = 41 # Spanwise points
Nz = 81 # Chordwise points

data = np.zeros((Nz,Nx))
with open('aerodynamicloadcrj700.dat', 'r') as f:
    for i,line in enumerate(f):
        numbers = line.split(',')
        for j,number in enumerate(numbers):
            data[i,j] = number
            
data = data*10**3

def meshCreator(Ncw,Nsw,chord,span,height):
    zmesh = np.zeros((Ncw,Nsw))
    xmesh = np.zeros((Ncw,Nsw))
    theta_z = []
    theta_x = []
    
    for i in range(1,Ncw+2):
        theta_z.append((i-1)*np.pi/Ncw)  
        
    for i in range(1,Nsw+2):
        theta_x.append((i-1)*np.pi/Nsw)
    
    for i in range(1,Ncw+1):
        for j in range(Nsw):
            zmesh[i-1,j-1] = height/2 - 1/2 * (chord/2 * (1-np.cos(theta_z[i-1])) + chord/2 * (1-np.cos(theta_z[i])))
        
    for j in range(Ncw):
        for i in range(1,Nsw+1):
            xmesh[j-1,i-1] = 1/2 * (span/2 * (1-np.cos(theta_x[i-1])) + span/2 * (1-np.cos(theta_x[i])))
        
    return xmesh,zmesh



def zload_zpoint(Ncw,Nsw,zmesh):
    loadlist = []
    load = 0
    for j in range(Nsw):
        for i in range(0,Ncw-2,2): 
            load += (zmesh[i+2,j] - zmesh[i,j])/6 * (data[i,j] + 4*data[i+1,j] + data[i+2,j])
        loadlist.append(load)
        load = 0
        
    z_F = 0
    F = 0
    z_pos = []
    for j in range(Nsw):
        for i in range(Ncw):
            z_F += data[i,j]*zmesh[i,j]
            F += data[i,j]
        z_pos.append(z_F/F)
        z_F = 0
        F = 0
        
    return loadlist,z_pos

###############################################################################

#def lagrange_basis(xval,j,xarray):
#    L = 1
#    for i in range(len(xarray)):
#        if i != j:
#            y = (xval - xarray[i])/(xarray[j] - xarray[i])
#            L *=y
#    return L
#
#def polynomial(xval,xarray,yarray):
#    poly = np.zeros((len(xarray),1))
#    dot = 0
#    for j in range(len(xarray)):
#        poly[j] = lagrange_basis(xval,j,xarray)
#    for k in range(len(xarray)):
#        dot += poly[k]*yarray[k]
#    return dot


#test = np.linspace(x[0,0],x[0,-1])
#pr = np.array([polynomial(x[0,i],x[0,:],torque) for i in range(Nx)])
#plt.figure(1)
#ax = plt.axes(projection='3d')
#ax.scatter(z,x,data)
#
##pr = pr/np.linalg.norm(pr)
#plt.figure(2)
##plt.plot(pr,'o')
#plt.plot(x[0,:],torque,'r')
#
#
#plt.figure(3)
#plt.plot(data[:,0])

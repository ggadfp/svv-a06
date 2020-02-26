# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:03:57 2020

@author: Ihab
"""

import numpy as np
import matplotlib.pyplot as plt
import readAeroload as aero
import Parameters as c
from mpl_toolkits.mplot3d import Axes3D
from shapely import geometry
#%%

def Centroid(choice):
    
    if choice == 'y':
        y_tilde = 0
        return y_tilde
    
    if choice == 'z':
        
        sum_areadist = c.Ask_semi*4*c.h/(2*3*np.pi) + c.Ask_incl*2*(c.h/2-c.Ca)/2
        sum_area = c.Ask_semi + c.Ask_incl*2 + c.Asp
        for i in range(len(c.Z_locations)):
            sum_areadist += c.Z_locations[i]*c.Boom_area[i]
            sum_area += c.Boom_area[i]
        z_tilde = sum_areadist/sum_area
        return z_tilde

z_t =Centroid('z') # Centroid location [m]
        
# z_t = -0.19406263838748938+c.h/2

#%%

# Computes moment of inertia in the yy plane
def MoI_yy(z_tilde):
    

    Iyy_skin = c.lsk**3 * c.tsk * (np.cos(c.beta_rad))**2 / 12 + c.Ask_incl * ( (c.h/2-c.Ca) / 2 - z_tilde)**2
    Iyy_spar = c.Asp * z_tilde**2
#    Iyy_circle = np.pi * r**3 * t_skin / 2 + A_circle * (4 * r / (3 * np.pi) - z_tilde)**2
    Iyy_circle = np.pi * (c.h/2)**3 * c.tsk / 2 + c.Ask_semi * (4 * c.h / (2*3 * np.pi) - z_tilde)**2
    I_yy = 2*Iyy_skin + Iyy_spar + Iyy_circle
    for i in range(len(c.Z_locations)):
        I_yy+= c.Boom_area[i]*(c.Z_locations[i]-z_tilde)**2
    return I_yy

# Computes moment of inertia in the zz plane
def MoI_zz():
    
#    Izz_skin = l_skin**3 * t_skin * (np.sin(alpha))**2 / 12 + A_skin * (l_skin * np.sin(alpha) / 2)**2
    Izz_skin = (2*c.lsk)**3 * c.tsk * (np.sin(c.beta_rad))**2 / 12
    Izz_spar = (c.h)**3 * c.tsp / 12
    Izz_circle = np.pi * (c.h/2)**3 * c.tsk / 2 
    I_zz = Izz_skin + Izz_spar + Izz_circle
    for i in range(len(c.Y_locations)):
        I_zz+= c.Boom_area[i]*c.Y_locations[i]**2
    return I_zz

Iyy = MoI_yy(z_t) # [m^4]
# Iyy = 4.363276766019503e-05
Izz = MoI_zz() # [m^4]

#%%


x,z = aero.meshCreator()
load,pointapp = aero.zload_zpoint(z)

torque = np.zeros((c.Nx,1))
for i in range(c.Nx):
    torque[i] = load[i]*pointapp[i]

#%%

def matrix(xvec):
    A = np.zeros((c.Nx,c.Nx))
    A[:,0] = 1
    xvec_new = (2*np.pi)/(xvec[-1] - xvec[0]) * xvec
    for i in range(c.Nx):
        for j in range(1,c.Nx):
            if 1<=j<=20:
                A[i,j] = np.cos(j*xvec_new[i])
            else:
                A[i,j] = np.sin(j*xvec_new[i])
                
    return np.matrix(A)

L = matrix(x[0,:])

#coeff_torque = np.linalg.inv(L).dot(torque)
#load = np.array(load).reshape((Nx,1))
#coeff_load = np.linalg.inv(L).dot(load)
#
#load_test = L.dot(coeff_load)
#torque_test = L.dot(coeff_torque)
#
#plt.figure(1)
#plt.plot(x[0,:],load_test,x[0,:],load)
#
#r = np.linalg.norm(load-load_test)/np.linalg.norm(load) * 100 # FAA checken
#
#plt.figure(2)
#plt.plot(x[0,:],torque_test,x[0,:],torque_test)





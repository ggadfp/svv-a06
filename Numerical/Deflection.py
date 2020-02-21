# -*- coding: utf-8 -*-
"""
Created on Thursday 20 2020

@author: Gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import Simulation as sim
from mpl_toolkits.mplot3d import Axes3D
#%%

def int_spline(s_coeff,degree,xloc,xvec):
    
    #adding the beginning and end points
    xvec = np.vstack(([0],xvec.reshape((len(xvec),1)),[sim.la]))
    s_coeff = np.vstack(([s_coeff[0,0],0,0,0],s_coeff,[s_coeff[-1,0],0,0,0]))
    
    #finding the index of the last node before the xloc
    diff = np.array(xvec - xloc)
    if xloc == 0:
        idx_fullspl = 0
    elif xloc == sim.la:
        idx_fullspl = len(s_coeff)-1
    else:
        idx_fullspl = int(np.where(diff < 0)[0][-1])    
  
    full_sum = 0
    
    #this one is the value of the function at that point
    if degree == 0:
        full_sum = (s_coeff[idx_fullspl,0] +
                    s_coeff[idx_fullspl,1]*(xloc-xvec[idx_fullspl]) +
                    s_coeff[idx_fullspl,2]*(xloc-xvec[idx_fullspl])**2 +
                    s_coeff[idx_fullspl,3]*(xloc-xvec[idx_fullspl])**3)
            
    #this one is to get the value of the integral
    for k in range(idx_fullspl+1):      
        if k < idx_fullspl:                
            if degree == 1:
                full_sum += (s_coeff[k,0]*(xvec[k+1]-xvec[k]) +
                             s_coeff[k,1]/2*(xvec[k+1]-xvec[k])**2 +
                             s_coeff[k,2]/3*(xvec[k+1]-xvec[k])**3 +
                             s_coeff[k,3]/4*(xvec[k+1]-xvec[k])**4)
            elif degree == 2:
                full_sum += (s_coeff[k,0]/2*(xvec[k+1]-xvec[k])**2 +
                             s_coeff[k,1]/6*(xvec[k+1]-xvec[k])**3 +
                             s_coeff[k,2]/12*(xvec[k+1]-xvec[k])**4 +
                             s_coeff[k,3]/20*(xvec[k+1]-xvec[k])**5)
            elif degree == 4:
                full_sum += (s_coeff[k,0]/24*(xvec[k+1]-xvec[k])**4 +
                             s_coeff[k,1]/120*(xvec[k+1]-xvec[k])**5 +
                             s_coeff[k,2]/360*(xvec[k+1]-xvec[k])**6 +
                             s_coeff[k,3]/840*(xvec[k+1]-xvec[k])**7)
            # print(xvec[k+1],full_sum)
        elif k == idx_fullspl:
            if degree == 1:
                full_sum += (s_coeff[k,0]*(xloc-xvec[k]) +
                             s_coeff[k,1]/2*(xloc-xvec[k])**2 +
                             s_coeff[k,2]/3*(xloc-xvec[k])**3 +
                             s_coeff[k,3]/4*(xloc-xvec[k])**4)
            elif degree == 2:
                full_sum += (s_coeff[k,0]/2*(xloc-xvec[k])**2 +
                             s_coeff[k,1]/6*(xloc-xvec[k])**3 +
                             s_coeff[k,2]/12*(xloc-xvec[k])**4 +
                             s_coeff[k,3]/20*(xloc-xvec[k])**5)
            elif degree == 4:
                full_sum += (s_coeff[k,0]/24*(xloc-xvec[k])**4 +
                             s_coeff[k,1]/120*(xloc-xvec[k])**5 +
                             s_coeff[k,2]/360*(xloc-xvec[k])**6 +
                             s_coeff[k,3]/840*(xloc-xvec[k])**7)
             
    return float(full_sum)

int_torque = int_spline(sim.s_coeff_torque,0,0,sim.x[0])
torque_dist = []
print(type(torque_dist))
for i in range(2001):
    torque_dist.append(int_spline(sim.s_coeff_torque,0,i*sim.la/2000,sim.x[0]))
plt.plot(np.linspace(0,sim.la,num=2001),torque_dist)
plt.show()

#%%
#We now wnat to create the 12x12 matrix and the vector of BCs e to solve for the unknown forces thanks to our deflection insights
mat_defl = np.zeros([12,12])
e = np.zeros([12,1])


# From BC 1 (check Simulation Plan)
# e = 

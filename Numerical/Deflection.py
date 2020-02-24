# -*- coding: utf-8 -*-
"""
Created on Thursday 20 2020

@author: Gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import Simulation as sim
import time
#%%
start_time = time.time()
#[H1y,H2y,H3y,H1z,H2z,H3z,A1,C1,C2,C3,C4,C5]
#We now wnat to create the 12x12 matrix and the vector of BCs e to solve for the unknown forces thanks to our deflection insights
mat_defl = np.zeros([12,12])
e = np.zeros([12,1])

# From BC 1 (check Simulation Plan)
e[0] = (sim.doubleinteg_spline(sim.s_coeff_load, sim.la, sim.x[0], 1000)+
        sim.P*np.sin(sim.theta_rad)*(sim.la-sim.x2-sim.xa/2))
mat_defl[0] = [
    (sim.la-sim.x1),
    (sim.la-sim.x2),
    (sim.la-sim.x3),
    0,
    0,
    0,
    np.sin(sim.theta_rad)*(sim.la-sim.x2+sim.xa/2),
    0,0,0,0,0
    ]
print("BC1 took", time.time() - start_time, "to run")
# From BC 2 (check Simulation Plan)
e[1] = sim.P*np.cos(sim.theta_rad)*(sim.la-sim.x2-sim.xa/2)
mat_defl[1] = [
    0,
    0,
    0,
    (sim.la-sim.x1),
    (sim.la-sim.x2),
    (sim.la-sim.x3),
    np.cos(sim.theta_rad)*(sim.la-sim.x2+sim.xa/2),
    0,0,0,0,0
    ] 
print("BC1+2 took", time.time() - start_time, "to run")
# From BC 3 (check Simulation Plan)
e[2] = (sim.integ_spline(sim.s_coeff_torque,sim.la,sim.x[0])+
        sim.P*sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad)))
mat_defl[2] = [
    0,0,0,0,0,0,
    sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad)),
    0,0,0,0,0
    ]
print("BC1+2+3 took", time.time() - start_time, "to run")
# From BC 4 (check Simulation Plan)
e[3] = (sim.integ_spline(sim.s_coeff_load,sim.la,sim.x[0])+
        sim.P*np.sin(sim.theta_rad))
mat_defl[3] = [
    1,
    1,
    1,
    0,
    0,
    0,
    np.sin(sim.theta_rad),
    0,0,0,0,0
    ]  
print("BC1+2+3+4 took", time.time() - start_time, "to run")
# From BC 5 (check Simulation Plan)
e[4] = sim.P*np.cos(sim.theta_rad)
mat_defl[4] = [
    0,
    0,
    0,
    1,
    1,
    1,
    np.cos(sim.theta_rad),
    0,0,0,0,0
    ]
print("BC1+2+3+4+5 took", time.time() - start_time, "to run")
# From BC 6 (check Simulation Plan)
e[5] = sim.quadrupleinteg_spline(sim.s_coeff_load,sim.x2,sim.x[0],1000)
mat_defl[5] = [
    (sim.x2-sim.x1)**3/6,
    0,
    0,
    0,
    0,
    0,
    np.sin(sim.theta_rad)/6*(sim.xa/2)**3,
    sim.x2*sim.E*sim.Izz,
    sim.E*sim.Izz,
    0,0,0
    ]
print("BC1+2+3+4+5+6 took", time.time() - start_time, "to run")
# From BC 7 (check Simulation Plan)
e[6] = 0
mat_defl[6] = [
    0,
    0,
    0,
    (sim.x2-sim.x1)**3/6,
    0,
    0,
    np.cos(sim.theta_rad)/6*(sim.xa/2)**3,
    0,
    0,
    sim.x2*sim.E*sim.Iyy,
    sim.E*sim.Iyy,
    0
    ]
print("BC1+2+3+4+5+6+7 took", time.time() - start_time, "to run")
# From BC 8 (check Simulation Plan)
e[7] = (sim.quadrupleinteg_spline(sim.s_coeff_load,sim.x1,sim.x[0],1000)/(sim.E*sim.Izz)+
        sim.d1*np.cos(sim.theta_rad))
mat_defl[7] = [
    0,0,0,0,0,0,
    0,
    sim.x1,
    1,
    0,0,0
    ]
print("BC1+2+3+4+5+6+7+8 took", time.time() - start_time, "to run")
# From BC 9 (check Simulation Plan)
e[8] = sim.d1*np.sin(sim.theta_rad)
mat_defl[8] = [
    0,0,0,0,0,0,
    0,
    0,
    0,
    sim.x1,
    1,
    0
    ]
print("BC1+2+3+4+5+6+7+8+9 took", time.time() - start_time, "to run")
# From BC 10 (check Simulation Plan)
e[9] = (sim.d3*np.cos(sim.theta_rad)+
        sim.quadrupleinteg_spline(sim.s_coeff_load,sim.x3,sim.x[0],1000)+
        sim.P*np.sin(sim.theta_rad)*(sim.x3-sim.x2-sim.xa/2)**3/6)
mat_defl[9] = [
    (sim.x3-sim.x1)**3/6,
    (sim.x3-sim.x2)**3/6,
    0,
    0,
    0,
    0,
    np.sin(sim.theta_rad)*(sim.x3-sim.x2+sim.xa/2)**3/6,
    sim.x3*sim.E*sim.Izz,
    sim.E*sim.Izz,
    0,0,0
    ]
print("BC1+2+3+4+5+6+7+8+9+10 took", time.time() - start_time, "to run")
# From BC 11 (check Simulation Plan)
e[10] = (-sim.d3*np.sin(sim.theta_rad)*(sim.E*sim.Iyy)+
         sim.P*np.cos(sim.theta_rad)*(sim.x3-sim.x2-sim.xa/2)**3/6)
mat_defl[10] = [
    0,0,0,
    (sim.x3-sim.x1)**3/6,
    (sim.x3-sim.x2)**3/6,
    0,
    np.cos(sim.theta_rad)*(sim.x3-sim.x2+sim.xa/2)**3/6,
    0,
    0,
    sim.x3*sim.E*sim.Iyy,
    sim.E*sim.Iyy,
    0
    ]
print("BC1+2+3+4+5+6+7+8+9+10+11 took", time.time() - start_time, "to run")
# From BC 12 (check Simulation Plan)
e[11] = sim.doubleinteg_spline(sim.s_coeff_torque,sim.x2-sim.xa/2,sim.x[0],1000)
mat_defl[11] = [
    0,0,0,0,0,0,0,0,0,0,0,
    sim.G*(sim.Iyy+sim.Izz)
    ]
print("BC1+2+3+4+5+6+7+8+9+10+11+12 took", time.time() - start_time, "to run")

#SOLVING FOR ALL THE UNKNOWN FORCES
unknowns = np.linalg.inv(mat_defl).dot(e)

print("Deflection took", time.time() - start_time, "to run")

#%%
def macaulay(xloc,start,n):
    d = xloc-start
    if d <= 0:
        return 0
    if d >0:
        return d**n

#%%
# MOMENT IN Z SPANWISE
    
def M_z(xloc):
    
    M_z = (sim.doubleinteg_spline(sim.s_coeff_load, xloc, sim.x[0], 100)
     +sim.P*np.sin(sim.theta_rad)*(sim.la-sim.x2-sim.xa/2)*macaulay(xloc,(sim.x2+sim.xa/2),1)
     - unknowns[0]*macaulay(xloc,sim.x1,1)
     - unknowns[1]*macaulay(xloc,sim.x2,1)
     - unknowns[2]*macaulay(xloc,sim.x3,1)
     - unknowns[6]*np.sin(sim.theta_rad)*macaulay(xloc,(sim.x2-sim.xa/2),1))
    # print(sim.doubleinteg_spline(sim.s_coeff_load, xloc, sim.x[0], 1000))
    # print(sim.P*np.sin(sim.theta_rad)*(sim.la-sim.x2-sim.xa/2)*macaulay(xloc,(sim.x2+sim.xa/2),1))
    return M_z

# print(M_z(sim.la))

plt.plot(np.linspace(0,sim.la,101),[M_z(i) for i in np.linspace(0,sim.la,101)])

#%%
# SHEAR FORCE IN Y SPANWISE

def S_y(xloc):
    S_y = (sim.integ_spline(sim.s_coeff_load,xloc,sim.x[0])
           + sim.P*np.sin(sim.theta_rad)*macaulay(xloc,(sim.x2+sim.xa/2),0)
           - unknowns[0]*macaulay(xloc,sim.x1,0)
           - unknowns[1]*macaulay(xloc,sim.x2,0)
           - unknowns[2]*macaulay(xloc,sim.x3,0)
           - unknowns[6]*np.sin(sim.theta_rad)*macaulay(xloc,(sim.x2-sim.xa/2),0)
           )
    
    return S_y

plt.plot(np.linspace(0,sim.la,101),[S_y(i) for i in np.linspace(0,sim.la,101)])

#%%
# SLOPE IN Y SPANWISE

def dvy_dx(xloc):
    
    dvy_dx = -(sim.tripleinteg_spline(sim.s_coeff_load,xloc,sim.x[0],100)
              + sim.P*np.sin(sim.theta_rad)/2*macaulay(xloc,(sim.x2+sim.xa/2),2)
              - unknowns[0]/2*macaulay(xloc,sim.x1,2)
              - unknowns[1]/2*macaulay(xloc,sim.x2,2)
              - unknowns[1]/2*macaulay(xloc,sim.x3,2)
              - unknowns[6]*np.sin(sim.theta_rad)/2*macaulay(xloc,(sim.x2-sim.xa/2),2) 
              )/(sim.E*sim.Izz) +unknowns[7]
    
    return dvy_dx

plt.plot(np.linspace(0,sim.la,101),[dvy_dx(i) for i in np.linspace(0,sim.la,101)])

#%%
# DEFLECTION IN Y SPANWISE

def deflection_y(xloc):
    
    def_y = -(sim.quadrupleinteg_spline(sim.s_coeff_load,xloc,sim.x[0],100)
              + sim.P*np.sin(sim.theta_rad)/6*macaulay(xloc,(sim.x2+sim.xa/2),3)
              - unknowns[0]/6*macaulay(xloc,sim.x1,3)
              - unknowns[1]/6*macaulay(xloc,sim.x2,3)
              - unknowns[1]/6*macaulay(xloc,sim.x3,3)
              - unknowns[6]*np.sin(sim.theta_rad)/6*macaulay(xloc,(sim.x2-sim.xa/2),3) 
              )/(sim.E*sim.Izz) +unknowns[7]*xloc+unknowns[8]
    return def_y

plt.plot(np.linspace(0,sim.la,101),[deflection_y(i) for i in np.linspace(0,sim.la,101)])

#%%




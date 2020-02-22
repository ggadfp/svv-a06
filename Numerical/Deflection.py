# -*- coding: utf-8 -*-
"""
Created on Thursday 20 2020

@author: Gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import Simulation as sim
#%%
#[H1y,H2y,H3y,H1z,H2z,H3z,A1,C1,C2,C3,C4,C5]
#We now wnat to create the 12x12 matrix and the vector of BCs e to solve for the unknown forces thanks to our deflection insights
mat_defl = np.zeros([12,12])
e = np.zeros([12,1])

# From BC 1 (check Simulation Plan)
e[0] = (int_spline(sim.s_coeff_load,2,sim.la,sim.x[0])+
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
# From BC 3 (check Simulation Plan)
e[2] = (int_spline(sim.s_coeff_torque,1,sim.la,sim.x[0])+
        sim.P*sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad)))
mat_defl[2] = [
    0,0,0,0,0,0,
    sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad)),
    0,0,0,0,0
    ]
# From BC 4 (check Simulation Plan)
e[3] = (int_spline(sim.s_coeff_load,1,sim.la,sim.x[0])+
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
# From BC 6 (check Simulation Plan)
e[5] = int_spline(sim.s_coeff_load,4,sim.x2,sim.x[0])
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
# From BC 8 (check Simulation Plan)
e[7] = (int_spline(sim.s_coeff_load,4,sim.x1,sim.x[0])/(sim.E*sim.Izz)+
        sim.d1*np.cos(sim.theta_rad))
mat_defl[7] = [
    0,0,0,0,0,0,
    0,
    sim.x1,
    1,
    0,0,0
    ]
# From BC 9 (check Simulation Plan)
e[8] = -sim.d1*np.sin(sim.theta_rad)
mat_defl[8] = [
    0,0,0,0,0,0,
    0,
    0,
    0,
    sim.x1,
    1,
    0
    ]
# From BC 10 (check Simulation Plan)
e[9] = (sim.d3*np.cos(sim.theta_rad)+
        int_spline(sim.s_coeff_load,4,sim.x3,sim.x[0])+
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
# From BC 12 (check Simulation Plan)
e[11] = int_spline(sim.s_coeff_torque,2,sim.x2-sim.xa/2,sim.x[0])
mat_defl[11] = [
    0,0,0,0,0,0,0,0,0,0,0,
    sim.G*(sim.Iyy+sim.Izz)
    ]

#SOLVING FOR ALL THE UNKNOWN FORCES
unknowns = np.linalg.inv(mat_defl).dot(e)
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:48:37 2020

@author: Ihab
"""


import numpy as np
import readAeroload as aero
import Simulation as sim
import matplotlib.pyplot as plt
import Deflection as d
import Parameters as c
import ShearDist as sh
import J as tor
#%%

def bending_stress(xloc,yloc,zloc):
    return d.M_z(xloc)*yloc/sim.Izz + d.M_y(xloc)*zloc/sim.Iyy



def shear_stress(xloc):
    def shear_z():
    
       q_Sz = d.S_z(xloc)*sh.base_shear_z(sim.Iyy)
       return q_Sz
    
    def shear_y():
    
        q_Sy = d.S_y(xloc)*sh.base_shear_y(sim.Izz)
        return q_Sy

    q01,q02,J = tor.torsion_shear(xloc)
    tau_yz_circle = (shear_z()[0] + shear_y()[0]+q01)/c.tsk
    tau_yz_circle_rev = (shear_z()[1] + shear_y()[1]+q01)/c.tsk
    tau_yz_incl_1 = (shear_z()[2] + shear_y()[2]+q02)/c.tsk
    tau_yz_incl_2 = (shear_z()[3] + shear_y()[3]+q02)/c.tsk
    tau_yz_spar = (shear_z()[4] + shear_y()[4]+q02-q01)/c.tsp
    tau_yz_spar_rev = (shear_z()[5] + shear_y()[5]+q02-q01)/c.tsp
    
    
    return tau_yz_circle,tau_yz_circle_rev,tau_yz_incl_1,tau_yz_incl_2,tau_yz_spar,tau_yz_spar_rev


test = shear_stress(0.5)
    
plt.plot(test[0])
    
    
    
    

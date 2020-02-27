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
#%%

def bending_stress(xloc,yloc,zloc):
    return d.M_z(xloc)*yloc/sim.Izz + d.M_y(xloc)*zloc/sim.Iyy

test = bending_stress(0.25,0,-sim.z_t-(c.Ca-c.h/2))
print(test*10**(-6))
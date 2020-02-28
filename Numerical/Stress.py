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
import TorsionShear as tor
from mpl_toolkits.mplot3d import Axes3D
#%%

z_1 = np.linspace(-sim.z_t,-sim.z_t+c.h/2,5000)
z_2 = np.linspace(-sim.z_t,-c.Ca+c.h/2-sim.z_t,5000)
z_3 = np.linspace(-c.Ca+c.h/2-sim.z_t,-sim.z_t,5000)
z_4 = -sim.z_t*np.ones(5000)

def circle(z1):
    return np.sqrt(np.abs((c.h/2)**2 - (z1+sim.z_t)**2))
def circle_rev(z1):
    return -np.sqrt(np.abs((c.h/2)**2 - (z1+sim.z_t)**2))
def len_incl_1(z1):
    return z1*np.tan(c.beta_rad) + sim.z_t*np.tan(c.beta_rad)+c.h/2
def len_incl_2(z1):
    return -(z1+sim.z_t)*np.tan(c.beta_rad) - c.h/2
def spar_h(z1):
    return np.linspace(0,-c.h/2,5000)
def spar_h_rev(z1):
    return np.linspace(0,c.h/2,5000)

def bending_stress(xloc):
    
    
    bending_circle = d.M_z(xloc)*circle(z_1)/sim.Izz + d.M_y(xloc)*z_1/sim.Iyy
    bending_circle_rev = d.M_z(xloc)*circle_rev(z_1)/sim.Izz + d.M_y(xloc)*z_1/sim.Iyy
    bending_incl_1 = d.M_z(xloc)*len_incl_1(z_2)/sim.Izz + d.M_y(xloc)*z_2/sim.Iyy
    bending_incl_2 = d.M_z(xloc)*len_incl_2(z_3)/sim.Izz + d.M_y(xloc)*z_3/sim.Iyy
    bending_spar = d.M_z(xloc)*spar_h(z_4)/sim.Izz + d.M_y(xloc)*z_4/sim.Iyy
    bending_spar_rev = d.M_z(xloc)*spar_h_rev(z_4)/sim.Izz + d.M_y(xloc)*z_4/sim.Iyy
    
    return bending_circle,bending_circle_rev,bending_incl_1,bending_incl_2,bending_spar,bending_spar_rev


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



def Von_Mises(xloc):
    
    vm_circle = np.sqrt(bending_stress(xloc)[0]**2 + 3*shear_stress(xloc)[0]**2)
    vm_circle_rev = np.sqrt(bending_stress(xloc)[1]**2 + 3*shear_stress(xloc)[1]**2)
    vm_incl_1 = np.sqrt(bending_stress(xloc)[2]**2 + 3*shear_stress(xloc)[2]**2)
    vm_incl_2 = np.sqrt(bending_stress(xloc)[3]**2 + 3*shear_stress(xloc)[3]**2)
    vm_spar = np.sqrt(bending_stress(xloc)[4]**2 + 3*shear_stress(xloc)[4]**2)
    vm_spar_rev = np.sqrt(bending_stress(xloc)[5]**2 + 3*shear_stress(xloc)[5]**2)
    
    return vm_circle,vm_circle_rev,vm_incl_1,vm_incl_2,vm_spar,vm_spar_rev

#%%
    
plt.plot(z_1,circle(z_1))
plt.plot(z_1,circle_rev(z_1))
plt.plot(z_2,len_incl_1(z_2))
plt.plot(z_3,len_incl_2(z_3))
plt.plot(z_4,spar_h(z_4))
plt.plot(z_4,spar_h_rev(z_4))
    

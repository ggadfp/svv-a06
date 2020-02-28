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
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
#%%

z_1 = np.linspace(-sim.z_t,-sim.z_t+c.h/2,5000)
z_2 = np.linspace(-sim.z_t,-c.Ca+c.h/2-sim.z_t,5000)
z_3 = np.linspace(-c.Ca+c.h/2-sim.z_t,-sim.z_t,5000)
z_4 = -sim.z_t*np.ones(5000)
z_5 = np.linspace(-sim.z_t+c.h/2,-sim.z_t,5000)
def circle(z1):
    return -np.sqrt(np.abs((c.h/2)**2 - (z1+sim.z_t)**2))
def circle_rev(z1):
    return np.sqrt(np.abs((c.h/2)**2 - (z1+sim.z_t)**2))
def len_incl_1(z1):
    return z1*np.tan(c.beta_rad) + sim.z_t*np.tan(c.beta_rad) +c.h/2
def len_incl_2(z1):
    return -(z1+sim.z_t)*np.tan(c.beta_rad) -c.h/2
def spar_h(z1):
    return np.linspace(0,-c.h/2,5000)
def spar_h_rev(z1):
    return np.linspace(0,c.h/2,5000)

y_1 = circle(z_5)
y_2 = circle_rev(z_1)
y_3 = len_incl_1(z_2)
y_4 = len_incl_2(z_3)
y_5 = spar_h(z_4)
y_6 = spar_h_rev(z_4)



# plt.plot(z_2,len_incl_1(z_2))
# plt.plot(z_3,len_incl_2(z_3))
def bending_stress(xloc):
    
    
    bending_circle = d.M_z(xloc)*y_1/sim.Izz + d.M_y(xloc)*z_5/sim.Iyy
    bending_circle_rev = d.M_z(xloc)*y_2/sim.Izz + d.M_y(xloc)*z_1/sim.Iyy
    bending_incl_1 = d.M_z(xloc)*y_3/sim.Izz + d.M_y(xloc)*z_2/sim.Iyy
    bending_incl_2 = d.M_z(xloc)*y_4/sim.Izz + d.M_y(xloc)*z_3/sim.Iyy
    bending_spar = d.M_z(xloc)*y_5/sim.Izz + d.M_y(xloc)*z_4/sim.Iyy
    bending_spar_rev = d.M_z(xloc)*y_6/sim.Izz + d.M_y(xloc)*z_4/sim.Iyy
    
    return bending_circle,bending_circle_rev,bending_incl_1,bending_incl_2,bending_spar,bending_spar_rev

def shear_flow(xloc):
    def shear_z():
    
       q_Sz = d.S_z(xloc)*sh.base_shear_z(sim.Iyy)
       return q_Sz
    
    def shear_y():
    
        q_Sy = d.S_y(xloc)*sh.base_shear_y(sim.Izz)
        return q_Sy

    q01,q02,J = tor.torsion_shear(xloc)
    q_yz_circle = ((shear_z()[0] + shear_y()[0])+q01)
    q_yz_circle_rev = ((shear_z()[1] + shear_y()[1])+q01)
    q_yz_incl_1 = ((shear_z()[2] + shear_y()[2])+q02)
    q_yz_incl_2 = ((shear_z()[3] + shear_y()[3])+q02)
    q_yz_spar = ((shear_z()[4] + shear_y()[4])+q02-q01)
    q_yz_spar_rev = ((shear_z()[5] + shear_y()[5])+q02-q01)
    
    
    return q_yz_circle,q_yz_circle_rev,q_yz_incl_1,q_yz_incl_2,q_yz_spar,q_yz_spar_rev

def shear_stress(xloc):
    def shear_z():
    
       q_Sz = d.S_z(xloc)*sh.base_shear_z(sim.Iyy)
       return q_Sz
    
    def shear_y():
    
        q_Sy = d.S_y(xloc)*sh.base_shear_y(sim.Izz)
        return q_Sy

    q01,q02,J = tor.torsion_shear(xloc)
    tau_yz_circle = ((shear_z()[0] + shear_y()[0] +q01) /c.tsk)
    tau_yz_circle_rev = ((shear_z()[1] + shear_y()[1] +q01)/c.tsk)
    tau_yz_incl_1 = ((shear_z()[2] + shear_y()[2] +q02)/c.tsk)
    tau_yz_incl_2 = ((shear_z()[3] + shear_y()[3] +q02)/c.tsk)
    tau_yz_spar = ((shear_z()[4] + shear_y()[4] +q02-q01)/c.tsp)
    tau_yz_spar_rev = ((shear_z()[5] + shear_y()[5] +q02-q01)/c.tsp)
    
    
    return tau_yz_circle,tau_yz_circle_rev,tau_yz_incl_1,tau_yz_incl_2,tau_yz_spar,tau_yz_spar_rev



def Von_Mises(xloc):
    
    r = bending_stress(xloc)
    s = shear_stress(xloc)
    vm_circle = np.sqrt(r[0]**2 + 3*s[0]**2)
    vm_circle_rev = np.sqrt(r[1]**2 + 3*s[1]**2)
    vm_incl_1 = np.sqrt(r[2]**2 + 3*s[2]**2)
    vm_incl_2 = np.sqrt(r[3]**2 + 3*s[3]**2)
    vm_spar = np.sqrt(r[4]**2 + 3*s[4]**2)
    vm_spar_rev = np.sqrt(r[5]**2 + 3*s[5]**2)
    
    return vm_circle,vm_circle_rev,vm_incl_1,vm_incl_2,vm_spar,vm_spar_rev

#%%
    
def forLooper(var1,var2,var3,var4,var5,var6):
    var = []
    for j in range(5000):
        var.append(var1[j])
    for j in range(5000):
        var.append(var2[j])
    for j in range(5000):
        var.append(var3[j])
    for j in range(5000):
        var.append(var4[j])
    for j in range(5000):
        var.append(var5[j])
    for j in range(5000):
        var.append(var6[j])
        
    return var

shearstress = shear_stress(0.5)
shearflow = shear_flow(0.554)
bendingstress = bending_stress(0.554)
vmstress = Von_Mises(0.554)
z_tot = forLooper(z_5,z_1,z_2,z_3,z_4,z_4)
y_tot = forLooper(y_1,y_2,y_3,y_4,y_5,y_6)
bending_tot = forLooper(bendingstress[0],bendingstress[1],bendingstress[2],bendingstress[3],bendingstress[4],bendingstress[5])
shearf_tot = forLooper(shearflow[0],shearflow[1],shearflow[2],shearflow[3],shearflow[4],shearflow[5])
vm_tot = forLooper(vmstress[0],vmstress[1],vmstress[2],vmstress[3],vmstress[4],vmstress[5])

# plt.figure(0)
# plt.title("Bending stress at span-location x = 0.554")
# plt.scatter(z_tot,y_tot,s=0.1,c=bending_tot,cmap=plt.cm.jet)
# plt.plot(0,0,'o',c='r',markersize=6,markeredgecolor = 'k')
# plt.xlabel('Z [m]')
# plt.ylabel('Y [m]')
# plt.xlim(-sim.z_t,-c.Ca+c.h/2-sim.z_t)
# plt.colorbar(orientation = 'horizontal', fraction = 0.05)
# plt.axis('equal')
# plt.savefig("Bendingstress_0.554.png")

# plt.figure(1)
# plt.title("Shear flow distribution at span-location x = 0.554")
# plt.scatter(z_tot,y_tot,s=0.1,c=shearf_tot,cmap=plt.cm.jet)
# plt.plot(0,0,'o',c='r',markersize=6,markeredgecolor = 'k')
# plt.xlabel('Z [m]')
# plt.ylabel('Y [m]')
# plt.xlim(-sim.z_t,-c.Ca+c.h/2-sim.z_t)
# plt.colorbar(orientation = 'horizontal', fraction = 0.05)
# plt.axis('equal')
# plt.savefig("Shearflow_0.554.png")

# plt.figure(2)
# plt.title("Von Mises stress distribution at span-location x = 0.554")
# plt.scatter(z_tot,y_tot,s=0.1,c=vm_tot,cmap=plt.cm.jet)
# plt.plot(0,0,'o',c='r',markersize=6,markeredgecolor = 'k')
# plt.xlabel('Z [m]')
# plt.ylabel('Y [m]')
# plt.xlim(-sim.z_t,-c.Ca+c.h/2-sim.z_t)
# plt.colorbar(orientation = 'horizontal', fraction = 0.05)
# plt.axis('equal')
# plt.savefig("Vonmises_0.554")




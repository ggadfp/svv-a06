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

z_t = Centroid('z') # Centroid location [m]

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


#%%

# Calculates the amount of point for each part of the structure
def node_maker():
    nodes_skin = np.linspace(0,c.lsk,10001)
    nodes_circle = np.linspace(0,np.pi/2,10001)
    nodes_circle_rev = np.linspace(-np.pi/2,0,10001)
    nodes_spar = np.linspace(0,c.h/2,10001)
    nodes_spar_rev = np.linspace(0,-c.h/2,10001)
    return nodes_skin,nodes_circle,nodes_circle_rev,nodes_spar,nodes_spar_rev


nodes_skin,nodes_circle,nodes_circle_rev,nodes_spar,nodes_spar_rev = node_maker()

# Function for the base shear flow
def base_shear(I_zz):
    
    # Function for the semi-circular part integrating the base shear flow from 0 to pi/2
    def circle():
        q_c = []
        
        I_c = -c.tsk*c.h*c.h/(2*2*I_zz)
        qval = 0
        count = 0
        for i in range(0,len(nodes_circle)-2,2):
            qval += I_c*(nodes_circle[i+2] - nodes_circle[i])/6 * (-np.sin(nodes_circle[i+2]) - 4*np.sin(nodes_circle[i+1]) - np.sin(nodes_circle[i]))
            if nodes_circle[i] >= c.angle_rad and count == 0:
                qval += -c.Boom_area[1]*c.Y_locations[1]/I_zz
                count += 1
            
            (q_c).append(qval)
        return np.array(q_c)
    
    # Function for the semi-circular part integrating the base shear flow from -pi/2 to 0
    def circle_rev():
        q_c = []
        I_c = -c.tsk*c.h*c.h/(2*2*I_zz) 
        qval = inclined_2()[-1] - spar_rev()[-1]
        count = 0
        for i in range(0,len(nodes_circle_rev)-2,2):
            qval += I_c*(nodes_circle_rev[i+2] - nodes_circle_rev[i])/6 * (-np.sin(nodes_circle_rev[i+2]) - 4*np.sin(nodes_circle_rev[i+1]) - np.sin(nodes_circle_rev[i]))
            if nodes_circle_rev[i] >= -c.angle_rad and count == 0:
                qval += -c.Boom_area[12]*c.Y_locations[12]/I_zz
                count += 1
            (q_c).append(qval)
        return np.array(q_c)
    
    # Function for the inclined section integrating the base shear flow from 0 to length of this section lsk
    def inclined_1():
        q_c = []
        I_c = -c.tsk/(I_zz)
        qval = circle()[-1] + spar()[-1]
        count = 0
        for i in range(0,len(nodes_skin)-2,2):
            qval += I_c*((nodes_skin[i+2] - nodes_skin[i])/6 * ( (-c.h/2 + nodes_skin[i+2]*np.sin(c.beta_rad)) + 4*(-c.h/2+nodes_skin[i+1]*np.sin(c.beta_rad)) + (-c.h/2 + nodes_skin[i]*np.sin(c.beta_rad))))
            if nodes_skin[i] >= np.abs(c.Z_locations[2]/np.cos(c.beta_rad)) and count == 0:
                qval += -c.Boom_area[2]*c.Y_locations[2]/I_zz
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[3]/np.cos(c.beta_rad)) and count == 1:
                qval += -c.Boom_area[3]*c.Y_locations[3]/I_zz
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[4]/np.cos(c.beta_rad)) and count == 2:
                qval += -c.Boom_area[4]*c.Y_locations[4]/I_zz
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[5]/np.cos(c.beta_rad)) and count == 3:
                qval += -c.Boom_area[5]*c.Y_locations[5]/I_zz
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[6]/np.cos(c.beta_rad)) and count == 4:
                qval += -c.Boom_area[6]*c.Y_locations[6]/I_zz
                count += 1
            (q_c).append(qval)
        return np.array(q_c)

    # Function for the inclined section integrating the base shear flow from 0 to length of this section lsk
    def inclined_2():
        q_c = []
        I_c = -c.tsk*np.sin(c.beta_rad)/(I_zz)
        qval = inclined_1()[-1]
        count = 0
        for i in range(0,len(nodes_skin)-2,2):
            qval += I_c*(nodes_skin[i+2] - nodes_skin[i])/6 * (nodes_skin[i+2] + 4*nodes_skin[i+1] + nodes_skin[i])
            if nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[7]/np.cos(c.beta_rad)) and count == 0:
                qval += -c.Boom_area[7]*c.Y_locations[7]/I_zz
                count+=1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[8]/np.cos(c.beta_rad)) and count == 1:
                qval += -c.Boom_area[8]*c.Y_locations[8]/I_zz
                count += 1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[9]/np.cos(c.beta_rad)) and count == 2:
                qval += -c.Boom_area[9]*c.Y_locations[9]/I_zz
                count += 1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[10]/np.cos(c.beta_rad)) and count == 3:
                qval += -c.Boom_area[10]*c.Y_locations[10]/I_zz
                count += 1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[11]/np.cos(c.beta_rad)) and count == 4:
                qval += -c.Boom_area[11]*c.Y_locations[11]/I_zz
                count += 1
            (q_c).append(qval)
        
        return np.array(q_c)
    # Function for the spar section integrating the base shear flow from 0 to h/2
    def spar():
        q_c = []
        I_c = c.tsp/I_zz
        qval = 0
        for i in range(0,len(nodes_spar)-2,2):
            qval += I_c*(nodes_spar[i+2] - nodes_spar[i])/6 * (nodes_spar[i+2] + 4*nodes_spar[i+1] + nodes_spar[i])
            (q_c).append(qval)
        return np.array(q_c)
    
    # Function for the spar section integrating the base shear flow from 0 to -h/2
    def spar_rev():
        q_c = []
        I_c = c.tsp/I_zz
        qval = 0
        for i in range(0,len(nodes_spar_rev)-2,2):
            qval += I_c*(nodes_spar_rev[i+2] - nodes_spar_rev[i])/6 * (nodes_spar_rev[i+2] + 4*nodes_spar_rev[i+1] + nodes_spar_rev[i])
            (q_c).append(qval)
        return np.array(q_c)
    
    
    return circle(),circle_rev(),inclined_1(),inclined_2(),spar(),spar_rev()


test1 = base_shear(Izz)
plt.plot(test1[0],nodes_circle[0:-1:2])
plt.plot(test1[1],nodes_circle_rev[0:-1:2])
def redundant_shear(I_zz):
    
    qb_circle,qb_circle_rev,qb_incl_1,qb_incl_2,qb_spar,qb_spar_rev = base_shear(I_zz)
    
    # Cell 1
    def circle_red():
        q_red = []
        I_c = c.h/(2*c.G*c.tsk)
        half_circle = nodes_circle[0:-1:2]
        q_int = 0
        for i in range(0,len(half_circle)-2,2):
            q_int = (half_circle[i+2] - half_circle[i])/6 * (qb_circle[i+2] + 4*qb_circle[i+1] + qb_circle[i])
            (q_red).append(q_int)
        return np.sum(q_red)*I_c
    
    def circle_red_rev():
        q_red = []
        I_c = c.h/(2*c.G*c.tsk)
        half_circle_rev = nodes_circle_rev[0:-1:2]
        q_int = 0
        for i in range(0,len(half_circle_rev)-2,2):
            q_int = (half_circle_rev[i+2] - half_circle_rev[i])/6 * (qb_circle_rev[i+2] + 4*qb_circle_rev[i+1] + qb_circle_rev[i])
            (q_red).append(q_int)
        return np.sum(q_red)*I_c
    
    def spar_red():
        q_red = []
        I_c = 1/(c.G*c.tsp)
        half_spar = nodes_spar[0:-1:2]
        q_int = 0
        for i in range(0,len(half_spar)-2,2):
            q_int = (half_spar[i+2] - half_spar[i])/6 * (qb_spar[i+2] + 4*qb_spar[i+1] + qb_spar[i])
            (q_red).append(q_int)
        return np.sum(q_red)*I_c
    
    def spar_red_rev():
        q_red = []
        I_c = 1/(c.G*c.tsp)
        half_spar_rev = nodes_spar_rev[0:-1:2]
        q_int = 0
        for i in range(0,len(half_spar_rev)-2,2):
            q_int = (half_spar_rev[i+2] - half_spar_rev[i])/6 * (qb_spar_rev[i+2] + 4*qb_spar_rev[i+1] + qb_spar_rev[i])
            (q_red).append(q_int)
        return np.sum(q_red)*I_c
    
    def incl_red():
        q_red_1 = []
        q_red_2 = []
        I_c = 1/(c.G*c.tsk)
        half_incl = nodes_skin[0:-1:2]
        q_int_1 = 0
        q_int_2 = 0
        for i in range(0,len(half_incl)-2,2):
            q_int_1 = (half_incl[i+2]-half_incl[i])/6 *(qb_incl_1[i+2] + 4*qb_incl_1[i+1] + qb_incl_1[i])
            q_int_2 = (half_incl[i+2]-half_incl[i])/6 *(qb_incl_2[i+2] + 4*qb_incl_2[i+1] + qb_incl_2[i])
            (q_red_1).append(q_int_1)
            (q_red_2).append(q_int_2)
            
        return np.sum(q_red_1)*I_c,np.sum(q_red_2)*I_c
    
    def denom_cell_1():
        return np.pi*c.h/(2*c.tsk) + c.h/c.tsp
    def denom_cell_2():
        return 2*c.lsk/c.tsk + c.h/c.tsp 
    
    D1 = denom_cell_1()
    D2 = denom_cell_2()
    Cr = circle_red()
    Cr_rev = circle_red_rev()
    Sr = spar_red()
    Sr_rev = spar_red_rev()
    Ir_1,Ir_2 = incl_red()
    A = np.matrix([[D1/c.G,-c.h/(c.G*c.tsp)],[-c.h/(c.G*c.tsp),D2/c.G]])
    b = np.array([[-Cr - Cr_rev - Sr - Sr_rev],[-Ir_1 - Ir_2 - Sr - Sr_rev]])
    v = np.linalg.solve(A,b)
        
    
    
    return v

    

test = redundant_shear(Izz) 
    
def Moment_Eq(I_zz):
    qb_circle,qb_circle_rev,qb_incl_1,qb_incl_2,qb_spar,qb_spar_rev = base_shear(I_zz)
    redundant_shear_flow = redundant_shear(I_zz)
    q_red_1 = redundant_shear_flow[0]
    q_red_2 = redundant_shear_flow[1]
    
    qt_circle = qb_circle + q_red_1
    qt_circle_rev = qb_circle_rev + q_red_1
    qt_incl_1 = qb_incl_1 + q_red_2
    qt_incl_2 = qb_incl_2 + q_red_2
    qt_spar = qb_spar + q_red_2 - q_red_1
    
#    plt.plot(qt_circle,nodes_circle[0:-1:2],'g')
#    plt.plot(qt_circle_rev,nodes_circle_rev[0:-1:2],'r')
    plt.figure(1)
    plt.plot(nodes_skin[0:-1:2],qt_incl_1)
    plt.plot(nodes_skin[0:-1:2],qt_incl_2)
    
    half_circle = nodes_circle[0:-1:2]
    half_circle_rev = nodes_circle_rev[0:-1:2]
    half_incl = nodes_skin[0:-1:2]
    
    # def force_circle():
    #     q_tot = []
    #     q_int = 0
    #     for i in range(0,len(half_circle)-2,2):
    #         q_int = (half_circle[i+2] - half_circle[i])/6 * (qt_circle[i+2] + 4*qt_circle[i+1] + qt_circle[i])
    #         (q_tot).append(q_int)
    #     return np.array(q_tot)
    
    # def force_circle_rev():
    #     q_tot = []
    #     q_int = 0
    #     for i in range(0,len(half_circle_rev)-2,2):
    #         q_int = (half_circle_rev[i+2] - half_circle_rev[i])/6 * (qt_circle_rev[i+2] + 4*qt_circle_rev[i+1] + qt_circle_rev[i])
    #         (q_tot).append(q_int)
    #     return np.array(q_tot)

    # def force_incl():
    #     q_tot_1 = []
    #     q_tot_2 = []
    #     q_int_1 = 0
    #     q_int_2 = 0
    #     for i in range(0,len(half_incl)-2,2):
    #         q_int_1 = (half_incl[i+2]-half_incl[i])/6 *(qt_incl_1[i+2] + 4*qt_incl_1[i+1] + qt_incl_1[i])
    #         q_int_2 = (half_incl[i+2]-half_incl[i])/6 *(qt_incl_2[i+2] + 4*qt_incl_2[i+1] + qt_incl_2[i])
    #         (q_tot_1).append(q_int_1)
    #         (q_tot_2).append(q_int_2)
    #     return np.array(q_tot_1), np.array(q_tot_2)
    
    arm_circ = c.h/2
    arm_incl = c.h/2 * np.sin(np.pi/2 - c.beta_rad)
    eta = 0
    eta += (np.sum(qt_circle)/len(qt_circle) + np.sum(qt_circle_rev)/len(qt_circle_rev))*arm_circ*c.h*np.pi/4
    eta+= (np.sum(qt_incl_1)/len(qt_incl_1) + np.sum(qt_incl_2)/len(qt_incl_2))*arm_incl*c.lsk

    # Fsum_circ = force_circle()
    # Fsum_circ_rev = force_circle_rev()
    # Fsum_incl_1,Fsum_incl_2 = force_incl()


    # eta = (np.sum(Fsum_circ + Fsum_circ_rev))*arm_circ + arm_incl*(np.sum(Fsum_incl_1 + Fsum_incl_2))
    
    return eta

shear_eta = Moment_Eq(Izz)

#%%

def camber_line()
    return np.where((x>=0)&(x<=(c.Ca)),


def dyc_over_dx( x, m, p, c ):
    return np.where((x>=0)&(x<=(c*p)),
                    ((2.0 * m) / np.power(p,2)) * (p - x / c),
                    ((2.0 * m ) / np.power(1-p,2)) * (p - x / c ))

def thickness( x, t, c ):
    term1 =  0.2969 * (np.sqrt(x/c))
    term2 = -0.1260 * (x/c)
    term3 = -0.3516 * np.power(x/c,2)
    term4 =  0.2843 * np.power(x/c,3)
    term5 = -0.1015 * np.power(x/c,4)
    return 5 * t * c * (term1 + term2 + term3 + term4 + term5)

def naca4(x, m, p, t, c=1):
    dyc_dx = dyc_over_dx(x, m, p, c)
    th = np.arctan(dyc_dx)
    yt = thickness(x, t, c)
    yc = camber_line(x, m, p, c)  
    return ((x - yt*np.sin(th), yc + yt*np.cos(th)), 
            (x + yt*np.sin(th), yc - yt*np.cos(th)))


#naca2412 
m = 0.02
p = 0.4
t = 0.12
c = 1.0

x = np.linspace(0,1,200)
for item in naca4(x, m, p, t, c):
    plt.plot(item[0], item[1], 'b')

plt.plot(x, camber_line(x, m, p, c), 'r')
plt.axis('equal')
plt.xlim((-0.05, 1.05))
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 16:58:55 2020

@author: halms
"""
import numpy as np
import matplotlib.pyplot as plt
import readAeroload as aero
import Parameters as c
import Simulation as sim
import matplotlib as mpl
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

#%%

# Function for the base shear flow
def base_shear_y(I_zz):
    
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


qb_circle,qb_circle_rev,qb_incl_1,qb_incl_2,qb_spar,qb_spar_rev = base_shear_y(sim.Izz)

def redundant_shear(I_zz):
    
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


    
def Scz(I_zz):
    qb_circle,qb_circle_rev,qb_incl_1,qb_incl_2,qb_spar,qb_spar_rev = base_shear_y(I_zz)
    redundant_shear_flow = redundant_shear(I_zz)
    q_red_1 = redundant_shear_flow[0]
    q_red_2 = redundant_shear_flow[1]
    
    qt_circle = qb_circle + q_red_1
    qt_circle_rev = qb_circle_rev + q_red_1
    qt_incl_1 = qb_incl_1 + q_red_2
    qt_incl_2 = qb_incl_2 + q_red_2
    qt_spar = qb_spar + q_red_2 - q_red_1

    
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

#%%

def base_shear_z(I_yy):
    
    # Function for the semi-circular part integrating the base shear flow from 0 to pi/2
    def circle():
        q_c = []
        I_c = -c.tsk*c.h/(2*I_yy)
        qval = -0.5*c.Boom_area[0]*(c.Z_locations[0]-sim.z_t)/I_yy

        count = 0
        for i in range(0,len(nodes_circle)-2,2):
            qval += I_c*(nodes_circle[i+2] - nodes_circle[i])/6 * (-c.h/2*np.cos(nodes_circle[i+2]) +
                                                                               4*(-c.h/2*np.cos(nodes_circle[i+1])) -c.h/2*np.cos(nodes_circle[i])+6*(c.h/2-sim.z_t))
            if nodes_circle[i] >= c.angle_rad and count == 0:
                qval += -c.Boom_area[1]*(c.Z_locations[1]-sim.z_t)/I_yy
                count += 1
            
            (q_c).append(qval)
        return np.array(q_c)
    
    # Function for the semi-circular part integrating the base shear flow from -pi/2 to 0
    def circle_rev():
        q_c = []
        I_c = -c.tsk*c.h/(2*I_yy) 
        qval = inclined_2()[-1] - spar_rev()[-1]
        count = 0
        for i in range(0,len(nodes_circle_rev)-2,2):
            qval += I_c*(nodes_circle_rev[i+2] - nodes_circle_rev[i])/6 * (-6*sim.z_t +c.h*0.5*np.cos(nodes_circle_rev[i+2]) + 4*c.h*0.5*np.cos(nodes_circle_rev[i+1]) + c.h*0.5*np.cos(nodes_circle_rev[i]))
            if nodes_circle_rev[i] >= -c.angle_rad and count == 0:
                qval += -c.Boom_area[12]*(c.Z_locations[12]-sim.z_t)/I_yy
                count += 1
            (q_c).append(qval)
            # q_c[-1] += -0.5*c.Boom_area[0]*(c.Z_locations[0]-sim.z_t)/I_yy
        return np.array(q_c)
    
    # Function for the inclined section integrating the base shear flow from 0 to length of this section lsk
    def inclined_1():
        q_c = []
        I_c = -c.tsk/(I_yy)
        qval = circle()[-1] + spar()[-1]
        # qval = 0
        count = 0
        dx = nodes_skin[1] - nodes_skin[0]
        for i in range(0,len(nodes_skin)-2,2):
            # qval += I_c*dx*(-sim.z_t - (nodes_skin[i])*np.cos(c.beta_rad))
            qval += I_c*(nodes_skin[i+2] - nodes_skin[i])/6 *(-6*sim.z_t -np.cos(c.beta_rad)*nodes_skin[i+2] + 4*(-np.cos(c.beta_rad)*nodes_skin[i+1])-np.cos(c.beta_rad)*nodes_skin[i])
            if nodes_skin[i] >= np.abs(c.Z_locations[2]/np.cos(c.beta_rad)) and count == 0:
                qval += -c.Boom_area[2]*(c.Z_locations[2]-sim.z_t)/I_yy
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[3]/np.cos(c.beta_rad)) and count == 1:
                qval += -c.Boom_area[3]*(c.Z_locations[3]-sim.z_t)/I_yy
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[4]/np.cos(c.beta_rad)) and count == 2:
                qval += -c.Boom_area[4]*(c.Z_locations[4]-sim.z_t)/I_yy
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[5]/np.cos(c.beta_rad)) and count == 3:
                qval += -c.Boom_area[5]*(c.Z_locations[5]-sim.z_t)/I_yy
                count += 1
            elif nodes_skin[i] >= np.abs(c.Z_locations[6]/np.cos(c.beta_rad)) and count == 4:
                qval += -c.Boom_area[6]*(c.Z_locations[6]-sim.z_t)/I_yy
                count += 1
            (q_c).append(qval)
        return np.array(q_c)

    # Function for the inclined section integrating the base shear flow from 0 to length of this section lsk
    def inclined_2():
        q_c = []
        I_c = -c.tsk/(I_yy)
        qval = inclined_1()[-1]
        count = 0
        for i in range(0,len(nodes_skin)-2,2):
            qval += I_c*(nodes_skin[i+2] - nodes_skin[i])/6 * (-(c.Ca - c.h/2)-sim.z_t + np.cos(c.beta_rad)*nodes_skin[i+2] + 4*(-(c.Ca - c.h/2)-sim.z_t + np.cos(c.beta_rad)*nodes_skin[i+1]) -(c.Ca - c.h/2)-sim.z_t + np.cos(c.beta_rad)*nodes_skin[i])
            if nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[7]/np.cos(c.beta_rad)) and count == 0:
                qval += -c.Boom_area[7]*(c.Z_locations[7]-sim.z_t)/I_yy
                count+=1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[8]/np.cos(c.beta_rad)) and count == 1:
                qval += -c.Boom_area[8]*(c.Z_locations[8]-sim.z_t)/I_yy
                count += 1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[9]/np.cos(c.beta_rad)) and count == 2:
                qval += -c.Boom_area[9]*(c.Z_locations[9]-sim.z_t)/I_yy
                count += 1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[10]/np.cos(c.beta_rad)) and count == 3:
                qval += -c.Boom_area[10]*(c.Z_locations[10]-sim.z_t)/I_yy
                count += 1
            elif nodes_skin[i] >= c.lsk - np.abs(c.Z_locations[11]/np.cos(c.beta_rad)) and count == 4:
                qval += -c.Boom_area[11]*(c.Z_locations[11]-sim.z_t)/I_yy
                count += 1
            (q_c).append(qval)
        
        return np.array(q_c)
    # Function for the spar section integrating the base shear flow from 0 to h/2
    def spar():
        q_c = []
        I_c = sim.z_t*c.tsp/I_yy
        qval = 0
        for i in range(0,len(nodes_spar)-2,2):
            qval += I_c*(nodes_spar[i+2] - nodes_spar[i])
            (q_c).append(qval)
        return np.array(q_c)
    
    # Function for the spar section integrating the base shear flow from 0 to -h/2
    def spar_rev():
        q_c = []
        I_c = sim.z_t*c.tsp/I_yy
        qval = 0
        for i in range(0,len(nodes_spar_rev)-2,2):
            qval += I_c*(nodes_spar_rev[i+2] - nodes_spar_rev[i])
            (q_c).append(qval)
        return np.array(q_c)
    
    
    return circle(),circle_rev(),inclined_1(),inclined_2(),spar(),spar_rev()

def Scy(I_yy):
    qb_circle,qb_circle_rev,qb_incl_1,qb_incl_2,qb_spar,qb_spar_rev = base_shear_z(I_yy)

    
    qt_circle = qb_circle 
    qt_circle_rev = qb_circle_rev
    qt_incl_1 = qb_incl_1
    qt_incl_2 = qb_incl_2 
    qt_spar = qb_spar

    plt.figure(1)
    plt.plot(nodes_skin[0:-1:2],qt_incl_1)
    plt.plot(nodes_skin[0:-1:2],qt_incl_2)
    
    half_circle = nodes_circle[0:-1:2]
    half_circle_rev = nodes_circle_rev[0:-1:2]
    half_incl = nodes_skin[0:-1:2]

    
    arm_circ = c.h/2
    arm_incl = c.h/2 * np.sin(np.pi/2 - c.beta_rad)
    xsi = 0
    xsi += (np.sum(qt_circle)/len(qt_circle) + np.sum(qt_circle_rev)/len(qt_circle_rev))*arm_circ*c.h*np.pi/4
    xsi+= (np.sum(qt_incl_1)/len(qt_incl_1) + np.sum(qt_incl_2)/len(qt_incl_2))*arm_incl*c.lsk

    
    return xsi

print(Scy(sim.Iyy))
print(Scz(sim.Izz))


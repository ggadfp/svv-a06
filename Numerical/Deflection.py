# -*- coding: utf-8 -*-
"""
Created on Thursday 20 2020

@author: Gabriel
"""

import numpy as np
import matplotlib.pyplot as plt
import Simulation as sim
import time
import CRJ700_VerificiationModel_copy_gabriel as verif
import J as tor
import ShearDist as sd
#%%
start_time = time.time()
#[H1y,H2y,H3y,H1z,H2z,H3z,A1,C1,C2,C3,C4,C5]
#We now wnat to create the 12x12 matrix and the vector of BCs e to solve for the unknown forces thanks to our deflection insights
mat_defl = np.zeros([12,12])
e = np.zeros([12,1])
nodes = 800
plot_points = 51
J = tor.torsion_unit()
scz = sd.Scz(sim.Izz)
# sim.s_coeff_load = sim.s_coeff_load * 0
# sim.s_coeff_torque = sim.s_coeff_torque * 0
# sim.P = 0
# sim.d1 = 0
# sim.d3 = 0

# From BC 1 (check Simulation Plan)
e[0] = (sim.doubleinteg_spline(sim.s_coeff_load, sim.la, sim.xfull, nodes)+
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
e[2] = (sim.integ_spline(sim.s_coeff_torque,sim.la,sim.xfull)+
        sim.P*sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad)))
mat_defl[2] = [
    0,0,0,0,0,0,
    sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad)),
    0,0,0,0,0
    ]
print("BC1+2+3 took", time.time() - start_time, "to run")
# From BC 4 (check Simulation Plan)
e[3] = (sim.integ_spline(sim.s_coeff_load,sim.la,sim.xfull)+
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
e[5] = sim.quadrupleinteg_spline(sim.s_coeff_load,sim.x2,sim.xfull,nodes)
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
e[7] = (sim.quadrupleinteg_spline(sim.s_coeff_load,sim.x1,sim.xfull,nodes)/(sim.E*sim.Izz)+
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
print("BC1+2+3+4+5+6+7+8+9 took", time.time() - start_time, "to run")
# From BC 10 (check Simulation Plan)
e[9] = (sim.d3*np.cos(sim.theta_rad)*(sim.E*sim.Izz)+
        sim.quadrupleinteg_spline(sim.s_coeff_load,sim.x3,sim.xfull,nodes)+
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
e[11] = sim.doubleinteg_spline(sim.s_coeff_torque,sim.x2-sim.xa/2,sim.xfull,nodes)
mat_defl[11] = [
    0,0,0,0,0,0,0,0,0,0,0,
    sim.G*(J)
    ]
print("BC1+2+3+4+5+6+7+8+9+10+11+12 took", time.time() - start_time, "to run")

#SOLVING FOR ALL THE UNKNOWN FORCES
unknowns = np.linalg.inv(mat_defl).dot(e)

#%%
def macaulay(xloc,start,n):
    d = xloc-start
    if d <= 0:
        return 0
    if d >0:
        return d**n

#%%
# SHEAR FORCE IN Y SPANWISE

def S_y(xloc):
    S_y = (sim.integ_spline(sim.s_coeff_load,xloc,sim.xfull)
           + sim.P*np.sin(sim.theta_rad)*macaulay(xloc,(sim.x2+sim.xa/2),0)
           - unknowns[0]*macaulay(xloc,sim.x1,0)
           - unknowns[1]*macaulay(xloc,sim.x2,0)
           - unknowns[2]*macaulay(xloc,sim.x3,0)
           - unknowns[6]*np.sin(sim.theta_rad)*macaulay(xloc,(sim.x2-sim.xa/2),0)
           )
    
    return S_y

Sy_dist = [S_y(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+Sy took", time.time() - start_time, "to run")

#%%
# MOMENT IN Z SPANWISE
    
def M_z(xloc):
    
    M_z = (sim.doubleinteg_spline(sim.s_coeff_load, xloc, sim.xfull, nodes)
     +sim.P*np.sin(sim.theta_rad)*macaulay(xloc,(sim.x2+sim.xa/2),1)
     - unknowns[0]*macaulay(xloc,sim.x1,1)
     - unknowns[1]*macaulay(xloc,sim.x2,1)
     - unknowns[2]*macaulay(xloc,sim.x3,1)
     - unknowns[6]*np.sin(sim.theta_rad)*macaulay(xloc,(sim.x2-sim.xa/2),1))
    
    return M_z

# print(M_z(sim.la))

Mz_dist = [M_z(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+Sy+Mz took", time.time() - start_time, "to run")

#%%
# SLOPE IN Y SPANWISE

def dvy_dx(xloc):
    
    dvy_dx = -(sim.tripleinteg_spline(sim.s_coeff_load,xloc,sim.xfull,nodes)
              + sim.P*np.sin(sim.theta_rad)/2*macaulay(xloc,(sim.x2+sim.xa/2),2)
              - unknowns[0]/2*macaulay(xloc,sim.x1,2)
              - unknowns[1]/2*macaulay(xloc,sim.x2,2)
              - unknowns[2]/2*macaulay(xloc,sim.x3,2)
              - unknowns[6]*np.sin(sim.theta_rad)/2*macaulay(xloc,(sim.x2-sim.xa/2),2) 
              )/(sim.E*sim.Izz) +unknowns[7]
    
    return dvy_dx

dvydx_dist = [dvy_dx(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+Sy+Mz+dvydx took", time.time() - start_time, "to run")

#%%
# DEFLECTION IN Y SPANWISE

def deflection_y(xloc):
    
    def_y = -(sim.quadrupleinteg_spline(sim.s_coeff_load,xloc,sim.xfull,nodes)
              + sim.P*np.sin(sim.theta_rad)/6*macaulay(xloc,(sim.x2+sim.xa/2),3)
              - unknowns[0]/6*macaulay(xloc,sim.x1,3)
              - unknowns[1]/6*macaulay(xloc,sim.x2,3)
              - unknowns[2]/6*macaulay(xloc,sim.x3,3)
              - unknowns[6]*np.sin(sim.theta_rad)/6*macaulay(xloc,(sim.x2-sim.xa/2),3) 
              )/(sim.E*sim.Izz) +unknowns[7]*xloc+unknowns[8]
    return def_y

# print(sim.d1*np.cos(sim.theta_rad),deflection_y(sim.x1))
# print(deflection_y(sim.x2))
# print(sim.d3*np.cos(sim.theta_rad),deflection_y(sim.x3))

deflectiony_dist = [deflection_y(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+Sy+Mz+dvydx+deflectiony took", time.time() - start_time, "to run")

#%%
#SHEAR FORCE IN Z SPANWISE

def S_z(xloc):
    Sz = (+ sim.P*np.cos(sim.theta_rad)*macaulay(xloc,(sim.x2+sim.xa/2),0)
          - unknowns[3]*macaulay(xloc,sim.x1,0)
          - unknowns[4]*macaulay(xloc,sim.x2,0)
          - unknowns[5]*macaulay(xloc,sim.x3,0)
          - unknowns[6]*np.cos(sim.theta_rad)*macaulay(xloc,(sim.x2-sim.xa/2),0)
          )
    return Sz

Sz_dist = [S_z(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+(all y)+Sz took", time.time() - start_time, "to run")

#%%
# MOMENT IN Y SPANWISE
    
def M_y(xloc):
    
    M_y = (+sim.P*np.cos(sim.theta_rad)*macaulay(xloc,(sim.x2+sim.xa/2),1)
           - unknowns[3]*macaulay(xloc,sim.x1,1)
           - unknowns[4]*macaulay(xloc,sim.x2,1)
           - unknowns[5]*macaulay(xloc,sim.x3,1)
           - unknowns[6]*np.cos(sim.theta_rad)*macaulay(xloc,(sim.x2-sim.xa/2),1))    
    return M_y

My_dist = [M_y(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+(all y)+Sz+My took", time.time() - start_time, "to run")

#%%
# SLOPE IN Z SPANWISE

def dvz_dx(xloc):
    
    dvz_dx = -(+ sim.P*np.cos(sim.theta_rad)/2*macaulay(xloc,(sim.x2+sim.xa/2),2)
              - unknowns[3]/2*macaulay(xloc,sim.x1,2)
              - unknowns[4]/2*macaulay(xloc,sim.x2,2)
              - unknowns[5]/2*macaulay(xloc,sim.x3,2)
              - unknowns[6]*np.cos(sim.theta_rad)/2*macaulay(xloc,(sim.x2-sim.xa/2),2) 
              )/(sim.E*sim.Iyy) +unknowns[9]
    
    return dvz_dx

dvzdx_dist = [dvz_dx(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+(all y)+Sz+My+dvzdx took", time.time() - start_time, "to run")

#%%
# DEFLECTION IN Z SPANWISE

def deflection_z(xloc):
    
    def_z = -(+ sim.P*np.cos(sim.theta_rad)/6*macaulay(xloc,(sim.x2+sim.xa/2),3)
              - unknowns[3]/6*macaulay(xloc,sim.x1,3)
              - unknowns[4]/6*macaulay(xloc,sim.x2,3)
              - unknowns[5]/6*macaulay(xloc,sim.x3,3)
              - unknowns[6]*np.cos(sim.theta_rad)/6*macaulay(xloc,(sim.x2-sim.xa/2),3) 
              )/(sim.E*sim.Iyy) +unknowns[9]*xloc+unknowns[10]
    return def_z

# print(sim.d1*np.sin(sim.theta_rad),deflection_z(sim.x1))
# print(deflection_z(sim.x2))
# print(sim.d3*np.sin(sim.theta_rad),deflection_z(sim.x3))

deflectionz_dist = [deflection_z(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+(all y)+Sz+My+dvzdx+deflectionz took", time.time() - start_time, "to run")

#%%
#TORQUE SPANWISE

def T_span(xloc):
    T = (sim.integ_spline(sim.s_coeff_torque,xloc,sim.xfull)
         +sim.P*sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad))*macaulay(xloc,(sim.x2+sim.xa/2),0)
         
         -unknowns[6]*sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad))*macaulay(xloc,(sim.x2-sim.xa/2),0)
         
         -unknowns[0]*scz*macaulay(xloc,sim.x1,0)*10
         -unknowns[1]*scz*macaulay(xloc,sim.x2,0)*10
         -unknowns[2]*scz*macaulay(xloc,sim.x3,0)*10)
    return T

T_dist = [T_span(i) for i in np.linspace(0,sim.la,plot_points)]
print("BCs+(all y)+(all z)+T took", time.time() - start_time, "to run")

#%%
#ANGLE OF TWIST SPANWISE

def twist(xloc):
    twist = (sim.doubleinteg_spline(sim.s_coeff_torque,xloc,sim.xfull,nodes)
             +sim.P*sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad))*macaulay(xloc,(sim.x2+sim.xa/2),1)
             -unknowns[6]*sim.h/2*(np.cos(sim.theta_rad)-np.sin(sim.theta_rad))*macaulay(xloc,(sim.x2-sim.xa/2),1)
             )/(sim.G*(J))+unknowns[11]
    return twist

twist_dist = [twist(i) for i in np.linspace(0,sim.la,plot_points)]

#%%
plt.close('all')
axis = np.linspace(0,sim.la,plot_points)
fig, axs = plt.subplots(2, 2)
plt.suptitle('Deflection in the y-axis',fontsize=20)

axs[0, 0].scatter(axis, Sy_dist,label='A06 Simulation')
axs[0, 0].plot(verif.x, verif.d3y, 'tab:pink',label='Verification code (N=20)')
axs[0, 0].set_title(r'$S_{y}(x)$')
axs[0, 0].set_ylabel('[N]')
axs[0, 0].legend(loc="upper right")
axs[0, 0].grid()

axs[0, 1].scatter(axis, Mz_dist,label='A06 Simulation')
axs[0, 1].plot(verif.x, verif.d2y, 'tab:orange',label='Verification code (N=20)')
axs[0, 1].set_title(r'$M_{z}(x)$')
axs[0, 1].set_ylabel('[Nm]')
axs[0, 1].grid()
axs[0, 1].legend(loc="lower right")

axs[1, 0].scatter(axis, dvydx_dist,label='A06 Simulation')
axs[1, 0].plot(verif.x, verif.d1y, 'tab:green',label='Verification code (N=20)')
axs[1, 0].set_title(r'$\frac{d\nu(x)}{dx}$')
axs[1, 0].set_ylabel('[-]')
axs[1, 0].grid()
axs[1, 0].legend(loc="upper right")
axs[1, 0].set_xlabel('Spanwise location [m]')


axs[1, 1].scatter(axis, deflectiony_dist,label='A06 Simulation')
axs[1, 1].plot(verif.x, verif.defy_v, 'tab:red',label='Verification code (N=20)')
axs[1, 1].set_title(r'$\nu(x)$')
axs[1, 1].set_ylabel('[m]')
axs[1, 1].grid()
axs[1, 1].legend(loc="lower right")
axs[1, 1].set_xlabel('Spanwise location [m]')

#%%

# plt.close('all')
axis = np.linspace(0,sim.la,plot_points)
fig2, axs2 = plt.subplots(2, 2)
plt.suptitle('Deflection in the z-axis',fontsize=20)

axs2[0, 0].scatter(axis, Sz_dist,label='A06 Simulation')
axs2[0, 0].plot(verif.x, verif.d3z, 'tab:pink',label='Verification code (N=20)')
axs2[0, 0].set_title(r'$S_{z}(x)$')
axs2[0, 0].set_ylabel('[N]')
axs2[0, 0].legend(loc="upper right")
axs2[0, 0].grid()

axs2[0, 1].scatter(axis, My_dist,label='A06 Simulation')
axs2[0, 1].plot(verif.x, verif.d2z, 'tab:orange',label='Verification code (N=20)')
axs2[0, 1].set_title(r'$M_{y}(x)$')
axs2[0, 1].set_ylabel('[Nm]')
axs2[0, 1].grid()
axs2[0, 1].legend(loc="upper right")

axs2[1, 0].scatter(axis, dvzdx_dist,label='A06 Simulation')
axs2[1, 0].plot(verif.x, verif.d1z, 'tab:green',label='Verification code (N=20)')
axs2[1, 0].set_title(r'$\frac{d\eta(x)}{dx}$')
axs2[1, 0].set_ylabel('[-]')
axs2[1, 0].grid()
axs2[1, 0].legend(loc="upper right")
axs2[1, 0].set_xlabel('Spanwise location [m]')


axs2[1, 1].scatter(axis, deflectionz_dist,label='A06 Simulation')
axs2[1, 1].plot(verif.x, verif.defz_v, 'tab:red',label='Verification code (N=20)')
axs2[1, 1].set_title(r'$\eta(x)$')
axs2[1, 1].set_ylabel('[m]')
axs2[1, 1].grid()
axs2[1, 1].legend(loc="upper right")
axs2[1, 1].set_xlabel('Spanwise location [m]')

#%%
axis = np.linspace(0,sim.la,plot_points)
fig3, axs3 = plt.subplots(2)
plt.suptitle('Twist around the x-axis',fontsize=20)

axs3[0].scatter(axis, T_dist, label='A06 Simulation')
axs3[0].plot(verif.x, verif.d2t, 'tab:orange',label='Verification code (N=20)')
axs3[0].set_ylabel('T(x) [Nm]')
axs3[0].grid()
axs3[0].legend(loc="lower right")
axs3[0].set_xlabel('Spanwise location [m]')

axs3[1].plot(axis, twist_dist, label='A06 Simulation')
axs3[1].plot(verif.x, verif.twist_v, 'tab:red',label='Verification code (N=20)')
axs3[1].set_ylabel(r'$\alpha(x)$''[rad]')
axs3[1].grid()
axs3[1].legend(loc="upper right")
axs3[1].set_xlabel('Spanwise location [m]')

#%%
# #PLOTTING THE ERROR OF THE DEFLECTION
# ey_P = np.zeros(len(verif.defz_v))
# ez_P = np.zeros(len(verif.defz_v))
# for i in range(len(verif.defz_v)):
#     ey_P[i] = (-deflectiony_dist[i][0]-verif.defy_v[i])
#     ez_P[i] = (deflectionz_dist[i][0]-verif.defz_v[i])
    
# plt.plot(axis,ey_q,label='Case q')
# plt.plot(axis,ey_P,label='Case P')
# plt.grid()
# plt.legend(loc="lower right")
# plt.xlabel('Spanwise location [m]')
# plt.ylabel('Discrepancies in 'r'$\nu(x)$' '[m]')
# plt.show()

# # print(max(abs(ey))/max(abs(verif.defy_v))*100)
# # print(max(abs(ez))/max(abs(verif.defz_v))*100)

# #%%

# plt.close('all')

# plt.plot(axis, -np.array(deflectiony_dist),label='A06 Simulation')
# plt.plot(verif.x, verif.defy_v, 'tab:red',label='Verification code (N=20)')
# plt.title(r'$\nu(x)$')
# plt.ylabel('[m]')
# plt.grid()
# plt.legend(loc="lower left")
# plt.xlabel('Spanwise location [m]')

print("Deflection took", time.time() - start_time, "to run")
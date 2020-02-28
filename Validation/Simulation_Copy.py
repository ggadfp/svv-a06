# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:03:57 2020

@author: Ihab
"""

import numpy as np
import matplotlib.pyplot as plt
import readAeroload_Copy
from mpl_toolkits.mplot3d import Axes3D
import time
#%%

# Number of nodes
Nx = 41 # Spanwise Nodes
Nz = 81 # Chordwise Nodes

# Geometry of Aileron in [m]
Ca = 0.605  # Chord length
la = 2.661 # Span
x1 = 0.172 # x-location of hinge 1
x2 = 1.211 # x-location of hinge 2
x3 = 2.591 # x-location of hinge 3
xa = 35e-2 # Distance between actuator 1 and 2
h = 20.5e-2 # Aileron height
tsk = 1.1e-3 # Thickness of the skin
tsp = 2.8e-3 # Thickness of the spar
lsk = np.sqrt((Ca-h/2)**2 + (h/2)**2) # Length of the inclined skin
Ask_incl = lsk*tsk # Area of inclined section
Ask_semi = np.pi*tsk*h/2 # Area of semi-circular section
Asp = h*tsp # Area of spar

peri = lsk*2 + np.pi * h/2 # Perimeter of aileron

# Geometry of Stiffener in [m]
tst = 2.8e-3 # Thickness of the stiffener
hst = 1.6e-2 # Height of the stiffener
wst = 1.9e-2 # Width of the stiffener
nst = 15 # [-] -- Number of stiffeners
Ast = (hst + wst)*tst # Stiffener area
deltast = peri/nst # Stiffener spacing

# Given conditions
d1 = 1.154e-2 # [m] -- Vertical displacement of hinge 1
d3 = 1.840e-2 # [m] -- Vertical displacement of hinge 3
theta_deg = 28 # [deg] -- Maximum upward deflection
theta_rad = theta_deg*np.pi/180 # [rad] -- Maximum upward deflection
P = 97.4e+3 # [N] -- Load acting on actuator 2

# Necessary angles
angle_rad = deltast/(h/2) # Angle between boom at symmetry axis and the next boom
angle_deg = angle_rad *180/np.pi # Same angle, but in degrees
beta_rad = np.arctan2((h/2),(Ca-h/2)) # Angle of the inclined section
beta_deg = beta_rad*180/np.pi # Same angle, but in degrees

# Material properties
E = 72.9e09 # [Pa] Young's modulus
G = 27.1e09 # [Pa] Shear modulus

# Other parameters
bridge = (-(np.pi/2 - angle_rad)*h/2 + deltast) # Bridging distance between semi-circular and inclined section
Boom_area = nst*[Ast] # All boom areas, where index 0 indicates the boom at the symmetry axis, going ccw

# Z-coordinates of booms wrt Hinge line
Z_locations = [h/2, h*np.cos(angle_rad)/2, -bridge*np.cos(beta_rad), -(bridge + deltast)*np.cos(beta_rad), - (bridge + 2*deltast)*np.cos(beta_rad),
                  -(bridge + 3*deltast)*np.cos(beta_rad), -(bridge + 4*deltast)*np.cos(beta_rad), -(bridge + 4*deltast)*np.cos(beta_rad),
                  -(bridge + 3*deltast)*np.cos(beta_rad),-(bridge + 2*deltast)*np.cos(beta_rad),-(bridge + deltast)*np.cos(beta_rad),
                  -bridge*np.cos(beta_rad),h*np.cos(angle_rad)/2] 

# Y-coordinates of booms wrt Hinge line
Y_locations = [0,-h*np.sin(angle_rad)/2,-h/2+bridge*np.sin(beta_rad),-h/2+(bridge+deltast)*np.sin(beta_rad),
               -h/2+(bridge+2*deltast)*np.sin(beta_rad),-h/2+(bridge+3*deltast)*np.sin(beta_rad),-h/2+(bridge+4*deltast)*np.sin(beta_rad),
               h/2-(bridge+4*deltast)*np.sin(beta_rad),h/2-(bridge+3*deltast)*np.sin(beta_rad),h/2-(bridge+2*deltast)*np.sin(beta_rad),
               h/2-(bridge+deltast)*np.sin(beta_rad),h/2-bridge*np.sin(beta_rad),h*np.sin(angle_rad)/2]
    
#plt.plot(-np.array(Z_locations),Y_locations)

#%%

def Centroid(z_loc,z_area,chord,r,A_circle,A_length,A_spar,choice):
    
    if choice == 'y':
        y_tilde = 0
        return y_tilde
    if choice == 'z':
        sum_areadist = A_circle*4*r/(3*np.pi) + A_length*2*(r-chord)/2
        sum_area = A_circle + A_length*2 + A_spar
        for i in range(len(z_loc)):
            sum_areadist += z_loc[i]*z_area[i]
            sum_area += z_area[i]
        z_tilde = sum_areadist/sum_area
        return z_tilde

z_t = Centroid(Z_locations,Boom_area,Ca,h/2,Ask_semi,Ask_incl,Asp,'z') # Centroid location [m]

#%%

# Computes moment of inertia in the yy plane
def MoI_yy(z_tilde,z_loc,boom_area,A_spar,A_circle,A_skin,t_skin,l_skin,chord,r,alpha):
    

    Iyy_skin = l_skin**3 * t_skin * (np.cos(alpha))**2 / 12 + A_skin * ( (r-chord) / 2 - z_tilde)**2
    Iyy_spar = A_spar * z_tilde**2
    Iyy_circle = (np.pi * r**3 * t_skin) / 2 + A_circle * (4 * r / (3 * np.pi) - z_tilde)**2
    I_yy = 2*Iyy_skin + Iyy_spar + Iyy_circle
    for i in range(len(z_loc)):
        I_yy+= boom_area[i]*(z_loc[i]-z_tilde)**2
    return I_yy

# Computes moment of inertia in the zz plane
def MoI_zz(y_loc,boom_area,A_spar,t_spar,A_circle,A_skin,t_skin,l_skin,chord,r,alpha):
    
#    Izz_skin = l_skin**3 * t_skin * (np.sin(alpha))**2 / 12 + A_skin * (l_skin * np.sin(alpha) / 2)**2
    Izz_skin = (2*l_skin)**3 * t_skin * (np.sin(alpha))**2 / 12
    Izz_spar = (2 * r)**3 * t_spar / 12
    Izz_circle = np.pi * r**3 * t_skin / 2 
    I_zz = Izz_skin + Izz_spar + Izz_circle
    for i in range(len(y_loc)):
        I_zz+= boom_area[i]*y_loc[i]**2
    return I_zz

Iyy = MoI_yy(z_t,Z_locations,Boom_area,Asp,Ask_semi,Ask_incl,tsk,lsk,Ca,h/2,beta_rad) # [m^4]
Izz = MoI_zz(Y_locations,Boom_area,Asp,tsp,Ask_semi,Ask_incl,tsk,lsk,Ca,h/2,beta_rad) # [m^4]

#%%


x,z = readAeroload_Copy.meshCreator(Nz,Nx,Ca,la,h)

# load,pointapp = readAeroload.zload_zpoint(Nz,Nx,z)
# torque = np.zeros((Nx,1))
# for i in range(Nx):
#     torque[i] = load[i]*pointapp[i]
    
load = np.ones(len(x[0,:]))*5.54e3
torque = load*(0.25*Ca-h/2)

#%%
#  WE'RE NOT USING THIS
    
# def matrix(Nsw,xvec):
#     A = np.zeros((Nsw,Nsw))
#     A[:,0] = 1
#     xvec_new = (2*np.pi)/(xvec[-1] - xvec[0]) * xvec
#     for i in range(Nsw):
#         for j in range(1,Nsw):
#             if 1<=j<=20:
#                 A[i,j] = np.cos(j*xvec_new[i])
#             else:
#                 A[i,j] = np.sin(j*xvec_new[i])
                
#     return np.matrix(A)

# L = matrix(Nx,x[0,:])

# coeff_torque = np.linalg.inv(L).dot(torque)
# load = np.array(load).reshape((Nx,1))
# coeff_load = np.linalg.inv(L).dot(load)

# load_test = L.dot(coeff_load)
# torque_test = L.dot(coeff_torque)

# plt.figure(1)
# plt.plot(x[0,:],load_test,x[0,:],load)

# r = np.linalg.norm(load-load_test)/np.linalg.norm(load) * 100 # FAA checken

# plt.figure(2)
# plt.plot(x[0,:],torque_test,x[0,:],torque_test)

#%%
def cubicspline(Nsw,xvec,data):
    A = np.zeros((Nsw-2,Nsw-2))
    f = np.zeros((Nsw-2,1))
    for i in range(1,Nsw-1):
        him1 = xvec[i]-xvec[i-1]
        hi = xvec[i+1]-xvec[i]
        
        #Am=f with m a vector conteinin (M1,M2,M3,...,MNsw-1)
        #entries to the matrix A, M0 and Mnsw-1
        if i == 1:
            A[i-1,i-1] = (him1+hi)/3
            A[i-1,i] = (hi)/6
        elif 1<i<Nsw-2:
            A[i-1,i-2] = (him1)/6
            A[i-1,i-1] = (him1+hi)/3
            A[i-1,i] = (hi)/6
        elif i == Nsw-2:
            A[i-1,i-2] = (him1)/6
            A[i-1,i-1] = (him1+hi)/3
            
        #entries of the solution vector f
        f[i-1,0] = (data[i+1]-data[i])/hi - (data[i]-data[i-1])/him1
    
    # m = A^-1 f
    m = np.vstack(([0],np.linalg.inv(A).dot(f),[0]))
    
    s_coeff = np.zeros([Nsw-1,4])
    for k in range(Nsw-1):
        hk = xvec[k+1]-xvec[k] #spacing between point k and the next one
        
        s_coeff[k,0] = data[k] #constant term
        s_coeff[k,1] = (data[k+1]-data[k])/hk-hk/3*m[k]-hk/6*m[k+1] #linear term (is multiplied by (x-xk))
        s_coeff[k,2] = m[k]/2 #quadratic term (is multiplied by (x-xk)^2)
        s_coeff[k,3] = (m[k+1]-m[k])/(6*hk) #cubic term (is multiplied by (x-xk)^3)
    
    return s_coeff
           

 #%%
def linearspline(Nsw,xvec,data):
    
    s_coeff = np.zeros([Nsw-1,4])
    for k in range(Nsw-1):
        
        s_coeff[k,0] = data[k] #constant term
        s_coeff[k,1] = (data[k+1]-data[k])/(xvec[k+1]-xvec[k]) #linear term (is multiplied by (x-xk))
        s_coeff[k,2] = 0 #quadratic term (is multiplied by (x-xk)^2)
        s_coeff[k,3] = 0 #cubic term (is multiplied by (x-xk)^3)
        
    last_value = (s_coeff[-1,0] +
                  s_coeff[-1,1]*(xvec[-1]-xvec[-2]) +
                  s_coeff[-1,2]*(xvec[-1]-xvec[-2])**2 +
                  s_coeff[-1,3]*(xvec[-1]-xvec[-2])**3)
    s_coeff = np.vstack(([s_coeff[0,0],0,0,0],s_coeff,[last_value,0,0,0]))
    
    return s_coeff
    
#%%   
s_coeff_torque =  linearspline(Nx,x[0,:],torque)
s_coeff_load = linearspline(Nx,x[0,:],load)

#checking if the splines have points and derivatives that connect
r_1 = [0]*(Nx-1)
for i in range(Nx-2):
    si = (s_coeff_load[i,0] +
          s_coeff_load[i,1]*(x[0,i+1]-x[0,i]))
    r_1[i] = abs(si-s_coeff_load[i+1,0])    

#%%
#adding the beginning and end points
xfull = np.vstack(([0],x[0].reshape((len(x[0]),1)),[la]))
   
#%%   
def function_value(s_coeff,xloc,xvec):
    
    #finding the index of the last node before the xloc
    diff = np.array(xvec - xloc)
    if xloc == 0:
        idx = 0
    elif xloc == la:
        idx = len(xvec)-2
    else:
        idx = int(np.where(diff < 0)[0][-1])
    #this one is the value of the function at that point
    value = (s_coeff[idx,0] +
            s_coeff[idx,1]*(xloc-xvec[idx]) +
            s_coeff[idx,2]*(xloc-xvec[idx])**2 +
            s_coeff[idx,3]*(xloc-xvec[idx])**3)
    
    return value

#%%
def integ_spline(s_coeff,xloc,xvec):       
    #finding the index of the last node before the xloc
    diff = np.array(xvec - xloc)
    if xloc == 0:
        idx = 0
    elif xloc == la:
        idx = len(s_coeff)-2
    else:
        idx = int(np.where(diff < 0)[0][-1])    
  
    full_sum = 0            
    #this one is to get the value of the different integrals
    for k in range(idx+1):     
        if k < idx:                
            full_sum += (s_coeff[k,0]*(xvec[k+1]-xvec[k]) +
                         (s_coeff[k,1]*(xvec[k+1]-xvec[k])**2)/2 +
                         s_coeff[k,2]/3*(xvec[k+1]-xvec[k])**3 +
                         s_coeff[k,3]/4*(xvec[k+1]-xvec[k])**4)
        elif k == idx:
            full_sum += (s_coeff[k,0]*(xloc-xvec[k]) +
                         s_coeff[k,1]/2*(xloc-xvec[k])**2 +
                         s_coeff[k,2]/3*(xloc-xvec[k])**3 +
                         s_coeff[k,3]/4*(xloc-xvec[k])**4)
    return float(full_sum)
#%%
def doubleinteg_spline(s_coeff,xloc,xvec,n):
    first_integral_full = [integ_spline(s_coeff,i,xvec) for i in np.linspace(0,la,n+1)]
    coeff_f_integral = linearspline(len(np.linspace(0,la,n+1)),np.linspace(0,la,n+1),first_integral_full)
    second_integral = integ_spline(coeff_f_integral,xloc,np.linspace(0,la,n+1))
    return second_integral
def doubleinteg_spline_list(s_coeff,xvec,n):
    first_integral_full = [integ_spline(s_coeff,i,xvec) for i in np.linspace(0,la,n+1)]
    coeff_f_integral = linearspline(len(np.linspace(0,la,n+1)),np.linspace(0,la,n+1),first_integral_full)
    second_integral = [integ_spline(coeff_f_integral,i,np.linspace(0,la,n+1)) for i in np.linspace(0,la,n+1)]
    return second_integral
#%%
def tripleinteg_spline(s_coeff,xloc,xvec,n):
    second_integral_full = doubleinteg_spline_list(s_coeff, xvec, n)
    coeff_s_integral = linearspline(len(np.linspace(0,la,n+1)),np.linspace(0,la,n+1),second_integral_full)
    third_integral = integ_spline(coeff_s_integral,xloc,np.linspace(0,la,n+1))
    return third_integral
def tripleinteg_spline_list(s_coeff,xvec,n):
    second_integral_full = doubleinteg_spline_list(s_coeff,xvec,n)
    coeff_s_integral = linearspline(len(np.linspace(0,la,n+1)),np.linspace(0,la,n+1),second_integral_full)
    third_integral = [integ_spline(coeff_s_integral,i,np.linspace(0,la,n+1)) for i in np.linspace(0,la,n+1)]
    return third_integral
#%%
def quadrupleinteg_spline(s_coeff,xloc,xvec,n):
    third_integral_full = tripleinteg_spline_list(s_coeff, xvec, n)
    coeff_t_integral = linearspline(len(np.linspace(0,la,n+1)),np.linspace(0,la,n+1),third_integral_full)
    fourth_integral = integ_spline(coeff_t_integral,xloc,np.linspace(0,la,n+1))
    return fourth_integral
def quadrupleinteg_spline_list(s_coeff,xvec,n):
    third_integral_full = tripleinteg_spline_list(s_coeff, xvec, n)
    coeff_t_integral = linearspline(len(np.linspace(0,la,n+1)),np.linspace(0,la,n+1),third_integral_full)
    fourth_integral = [integ_spline(coeff_t_integral,i,np.linspace(0,la,n+1)) for i in np.linspace(0,la,n+1)]
    return fourth_integral 
#%%
# start_time = time.time()

# n = 100 #number of interpolating splines

# test_data = x[0]
# coeff_testing = linearspline(Nx,x[0],test_data)
# test_returned_values = [function_value(coeff_testing,i,xfull) for i in np.linspace(0,la,n+1)]
# integral_test = [integ_spline(coeff_testing,i,xfull) for i in np.linspace(0,la,n+1)]
# print("Integral 1 took", time.time() - start_time, "to run")
# integral_test_1 = doubleinteg_spline_list(coeff_testing,xfull,n)
# print("Integral 1+2 took", time.time() - start_time, "to run")
# integral_test_2 = tripleinteg_spline_list(coeff_testing, xfull,n)
# print("Integral 1+2+3 took", time.time() - start_time, "to run")
# integral_test_3 = quadrupleinteg_spline_list(coeff_testing, xfull,n)

# plt.plot(x[0],test_data)
# plt.scatter(np.linspace(0,la,len(test_returned_values)),test_returned_values)
# plt.plot(test_data,np.array(test_data)**2/2)
# plt.scatter(np.linspace(0,la,len(integral_test)),integral_test)
# plt.plot(test_data,np.array(test_data)**3/6)
# plt.scatter(np.linspace(0,la,len(integral_test_1)),integral_test_1)
# plt.plot(test_data,np.array(test_data)**4/24)
# plt.scatter(np.linspace(0,la,len(integral_test_2)),integral_test_2)
# plt.plot(test_data,np.array(test_data)**5/120)
# plt.scatter(np.linspace(0,la,len(integral_test_3)),integral_test_3)
# plt.show()

# print("Simulation.py took", time.time() - start_time, "to run")


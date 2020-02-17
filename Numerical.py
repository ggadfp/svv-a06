# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:03:57 2020

@author: halms
"""

import numpy as np
import matplotlib.pyplot as plt
import readAeroload
from mpl_toolkits.mplot3d import Axes3D
#%%

# Number of nodes
Nx = 41 # Spanwise Nodes
Nz = 81 # Chordwise Nodes

# Geometry of Aileron in [m]
Ca = 0.484  # Chord length
la = 1.691 # Span
x1 = 0.149 # x-location of hinge 1
x2 = 0.554 # x-location of hinge 2
x3 = 1.541 # x-location of hinge 3
xa = 27.2e-2 # Distance between actuator 1 and 2
h = 17.3e-2 # Aileron height
tsk = 1.1e-3 # Thickness of the skin
tsp = 2.5e-3 # Thickness of the spar
lsk = np.sqrt((Ca-h/2)**2 + (h/2)**2) # Length of the inclined skin
Ask_incl = lsk*tsk # Area of inclined section
Ask_semi = np.pi*tsk*h/2 # Area of semi-circular section
Asp = h*tsp # Area of spar

peri = lsk*2 + np.pi * h/2 # Perimeter of aileron

# Geometry of Stiffener in [m]
tst = 1.2e-3 # Thickness of the stiffener
hst = 1.4e-2 # Height of the stiffener
wst = 1.8e-2 # Width of the stiffener
nst = 13 # [-] -- Number of stiffeners
Ast = (hst + wst)*tst # Stiffener area
deltast = peri/nst # Stiffener spacing

# Given conditions
d1 = 0.681e-2 # [m] -- Vertical displacement of hinge 1
d3 = 2.03e-2 # [m] -- Vertical displacement of hinge 3
theta_deg = 26 # [deg] -- Maximum upward deflection
theta_rad = theta_deg*np.pi/180 # [rad] -- Maximum upward deflection
P = 37.9e+3 # [N] -- Load acting on actuator 2

# Necessary angles
angle_rad = deltast/(h/2) # Angle between boom at symmetry axis and the next boom
angle_deg = angle_rad *180/np.pi # Same angle, but in degrees
beta_rad = np.arctan2((h/2),(Ca-h/2)) # Angle of the inclined section
beta_deg = beta_rad*180/np.pi # Same angle, but in degrees

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
    Iyy_circle = np.pi * r**3 * t_skin / 2 + A_circle * (4 * r / (3 * np.pi) - z_tilde)**2
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


x,z = meshCreator(Nz,Nx,Ca,la,h)
load,pointapp = zload_zpoint(Nz,Nx,z)

torque = np.zeros((Nx,1))
for i in range(Nx):
    torque[i] = load[i]*pointapp[i]

def lagrange_basis(xval,j,xarray):
    L = 1
    for i in range(len(xarray)):
        if i != j:
            y = (xval - xarray[i])/(xarray[j] - xarray[i])
            L *=y
    return L

def polynomial(xval,xarray,yarray):
    poly = np.zeros((len(xarray),1))
    dot = 0
    for j in range(len(xarray)):
        poly[j] = lagrange_basis(xval,j,xarray)
    for k in range(len(xarray)):
        dot += poly[k]*yarray[k]
    return dot


#test = np.linspace(x[0,0],x[0,-1])
pr = np.array([polynomial(x[0,i],x[0,:],torque) for i in range(Nx)])
plt.figure(1)
ax = plt.axes(projection='3d')
ax.scatter(z,x,data)

#pr = pr/np.linalg.norm(pr)
plt.figure(2)
#plt.plot(pr,'o')
plt.plot(x[0,:],torque,'r')


plt.figure(3)
plt.plot(data[:,0])
#%%


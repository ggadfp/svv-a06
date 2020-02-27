# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 15:33:26 2020

@author: Ihab
"""

import numpy as np
import matplotlib.pyplot as plt

#%%
# Material properties
E = 73.1e9 # [Pa] -- E-modulus
G = 28e9 # [Pa] -- Shear modulus

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
    
# plt.scatter(-np.array(Z_locations),Y_locations)

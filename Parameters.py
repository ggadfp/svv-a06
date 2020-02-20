# -*- coding: utf-8 -*-
"""
Created on Sun Feb 16 13:03:57 2020

@author: Sabri
"""


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

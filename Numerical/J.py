# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 09:48:37 2020

@author: Ihab
"""


import numpy as np
import matplotlib.pyplot as plt
import Deflection as d
import Parameters as c

#%%

def torsion_shear(xloc):
    A1 = 1/2*np.pi*((c.h/2)**2)
    A2 = 1/2*c.h*(c.Ca-c.h/2)

    l2 = np.pi * c.h/2

    T = d.T_span(xloc)[0]
    c1 = 1/(2*A1*c.G)*(l2/c.tsk + c.h/c.tsp) + 1/(2*A2*c.G)*(c.h/c.tsp)
    c2 = -(1/(2*A1*c.G)*c.h/c.tsp + 1/(2*A2*c.G)*(2*c.lsk/c.tsk + c.h/c.tsp))
    c3 = -1/(2*A1*c.G)*(l2/c.tsk + c.h/c.tsp)
    c4 = 1/(2*A1*c.G)*c.h/c.tsp

    # c11 = 1/(2*A1)*(l2/tsk + l1/tsp)
    # c12 = -1/(2*A1)*l1/tsp
    # c21 = -l1/tsp*1/(2*A2)
    # c22 = 1/A2*l3/tsk + l1/tsp*1/(2*A2)


    mat = np.array([[2*A1, 2*A2, 0], [c1, c2, 0], [c3, c4, 1]])
    sol = np.array([[T],[0],[0]])

# mat = np.array([[c11, c12, -1], [c21, c22, -1],[2*A1, 2*A2, 0]])
# sol = np.array([[0],[0],[T]])

    ans = np.matmul(np.linalg.inv(mat), sol)
# ans = np.linalg.solve(mat,sol)

    J = T/(c.G*ans[2])
    
    return ans[0],ans[1],J

def torsion_const():
    return torsion_shear(0.5)[2]




#Calculation of the Torsional stiffness J
import numpy as np

# J = T/(G*(dalpha/dx))
G = 27.1e09 # [Pa] Shear modulus
tsk = 1.1e-3 # Thickness of the skin [m]
tsp = 2.5e-3 # Thickness of the spar [m]
Ca = 0.484  # Chord length [m]
h = 17.3e-2 # Aileron height [m]


A1 = 1/2*np.pi*((h/2)**2)
A2 = 1/2*h*(Ca-h/2)

l1 = h
l2 = np.pi * h/2
l3 = np.sqrt((Ca-h/2)**2 + (h/2)**2)

T = 1 #unit torque

c1 = 1/(2*A1*G)*(l2/tsk + l1/tsp) + 1/(2*A2*G)*(l1/tsp)
c2 = -(1/(2*A1*G)*l1/tsp + 1/(2*A2*G)*(2*l3/tsk + l1/tsp))
c3 = -1/(2*A1*G)*(l2/tsk + l1/tsp)
c4 = 1/(2*A1*G)*l1/tsp

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

J = T/(G*ans[2])
actual_j = 1.8e-4
print(J)

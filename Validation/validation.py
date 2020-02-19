# Validation
import numpy as np


# nodes
a = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
a = a.astype(np.float)

# output file
b = np.genfromtxt("B737.inp", dtype=str, skip_header=14146, skip_footer=(14594 - 14178), delimiter=",", comments="*")
b = b.astype(np.float)

# elements
e = np.genfromtxt("B737.inp", dtype=str, skip_header=6598, skip_footer=(14594 - 13233), delimiter=",")
e = e.astype(np.float)

f = open("B737.inp", "r")
lines = f.readlines()
f.close()

#print(a)
#print(e)

# plotting aileron

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

x = a[:,1]
y = a[:,2]
z = a[:,3]

plt.scatter(z,y)
#plt.show()

# nodes in elements
n_in_e = e[:,[1,2,3,4]]
#print(n_in_e)
empt = np.array([[]])
print(empt)
for i in range(len(n_in_e)):

    e_of_n = np.argwhere(n_in_e==i)
    np.vstack([empt, n_in_e])

print(empt)
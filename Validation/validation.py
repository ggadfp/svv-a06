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

# bending stresses
bendstr1 = np.genfromtxt("B737.rpt", dtype=str, skip_header=20, skip_footer=(59956 - 5964), delimiter="")
bendstr1 = bendstr1.astype(np.float)
bendstr2 = np.genfromtxt("B737.rpt", dtype=str, skip_header=5816, skip_footer=(59956 - 6673 - 158), delimiter="")
bendstr2 = bendstr2.astype(np.float)
print(bendstr2)

print(len(bendstr1))
print(bendstr1)

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

test = np.argwhere(n_in_e==0)[:,0]

print(test)
avvallis = []


for i in range(len(n_in_e)):
    #print(i)
    val = []
    e_of_n = np.argwhere(n_in_e==i)[:,0] #corresponding elements where a node appears in
    

    for j in range(len(e_of_n)):
        #print(j)
        ind2 = e_of_n[j]
        ind1 = e_of_n[j]
        #print(ind, j)
        if 1732 <= ind2 <= 2139:
            ind2 = ind2 - 1732
            value = bendstr2[:,2][ind2]
            ind1 = 0        
        elif 3835 <= ind2 <= 3946:
            ind2 = ind2 - 1732 - 1695
            value = bendstr2[:,2][ind2]
            ind1 = 0
        elif 5375 <= ind2 <= 5430:
            ind2 = ind2 - 1732 - 1695 - 1428
            value = bendstr2[:,2][ind2]
            ind1 = 0
        elif 6229 <= ind2 <= 6508:
           ind2 = ind2 - 1732 - 1695 - 1428 - 798
           value = bendstr2[:,2][ind2]
           ind1 = 0
        else:
            ind2 = 0
        
        if 0 <= ind1 <= 1731:
            ind1 = ind1 
            #print(ind1, j)
        elif 2140 <= ind1 <= 3834:
            ind1 = ind1 - 408 #112
            #print(ind1, j)
        elif 3947 <= ind1 <= 5374:
            ind1 = ind1 - 408 - 112 
            #print(ind1, j)
        elif 5431 <= ind1 <= 6228:
            ind1 = ind1 - 408 - 112 - 56 
            #print(ind1, j)
        elif 6509 <= ind1 <= 6634:
            ind1 = ind1 - 408 - 112 - 56 - 280 
            #print(ind1, j)
        
        
        value = bendstr1[:,2][ind1]
        
        val.append(value)
        
        
    
    #print("donzo")
    avval = np.mean(val)
    avvallis.append(avval)
    #print(val)
    #print(avval)
    #print(avvallis)



avvallis = np.delete(avvallis, 0)
avvallis = np.resize(avvallis, (len(avvallis) - 45))
    

print(avvallis)
print(len(avvallis))
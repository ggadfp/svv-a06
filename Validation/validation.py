# Validation
import numpy as np
import itertools

## ----------------IMPORTING FROM B737 INP and RPT-------------------

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

# jammed bending stresses
jambendstr1 = np.genfromtxt("B737.rpt", dtype=str, skip_header=6705, skip_footer=(59956 - 12483 - 147), delimiter="")
jambendstr1 = jambendstr1.astype(np.float)
jambendstr2 = np.genfromtxt("B737.rpt", dtype=str, skip_header=12501, skip_footer=(59956 - 13357 - 140), delimiter="")
jambendstr2 = jambendstr2.astype(np.float)

# jammed stresses
jamstr1 = np.genfromtxt("B737.rpt", dtype=str, skip_header=13390, skip_footer=(59956 - 19168 - 128), delimiter="")
jamstr1 = jamstr1.astype(np.float)
jamstr2 = np.genfromtxt("B737.rpt", dtype=str, skip_header=19186, skip_footer=(59956 - 20042 - 121), delimiter="")
jamstr2 = jamstr2.astype(np.float)

# bending displacements
bendu = np.genfromtxt("B737.rpt", dtype=str, skip_header=20074, skip_footer=(59956 - 26662 - 109), delimiter="")
bendu = bendu.astype(np.float)/1000

# jammed bending displacements
jambendu = np.genfromtxt("B737.rpt", dtype=str, skip_header=26724, skip_footer=(59956 - 33312 - 90), delimiter="")
jambendu = jambendu.astype(np.float)/1000

# jammed displacements
jamu = np.genfromtxt("B737.rpt", dtype=str, skip_header=33374, skip_footer=(59956 - 39962 - 71), delimiter="")
jamu = jamu.astype(np.float)/1000

# xyz coordinates of nodes
xyzn = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
xyzn = xyzn.astype(np.float)
xyzn[:,1] = xyzn[:,1]/1000
xyzn[:,2] = xyzn[:,2]/1000
xyzn[:,3] = xyzn[:,3]/1000


# reaction forces of assembly
reac = np.genfromtxt("B737.rpt", dtype=str, skip_header=59928, skip_footer=(5), delimiter="")
reac = reac.astype(np.float)

# nodes and xyz of assembly
assemb = []
for q in range(16):
    
    assem = np.genfromtxt("B737.inp", dtype=str, skip_header=(14146 + 2*q), skip_footer=(14594 - 14148 - 2*q), delimiter=",")
    assem = assem.astype(np.float)
    assem = list(assem)
    assemb.append(assem) #ff fixen

assemb = np.array(assemb)


## ------------------------------Extrapolating stresses from the data and putting them in coherent lists sorted in node number--------------------------------

# nodes in elements
n_in_e = e[:,[1,2,3,4]]
#print(n_in_e)


# von mises loc1
m1avvallis_b = []
m1avvallis_jb = []
m1avvallis_j = []

# von mises loc2
m2avvallis_b = []
m2avvallis_jb = []
m2avvallis_j = []

# S12 loc1
s1avvallis_b = []
s1avvallis_jb = []
s1avvallis_j = []

# S12 loc2
s2avvallis_b = []
s2avvallis_jb = []
s2avvallis_j = []


for i in range(len(n_in_e)):
    
    
    # Von Mises loc1
    m1val_b = []
    m1val_jb = []
    m1val_j = []

    # Von Mises loc2
    m2val_b = []
    m2val_jb = []
    m2val_j = []

    # S12 loc1
    s1val_b = []
    s1val_jb = []
    s1val_j = []

    # S12 loc2
    s2val_b = []
    s2val_jb = []
    s2val_j = []


    e_of_n = np.argwhere(n_in_e==i)[:,0] #corresponding elements where a node appears in
    #print(np.argwhere(n_in_e==i))
    

    for j in range(len(e_of_n)):
        #print(j)
        ind2 = e_of_n[j]
        ind1 = e_of_n[j]
        #print(ind, j)
        if 1732 <= ind2 <= 2139:
            ind2 = ind2 - 1732
            
            # Von Mises loc1
            m1value_b = bendstr2[:,2][ind2]
            m1value_jb = jambendstr2[:,2][ind2]
            m1value_j = jamstr2[:,2][ind2]

            # Von Mises loc2
            m2value_b = bendstr2[:,3][ind2]
            m2value_jb = jambendstr2[:,3][ind2]
            m2value_j = jamstr2[:,3][ind2]

            # S12 loc1
            s1value_b = bendstr2[:,4][ind2]
            s1value_jb = jambendstr2[:,4][ind2]
            s1value_j = jamstr2[:,4][ind2]

            # S12 loc2
            s2value_b = bendstr2[:,5][ind2]
            s2value_jb = jambendstr2[:,5][ind2]
            s2value_j = jamstr2[:,5][ind2]


            ind1 = 0        
        elif 3835 <= ind2 <= 3946:
            ind2 = ind2 - 1732 - 1695
            
            # Von Mises loc1
            value_b = bendstr2[:,2][ind2]
            value_jb = jambendstr2[:,2][ind2]
            value_j = jamstr2[:,2][ind2]

            # Von Mises loc2
            m2value_b = bendstr2[:,3][ind2]
            m2value_jb = jambendstr2[:,3][ind2]
            m2value_j = jamstr2[:,3][ind2]

            # S12 loc1
            s1value_b = bendstr2[:,4][ind2]
            s1value_jb = jambendstr2[:,4][ind2]
            s1value_j = jamstr2[:,4][ind2]

            # S12 loc2
            s2value_b = bendstr2[:,5][ind2]
            s2value_jb = jambendstr2[:,5][ind2]
            s2value_j = jamstr2[:,5][ind2]

            ind1 = 0
        elif 5375 <= ind2 <= 5430:
            ind2 = ind2 - 1732 - 1695 - 1428
            
            # Von Mises loc1
            value_b = bendstr2[:,2][ind2]
            value_jb = jambendstr2[:,2][ind2]
            value_j = jamstr2[:,2][ind2]

            # Von Mises loc2
            m2value_b = bendstr2[:,3][ind2]
            m2value_jb = jambendstr2[:,3][ind2]
            m2value_j = jamstr2[:,3][ind2]

            # S12 loc1
            s1value_b = bendstr2[:,4][ind2]
            s1value_jb = jambendstr2[:,4][ind2]
            s1value_j = jamstr2[:,4][ind2]

            # S12 loc2
            s2value_b = bendstr2[:,5][ind2]
            s2value_jb = jambendstr2[:,5][ind2]
            s2value_j = jamstr2[:,5][ind2]

            ind1 = 0
        elif 6229 <= ind2 <= 6508:
            ind2 = ind2 - 1732 - 1695 - 1428 - 798
            # Von Mises loc1
            value_b = bendstr2[:,2][ind2]
            value_jb = jambendstr2[:,2][ind2]
            value_j = jamstr2[:,2][ind2]

            # Von Mises loc2
            m2value_b = bendstr2[:,3][ind2]
            m2value_jb = jambendstr2[:,3][ind2]
            m2value_j = jamstr2[:,3][ind2]

            # S12 loc1
            s1value_b = bendstr2[:,4][ind2]
            s1value_jb = jambendstr2[:,4][ind2]
            s1value_j = jamstr2[:,4][ind2]

            # S12 loc2
            s2value_b = bendstr2[:,5][ind2]
            s2value_jb = jambendstr2[:,5][ind2]
            s2value_j = jamstr2[:,5][ind2]

            ind1 = 0
        else:
            ind2 = 0
        
        if 0 <= ind1 <= 1731:
            ind1 = ind1 
            
            # Von Mises loc1
            m1value_b = bendstr1[:,2][ind1]
            m1value_jb = jambendstr1[:,2][ind1]
            m1value_j = jamstr1[:,2][ind1]

            # Von Mises loc2
            m2value_b = bendstr1[:,3][ind1]
            m2value_jb = jambendstr1[:,3][ind1]
            m2value_j = jamstr1[:,3][ind1]

            # S12 loc1
            s1value_b = bendstr1[:,4][ind1]
            s1value_jb = jambendstr1[:,4][ind1]
            s1value_j = jamstr1[:,4][ind1]

            # S12 loc2
            s2value_b = bendstr1[:,5][ind1]
            s2value_jb = jambendstr1[:,5][ind1]
            s2value_j = jamstr1[:,5][ind1]

            #print(ind1, j)
        elif 2140 <= ind1 <= 3834:
            ind1 = ind1 - 408 #112
            
            # Von Mises loc1
            m1value_b = bendstr1[:,2][ind1]
            m1value_jb = jambendstr1[:,2][ind1]
            m1value_j = jamstr1[:,2][ind1]

            # Von Mises loc2
            m2value_b = bendstr1[:,3][ind1]
            m2value_jb = jambendstr1[:,3][ind1]
            m2value_j = jamstr1[:,3][ind1]

            # S12 loc1
            s1value_b = bendstr1[:,4][ind1]
            s1value_jb = jambendstr1[:,4][ind1]
            s1value_j = jamstr1[:,4][ind1]

            # S12 loc2
            s2value_b = bendstr1[:,5][ind1]
            s2value_jb = jambendstr1[:,5][ind1]
            s2value_j = jamstr1[:,5][ind1]

            #print(ind1, j)
        elif 3947 <= ind1 <= 5374:
            ind1 = ind1 - 408 - 112 
            
            # Von Mises loc1
            m1value_b = bendstr1[:,2][ind1]
            m1value_jb = jambendstr1[:,2][ind1]
            m1value_j = jamstr1[:,2][ind1]

            # Von Mises loc2
            m2value_b = bendstr1[:,3][ind1]
            m2value_jb = jambendstr1[:,3][ind1]
            m2value_j = jamstr1[:,3][ind1]

            # S12 loc1
            s1value_b = bendstr1[:,4][ind1]
            s1value_jb = jambendstr1[:,4][ind1]
            s1value_j = jamstr1[:,4][ind1]

            # S12 loc2
            s2value_b = bendstr1[:,5][ind1]
            s2value_jb = jambendstr1[:,5][ind1]
            s2value_j = jamstr1[:,5][ind1]

            #print(ind1, j)
        elif 5431 <= ind1 <= 6228:
            ind1 = ind1 - 408 - 112 - 56 
            
            # Von Mises loc1
            m1value_b = bendstr1[:,2][ind1]
            m1value_jb = jambendstr1[:,2][ind1]
            m1value_j = jamstr1[:,2][ind1]

            # Von Mises loc2
            m2value_b = bendstr1[:,3][ind1]
            m2value_jb = jambendstr1[:,3][ind1]
            m2value_j = jamstr1[:,3][ind1]

            # S12 loc1
            s1value_b = bendstr1[:,4][ind1]
            s1value_jb = jambendstr1[:,4][ind1]
            s1value_j = jamstr1[:,4][ind1]

            # S12 loc2
            s2value_b = bendstr1[:,5][ind1]
            s2value_jb = jambendstr1[:,5][ind1]
            s2value_j = jamstr1[:,5][ind1]

            #print(ind1, j)
        elif 6509 <= ind1 <= 6634:
            ind1 = ind1 - 408 - 112 - 56 - 280 
            
            # Von Mises loc1
            m1value_b = bendstr1[:,2][ind1]
            m1value_jb = jambendstr1[:,2][ind1]
            m1value_j = jamstr1[:,2][ind1]

            # Von Mises loc2
            m2value_b = bendstr1[:,3][ind1]
            m2value_jb = jambendstr1[:,3][ind1]
            m2value_j = jamstr1[:,3][ind1]

            # S12 loc1
            s1value_b = bendstr1[:,4][ind1]
            s1value_jb = jambendstr1[:,4][ind1]
            s1value_j = jamstr1[:,4][ind1]

            # S12 loc2
            s2value_b = bendstr1[:,5][ind1]
            s2value_jb = jambendstr1[:,5][ind1]
            s2value_j = jamstr1[:,5][ind1]

            #print(ind1, j)
        else:
            ind1 = 0
        
        
        
        # Von Mises loc1
        m1val_b.append(m1value_b)
        m1val_jb.append(m1value_jb)
        m1val_j.append(m1value_j)

        # Von Mises loc2
        m2val_b.append(m2value_b)
        m2val_jb.append(m2value_jb)
        m2val_j.append(m2value_j)

        # S12 loc1
        s1val_b.append(s1value_b)
        s1val_jb.append(s1value_jb)
        s1val_j.append(s1value_j)

        # S12 loc2
        s2val_b.append(s2value_b)
        s2val_jb.append(s2value_jb)
        s2val_j.append(s2value_j)
        
        
    
    
    ## Taking the mean of the stresses of the elements where the node appears in and take that value as the stresses in the node
    # Von Mises loc1
    m1avval_b = np.mean(m1val_b)
    m1avvallis_b.append(m1avval_b)
    m1avval_jb = np.mean(m1val_jb)
    m1avvallis_jb.append(m1avval_jb)
    m1avval_j = np.mean(m1val_j)
    m1avvallis_j.append(m1avval_j)

    # Von Mises loc2
    m2avval_b = np.mean(m2val_b)
    m2avvallis_b.append(m2avval_b)
    m2avval_jb = np.mean(m2val_jb)
    m2avvallis_jb.append(m2avval_jb)
    m2avval_j = np.mean(m2val_j)
    m2avvallis_j.append(m2avval_j)

    # S12 loc1
    s1avval_b = np.mean(s1val_b)
    s1avvallis_b.append(s1avval_b)
    s1avval_jb = np.mean(s1val_jb)
    s1avvallis_jb.append(s1avval_jb)
    s1avval_j = np.mean(s1val_j)
    s1avvallis_j.append(s1avval_j)

    # S12 loc2
    s2avval_b = np.mean(s2val_b)
    s2avvallis_b.append(s2avval_b)
    s2avval_jb = np.mean(s2val_jb)
    s2avvallis_jb.append(s2avval_jb)
    s2avval_j = np.mean(s2val_j)
    s2avvallis_j.append(s2avval_j)

    


# Von Mises loc1
m1avvallis_b = np.delete(m1avvallis_b, 0)
m1avvallis_b = np.resize(m1avvallis_b, (len(m1avvallis_b) - 45))
m1avvallis_jb = np.delete(m1avvallis_jb, 0)
m1avvallis_jb = np.resize(m1avvallis_jb, (len(m1avvallis_jb) - 45))
m1avvallis_j = np.delete(m1avvallis_j, 0)
m1avvallis_j = np.resize(m1avvallis_j, (len(m1avvallis_j) - 45))

# Von Mises loc2
m2avvallis_b = np.delete(m2avvallis_b, 0)
m2avvallis_b = np.resize(m2avvallis_b, (len(m2avvallis_b) - 45))
m2avvallis_jb = np.delete(m2avvallis_jb, 0)
m2avvallis_jb = np.resize(m2avvallis_jb, (len(m2avvallis_jb) - 45))
m2avvallis_j = np.delete(m2avvallis_j, 0)
m2avvallis_j = np.resize(m2avvallis_j, (len(m2avvallis_j) - 45))

# S12 loc1
s1avvallis_b = np.delete(s1avvallis_b, 0)
s1avvallis_b = np.resize(s1avvallis_b, (len(s1avvallis_b) - 45))
s1avvallis_jb = np.delete(s1avvallis_jb, 0)
s1avvallis_jb = np.resize(s1avvallis_jb, (len(s1avvallis_jb) - 45))
s1avvallis_j = np.delete(s1avvallis_j, 0)
s1avvallis_j = np.resize(s1avvallis_j, (len(s1avvallis_j) - 45))

# S12 loc2
s2avvallis_b = np.delete(s2avvallis_b, 0)
s2avvallis_b = np.resize(s2avvallis_b, (len(s2avvallis_b) - 45))
s2avvallis_jb = np.delete(s2avvallis_jb, 0)
s2avvallis_jb = np.resize(s2avvallis_jb, (len(s2avvallis_jb) - 45))
s2avvallis_j = np.delete(s2avvallis_j, 0)
s2avvallis_j = np.resize(s2avvallis_j, (len(s2avvallis_j) - 45))


# Von Mises loc1
m1xyzn_b = np.c_[xyzn,m1avvallis_b]
m1xyzn_jb = np.c_[xyzn,m1avvallis_jb]
m1xyzn_j = np.c_[xyzn,m1avvallis_j]

# Von Mises loc2
m2xyzn_b = np.c_[xyzn,m2avvallis_b]
m2xyzn_jb = np.c_[xyzn,m2avvallis_jb]
m2xyzn_j = np.c_[xyzn,m2avvallis_j]


# S12 loc1
s1xyzn_b = np.c_[xyzn,s1avvallis_b]
s1xyzn_jb = np.c_[xyzn,s1avvallis_jb]
s1xyzn_j = np.c_[xyzn,s1avvallis_j]

# S12 loc2
s2xyzn_b = np.c_[xyzn,s2avvallis_b]
s2xyzn_jb = np.c_[xyzn,s2avvallis_jb]
s2xyzn_j = np.c_[xyzn,s2avvallis_j]



## ---------------------------Extrapolating displacements of the hingeline from the data and putting them in coherent lists sorted in x-----------------

# extracting all points on the hinge line from xyzn

nn = xyzn[:,0]
xn = xyzn[:,1]
yn = xyzn[:,2]
zn = xyzn[:,3]
nhinge = []
xhinge = []


for u in range(len(nn)):
    
    if yn[u] == 0 and zn[u] == 0:
        nhinge.append(nn[u])
        xhinge.append(xn[u])



# getting u1-3 for hingeline

bu1list = []
bu2list = []
bu3list = []
jbu1list = []
jbu2list = []
jbu3list = []
ju1list = []
ju2list = []
ju3list = []

bendu1 = bendu[:, 2]
bendu2 = bendu[:, 3]
bendu3 = bendu[:, 4]
jambendu1 = jambendu[:, 2]
jambendu2 = jambendu[:, 3]
jambendu3 = jambendu[:, 4]
jamu1 = jamu[:, 2]
jamu2 = jamu[:, 3]
jamu3 = jamu[:, 4]


for q in range(len(nhinge)):

    node = int(nhinge[q] - 1)
    bu1list.append(bendu1[node])
    bu2list.append(bendu2[node])
    bu3list.append(bendu3[node])
    jbu1list.append(jambendu1[node])
    jbu2list.append(jambendu2[node])
    jbu3list.append(jambendu3[node])
    ju1list.append(jamu1[node])
    ju2list.append(jamu2[node])
    ju3list.append(jamu3[node])
    

nxn = np.vstack([nhinge,xhinge])


#1

bnxu1 = np.vstack([nxn,bu1list])
bnxu12 = np.vstack([bnxu1,bu2list])
bnxu123 = np.vstack([bnxu12,bu3list])
bnxu123 = np.transpose(bnxu123)

bnxu123 = bnxu123[np.argsort(bnxu123[:, 1])]

#2

jbnxu1 = np.vstack([nxn,jbu1list])
jbnxu12 = np.vstack([jbnxu1,jbu2list])
jbnxu123 = np.vstack([jbnxu12,jbu3list])
jbnxu123 = np.transpose(jbnxu123)

jbnxu123 = jbnxu123[np.argsort(jbnxu123[:, 1])]

#3

jnxu1 = np.vstack([nxn,ju1list])
jnxu12 = np.vstack([jnxu1,ju2list])
jnxu123 = np.vstack([jnxu12,ju3list])
jnxu123 = np.transpose(jnxu123)

jnxu123 = jnxu123[np.argsort(jnxu123[:, 1])]

#------------------Displacements validation---------
import Deflection as defl
import Simulation as sim
# xloc = np.array([jbnxu123[:,1]])
# xloc = np.transpose(xloc)

Usim2l = []
Usim3l = []

for w in range(len(jbnxu123[:,1])):
    xloc = jbnxu123[w,1]
    Usim2 = defl.deflection_y(xloc)
    Usim3 = defl.deflection_z(xloc)
    Usim2 = float(Usim2)
    Usim3 = float(Usim3)
    Usim2l.append(Usim2)
    Usim3l.append(Usim3)
    

#print(Usim2l)
#print(Usim3l)
#print(jbnxu123[:,3])
#print(jbnxu123[:,4])

Usim2l = np.array([Usim2l])
Usim3l = np.array([Usim3l])
Uval2 = np.array([jbnxu123[:,3]])
Uval3 = np.array([jbnxu123[:,4]])

Usim2l = np.transpose(Usim2l)
Usim3l = np.transpose(Usim3l)
Uval2 = np.transpose(Uval2)
Uval3 = np.transpose(Uval3)

diff2 = Usim2l - Uval2
diff3 = Usim3l - Uval3
fac2 = Usim2l / Uval2
fac3 = Usim3l / Uval3

rmse2 = np.sqrt((np.square(Usim2l-Uval2)).mean())
rmse3 = np.sqrt((np.square(Usim3l-Uval3)).mean())

print(rmse2)
print(rmse3)

## -------------------Plotting displacements against x location--------------------

import matplotlib.pyplot as plt

axis = np.linspace(0,sim.la,defl.plot_points)

plt.plot(jbnxu123[:,1], jbnxu123[:,2])
plt.show()

plt.plot(jbnxu123[:,1], jbnxu123[:,3], 'r', label='Validation')
plt.scatter(axis, defl.deflectiony_dist,label='A06 Simulation')
#plt.plot(verif.x, verif.defy_v, 'tab:red',label='Verification code (N=20)')
plt.title(r'$\nu(x)$')
plt.ylabel('[m]')
plt.grid()
plt.legend(loc="lower right")
plt.xlabel('Spanwise location [m]')
plt.show()

plt.plot(jbnxu123[:,1], jbnxu123[:,4], 'r', label='Validation')
plt.scatter(axis, defl.deflectionz_dist,label='A06 Simulation')
#plt.plot(verif.x, verif.defz_v, 'tab:red',label='Verification code (N=20)')
plt.title(r'$\eta(x)$')
plt.ylabel('[m]')
plt.grid()
plt.legend(loc="upper right")
plt.xlabel('Spanwise location [m]')
plt.show()


#------stresses validation--------------
from Stress import *

rmsc = []
rmscr = []
rmsi1 = []
rmsi2 = []
rmssp = []
rmsspr = []

for t in range(len(jbnxu123[:,1])):
    x_cind = np.argwhere(xyzn[:,1] == jbnxu123[t,1])[:,0]
    #print(np.argwhere(xyzn[:,1] == jbnxu123[t,1])[:,0])
    xlo = float(jbnxu123[t,1])
    vm_circle,vm_circle_rev,vm_incl_1,vm_incl_2,vm_spar,vm_spar_rev = Von_Mises(xlo)

    z1 = []
    z5 = []
    z2 = []
    z3 = []
    z4 = []

    vm_val_circ = []
    vm_val_circ_rev = []
    vm_val_incl_1 = []
    vm_val_incl_2 = []
    vm_val_spar = []
    vm_val_spar_rev = []

    for d in range(len(x_cind)):

        if xyzn[x_cind[d],3] < 0 and xyzn[x_cind[d],2] >= 0:
            z1e = float(xyzn[d,3])
            z1.append(z1e)
            vm_val_circ_rev_e = float(m1xyzn_jb[d,4])
            vm_val_circ_rev.append(vm_val_circ_rev_e)
                
        elif xyzn[x_cind[d],3] < 0 and xyzn[x_cind[d],2] <= 0:
            z5e = float(xyzn[d,3])
            z5.append(z5e)
            vm_val_circ_e = float(m1xyzn_jb[d,4])
            vm_val_circ.append(vm_val_circ_e)
                
        elif xyzn[x_cind[d],3] > 0 and xyzn[x_cind[d],2] <= 0:
            z2e = float(xyzn[d,3])
            z2.append(z2e)
            vm_val_incl_1_e = float(m1xyzn_jb[d,4])
            vm_val_incl_1.append(vm_val_incl_1_e)
                
        elif xyzn[x_cind[d],3] > 0 and xyzn[x_cind[d],2] >= 0:
            z3e = float(xyzn[d,3])
            z3.append(z3e)
            vm_val_incl_2_e = float(m1xyzn_jb[d,4])
            vm_val_incl_2.append(vm_val_incl_2_e)
                
        elif xyzn[x_cind[d],3] == 0 and xyzn[x_cind[d],2] <= 0:
            z4e = float(xyzn[d,3])
            z4.append(z4e)
            vm_val_spar_e = float(m1xyzn_jb[d,4])
            vm_val_spar.append(vm_val_spar_e)

        elif xyzn[x_cind[d],3] == 0 and xyzn[x_cind[d],2] <= 0:
            z4e = float(xyzn[d,3])
            z4.append(z4e)
            vm_val_spar_rev_e = float(m1xyzn_jb[d,4])
            vm_val_spar_rev.append(vm_val_spar_rev_e)
                
        else:
            print('bad coding', d)

    vm_val_circ = np.array([vm_val_circ])
    vm_val_circ_rev = np.array([vm_val_circ_rev])
    vm_val_incl_1 = np.array([vm_val_incl_1])
    vm_val_incl_2 = np.array([vm_val_incl_2])
    vm_val_spar = np.array([vm_val_spar])
    vm_val_spar_rev = np.array([vm_val_spar_rev])

    rmsecirc = np.sqrt((np.square(vm_val_circ-vm_circle)).mean())
    rmsecircr = np.sqrt((np.square(vm_val_circ_rev-vm_circle_rev)).mean())
    rmsein1 = np.sqrt((np.square(vm_val_incl_1-vm_incl_1)).mean())
    rmsein2 = np.sqrt((np.square(vm_val_incl_2-vm_incl_2)).mean())
    rmsesp = np.sqrt((np.square(vm_val_spar-vm_spar)).mean())
    rmsespr = np.sqrt((np.square(vm_val_spar_rev-vm_spar_rev)).mean())
        
    rmsc.append(rmsecirc)
    rmscr.append(rmsecircr)
    rmsi1.append(rmsein1)
    rmsi2.append(rmsein2)
    rmssp.append(rmsesp)
    rmsspr.append(rmsespr)

print(rmsc)
    

    
    
    

    #print(z1)
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
bendu = bendu.astype(np.float)

# jammed bending displacements
jambendu = np.genfromtxt("B737.rpt", dtype=str, skip_header=26724, skip_footer=(59956 - 33312 - 90), delimiter="")
jambendu = jambendu.astype(np.float)

# jammed displacements
jamu = np.genfromtxt("B737.rpt", dtype=str, skip_header=33374, skip_footer=(59956 - 39962 - 71), delimiter="")
jamu = jamu.astype(np.float)

# xyz coordinates of nodes
xyzn = np.genfromtxt("B737.inp", dtype=str, skip_header=9, skip_footer=(14594 - 6598), delimiter=",")
xyzn = xyzn.astype(np.float)

#print(xyzn)


f = open("B737.inp", "r")
lines = f.readlines()
f.close()



# # plotting aileron

# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# x = a[:,1]
# y = a[:,2]
# z = a[:,3]

# plt.scatter(z,y)
# #plt.show()

#### stresses

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
    #print(i)
    
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

#print(np.mean(avvallis_b))
#print(len(avvallis_b))


### displacements

#extracting all points on the hinge line from xyzn

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

#print(nhinge)
#print(len(nhinge))

#getting u1-3 for hingeline

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
    
# u1list = np.array(u1list)
# u2list = np.array(u2list)
# u3list = np.array(u3list)

#print(u1list)



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

# print(bnxu123)
# print(jbnxu123)
# print(jnxu123)
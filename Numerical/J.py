#Calculation of the Torsional stiffness J

# J = T/(G*(dalpha/dx))
G = 27.1e09 # [Pa] Shear modulus
tsk = 1.1e-3 # Thickness of the skin
tst = 1.2e-3 # Thickness of the stiffener

def integral(): #cross sectional integral of (qds)/(Gt)
    q = ...


dalpha_dx = 1/(2*Ar)*integral()

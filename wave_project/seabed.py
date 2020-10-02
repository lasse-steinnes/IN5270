import numpy as np
# X, Y = np.meshgrid(x,y, indexing = 'ij') <- Only for plotting!

# H0 and z is height of stillwater
# B is how high the ground level is
# g er hastighet

def g_floor(opt,x,y, Lx, Ly): # must be a n,n mesh Only make this once for each plot to save time
    """
    Calculates the shape of sea
    """
    B0 = 1; Ba = B0*2; Bs = Lx*0.3; Bmx = 0.5*Lx; Bmy = 0.5*Ly; b = 1.1;
    H0 = 3*B0; I0 = 0; Ia = 6*B0; g_a = 9.81;
    z = np.full((len(x),len(y)),H0)
    # circular
    if opt == 1:
        B = np.full((len(x),len(y)),I0)
        for i in range(len(x)):
            for j in range(len(y)):
                B[i,j] = B[i,j] + Ia*np.exp(-((x[i]-Bmx)/Bs)**2 -(y[j]-Bmy/(b*Bs))**2)

    # less smooth
    elif opt == 2:
        B = np.full((len(x),len(y)),B0)
        for i in range(len(x)):
            for j in range(len(y)):
                if 0 <= np.sqrt((x[i]-Bmx)**2 + (y[j]-Bmy)**2) <= Bs:
                    B[i,j] = B[i,j] + Ba*np.cos((np.pi*x[i]-Bmx)/(2*Bs))*np.cos((np.pi*y[j]-Bmy)/(2*Bs))
    # rectangle
    elif opt == 3:
        B = np.full((len(x),len(y)),B0)
        for i in range(len(x)):
            for j in range(len(y)):
                if Bmx - Bs <= x[i] <= Bmx + Bs  and Bmy - b*Bs <= y[j] <= Bmy + b*Bs:
                    B[i,j] = B[i,j] + Ba
    return z,B,g_a*(H0-B)

# gravitational constant m/s^2
## Calculate velocities at each point

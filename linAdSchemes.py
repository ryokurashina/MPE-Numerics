"""""
Contains various numerical schemes for the linear advection equation.

List of numerical schemes:
* FTCS
* CTCS
* BTCS
"""""

# Import packages and libraries
import numpy as np

# FTCS scheme with periodic boundary conditions
def FTCS(phi, ntime, c):
    N = len(phi)
    # Initialise for loop
    phi_old = np.array(phi)
    phi_new = np.zeros(N)
    for i in range(ntime):
        phi_new[1:N-1] = phi_old[1:N-1]-c/2*(phi_old[2:N]-phi_old[0:N-2])
        # Periodic B.C.
        phi_new[0] = phi_old[0]-c/2*(phi_old[1]-phi_old[N-2])
        phi_new[N-1] = phi_new[0]
        phi_old = phi_new.copy()
    return phi_new

def CTCS(phi, ntime, c):
    N = len(phi)
    # Initialise for loop
    phi_old = FTCS(phi, 1, c)
    phi_v_old = np.array(phi)
    phi_new = np.zeros(N)
    for i in range(ntime-1):
        ##### THIS SHOULD BE C NOT C/2... #####
        phi_new[1:N-1] = phi_v_old[1:N-1]-c/2*(phi_old[2:N]-phi_old[0:N-2])
        # Periodic B.C.
        phi_new[0] = phi_v_old[0]-c/2*(phi_old[1]-phi_old[N-2])
        phi_new[N-1] = phi_new[0]
        phi_old = phi_new.copy()
        phi_v_old = phi_old.copy()
    return phi_new

def BTCS(phi, ntime, c):
    N = len(phi)
    # Initialise for loop
    phi_old = np.array(phi)
    phi_new = np.zeros(N)
    # Construct matrix to invert
    A = np.zeros((N,N))
    stencil = np.array([-c/2, 1, c/2])
    for i in range(1,N-1):
        A[i,i-1:i+2] = stencil
    A[0, 0] = 1
    A[0, 1] = c/2
    A[0, N-1] = -c/2
    A[N-1, N-1] = 1
    A[N-1, 0] = c/2
    A[N-1, N-2] = -c/2
    for i in range(ntime):
        phi_new = np.linalg.solve(A, phi_old)
        phi_old = phi_new.copy()
    return phi_new

"""""
Contains various numerical schemes for the linearised 1D-SW equations with f = 0.

List of numerical schemes:
* Unstaggered FW (U_FW)
"""""

def UFW(phi, h, ntime, c):
    # Gravitational aceleration
    g = 9.81
    H = c**2/g
    # print("H = %f" %(H))
    N = len(phi)
    # Initialise for loop
    phi_old = np.array(phi)
    phi_new = np.zeros(N)
    h_old = np.array(h)
    h_new = np.zeros(N)
    for i in range(ntime):
        phi_new[1:N-1] = phi_old[1:N-1]-c/2*np.sqrt(g/H)*(h_old[2:N]-h_old[0:N-2])
        phi_new[0] = phi_old[0]-c/2*np.sqrt(g/H)*(h_old[1]-h_old[N-2])
        phi_new[N-1] = phi_new[0]
        h_new[1:N-1] = h_old[1:N-1]-c/2*np.sqrt(H/g)*(phi_new[2:N]-phi_new[0:N-2])
        h_new[0] = h_old[0]-c/2*np.sqrt(H/g)*(phi_new[1]-phi_new[N-2])
        h_new[N-1] = h_new[0]
        phi_old = phi_new.copy()
        h_old = h_new.copy()
    return phi_new, h

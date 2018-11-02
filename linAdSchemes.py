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
        phi_new = phi_old-.5*c*(np.roll(phi_old,1)-np.roll(phi_old,-1))
        phi_old = phi_new.copy()
    return phi_new

#### THIS DOESN'T WORK FOR SOME REASON #####
def CTCS(phi, ntime, c):
    N = len(phi)
    # Initialise for loop, do one FTCS so we have data for 2 time-steps
    phi_old = FTCS(phi, 1, c)
    phi_v_old = np.array(phi)
    phi_new = np.zeros(N)
    for i in range(ntime-1):
        phi_new = phi_v_old-c*(np.roll(phi_old,1)-np.roll(phi_old,-1))
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

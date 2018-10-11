"""""
Contains various numerical schemes for the linear advection equation.

List of numerical schemes:
* FTCS
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
        phi_new[0] = phi_old[0]-c/2*(phi_old[1]-phi_old[N-1])
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
        phi_new[1:N-1] = phi_v_old[1:N-1]-c/2*(phi_old[2:N]-phi_old[0:N-2])
        # Periodic B.C.
        phi_new[0] = phi_v_old[0]-c/2*(phi_old[1]-phi_old[N-1])
        phi_new[N-1] = phi_new[0]
        phi_old = phi_new.copy()
        phi_v_old = phi_old.copy()
    return phi_new

"""""
Script for post processing of data
"""""

# Import packages and libraries
import numpy as np
import matplotlib.pyplot as plt
from linAdSchemes import *
from initialConditions import *

def main():
    # Set resolution
    N = 1001
    x = np.linspace(0, 1, N)
    a = .25
    b = .5
    # Initial condition
    phi_init = squareWave(x, a, b)
    # Initialise timestep, grid-size and advection velocity
    dx = 1/(N-1)
    dt = 1e-4
    n_steps = 1000
    u = .1
    T = n_steps*dt
    print((a+u*T) % 1)
    print((b+u*T) % 1)
    phi_exact = squareWave(x, (a+u*T) % 1, (b+u*T) % 1)
    # Courant number
    c = u*dt/dx
    print("Courant Number (c): %f" %(c))

    """""
    Pick the numerical scheme here:
    * FTCS (Linear Advection)
    * CTCS (Linear Advection)
    """""

    # FTCS
    phi1 = FTCS(phi_init,n_steps,c)
    # CTCS
    phi2 = CTCS(phi_init,n_steps,c)
    print(len(phi1))
    print(len(phi2))
    plt.plot(x,phi_init)
    plt.plot(x,phi1)
    # plt.plot(x,phi2)
    plt.plot(x,phi_exact)
    plt.show()
main()

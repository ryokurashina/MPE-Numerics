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
    b = .75
    # Initial condition
    # phi_init = squareWave(x, a, b)
    phi_init = cosBell(x, alpha=a, beta=b)
    h_init = cosBell(x, alpha=.25, beta=.75)
    # Initialise timestep, grid-size and advection velocity
    dx = 1/(N-1)
    dt = 1e-3
    n_steps = 100
    u = .1
    T = n_steps*dt
    print((a+u*T) % 1)
    print((b+u*T) % 1)
    # phi_exact = squareWave(x, (a+u*T) % 1, (b+u*T) % 1)
    phi_exact = cosBell(x, alpha = (a+u*T) % 1, beta = (b+u*T) % 1)
    # Courant number
    c = u*dt/dx
    print("Courant Number (c): %f" %(c))

    """""
    Pick the numerical scheme here:
    * FTCS (Linear Advection)
    * CTCS (Linear Advection)
    """""

    # FTCS
    # phi = FTCS(phi_init, n_steps, c)
    # CTCS
    phi = CTCS(phi_init, n_steps, c)
    # BTCS
    # phi = BTCS(phi_init, n_steps, c)

    # SW
    # phi, h = UFW(phi_init, h_init, 100, 0.5)
    plt.plot(x,phi_init)
    plt.plot(x,phi)
    # plt.plot(x,h)

    plt.plot(x,phi_exact)
    plt.show()

main()

"""""
Main code to run numerical schemes and plot results
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
    n_steps = 250
    u = .2
    T = n_steps*dt
    print((a+u*T) % 1)
    print((b+u*T) % 1)
    # phi_exact = squareWave(x, (a+u*T) % 1, (b+u*T) % 1)
    phi_exact = cosBell(x, alpha = (a-u*T) % 1, beta = (b-u*T) % 1)
    # Courant number
    c = u*dt/dx
    print("Courant Number (c): %f" %(c))

    """""
    Pick the numerical scheme here:
    * FTCS (Linear Advection)
    * CTCS (Linear Advection)
    * BTCS (Linear Advection)
    * Unstaggered SW
    """""

    # FTCS
    phi1 = FTCS(phi_init, n_steps, c)

    # CTCS
    # phi2 = CTCS(phi_init, n_steps, c)

    # BTCS
    # phi3 = BTCS(phi_init, n_steps, c)

    # SW
    # phi, h = USW(phi_init, h_init, 100, 0.5)

    # Plot results
    plt.plot(x, phi_init, label="Initial")
    plt.plot(x, phi1, label="FTCS")
    # plt.plot(x, phi2, label="CTCS")
    # plt.plot(x,h)
    # plt.plot(x,phi_exact, label="Exact")
    plt.legend()
    plt.show()

main()

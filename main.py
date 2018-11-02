"""""
Main code to run numerical schemes and plot results
"""""

# Import packages and libraries
import numpy as np
import matplotlib.pyplot as plt
from linAdSchemes import *
from initialConditions import *
from sWSchemes import *
from math import pi

def main():
    # Set resolution
    N = 1001
    x = np.linspace(0, 1, N)

    a = 0.
    b = 1.
    g = 9.81
    H = 1
    # Initial condition
    # u_init = squareWave(x, a, b)
    # h_init = cosBell(x, alpha=a, beta=b)
    # u_init = np.sqrt(g/H)*cosBell(x, alpha=a, beta=b)
    k = 5
    h_init = np.sin(k*pi*x)
    f1 = np.zeros(N)
    f2 = np.zeros(N)
    print(f)

    # Initialise timestep, grid-size and advection velocity
    dx = 1/(N-1)
    dt = 1e-3
    n_steps = 100
    u = .1
    T = n_steps*dt

    # Courant number
    c = u*dt/dx
    print("Courant Number (c): %f" %(c))

    """""
    Pick the numerical scheme here:
    * Unstaggered SW
    """""

    # SW
    u, h = USW(u_init, h_init, f1, f2, dt, n_steps, c)

    # Plot results
    plt.plot(x,u, label="u")
    plt.plot(x,h, label="h")
    # plt.plot(x,u_exact, label="Exact")
    plt.legend()
    plt.show()

main()

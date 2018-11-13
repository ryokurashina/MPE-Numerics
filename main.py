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
    # Set resolution and parameters
    N = 64
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81
    H = 1

    # Initial condition
    k = 4
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)
    f1 = np.zeros(N)
    f2 = np.zeros(N)

    # Initialise timestep, grid-size and advection velocity
    dx = 1/(N-1)
    dt = 1e-3
    n_steps = 200
    T = n_steps*dt

    # Courant number
    c = np.sqrt(g*H)*dt/dx
    print("Courant Number (c): %f" %(c))

    """""
    Pick the numerical scheme here:
    * Unstaggered SW
    """""

    # SW
    u, h = SFB(u_init, h_init, n_steps, c, H)
    u_exact, h_exact = trav_wave(x, T, k, H, 1)
    # Plot results
    plt.figure(1)
    plt.plot(x,u, label="Numerical Solution")
    plt.plot(x,u_exact, label="Exact Solution")
    plt.plot(x,u_init, '--', label="Initial Condition")
    plt.xlabel('x')
    plt.ylabel('u')
    plt.legend()
    plt.show()

    plt.figure(2)
    plt.plot(x,h, label="Numerical Solution")
    plt.plot(x,h_exact, label="Exact Solution")
    plt.plot(x,h_init, '--',label="Initial Condition")
    plt.xlabel('x')
    plt.ylabel('h')
    plt.legend()
    plt.show()

main()

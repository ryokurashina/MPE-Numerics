# Import packages and libraries
import numpy as np
import matplotlib.pyplot as plt

from initialConditions import *
from linAdSchemes import *
from sWSchemes import *
from analyse import *
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
    n_steps = 500
    T = n_steps*dt

    x_shift = x + 0.5*dx*np.ones_like(x)

    # Courant number
    c = np.sqrt(g*H)*dt/dx
    print("Courant Number (c): %f" %(c))

    # Array to store L2-errors
    e1_u = np.zeros(n_steps)
    e2_u = np.zeros(n_steps)
    e1_h = np.zeros(n_steps)
    e2_h = np.zeros(n_steps)
    """""
    Pick the numerical scheme here:
    * Unstaggered F-B SW (UFB)
    * Staggered F-B SW (SFB)
    """""

    u_old1 = u_init
    h_old1 = h_init
    u_old2 = u_init
    h_old2 = h_init

    for i in range(n_steps):
        t = (i+1)*dt
        # SW
        u_new1, h_new1 = UFB(u_old1, h_old1, 1, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, 1, c, H , x, x_shift)
        u_exact, h_exact = trav_wave(x, t, k, H, 1)

        e1_u[i] = L2_error(u_new1, u_exact)
        e2_u[i] = L2_error(u_new2, u_exact)

        e1_h[i] = L2_error(h_new1, h_exact)
        e2_h[i] = L2_error(h_new2, h_exact)

        u_old1 = u_new1.copy()
        h_old1= h_new1.copy()

        u_old2 = u_new1.copy()
        h_old2= h_new1.copy()

    # Plot results
    plt.figure(1)
    plt.plot(range(1,n_steps+1), e1_u, label="UFB")
    plt.plot(range(1,n_steps+1), e2_u, label="SFB")
    plt.title('L2-error norms for u')
    plt.xlabel('x')
    plt.ylabel('L2_error')
    plt.legend()
    plt.show()

    # Plot results
    plt.figure(2)
    plt.plot(range(1,n_steps+1), e1_h, label="UFB")
    plt.plot(range(1,n_steps+1), e2_h, label="SFB")
    plt.title('L2-error norms for h')
    plt.xlabel('Timestep')
    plt.ylabel('L2_error')
    plt.legend()
    plt.show()

main()

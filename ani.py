"""
Animation code to see how numerical scheme develops in real time.

It has been set up to show the strange "oscillation" of the solution about the
exact solution.
"""

import os
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt

from initialConditions import *
from linAdSchemes import *
from sWSchemes import *
from analyse import *


def animation(N, k, dt, n_steps, H):
    """Animation function"""

    # Parameters
    x_ = np.linspace(0, 1, N+1)
    x = x_[0:N]
    g = 9.81
    H = 1
    dx = 1/(N-1)
    x_shift = x + 0.5*dx*np.ones_like(x)
    dt = 1e-3

    # Initial condition
    k = 1
    t = 0
    u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
    h_init = np.cos(2*pi*k*x)

    # Courant number
    c = np.sqrt(g*H)*dt/dx
    print("Courant Number (c): %f" %(c))

    plt.ion()
    plt.pause(1)

    u_old1 = u_init
    h_old1 = h_init

    u_old2 = u_init
    h_old2 = h_init

    t = 0

    # Loop through 1 time-step at a time nd plot for animation
    for i in range(n_steps):
        # Compute the solution for one time-step
        u_new1, h_new1 = UFB(u_old1, h_old1, 1, c, H)
        u_new2, h_new2 = SFB(u_old2, h_old2, 1, c, H, x, x_shift)
        u_exact, h_exact = trav_wave(x, t, k, H, 1)
        # Plot results
        plt.title('Animation')
        plt.plot(x,u_new2,label='u (SFB)')
        plt.plot(x,u_exact,label='u_exact')
        plt.xlabel('x')
        plt.legend()
        plt.draw()
        plt.pause(0.1)
        plt.clf()
        u_old1 = u_new1.copy()
        h_old1 = h_new1.copy()
        u_old2 = u_new2.copy()
        h_old2 = h_new2.copy()
        # Update time
        t = (i+1)*dt
    plt.ioff()

animation(16, 1, 1e-3, 500, 1)

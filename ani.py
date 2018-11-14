"""""
Animation code to see how numerical scheme develops in real time. This isn't
so good for post-processing or analysis of results but more as a visual aid for
the user
"""""

import os
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt

from initialConditions import *
from linAdSchemes import *
from sWSchemes import *
from analyse import *

# Parameters
N = 64
x_ = np.linspace(0, 1, N+1)
x = x_[0:N]
g = 9.81
H = 1
dx = 1/(N-1)
x_shift = x + 0.5*dx*np.ones_like(x)
dt = 1e-3
n_steps = 500

# Initial condition
k = 4
t = 0
u_init = -np.sqrt(g/H)*np.cos(2*pi*k*x)
h_init = np.cos(2*pi*k*x)

# Courant number
c = np.sqrt(g*H)*dt/dx
print("Courant Number (c): %f" %(c))

plt.ion()
plt.pause(1)

u_old = u_init
h_old = h_init

t = 0

# Loop through 1 time-step at a time and plot for animation
for i in range(n_steps):
    # Update f1 and f2
    # Compute the solution for one time-step
    u_new, h_new = SFB(u_old, h_old, 1, c, H, x, x_shift)
    u_exact, h_exact = trav_wave(x, t, k, H, 1)
    # Plot results
    plt.title('Animation')
    plt.plot(x,h_new,label='h')
    plt.plot(x,u_new,label='u')
    plt.plot(x,h_exact,label='h_exact')
    plt.plot(x,u_exact,label='u_exact')
    # plt.plot(x,u_exact,label='u_exact')
    plt.xlabel('x')
    plt.legend()
    plt.draw()
    plt.pause(0.1)
    plt.clf()
    u_old = u_new.copy()
    h_old = h_new.copy()
    # Update time
    t = (i+1)*dt
plt.ioff()

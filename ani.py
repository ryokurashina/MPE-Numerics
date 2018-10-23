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

N = 1001
x = np.linspace(0, 1, N)
a = 0
b = .5
# Initial condition
phi_init = cosBell(x,alpha=0.,beta=1.)
h_init = cosBell(x, alpha=0., beta=1.)
# Initialise timestep, grid-size and advection velocity
dx = 1/(N-1)
dt = 1e-3
n_steps = 200
u = .1
T = n_steps*dt

# Courant number
c = u*dt/dx
print("Courant Number (c): %f" %(c))
plt.ion()
plt.pause(5)

phi_old = phi_init
h_old = h_init

# Loop through 1 time-step at a time and plot for animation
for i in range(n_steps):
    phi_new, h_new = USW(phi_old, h_old, 1, 0.5)
    plt.title('Animation')
    plt.plot(x,phi_new,label='u')
    plt.plot(x,h_new,label='h')
    plt.xlabel('x')
    plt.legend()
    plt.draw()
    plt.pause(0.2)
    plt.clf()
    phi_old = phi_new.copy()
    h_old = h_new.copy()

plt.ioff()

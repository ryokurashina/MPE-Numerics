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
from postProcess import *

N = 251
x_ = np.linspace(0, 1, N+1)
x = x_[0:N]
g = 9.81
H = 1
# Initial condition
# u_init = np.sqrt(g/H)*cosBell(x,alpha=0.,beta=1.)
# u_init = np.zeros(N)
# k = 2
# h_init = np.sin(k*2*pi*x)
k = 2
u_init = exact_sol(k,x,0)
h_init = exact_sol(k,x,0)

f1 = np.zeros(N)
f2 = np.zeros(N)

# Initialise timestep, grid-size and advection velocity
dx = 1/(N-1)
dt = 1e-3
n_steps = 500

# Courant number
c = np.sqrt(g*H)*dt/dx
print("Courant Number (c): %f" %(c))

# f1, f2 = source_f(x, 0, k, c)

plt.ion()
plt.pause(1)

u_old = u_init
h_old = h_init

t = 0

# Loop through 1 time-step at a time and plot for animation
for i in range(n_steps):
    f1, f2 = source_f(x, t, k, c)
    u_new, h_new = USW(u_old, h_old, f1, f2, dt, 1, c)
    u_exact = exact_sol(k ,x, t)
    plt.title('Animation')
    plt.plot(x,h_new,label='u')
    # plt.plot(x,h_new,label='h')
    plt.plot(x,u_exact,label='u_exact')
    plt.xlabel('x')
    plt.legend()
    plt.draw()
    plt.pause(0.2)
    plt.clf()
    u_old = u_new.copy()
    h_old = h_new.copy()
    # Update time
    t = (i+1)*dt
plt.ioff()

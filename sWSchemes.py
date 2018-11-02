"""""
Contains various numerical schemes for the linearised 1D-SW equations with f = 0.

List of numerical schemes:
* Unstaggered SW (USW)
"""""

import numpy as np
import matplotlib.pyplot as plt
from postProcess import *
from math import pi

def exact_sol(k,x,t):
    # Source functions have been designed such that both u and h' are the same
    f = np.sin(2*k*pi*x)*np.cos(2*pi*t)
    return f

def source_f(x, t, k, c):
    # Gravitational acceleration
    g = 9.81
    H = c**2/g
    f1 = 2*k*pi*g*np.cos(2*k*pi*x)*np.sin(2*pi*t)-2*pi*np.sin(2*k*pi*x)*np.sin(2*pi*t)
    f2 = 2*k*pi*H*np.cos(2*k*pi*x)*np.sin(2*pi*t)-2*pi*np.sin(2*k*pi*x)*np.sin(2*pi*t)
    return f1, f2

def USW(u, h, f1, f2, dt, ntime, c):
    # Gravitational aceleration
    g = 9.81
    H = c**2/g
    # print("H = %f" %(H))
    N = len(u)
    # Initialise for loop
    u_old = np.array(u)
    u_new = np.zeros(N)
    h_old = np.array(h)
    h_new = np.zeros(N)
    for i in range(ntime):
        u_new = u_old-c/2*np.sqrt(g/H)*(np.roll(h_old,1)-np.roll(h_old,-1))+f1*dt
        h_new = h_old-c/2*np.sqrt(H/g)*(np.roll(u_new,1)-np.roll(u_new,-1))+f2*dt
        print(L2_error(h_old,h_new))
        u_old = u_new.copy()
        h_old = h_new.copy()
    return u_new, h_new

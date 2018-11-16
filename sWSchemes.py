"""""
Contains various numerical schemes for the linearised 1D-SW equations with f = 0.

List of numerical schemes:
* Unstaggered SW (USW)
"""""

import numpy as np
import matplotlib.pyplot as plt

from analyse import *
from math import pi
from scipy.interpolate import interp1d

def trav_wave(x, t, k, H, direction):
    g = 9.81
    t_vec = t*np.ones_like(x)
    # Take left-travelling wave if 0
    if direction == 0:
        u = np.sqrt(g/H)*np.cos(2*pi*k*(x+np.sqrt(g*H)*t))
        h = np.cos(2*pi*k*(x+np.sqrt(g*H)*t))
    elif direction == 1:
        u = -np.sqrt(g/H)*np.cos(2*pi*k*(x-np.sqrt(g*H)*t))
        h = np.cos(2*pi*k*(x-np.sqrt(g*H)*t))
    return u, h

def UFB(u, h, ntime, c, H):
    # Gravitational aceleration
    g = 9.81
    N = len(u)
    # Initialise for loop
    u_old = np.array(u)
    u_new = np.zeros(N)
    h_old = np.array(h)
    h_new = np.zeros(N)
    for i in range(ntime):
        u_new = u_old-c/2*np.sqrt(g/H)*(np.roll(h_old,1)-np.roll(h_old,-1))
        h_new = h_old-c/2*np.sqrt(H/g)*(np.roll(u_new,1)-np.roll(u_new,-1))
        u_old = u_new.copy()
        h_old = h_new.copy()
    return u_new, h_new

def SFB(u, h, ntime, c, H, x, x_shift):
    # Gravitational aceleration
    g = 9.81
    N = len(u)
    # Initialise for loop
    u_old = np.array(u)
    u_new = np.zeros(N)
    h_old = np.array(h)
    h_new = np.zeros(N)
    # Stagger u to the right
    x_ = np.concatenate((x, np.array([1])))
    u_end = np.array([u_old[0]])
    u = np.concatenate((u, u_end))
    # Create functions for cubically interpolated values
    f = interp1d(x_, u, kind='cubic')
    u_stag_old = f(x_shift)
    for i in range(ntime):
        u_stag_new = u_stag_old-c*np.sqrt(g/H)*(np.roll(h_old,1)-h_old)
        h_new = h_old-c*np.sqrt(H/g)*(u_stag_new-np.roll(u_stag_new,-1))
        u_stag_old = u_stag_new.copy()
        h_old = h_new.copy()
    # Stagger u back to the left again
    x_shift_ = np.concatenate((np.array([-x_shift[0]]),x_shift))
    u_start = np.array([u_stag_new[N-1]])
    u = np.concatenate((u_start, u_stag_new))
    f_ = interp1d(x_shift_, u, kind='cubic')
    u_new = f_(x)
    return u_new, h_new

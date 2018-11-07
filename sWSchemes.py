"""""
Contains various numerical schemes for the linearised 1D-SW equations with f = 0.

List of numerical schemes:
* Unstaggered SW (USW)
"""""

import numpy as np
import matplotlib.pyplot as plt

from postProcess import *
from math import pi
from scipy.interpolate import interp1d

def trav_wave(x, t, k, H, direction):
    g = 9.81
    t_vec = t*np.ones_like(x)
    # Take left travelling wave if 0
    if direction == 0:
        u = np.sqrt(g/H)*np.cos(2*pi*k*(x+np.sqrt(g*H)*t))
        h = np.cos(2*pi*k*(x+np.sqrt(g*H)*t))
    elif direction == 1:
        u = -np.sqrt(g/H)*np.cos(2*pi*k*(x-np.sqrt(g*H)*t))
        h = np.cos(2*pi*k*(x-np.sqrt(g*H)*t))
    return u, h

def USW(u, h, ntime, c, H):
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

# Linear interpolator for two arrays, takes a left and right array to calculate
# middle value.
def lin_interp(l, r):
    m = 0.5*(l+r)
    return m

def SSW(u, h, ntime, c, H):
    # Gravitational aceleration
    g = 9.81
    N = len(u)
    # Initialise for loop
    u_old = np.array(u)
    u_new = np.zeros(N)
    h_old = np.array(h)
    h_new = np.zeros(N)
    # Stagger u to the right
    u_stag_old = 0.5*(u_old+np.roll(u_old, 1))
    for i in range(ntime):
        u_stag_new = u_stag_old-c/2*np.sqrt(g/H)*(np.roll(h_old,1)-np.roll(h_old,-1))
        h_new = h_old-c/2*np.sqrt(H/g)*(np.roll(u_stag_new,1)-np.roll(u_stag_new,-1))
        u_stag_old = u_stag_new.copy()
        h_old = h_new.copy()
    # Stagger u back to the left again
    u_new = 0.5*(np.roll(u_stag_new, -1)+u_stag_new)
    return u_new, h_new

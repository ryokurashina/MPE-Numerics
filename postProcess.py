"""""
Script for post processing of data
"""""

import numpy as np
import matplotlib.pyplot as plt
from linAdSchemes import *
from initialConditions import *

def main():
    N = 1001
    x = np.linspace(0, 1, N )
    phi_init = squareWave(x, 0, .25)
    # Initialise timesteps and advection velocity to calculate Courant no.
    dx = 1/(N-1)
    dt = 1e-3
    u = .1
    c = u*dt/dx
    print("Courant Number (c): %f" %(c))

    phi = FTCS(phi_init,100,c)

    plt.plot(x,phi)
    plt.show()
main()

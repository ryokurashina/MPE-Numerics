"""""
File that contains functions used for post-processing of data/schemes.
"""""
# Import packages and libraries
import numpy as np

# Finds the total amount of a quantity phi
def total_phi(phi):
    phi = np.array(phi)
    total_phi = np.sum(phi)
    return total_phi

# Calculates the L2-norm error for two arrays for functions on the unit interval
def L2_error(phi1, phi2):
    N = len(phi1)
    errors = np.absolute(np.array(phi1-phi2))
    L2 = sum(np.square(errors))/(N-1)
    

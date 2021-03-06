"""""
File that contains functions used for analysis of data/schemes.
"""""
# Import packages and libraries
import numpy as np


def total(phi):
    """ Finds the total amount of a quantity phi """
    N = len(phi)
    phi = np.array(phi)
    total_phi = np.sum(phi)/(N-1)
    return total_phi

def total_abs(phi):
    """ Calculates the L1-norm of a quantity phi """
    N = len(phi)
    phi = np.array(phi)
    total_phi = np.sum(np.absolute(phi))/(N-1)
    return total_phi

def expect(phi):
    """  Calculates the average value of a quantity phi """
    N = len(phi)
    # Sum and divide by total number of points
    phi_expect = np.sum(phi)/N
    return phi_expect


def p_moment(phi, p):
    """  Find the centred p-th moment of a quantity phi """
    N = len(phi)
    phi = np.array(phi)
    # Expectation of phi
    phi_expect = expect(phi)
    # Center phi
    phi_center = phi-phi_expect*np.ones_like(phi)
    mom_phi = np.sum(np.power(phi_center, p))/(N-1)
    return mom_phi


def L2_error(phi1, phi2):
    """
    Calculates the L2-norm error for two arrays for functions on the unit
    interval.
    """
    N = len(phi1)
    errors = np.absolute(np.array(phi1-phi2))
    L2 = sum(np.square(errors))/(N-1)
    return L2

"""""
Main code to run numerical schemes and plot results.
"""""

# Import packages and libraries
import numpy as np
import matplotlib.pyplot as plt
import os

from sWSchemes import *
from postProcess import *
from math import pi

# Run these functions to produce all plots required for report
dispersion(8, 1, 1e-3, 100, 1)
L1_growth(8, 1, 1e-3, 500, 1)
L2_convergence_x(1, 1e-4, 200, 1)
calc_moments(8, 1, 1e-3, 1000, 1, 2)

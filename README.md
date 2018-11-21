# MPE-Numerics

Files and descriptions of files are: 

* main.py
Main file to be run to produce all plots. Plots are saved in the /plots/ directory and all plots will be ignored by git because of the .gitignore file.

* ani.py
File that produces an animation of the SFB code to show the strange oscillatory behaviour of the scheme (details will be in final report).

* analyse.py 
File that contains functions for post-processing analysis of data.

* postProcess.py
File that contains functions that produces plots to show various results such as spatial convergence of schemes and dispersion properties of schemes.

* sWSchemes.py
File that contains the numerical schemes for solving the 1-D linearised shallow water equations.

N.B. The scipy package may need to be installed in order to use its 1-D interpolation function.

 ### Run main.py to get all the plots for the report which will be saved into the /plots/ directory.

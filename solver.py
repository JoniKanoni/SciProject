"""
Module that Creates a Hamiltonian and contains the corresponding solver
"""
import numpy as np

def hamilton_operator(potential, mass, xnum):
    """
    Function that creates the hamiltonian based on a potential grid and mass
    of the particle
    
    Args:
        potential:  2D array that contains the potential in the second column
                    and the corresponding x-values in the first column
        mass:       Mass of the particle
        xnum:       Number of points in the grid
        
    returns:
        ham:        Matrix of the hamiltonian on the spacegrid of the potential
    """

    ham = np.zeros((xnum-2, xnum-2))
    a = 1/(mass * (potential[1,0] - potential[0,0]))
    for ii in range(0,xnum-2):
        if ii - 1 >= 0:
            ham[ii, ii - 1] = - 1 / (2 * a)
        if ii + 1 < xnum-2:
            ham[ii, ii + 1] = - 1 / (2 * a)
        ham[ii, ii] = potential[ii+1, 1]
    return ham

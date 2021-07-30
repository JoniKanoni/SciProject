"""
Module that Creates a Hamiltonian and contains the corresponding solver
"""
import numpy as np
import scipy.linalg as la

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

def diag_solver(mass, potential, xnum, first_val, last_val):
    """
    Function that diagonalizes the hamiltonian
    
    Args:
        mass:       Mass of the particle
        potential:  Potential grid with corresponding x-values.
                    (x-values in first column)
                    The x-values have to be equidistant
        xnum:       Number of x-values in the grid
        first_val:  Number of first desired eigenvalue
        last_val:   Number of last desired eigenvalue
        
    Returns:
        eigvals:    1D Array that contains the n eigenvalues
        wavefuncs:  n-D Array that contains the corresponding wavefuncs
    """

    a = 1/(mass * (potential[1,0] - potential[0,0])**2)
    diagonal = potential[1:xnum,1] + a
    neben_diagonal = np.zeros((xnum-1,1))
    neben_diagonal[:] = - a / 2
    eigvals, wavefuncs = la.eigh_tridiagonal(diagonal, neben_diagonal, (first_val, last_val))
    return eigvals, wavefuncs
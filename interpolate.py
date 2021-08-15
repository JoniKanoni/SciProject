'''
Module with multiple options to interpolate a given function
'''

import numpy as np
from scipy.interpolate import CubicSpline, interp1d, lagrange

def interpol(xvalues, yvalues, inttype):
    '''
    Function that creates an interpolation function based on Input

    Args:
        xvalues:    known x-values
        yvalues:    known y-values
        inttype:    Type of desired kind of interpolation

    returns:
        int_fct:    function that returns new y-value for arbitrary x-value,
                    based on the interpolation
    '''

    if inttype == 'linear':
        int_fct = interp1d(xvalues, yvalues, kind='linear')
    elif inttype == 'cspline':
        int_fct = CubicSpline(xvalues, yvalues)
    elif inttype == 'polynomial':
        int_fct = lagrange(xvalues, yvalues)
    return int_fct

def potential_grid(xmin, xmax, xnum, xvalues, potvalues, inttype):
    """
    Function that creates a grid (2D array) of Potential values and
    their corresponding x-values
    Args:
        xmin:       smallest considered x-value
        xmax:       largest cosidered x-value
        xnum:       number of desired x-values
        int_fct:    function used to create Potentialvalues with interpolation
        
    Returns:
        potential:  2D Array with the x-values in the first column and
                    the potential in the second column
    """

    potential = np.zeros((xnum,2))
    potential[:,0] = np.linspace(xmin, xmax, xnum)
    potential[:,1] = interpol(xvalues, potvalues, inttype)(potential[:,0])
    return potential

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
        hamiltonian:        Matrix of the hamiltonian on the spacegrid of the potential
    """

    hamiltonian = np.zeros((xnum-2, xnum-2))
    a = 1/(mass * (potential[1,0] - potential[0,0])**2)
    for ii in range(0,xnum-2):
        if ii - 1 >= 0:
            hamiltonian[ii, ii - 1] = - a / 2
        if ii + 1 < xnum-2:
            hamiltonian[ii, ii + 1] = - a / 2
        hamiltonian[ii, ii] = potential[ii+1, 1] + a
    return hamiltonian

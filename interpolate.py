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
        int_fct = CubicSpline(xvalues, yvalues, bc_type='natural')
    elif inttype == 'polynomial':
        int_fct = lagrange(xvalues, yvalues)
    return int_fct

def potential_grid(xmin, xmax, xnum, int_fct):
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
    potential[:,1] = int_fct(potential[:,0])
    return potential

'''
Module with multiple options to interpolate a given function
'''

from scipy.interpolate import CubicSpline, interp1d, lagrange

def interpol(xx, yy, inttype):
    '''
    Function that creates an interpolation function based on Input
    
    Args:
        xx:         known x-values
        yy:         known y-values
        inttype:    Type of desired kind of interpolation
        
    returns:
        int_fct:    function that returns new y-value for arbitrary x-value,
                    based on the interpolation
    '''

    if inttype == 'linear':
        int_fct = interp1d(xx,yy, kind='linear')
    elif inttype == 'cspline':
        int_fct = CubicSpline(xx,yy, bc_type='natural')
    elif inttype == 'polynomial':
        int_fct = lagrange(xx,yy)
    return int_fct 

    
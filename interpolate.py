'''
Module with multiple options to interpolate a given function
'''

from scipy.interpolate import CubicSpline, interp1d, lagrange

def interpolatierer(xx, yy, inttype):
    '''
    dgiuohdshodgiho
    '''

    if inttype == 'linear':
        resultpot = interp1d(xx,yy, kind='linear')
    elif inttype == 'cspline':
        resultpot = CubicSpline(xx,yy, bc_type='natural')
    elif inttype == 'polynomial':
        resultpot = lagrange(xx,yy)
    return resultpot

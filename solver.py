"""
Module that Creates a Hamiltonian and contains the corresponding solver
"""
import numpy as np
import scipy.linalg as la



def Integrate(fct,xvals):
    """
    function that numerically integrates a 1D function in the form of an array over a given
    xvalue grid
    
    Args:
        fct:        1D array that contains the values of the function
        xvals:      1D array of corresponding xvalues (have to be equidistant)
                    (Integration considers order by comparing the first two elements
                    ---> Switching these elements results in sign change)
        
    Returns:
        Integral:   Value of the Integral (float)
    """
    deltax = (xvals[1] - xvals[0])
    Integral = np.sum(fct[1:]) * deltax
    return Integral



def QM_Norming(fct, xvals):
    """
    Function that norms a L2-function (standard norm of QM) based on given xvalue grid.
    The function is assumed to be 0 outside the grid.
    
    Args:
        fct:        1D array that contains the values of the function
        xvals:      1D array of corresponding xvalues (have to be equidistant)
                    (Integration goes from first xvalue to last xvalue
                    ---> Order of values matters!!!)
    
    Returns:
        fct_normed: 1D Array that contains the values of the normed function
    """
    fct_normed = fct / np.sqrt(Integrate(np.square(np.abs(fct)), xvals))
    return fct_normed



def QM_position_info(wavefuncs, xvalues = np.array([])):
    """
    Function that computes the expectation value for the position of a particle, as well
    as the uncertainty of the position based on given 1D-wavefunction and xvalue grid.
    The function is assumed to be 0 outside the grid.
    
    Args:
        wavefunc:   1D array of values of a wavefunction
        xvals:      1D array of the corresponding x values
    
    Returns:
        exp_x:      Expectationvalue for the position (float)
        unc_x:      Uncertainty of position (float)
    """
    eigcount = np.shape(wavefuncs)[1]
    if xvalues == np.array([]):
        expvalues = np.zeros((eigcount-1, 2))
        for ii in range (0,eigcount-1):
            expvalues[ii,0] = Integrate(np.square(np.abs(wavefuncs[:,ii+1]))*wavefuncs[:, 0],wavefuncs[:, 0])
            exp_x2 = Integrate(np.square(np.abs(wavefuncs[:,ii+1]))*np.square(wavefuncs[:, 0]),wavefuncs[:, 0])
            expvalues[ii,1] = np.sqrt(exp_x2 - expvalues[ii,0])
    else:
        expvalues = np.zeros((eigcount, 2))
        for ii in range (0,eigcount):
            expvalues[ii,0] = Integrate(np.square(np.abs(wavefuncs[:,ii]))*xvalues,xvalues)
            exp_x2 = Integrate(np.square(np.abs(wavefuncs[:,ii]))*np.square(xvalues),xvalues)
            expvalues[ii,1] = np.sqrt(exp_x2 - np.square(expvalues[ii,0]))
    return expvalues

def QM_wavefct(hamiltonian, xnum, first_val, last_val, xvalues):
    """
    Function that computes the wavefunction based on a hamiltonian in position basis
    """
    energies, eigenvec = la.eigh(hamiltonian, eigvals = (first_val-1,last_val-1))
    wavefuncs = np.zeros((xnum-2, last_val - first_val + 2))
    wavefuncs[:,0] = xvalues
    for ii in range (1,last_val - first_val + 2):
        wavefuncs[:,ii] = QM_Norming(eigenvec[:,ii-1], xvalues)
    return energies, wavefuncs
    
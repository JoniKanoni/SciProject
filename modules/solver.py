"""
Module that calculates qunatummechanical quantities based on a given hamiltonian
"""
import numpy as np
import scipy.linalg as la



def _integrate(fct,xvalues):
    """
    Function that numerically integrates a 1D function
    in the form of an array over a given xvalue grid

    Args:
        fct:        1D array that contains the values of the function
        xvalues:      1D array of corresponding xvalues (have to be equidistant)
                    (Integration considers order by comparing the first two elements
                    ---> Switching these elements results in sign change)

    Returns:
        integral:   Value of the Integral (float)
    """

    delta_x = (xvalues[1] - xvalues[0])
    integral = np.sum(fct[1:]) * delta_x
    return integral



def _qm_norming(fct, xvalues):
    """
    Function that norms a L2-function based on given xvalue grid.
    The function is assumed to be 0 outside the grid.

    Args:
        fct:        1D array that contains the values of the function
        xvalues:    1D array of corresponding xvalues (have to be equidistant)
                    (Integration only uses first and second values
                    ---> Order of these values matters!!!)
    Returns:
        fct_normed: 1D Array that contains the values of the normed function
    """
    fct_normed = fct / np.sqrt(_integrate(np.square(np.abs(fct)), xvalues))
    return fct_normed



def qm_position_info(wavefuncs, xvalues):
    """
    Function that computes the expectation value for the position of a particle, as well
    as the variance of the position based on given 1D-wavefunction and xvalue grid.
    The function is assumed to be 0 outside the grid.

    Args:
        wavefunc:   1D array of values of a wavefunction
        xvalues:    1D array of the corresponding x values

    Returns:
        expvalues:  2D array that contains the expectation value and variance of position
    """
    eig_count = np.shape(wavefuncs)[1]
    expvalues = np.zeros((eig_count, 2))
    for ii in range (0,eig_count):
        expvalues[ii,0] = _integrate(np.square(np.abs(wavefuncs[:,ii]))*xvalues,xvalues)
        exp_x2 = _integrate(np.square(np.abs(wavefuncs[:,ii]))*np.square(xvalues),xvalues)
        expvalues[ii,1] = np.sqrt(exp_x2 - np.square(expvalues[ii,0]))
    return expvalues

def qm_wavefct(hamiltonian, xnum, first_val, last_val, xvalues):
    """
    Function that computes the wavefunction based on a hamiltonian in position basis.
    The position has to be 1-dimensional

    Args:
        hamiltonian:    The hamiltonian in position basis (ndarray)
        xnum:           Number of considered xvalues
        first_val:      first desired eigenvalue as int
        last_val:       last desired eigenvalue as int
        xvalues:        1d array that contains all position values

    Returns:
        energies:       Eigenvalues of the given hamiltonian in given range (1darray)
        wavefuncs:      Corresponding wavefunctions on grid of given xvalues
    """
    energies, eigen_vec = la.eigh(hamiltonian, eigvals = (first_val-1,last_val-1))
    wavefuncs = np.zeros((xnum-2, last_val - first_val + 2))
    wavefuncs[:,0] = xvalues
    for ii in range (1,last_val - first_val + 2):
        wavefuncs[:,ii] = _qm_norming(eigen_vec[:,ii-1], xvalues)
    return energies, wavefuncs
    
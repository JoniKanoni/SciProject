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
    a = 1/(mass * (potential[1,0] - potential[0,0])**2)
    for ii in range(0,xnum-2):
        if ii - 1 >= 0:
            ham[ii, ii - 1] = - a / 2
        if ii + 1 < xnum-2:
            ham[ii, ii + 1] = - a / 2
        ham[ii, ii] = potential[ii+1, 1] + a
    return ham

def diag_solver(hamiltonian, xnum, first_val, last_val):
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
    eigvals, wavefuncs = la.eigh(hamiltonian, eigvals = (0,5))
    """
    wavefuncs = np.zeros((xnum+2,last_val-first_val+2))
    print(np.shape(h_wavefuncs))
    print(np.shape(wavefuncs))
    wavefuncs[:,1:xnum+1] = h_wavefuncs[:,:]
    """
    return eigvals, wavefuncs

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
    Integral = np.sum(fct) * (xvals[1] - xvals[0])
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



def QM_position_info(wavefunc, xvals):
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
    exp_x = Integrate(np.square(np.abs(wavefunc))*xvals,xvals)
    exp_x2 = Integrate(np.square(np.abs(wavefunc))*np.square(xvals),xvals)
    unc_x = np.sqrt(exp_x2 - exp_x**2)
    return exp_x, unc_x
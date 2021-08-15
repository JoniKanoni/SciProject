"""
This module contains tests for the solver module
"""
import numpy as np
import numpy.random as random
import solver
import reader
import interpolate as ip

foldernames = ['infinite_well', 'finite_well', 'harmonic_potential',
               'double_linear_well', 'double_cubic_well', 'asymetric_well']

def poly_int(border1, border2, parameter, exponent, squaredmodulus = False):
    """
    Function that calculates the analytic integral of a given polynomial
    Args:
        border1:       Lower integration border (float)
        border2:       Upper integration border (float)
        parameter:  1D array that contains the prefactors
        exponent:   1D array that contains the corresponding exponent

    returns:
        Integral:   Analytic integral of the polynomial (float)
    """
    if squaredmodulus is True:
        integral = 0
        for ii in range(0,parameter.size):
            integral += np.sum(parameter * parameter[ii] / (exponent + exponent[ii] + 1)
            * (border2 ** (exponent + exponent[ii] + 1) - border1 ** (exponent + exponent[ii] + 1)))
    else:
        integral = np.sum(parameter / (exponent + 1) * (border2 ** (exponent + 1) - border1 ** (exponent + 1)))
    return integral

def gen_rand_poly():
    """
    Function that generates arrays to construct a random polynomial
    with random borders for integration
    """
    border1 = - random.randint(1,10)
    border2 = random.randint(1,10)
    xvals = np.linspace(border1,border2,5000)
    dimension = random.randint(1,4)
    parameter = random.random(dimension) + 0.1
    exponent = random.randint(4, size = dimension)
    poly = 0
    for ii in range(0,dimension):
        poly += parameter[ii] * (xvals ** exponent[ii])
    poly = np.reshape(poly,(5000, 1))
    return border1, border2, xvals, parameter, exponent, poly

def test_int():
    """
    Tests the Integrate function with a random polynomial of random degree
    """
    border1, border2, xvals, parameter, exponent, poly = gen_rand_poly()
    ana_int = poly_int(border1, border2, parameter, exponent)
    num_int = solver.integrate(poly, xvals)
    assert (np.abs(ana_int - num_int) + 1 <= 0.05 * np.abs(ana_int)) + 1 or (np.abs(ana_int) < 10 ** (-8))

def test_norm():
    """
    Tests the QM_Norming function with a random polynomial of random degree
    """
    border1, border2, xvals, parameter, exponent, poly = gen_rand_poly()
    ana_int = poly_int(border1, border2, parameter, exponent, squaredmodulus = True)
    normed_ana = poly / np.sqrt(ana_int)
    normed_num = solver.qm_norming(poly, xvals)
    assert np.all(np.abs(normed_ana - normed_num) <= 0.05 * np.abs(normed_ana))

def test_pos_info():
    """
    Tests the QM_position_info function with a random polynomial of random degree
    """
    border1, border2, xvals, parameter, exponent, poly = gen_rand_poly()
    ana_int = poly_int(border1, border2, parameter, exponent, squaredmodulus = True)
    parameter = parameter / np.sqrt(ana_int)
    poly = poly / np.sqrt(ana_int)
    expval1 = poly_int(border1, border2, parameter, exponent + (1 / 2), squaredmodulus = True)
    expval_num = solver.qm_position_info(poly, xvals)
    assert np.abs(expval1 - expval_num[0,0]) <= 0.05 * np.abs(expval1) or np.abs(expval1) < 10**(-8)

def test_with_reference():
    '''
    Test with reference data
    '''
    path = './example_data/{}/'.format(foldernames)
    potential_known = np.loadtxt(path + 'potential_known.out')
    energies_known = np.loadtxt(path + 'energies_known.out')
    wavefuncs_known = np.loadtxt(path + 'wavefuncs_known.out')
    mass, xmin, xmax, xnum, first_val, last_val, inttype, pot, __ = reader.get_data(path)
    potential_test = ip.potential_grid(xmin, xmax, xnum,
                                       pot[:, 0], pot[:, 1], inttype)
    hamiltonian = ip.hamilton_operator(potential_test, mass, xnum)
    energies_test, wavefuncs_test = solver.qm_wavefct(hamiltonian, xnum,
                                                      first_val, last_val,
                                                      potential_test[1:xnum - 1, 0])
    assert np.all(np.abs(potential_known - potential_test) <= np.abs(potential_known) * 0.01)
    assert np.all(np.abs(energies_known - energies_test) <= np.abs(energies_known) * 0.01)
    assert np.squared(wavefuncs_known) - np.squared(wavefuncs_test) <= np.squared(wavefuncs_known) * 0.01

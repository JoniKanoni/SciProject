"""
This module contains tests for the solver module
"""
import numpy as np
import numpy.random as random
import solver

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
    if squaredmodulus == True:
        Integral = 0
        for ii in range(0,parameter.size):
            Integral += np.sum(parameter * parameter[ii] / (exponent + exponent[ii] + 1) * (border2 ** (exponent + exponent[ii] + 1) - border1 ** (exponent + exponent[ii] + 1)))
    else:
        Integral = np.sum(parameter / (exponent + 1) * (border2 ** (exponent + 1) - border1 ** (exponent + 1)))
    return Integral

def test_Int():
    """
    Tests the Integrate function with a random polynomial of random degree
    """
    border1 = - random.randint(1,10)
    border2 = random.randint(1,10)
    xvals = np.linspace(border1,border2,2000)
    dimension = random.randint(1,10)
    parameter = random.random(dimension)
    exponent = random.randint(10, size = dimension)
    ana_int = poly_int(border1, border2, parameter, exponent)
    poly = 0
    for ii in range(0,dimension):
        poly += parameter[ii] * (xvals ** exponent[ii])
    num_int = solver.Integrate(poly, xvals)
    print(ana_int, num_int)
    assert np.abs((ana_int - num_int) / ana_int) < 0.01

def test_Norm():
    """
    Tests the QM_Norming function with a random polynomial of random degree
    """
    border1 = - random.randint(1,10)
    border2 = random.randint(1,10)
    xvals = np.linspace(border1,border2,2000)
    dimension = random.randint(1,10)
    parameter = random.random(dimension)
    exponent = random.randint(10, size = dimension)
    ana_int = poly_int(border1, border2, parameter, exponent, squaredmodulus = True)
    poly = 0
    for ii in range(0,dimension):
        poly += parameter[ii] * (xvals ** exponent[ii])
    Normed_ana = poly / np.sqrt(ana_int)
    Normed_num = solver.QM_Norming(poly, xvals)
    assert np.all(np.abs((Normed_ana - Normed_num) / Normed_ana) < 0.01)

def test_Pos_info():
    """
    Tests the QM_position_info function with a random polynomial of random degree
    """
    border1 = - random.randint(1,10)
    border2 = random.randint(1,10)
    xnum = 10000000
    xvals = np.linspace(border1,border2,xnum)
    dimension = random.randint(1,10)
    parameter = random.randint(10, size = dimension)
    exponent = random.randint(10, size = dimension)
    ana_int = poly_int(border1, border2, parameter, exponent, squaredmodulus = True)
    parameter = parameter / ana_int
    expval1 = poly_int(border1, border2, parameter, exponent + (1 / 2), squaredmodulus = True)
    expval2 = poly_int(border1, border2, parameter, exponent + 1, squaredmodulus = True)
    poly = 0
    for ii in range(0,dimension):
        poly += parameter[ii] * (xvals ** exponent[ii])
    poly = np.reshape(poly,(xnum, 1))
    expval1_num = solver.QM_position_info(poly, xvals)
    variance = np.sqrt(expval2 - np.square(expval1))
    print('!!!!!!!!!!!EXPVAL2:', expval2, '!!!!!!!!')
    assert np.abs((expval1 - expval1_num[0,0]) / expval1) < 0.01 and np.abs((variance - expval1_num[0,1]) / variance) < 0.01

"""
This module contains tests for the solver module
"""
import numpy as np
import pytest
from modules import solver
from modules import reader
from modules import interpolate as ip

foldernames = ['infinite_well', 'finite_well', 'harmonic',
               'double_well_linear', 'double_well_kubic', 'asym_well']


@pytest.mark.parametrize("foldernames", foldernames)
def test_with_reference_potential(foldernames):
    '''
    Tests the interpolation of the potential with reference data
    '''
    path = './results/{}/'.format(foldernames)
    potential_known = np.loadtxt(path + 'potential.dat')
    mass, xmin, xmax, xnum, first_val, last_val, inttype, pot, __ = reader.get_data(path)
    potential_test = ip.potential_grid(xmin, xmax, xnum,
                                       pot[:, 0], pot[:, 1], inttype)
    assert np.all(np.abs(potential_known - potential_test) <= np.abs(potential_known) * 0.01)

@pytest.mark.parametrize("foldernames", foldernames)
def test_with_reference_energies_wavefuncs(foldernames):
    '''
    Tests the eigenvalues and wavefunctions with reference data
    '''
    path = './results/{}/'.format(foldernames)
    energies_known = np.loadtxt(path + 'energies.dat')
    wavefuncs_known = np.loadtxt(path + 'wavefuncs.dat')
    potential_known = np.loadtxt(path + 'potential.dat')
    mass, xmin, xmax, xnum, first_val, last_val, inttype, pot, __ = reader.get_data(path)
    hamiltonian = ip.hamilton_operator(potential_known, mass, xnum)
    energies_test, wavefuncs_test = solver.qm_wavefct(hamiltonian, xnum,
                                                      first_val, last_val,
                                                      potential_known[1:xnum - 1, 0])
    assert np.all(np.abs(energies_known) - np.abs(energies_test) <= np.abs(energies_known) * 0.01)
    assert np.all(np.square(wavefuncs_known) - np.square(wavefuncs_test) <= np.square(wavefuncs_known) * 0.1 + 10 ** -15)

@pytest.mark.parametrize("foldernames", foldernames)
def test_with_reference_expvals(foldernames):
    '''
    Tests the expectation value and variance with reference data
    '''
    path = './results/{}/'.format(foldernames)
    expvalues_known = np.loadtxt(path + 'expvalues.dat')
    wavefuncs_known = np.loadtxt(path + 'wavefuncs.dat')
    expvalues_test = solver.qm_position_info(wavefuncs_known[:,1:], wavefuncs_known[:,0])
    assert np.all(np.abs(expvalues_known - expvalues_test) <= np.abs(expvalues_known) * 0.01)
import os.path
import matplotlib.pyplot as plt
import numpy as np


def read_data():
    '''
    Reads data generated by the solver, returns data as seperate arrays


    Args:


    Returns:


    '''
    path = input('Pls give Ordner: ')

    energies_path = os.path.join(path, 'energies.dat')
    potential_path = os.path.join(path, 'potential.dat')
    wavefuncs_path = os.path.join(path, 'wavefuncs.dat') 
    expvalues_path = os.path.join(path, 'expvalues.dat')

    energies = np.loadtxt(energies_path)
    potential = np.loadtxt(potential_path)
    wavefuncs = np.loadtxt(wavefuncs_path)
    expvalues = np.loadtxt(expvalues_path)

    return energies, potential, wavefuncs, expvalues
    
a, b, c, d = read_data()



def ploty():
    '''

    '''


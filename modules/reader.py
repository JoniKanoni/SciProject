'''
Module containing functions for reading input data and saving ouput data
'''
import os.path
import numpy as np




def get_data(input_path = ''):
    '''
    Reads data from an input file

    Args:
        input_path:      String that contains the path to the data that is about
                        to be read

    Returns:
        masse:          mass of the particle as float
        xmin:           first considered x-value for interpolation as float
        xmax:           last considered x-value for interpolation as float
        xnum:           Number of Points used for interpolation as int
        first, last:    first and last desired eigenvalue as int
        int_type:       Desired interpolation type as string
                        ('linear', 'natural' or 'polynomial')
        pot_points:     2D Array of known points of the potential ([:,1])
                        and corresponding x-values([:,0])
        input_path:     String that contains the path to the data that is about
                        to be read
    '''
    try:
        input_data = os.path.join(input_path, 'schrodinger.inp')
        open(input_data)
    except:
        input_path = input("Please enter the path to the input file: ")
        input_data = os.path.join(input_path, 'schrodinger.inp')
    with open(input_data, "r") as data:
        masse = float(data.readline().split()[0])
        minmax = np.array(data.readline().split()[0:3]).astype(float)
        xmin, xmax, xnum = minmax[0], minmax[1], int(minmax[2])
        first_last = np.array(data.readline().split()[0:2]).astype(int)
        first, last = first_last[0], first_last[1]
        int_type = data.readline().split()[0]
        data.readline()
        pot_points = np.loadtxt(data.readlines())

    return masse, xmin, xmax, xnum, first, last, int_type, pot_points, input_path




def save_data(save_path, potential, energies, wavefuncs, expvalues):
    '''
    Function that saves computed data of a quantummechanical system
    Saves files: potential.dat, energies.dat, wavefuncs.dat, expvalues.dat 
    at savepath location.

    Args:
        save_path:           path to save files
        potential:          array of potential values
        energies:           array of energie values
        wavefuncs:          array of wavefunction values
        expvalues:          array of expected values

    Returns:
        None
    '''
    np.savetxt(os.path.join(save_path, 'potential.dat'), potential)

    np.savetxt(os.path.join(save_path, 'energies.dat'), energies)

    np.savetxt(os.path.join(save_path, 'wavefuncs.dat'), wavefuncs)

    np.savetxt(os.path.join(save_path, 'expvalues.dat'), expvalues)

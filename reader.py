'''
Module that reads data
'''
import os.path
import numpy as np




def getdata(inputpath = ''):
    '''
    Reads data from input file
    Args:
        None yet
    Returns:
        masse:           mass of the particle as float
        xmin:           first considered x-value for interpolation as float
        xmax:           last considered x-value for interpolation as float
        xnum:           Number of Points used for interpolation
        first, last:    first and last desired eigenvalue as int
        interpolation:  Desired interpolation type as string
        numinterpol:    Number of given x-values with corresponding potential
        potential:      2D Array of the Potential ([:,1]) and corresponding
                        x-values([:,0])
    '''
    try:
        inputdata = os.path.join(inputpath, 'schrodinger.inp')
        open(inputdata)
    except:
        inputpath = input("give data: pls give path danke schön wunder bar ")
        inputdata = os.path.join(inputpath, 'schrodinger.inp')
    with open(inputdata, "r") as data:
        masse = float(data.readline().split()[0])
        minmax = np.array(data.readline().split()[0:3]).astype(float)
        xmin, xmax, xnum = minmax[0], minmax[1], int(minmax[2])
        firstlast = np.array(data.readline().split()[0:2]).astype(int)
        first, last = firstlast[0], firstlast[1]
        interpolation = data.readline().split()[0]
        numinterpol = int(data.readline().split()[0])
        pot = np.loadtxt(data.readlines())

    return masse, xmin, xmax, xnum, first, last, interpolation, numinterpol, pot, inputpath




def savedata(savepath, potential, energies, wavefuncs, expvalues):
    '''
    Args:
        savepath:           path to save files
        potential:          array of potential values
        energies:           array of energie values
        wavefuncs:          array of wavefunction values
        expvalues:          array of expected values
    Returns:
        Saved files, potential.dat, energies.dat, wavefuncs.dat, expvalues.dat at savepath location.
    '''
    np.savetxt(os.path.join(savepath, 'potential.dat'), potential)

    np.savetxt(os.path.join(savepath, 'energies.dat'), energies)

    np.savetxt(os.path.join(savepath, 'wavefuncs.dat'), wavefuncs)

    np.savetxt(os.path.join(savepath, 'expvalues.dat'), expvalues)

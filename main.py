"""
Main routine for solving a stationary 1D Schroedinger equation
"""

import reader
import interpolate as ip
import solver
#import plotter
import plotter
import numpy as np
#Deleted a bunch of useless commands and maybe some usefull here but simple ones for merge
#and to make everything simple and understandable again


def main():
    #get data
    mass, xmin, xmax, xnum, first_val, last_val, inttype, numinterpol, pot, inputpath = reader.getdata()

    #interpolate data and create hamiltonian
    potential = ip.potential_grid(xmin, xmax, xnum, pot[:, 0], pot[:, 1], inttype)
    hamiltonian = ip.hamilton_operator(potential, mass, xnum)

    #find energies/eigenvalues and wavefunctions
    energies, wavefuncs = solver.qm_wavefct(hamiltonian, xnum, first_val, last_val, potential[1:xnum - 1, 0])
    
    #find expectation value for position as well as variance
    expvalues = solver.qm_position_info(wavefuncs[:,1:], wavefuncs[:,0])
    

    reader.savedata(inputpath, potential, energies, wavefuncs, expvalues)

if __name__ == '__main__':
    main()

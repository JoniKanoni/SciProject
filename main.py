#!/usr/bin/env python3
""" Main routine for solving a stationary 1D Schroedinger equation """

from modules import reader
from modules import interpolate as ip
from modules import solver
import argparse

_DESCRIPTION = "Solves a 1D Schrodinger equation for given potential."

def main():
    '''
    Main function for solving a stationary 1D Schroedinger equation
    '''
    
    #get data
    parser = argparse.ArgumentParser(description=_DESCRIPTION)
    msg = 'Path to input file (default: .)'
    parser.add_argument('-i', '--input', type=str, default='.', help=msg)
    args = parser.parse_args()

    mass, xmin, xmax, xnum, first_val, last_val, inttype, pot, inputpath = reader.get_data(args.input)

    #interpolate data and create hamiltonian
    potential = ip.potential_grid(xmin, xmax, xnum, pot[:, 0], pot[:, 1], inttype)
    hamiltonian = ip.hamilton_operator(potential, mass, xnum)

    #find energies/eigenvalues and wavefunctions
    energies, wavefuncs = solver.qm_wavefct(hamiltonian, xnum, first_val,
                                            last_val, potential[1:xnum - 1, 0])

    #find expectation value for position as well as variance
    expvalues = solver.qm_position_info(wavefuncs[:,1:], wavefuncs[:,0])


    reader.save_data(inputpath, potential, energies, wavefuncs, expvalues)

if __name__ == '__main__':
    main()

"""
Main routine for solving a stationary 1D Schroedinger equation
"""

import reader
import interpolate as ip
import solver
import numpy as np


def main():
    mass, xmin, xmax, xnum, first_val, last_val, inttype, numinterpol, pot = reader.getdata()
    int_fct = ip.interpol(pot[:, 0], pot[:, 1], inttype)
    potential = ip.potential_grid(xmin, xmax, xnum, int_fct)
    print(np.shape(potential), potential)
    hamiltonian = solver.hamilton_operator(potential, mass, xnum)
    print(hamiltonian)
    eigenvalue, wavefuncs = (solver.diag_solver(hamiltonian, xnum, first_val, last_val))
    print(eigenvalue)
    #reader.savedata('results/infinite_well/', potential, potential, wavefuncs, wavefuncs)
    


if __name__ == '__main__':
    main()

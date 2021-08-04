"""
Main routine for solving a stationary 1D Schroedinger equation
"""

import reader
import interpolate as ip
import solver
import plotter
import numpy as np
#Deleted a bunch of useless commands and maybe some usefull here but simple ones for merge
#and to make everything simple and understandable again


def main():
    mass, xmin, xmax, xnum, first_val, last_val, inttype, numinterpol, pot = reader.getdata()
    int_fct = ip.interpol(pot[:, 0], pot[:, 1], inttype)
    potential = ip.potential_grid(xmin, xmax, xnum, int_fct)
    eigenvalue, wavefuncs = (solver.diag_solver(mass, potential, xnum, first_val, last_val))

if __name__ == '__main__':
    main()

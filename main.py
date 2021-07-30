"""
Main routine for solving a stationary 1D Schroedinger equation
"""

import reader
import interpolate as ip
import solver



def main():
    mass, xmin, xmax, xnum, first_val, last_val, inttype, numinterpol, pot = reader.getdata()
    int_fct = ip.interpol(pot[:, 0], pot[:, 1], inttype)
    potential = ip.potential_grid(xmin, xmax, xnum, int_fct)
    print(solver.diag_solver(mass, potential, xnum, first_val, last_val))
    
    


if __name__ == '__main__':
    main()

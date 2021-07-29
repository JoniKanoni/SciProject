"""
Main routine for solving a stationary 1D Schroedinger equation
"""

import reader
import interpolate as ip
import solver



def main():
    mass, xmin, xmax, xnum, first_egival, last_eigval, inttype, numinterpol, pot = reader.getdata()
    int_fct = ip.interpol(pot[:, 0], pot[:, 1], inttype)
    potential = ip.potential_grid(xmin, xmax, xnum, int_fct)
    print(solver.hamilton_operator(potential, mass, xnum))
    
    


if __name__ == '__main__':
    main()

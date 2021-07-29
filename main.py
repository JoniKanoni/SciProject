"""
Main routine for solving a stationary 1D Schroedinger equation
"""

import reader
import interpolate as ip



def main():
    mass, xmin, xmax, xnum, first_egival, last_eigval, inttype, numinterpol, pot = reader.getdata()
    print(ip.interpol(pot[:, 0], pot[:, 1], inttype)(0))


if __name__ == '__main__':
    main()
